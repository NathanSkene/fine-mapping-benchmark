"""Actual learned-prior training and analytic posterior checks on tiny loci.

Opt in with FMB_NEURAL_SMOKE=1 in an isolated environment with the bundled
BEATRICE dependencies. The default dependency-free suite skips these tests.
"""
import math
import os
from pathlib import Path
import pickle
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]


@unittest.skipUnless(os.getenv("FMB_NEURAL_SMOKE") == "1", "optional Torch training smoke")
class NeuralPosteriorSmoke(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(ROOT / "BEATRICE_annot_sparse"))
        import numpy as np
        import torch
        from beatrice_annot import run_beatrice_annot
        from scripts.joint_trainer import LinearPrior, run_joint
        from scripts.trainer_annot import LassoNetPrior
        cls.np, cls.torch = np, torch
        cls.run_single = staticmethod(run_beatrice_annot)
        cls.run_joint = staticmethod(run_joint)
        cls.LinearPrior, cls.LassoNetPrior = LinearPrior, LassoNetPrior
        torch.set_num_threads(1)

    def fixture(self, directory, name, offset=0):
        np = self.np
        z = np.array([6, .2, -.3, 4, 0, .1, -.4, .3], dtype=float)
        if offset:
            z = np.roll(z, offset)
        v = np.column_stack([np.arange(len(z)) % 2, np.linspace(-1, 1, len(z))])
        target = Path(directory) / name
        target.mkdir()
        zp, ldp, ap = [target / n for n in ["data.z", "data.ld", "data.annot"]]
        names = ["v" + str(i) for i in range(len(z))]
        zp.write_text("\n".join(f"{n} {x}" for n, x in zip(names, z)) + "\n")
        np.savetxt(ldp, np.eye(len(z)))
        ap.write_text("\n".join(f"{n} {x[0]} {x[1]}" for n, x in zip(names, v)) + "\n")
        return {"z": str(zp), "LD": str(ldp), "annot": str(ap),
                "target": str(target), "N": 1000}, z, v

    def controlled_training(self, kind, n_samples=1, bad_sample=False):
        """Run the production training path with repeatable soft samples."""
        torch = self.torch
        from scripts import trainer_annot, joint_trainer

        class Posterior(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.theta = torch.nn.Parameter(torch.tensor(0.))

            def forward(self, temperature, samples):
                c = torch.tensor([[.8, .6, .1, .2]]).repeat(samples, 1) + self.theta * 0
                if bad_sample:
                    c = c * float("nan")
                imp = torch.sigmoid(self.theta).expand(4)
                return c, c[:1], c[:1], imp

        class Prior(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.theta = torch.nn.Parameter(torch.tensor(0.))

            def forward(self, v, temperature, samples):
                return None, None, None, torch.sigmoid(self.theta).expand(4)

            def proximal_l1(self, *args):
                pass

            def apply_hierarchy_constraint(self):
                pass

            def l1_penalty(self):
                return self.theta * 0

        model, prior = Posterior(), Prior()
        opt = torch.optim.SGD(list(model.parameters()) + list(prior.parameters()), lr=0)
        z, ld, v = torch.tensor([2., -1., .1, .4]), torch.eye(4), torch.ones(4, 1)
        temperature = torch.tensor(1.)
        if kind == "joint":
            run = lambda: float(joint_trainer._region_elbo(
                model, prior, z.unsqueeze(1), z, ld, v, temperature, n_samples,
                torch.tensor(.05), 10, 2, 0, 0, {})[0].detach())
        else:
            cls = getattr(trainer_annot, kind)
            if kind == "finemapper":
                mapper = cls(model, opt, None)
                run = lambda: float(mapper.train(z, ld, temperature, n_samples, .05,
                                                10, torch.full((4,), .5), 1, {}, 0, 2, 0)[0][0])
            else:
                args = [model, prior, opt, None] + ([0, 0] if kind == "finemapper_lassonet" else [0])
                mapper = cls(*args)
                run = lambda: float(mapper.train(z, ld, v, temperature, n_samples,
                                                .05, 10, 1, {}, 0, 2, 0)[0][0])
            mapper.S = z.unsqueeze(1)
        return run, opt

    def test_rejected_likelihoods_do_not_dilute_finite_sample_mean(self):
        torch = self.torch
        for kind in ["finemapper", "finemapper_annot", "finemapper_lassonet", "joint"]:
            with self.subTest(kind=kind):
                baseline, _ = self.controlled_training(kind)
                expected = baseline()
                mixed, _ = self.controlled_training(kind, n_samples=2)
                real_logdet = torch.logdet
                calls = 0

                def first_invalid(matrix):
                    nonlocal calls
                    calls += 1
                    return torch.tensor(float("nan")) if calls == 1 else real_logdet(matrix)

                with mock.patch.object(torch, "logdet", side_effect=first_invalid):
                    self.assertAlmostEqual(mixed(), expected, places=6)

    def test_no_finite_likelihood_fails_before_optimizer_step(self):
        for kind in ["finemapper", "finemapper_annot", "finemapper_lassonet", "joint"]:
            with self.subTest(kind=kind):
                run, opt = self.controlled_training(kind, n_samples=2)
                with mock.patch.object(self.torch, "logdet", return_value=self.torch.tensor(float("nan"))), \
                        mock.patch.object(opt, "step") as step:
                    with self.assertRaisesRegex(FloatingPointError, "no finite likelihood among 2 samples"):
                        run()
                    step.assert_not_called()

    def test_nonfinite_soft_samples_fail_before_matrix_calculation(self):
        for kind in ["finemapper", "finemapper_annot", "finemapper_lassonet", "joint"]:
            with self.subTest(kind=kind):
                run, _ = self.controlled_training(kind, n_samples=2, bad_sample=True)
                with mock.patch.object(self.torch, "inverse") as inverse:
                    with self.assertRaisesRegex(FloatingPointError, "no finite likelihood among 2 samples"):
                        run()
                    inverse.assert_not_called()

    def assert_posterior(self, region, z, v, prior):
        np, torch = self.np, self.torch
        target = Path(region["target"])
        with (target / "res").open("rb") as f:
            result = pickle.load(f)
        context = result["posterior_context"]
        self.assertEqual(result["memo_format"], "beatrice-log-bf-final-prior-v1")
        self.assertIn((), context["log_bf"])
        self.assertGreater(len(context["log_bf"]), 2)
        with torch.no_grad():
            prior.eval()
            _, _, _, expected_prior = prior(torch.tensor(v, dtype=torch.float32), torch.tensor(1.), 1)
        final_prior = np.array(context["final_prior"])
        np.testing.assert_allclose(final_prior, expected_prior.numpy(), rtol=1e-6, atol=1e-7)
        self.assertTrue(np.all((final_prior > 0) & (final_prior < 1)))
        self.assertGreater(np.ptp(final_prior), 1e-5)
        # Closed-form independent-LD Gaussian BF, followed by explicit Bernoulli
        # configuration probabilities. This does not use posterior_memo.
        n_w = region["N"] * .05
        terms = -.5 * math.log1p(n_w) + .5 * n_w / (1 + n_w) * z**2
        log_weights = {}
        for config, cached_bf in context["log_bf"].items():
            analytic_bf = math.fsum(float(terms[j]) for j in config)
            self.assertAlmostEqual(cached_bf, analytic_bf, delta=2e-5)
            log_prior = math.fsum(math.log(final_prior[j] if j in config else 1 - final_prior[j])
                                 for j in range(len(z)))
            log_weights[config] = analytic_bf + log_prior
        centre = max(log_weights.values())
        weights = {k: math.exp(x - centre) for k, x in log_weights.items()}
        normalizer = math.fsum(weights.values())
        weights = {k: x / normalizer for k, x in weights.items()}
        expected_pip = [math.fsum(x for k, x in weights.items() if j in k) for j in range(len(z))]
        for config, weight in weights.items():
            self.assertAlmostEqual(result["memo"][config], weight, delta=5e-6)
        self.assertAlmostEqual(math.fsum(result["memo"].values()), 1, places=12)
        np.testing.assert_allclose(result["pip"], expected_pip, rtol=1e-5, atol=1e-7)
        exported = np.genfromtxt(target / "pip.csv", delimiter=",", names=True,
                                 dtype=None, encoding="utf-8")
        np.testing.assert_allclose(exported["pip"], expected_pip, rtol=1e-5, atol=1e-7)
        lines = (target / "credible_set.txt").read_text().splitlines()
        sets = [tuple(sorted(map(int, line.split()))) for line in lines if line.strip()]
        self.assertEqual(len(sets), len(set(sets)))

    def test_single_lassonet_training_expansion_and_checkpoint_agree(self):
        with tempfile.TemporaryDirectory(prefix="fmb-neural-") as directory:
            region, z, v = self.fixture(directory, "single")
            self.run_single({**region, "max_iter": 4, "neural_network": [4, 8, 4],
                             "prior_neural_network": [4, 4], "sparse_concrete": 2,
                             "n_caus": 2, "plot_loss": False, "return_weights": True})
            checkpoint = self.torch.load(Path(region["target"]) / "prior_network_weights.pt",
                                         map_location="cpu", weights_only=False)
            prior = self.LassoNetPrior(m=2, hidden_dims=[4, 4], M=10)
            prior.load_state_dict(checkpoint["state_dict"])
            self.assert_posterior(region, z, v, prior)

    def test_joint_shared_linear_training_expansion_and_checkpoint_agree(self):
        with tempfile.TemporaryDirectory(prefix="fmb-neural-") as directory:
            a, za, va = self.fixture(directory, "joint_a")
            b, zb, vb = self.fixture(directory, "joint_b", offset=2)
            path = str(Path(directory) / "shared_prior.pt")
            self.run_joint({"regions": [a, b], "prior_head": "linear", "NN": [4, 8, 4],
                            "prior_neural_network": [4, 4], "max_iter": 4, "sigma_sq": .05,
                            "sparse_concrete": 2, "MCMC_samples": 1, "temp_lower_bound": .01,
                            "gamma": .1, "prior_regularisation": 1, "n_causal": 2,
                            "coverage_ths": .95, "selection_prob": .05, "key_thres": .2,
                            "allow_duplicates": False, "save_prior_weights": path})
            checkpoint = self.torch.load(path, map_location="cpu", weights_only=False)
            prior = self.LinearPrior(m=2)
            prior.load_state_dict(checkpoint["state_dict"])
            self.assert_posterior(a, za, va, prior)
            self.assert_posterior(b, zb, vb, prior)


if __name__ == "__main__":
    unittest.main()
