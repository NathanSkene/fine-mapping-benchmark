"""Hand-computed posterior fixtures; no Torch/NumPy/network is required.

Run: python3 -m unittest discover -s tests/python -p 'test_*.py'
"""
import ast
import importlib.util
import itertools
import math
from pathlib import Path
import pickle
import unittest

ROOT = Path(__file__).resolve().parents[2]
MODULE = ROOT / "BEATRICE_annot_sparse" / "scripts" / "posterior_memo.py"
spec = importlib.util.spec_from_file_location("posterior_memo", MODULE)
pm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pm)


def pips(memo, p):
    return [math.fsum(value for key, value in memo.items() if j in key) for j in range(p)]


class PosteriorMemoTests(unittest.TestCase):
    def assert_probabilities(self, actual, expected):
        self.assertEqual(set(actual), set(expected))
        for key, value in expected.items():
            self.assertAlmostEqual(actual[key], value, places=13)
        self.assertAlmostEqual(math.fsum(actual.values()), 1, places=14)

    def test_zero_evidence_recovers_final_bernoulli_prior_and_null(self):
        memo = {(): 0, (0,): 0, (1,): 0, (0, 1): 0}
        context = pm.freeze_memo(memo, [0.8, 0.2])
        self.assert_probabilities(memo, {(): 0.16, (0,): 0.64, (1,): 0.04, (0, 1): 0.16})
        self.assertEqual(context["log_bf"], {(): 0, (0,): 0, (1,): 0, (0, 1): 0})
        for actual, expected in zip(pips(memo, 2), [0.8, 0.2]):
            self.assertAlmostEqual(actual, expected)

    def test_discovery_order_does_not_choose_historical_priors(self):
        # BF(singleton)=1/sqrt(2) for z=0,n=1,W=1 under independent LD.
        # Updating prior .1 -> .8 must reweight existing configurations.
        log_bf = -0.5 * math.log(2)
        memo = {(0,): log_bf}
        context = pm.freeze_memo(memo, [0.8, 0.2])
        causal_mass = 0.8 * 0.8 / math.sqrt(2)
        null_mass = 0.2 * 0.8
        self.assertAlmostEqual(memo[(0,)], causal_mass / (causal_mass + null_mass))
        opposite_order = {(): 0, (0,): log_bf}
        frozen = pm.freeze_memo(opposite_order, [0.8, 0.2])
        self.assertEqual(memo, opposite_order)
        self.assertEqual(context, frozen)
        stale_causal_mass = 0.1 * 0.8 / math.sqrt(2)
        self.assertNotAlmostEqual(memo[(0,)], stale_causal_mass / (stale_causal_mass + null_mass))

    def test_expansion_uses_same_learned_prior_as_cached_models(self):
        memo = {(0,): math.log(2)}
        context = pm.freeze_memo(memo, [0.8, 0.2])
        pm.add_configuration(memo, context, (1,), math.log(3))
        self.assert_probabilities(memo, {(): 0.16 / 1.56, (0,): 1.28 / 1.56, (1,): 0.12 / 1.56})
        pm.add_configuration(memo, context, (0, 1), math.log(4))
        self.assert_probabilities(memo, {
            (): 0.16 / 2.2, (0,): 1.28 / 2.2, (1,): 0.12 / 2.2, (0, 1): 0.64 / 2.2})
        self.assertEqual(memo, pm.posterior_weights(context))
        self.assertEqual(context["final_prior"], [0.8, 0.2])

    def test_log_sum_exp_does_not_clip_distinct_strong_configurations_into_a_tie(self):
        # The old mean-centering plus 1e15 clipping tied the two top models.
        memo = {(0,): -200, (1,): 100, (2,): 200}
        pm.freeze_memo(memo, [0.1, 0.1, 0.1])
        self.assertAlmostEqual(memo[(2,)], 1)
        self.assertAlmostEqual(memo[(1,)] / memo[(2,)], math.exp(-100), delta=1e-56)
        self.assertGreater(memo[(2,)], memo[(1,)])
        extreme = {(0,): 10000, (1,): -10000}
        pm.freeze_memo(extreme, [0.5, 0.5])
        self.assert_probabilities(extreme, {(): 0, (0,): 1, (1,): 0})

    def test_boundary_priors_preserve_impossible_configuration_probabilities(self):
        context = pm.make_context({(0,): 0, (1,): 0, (0, 1): 0}, [1, 0])
        self.assert_probabilities(pm.posterior_weights(context), {(): 0, (0,): 1, (1,): 0, (0, 1): 0})
        with self.assertRaisesRegex(ValueError, "positive probability"):
            pm.posterior_weights(pm.make_context({(1,): 0}, [1, 0]))

    def test_versioned_archive_round_trip_and_legacy_rejection(self):
        memo = {(0,): math.log(2)}
        context = pm.freeze_memo(memo, [0.8, 0.2])
        archive = {"memo_format": pm.MEMO_FORMAT, "memo": memo, "posterior_context": context}
        loaded = pickle.loads(pickle.dumps(archive))
        self.assertEqual(pm.load_context(loaded), context)
        self.assertEqual(pm.posterior_weights(pm.load_context(loaded)), memo)
        with self.assertRaisesRegex(ValueError, "regenerate"):
            pm.load_context({"memo": memo, "mean_memo": 1})
        with self.assertRaisesRegex(ValueError, "regenerate"):
            pm.posterior_weights(1)

    def test_invalid_likelihoods_priors_and_keys_are_rejected(self):
        bad = [({(0,): float("nan")}, [.5]), ({(): 1}, [.5]),
               ({(1,): 0}, [.5]), ({(0, 0): 0}, [.5]),
               ({(1, 0): 0}, [.5, .5]), ({(0,): 0}, [float("inf")]),
               ({(0,): 0}, [-.1]), ({(0,): 0}, [1.1])]
        for memo, prior in bad:
            with self.subTest(memo=memo, prior=prior), self.assertRaises(ValueError):
                pm.make_context(memo, prior)

    def test_all_training_caches_store_prior_free_likelihoods(self):
        # Check each production cache boundary, not merely the new algebra helper.
        for name in ["trainer_annot.py", "joint_trainer.py"]:
            tree = ast.parse((MODULE.parent / name).read_text())
            funcs = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)
                     and n.name in ["abf", "_abf"]]
            self.assertEqual(len(funcs), 3 if name == "trainer_annot.py" else 1)
            for fn in funcs:
                res = [n for n in ast.walk(fn) if isinstance(n, ast.Assign)
                       and any(isinstance(t, ast.Name) and t.id == "res" for t in n.targets)]
                self.assertEqual(len(res), 1)
                names = {n.id for n in ast.walk(res[0].value) if isinstance(n, ast.Name)}
                self.assertNotIn("prior", names)
                self.assertNotIn("p0", names)

    def test_frozen_archives_and_expansion_use_the_versioned_context(self):
        text = (MODULE.parent / "generate_credible_sets.py").read_text()
        self.assertIn("mean_memo = load_context(result)", text)
        self.assertIn("p0 = gpu_t(np.asarray(mean_memo['final_prior']", text)
        self.assertNotIn("p0 = gpu_t(np.array([1/len(Z)]", text)
        self.assertIn("return add_configuration(memo, mean_memo", text)
        for name in ["trainer_annot.py", "joint_trainer.py"]:
            text = (MODULE.parent / name).read_text()
            self.assertIn("'memo_format': MEMO_FORMAT", text)
            self.assertIn("'posterior_context': mean_memo", text)


if __name__ == "__main__":
    unittest.main()
