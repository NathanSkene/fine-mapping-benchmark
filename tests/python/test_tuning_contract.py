"""Scientific contracts for tuning, without a neural fit or Optuna campaign."""
import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("fb_tuning", ROOT / "scripts/tuning/optuna_fb.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def args(**updates):
    values = dict(sim="unused.rds", objective_script="scripts/tuning/fb_objective.R",
                  regions=0, scenarios=2, rscript="Rscript", beatrice_dir="BEATRICE_annot_sparse",
                  python="python", max_iter=2000, rotate_regions=False, trial_timeout=None,
                  objective="multi", multi_second="mass", penalty=.5, mass_penalty=.5,
                  tune_gamma_coverage=False, abort_quantile=0, abort_after=3,
                  abort_min_trials=30, startup_trials=20, seed=17,
                  fix=dict(lambda_l1=.01, prior_regularisation=1., sigma_sq=.05,
                           hierarchy_M=10., n_caus=5, sparse_concrete=50))
    values.update(updates)
    return SimpleNamespace(**values)


class Study:
    study_name = "test"
    def __init__(self, trials=None):
        self.user_attrs = {}
        self.trials = trials or []
    def get_trials(self, deepcopy=False):
        return self.trials
    def set_user_attr(self, key, value):
        self.user_attrs[key] = value


class Trial:
    number = 0
    def __init__(self):
        self.user_attrs = {}
    def set_user_attr(self, key, value):
        self.user_attrs[key] = value
    def report(self, value, step):
        raise AssertionError("Pareto trial must not use a scalar pruner")


class TuningContractTests(unittest.TestCase):
    def test_runtime_diagnostics_use_tie_aware_parameter_ranks(self):
        self.assertEqual(module._rank([1, 1, 3]), [.5, .5, 2.])
        self.assertEqual(module._rank([2, 1, 2]), [1.5, 0., 1.5])

    def test_partial_fit_trial_is_not_silently_promoted(self):
        proc = SimpleNamespace(returncode=0, stderr="", stdout="status=ok\nap=0.99\nn_failed=1\n")
        with patch.object(module.subprocess, "run", return_value=proc):
            with self.assertRaisesRegex(RuntimeError, "partial-success"):
                module.run_one_scenario(args(), {}, 1)
        proc.stdout = "status=ok\nap=0.99\n"
        with patch.object(module.subprocess, "run", return_value=proc):
            with self.assertRaisesRegex(RuntimeError, "failure count"):
                module.run_one_scenario(args(), {}, 1)

    def test_subprocess_exit_failure_overrides_success_text(self):
        proc = SimpleNamespace(returncode=1, stderr="crashed", stdout="status=ok\nn_failed=0\nap=1\n")
        with patch.object(module.subprocess, "run", return_value=proc):
            with self.assertRaisesRegex(RuntimeError, "crashed"):
                module.run_one_scenario(args(), {}, 1)

    def test_preflight_requires_full_grid_and_balanced_panel(self):
        proc = SimpleNamespace(returncode=0, stderr="", stdout="status=ok\nn_scenarios=4\nn_regions=3\n")
        with patch.object(module.subprocess, "run", return_value=proc):
            with self.assertRaisesRegex(ValueError, "complete cached grid"):
                module.preflight(args())
            self.assertEqual(module.preflight(args(scenarios=4))["n_regions"], "3")
        with self.assertRaisesRegex(ValueError, "aliases"):
            module.preflight(args(rotate_regions=True))
        with self.assertRaisesRegex(ValueError, "tradeoffs"):
            module.preflight(args(abort_quantile=.25))

    def test_contract_rejects_changed_or_unfingerprinted_history(self):
        with tempfile.TemporaryDirectory() as directory:
            storage = str(Path(directory) / "journal")
            study = Study()
            module.bind_study_contract(study, {"data": "A"}, storage)
            module.bind_study_contract(study, {"data": "A"}, storage)
            with self.assertRaisesRegex(ValueError, "changed"):
                module.bind_study_contract(study, {"data": "B"}, storage)
            competing = Study()
            with self.assertRaisesRegex(ValueError, "sidecar"):
                module.bind_study_contract(competing, {"data": "B"}, storage)
            with self.assertRaisesRegex(ValueError, "without a benchmark"):
                module.bind_study_contract(Study([object()]), {"data": "A"}, storage)

    def test_dataset_and_source_changes_change_fingerprint(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            tuning = repo / "scripts/tuning"
            tuning.mkdir(parents=True)
            (repo / "R").mkdir()
            beatrice = repo / "BEATRICE"
            beatrice.mkdir()
            objective = tuning / "fb_objective.R"
            objective.write_text("initial")
            (tuning / "tuning_contract.R").write_text("helper")
            data = repo / "sim.rds"
            data.write_bytes(b"dataset A")
            settings = args(sim=str(data), objective_script=str(objective), beatrice_dir=str(beatrice))
            original = module.contract_digest(module.make_contract(settings))
            data.write_bytes(b"dataset B")
            changed_data = module.contract_digest(module.make_contract(settings))
            self.assertNotEqual(original, changed_data)
            objective.write_text("revised scoring")
            self.assertNotEqual(changed_data, module.contract_digest(module.make_contract(settings)))

    def test_mass_errors_do_not_cancel_and_pareto_trials_are_not_pruned(self):
        trial = Trial()
        metrics = [dict(ap=.8, total_mass_ratio=.5, max_fdr_violation_n20=0.,
                        hi_pip_reliab=1., ece_hi=.1, fdr_at_95=0., seconds=0.),
                   dict(ap=.6, total_mass_ratio=1.5, max_fdr_violation_n20=0.,
                        hi_pip_reliab=1., ece_hi=.1, fdr_at_95=0., seconds=0.)]
        with patch.object(module, "run_one_scenario", side_effect=metrics):
            score = module.build_objective(args())(trial)
        self.assertAlmostEqual(score[0], .7)
        self.assertAlmostEqual(score[1], .5)
        self.assertEqual(trial.user_attrs["mean_total_mass_ratio"], 1.)
        self.assertEqual(trial.user_attrs["mean_absolute_mass_ratio_error"], .5)


if __name__ == "__main__":
    unittest.main()
