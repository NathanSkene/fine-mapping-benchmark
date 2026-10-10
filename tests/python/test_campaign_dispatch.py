"""Scheduler observation is evidence, never an optimistic restart hint."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

MODULE_PATH = Path(__file__).resolve().parents[2] / "scripts/campaign/dispatch.py"
SPEC = importlib.util.spec_from_file_location("campaign_dispatch", MODULE_PATH)
DISPATCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DISPATCH)


class DispatchTests(unittest.TestCase):
    def test_dispatch_projection_and_plan_bytes_are_frozen(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plan = root / "plan.rds"; plan.write_bytes(b"original immutable R plan")
            names = ("tasks.csv", "datasets.csv", "generation.csv")
            for name in names: (root / name).write_text("declared inputs\n")
            value = {"schema": 1, "plan_sha256": "a" * 64,
                     "plan_rds_sha256": DISPATCH.sha256(plan),
                     "files": {name: DISPATCH.sha256(root / name) for name in names}}
            (root / "dispatch-manifest.json").write_text(json.dumps(value))
            self.assertEqual(len(DISPATCH.validate_dispatch_manifest(root, plan, "a" * 64)), 64)
            (root / "tasks.csv").write_text("altered task assignments\n")
            with self.assertRaisesRegex(RuntimeError, "projection changed"):
                DISPATCH.validate_dispatch_manifest(root, plan, "a" * 64)
            plan.write_bytes(b"altered plan")
            with self.assertRaisesRegex(RuntimeError, "commitment changed"):
                DISPATCH.validate_dispatch_manifest(root, plan, "a" * 64)

    def test_live_and_completed_jobs_are_distinguished(self):
        def response(state):
            return lambda *a, **k: subprocess.CompletedProcess(a, 0, json.dumps(
                {"Jobs": {"123.pbs": {"job_state": state, "Exit_status": 0}}}), "")
        self.assertEqual(DISPATCH.observe_job("123.pbs", response("R"))["state"], "live")
        self.assertEqual(DISPATCH.observe_job("123.pbs", response("F"))["state"], "terminal")

    def test_timeout_and_unknown_network_errors_do_not_authorize_restart(self):
        def timeout(*args, **kwargs):
            raise subprocess.TimeoutExpired("qstat", 30)
        self.assertEqual(DISPATCH.observe_job("123.pbs", timeout)["state"], "observation_error")
        response = lambda *a, **k: subprocess.CompletedProcess(a, 1, "", "server unavailable")
        self.assertEqual(DISPATCH.observe_job("123.pbs", response)["state"], "observation_error")

    def test_explicit_missing_handle_is_authoritative(self):
        response = lambda *a, **k: subprocess.CompletedProcess(a, 153, "", "qstat: Unknown Job Id 123.pbs")
        self.assertEqual(DISPATCH.observe_job("123.pbs", response)["state"], "missing")
        empty = lambda *a, **k: subprocess.CompletedProcess(a, 0, '{"Jobs":{}}', "")
        self.assertEqual(DISPATCH.observe_job("123.pbs", empty)["state"], "observation_error")

    def test_terminal_receipt_commits_checkpoint_bytes_and_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "completed.rds").write_bytes(b"actual frozen checkpoint")
            value = {"schema": 1, "terminal": True, "id": "task_1", "plan_sha256": "a" * 64,
                     "checkpoint_sha256": DISPATCH.sha256(root / "completed.rds")}
            (root / "terminal.json").write_text(json.dumps(value))
            self.assertTrue(DISPATCH.terminal_receipt(root, "task_1", "a" * 64))
            with self.assertRaisesRegex(RuntimeError, "Foreign"):
                DISPATCH.terminal_receipt(root, "task_2", "a" * 64)
            (root / "completed.rds").write_bytes(b"changed checkpoint")
            with self.assertRaisesRegex(RuntimeError, "changed"):
                DISPATCH.terminal_receipt(root, "task_1", "a" * 64)

    def test_restart_revalidates_previously_cached_completed_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); plan = root / "plan.rds"
            plan.write_bytes(b"frozen RDS bytes")
            (root / "plan.sha256").write_text("a" * 64)
            bootstrap = root / "bootstrap.sh"; bootstrap.write_text("# frozen environment\n")
            (root / "generation.csv").write_text("generation_id\n1\n")
            (root / "tasks.csv").write_text(
                "index,task_id,stage,method,dataset_ids\n1,task_1,challenge,exact,dataset_1\n")
            (root / "datasets.csv").write_text(
                f"dataset_id,sim_sha256,sim_file\ndataset_1,generated,{root}/datasets/generation_00001/panel.rds\n")
            names = ("tasks.csv", "datasets.csv", "generation.csv")
            (root / "dispatch-manifest.json").write_text(json.dumps({
                "schema": 1, "plan_sha256": "a" * 64, "plan_rds_sha256": DISPATCH.sha256(plan),
                "files": {name: DISPATCH.sha256(root / name) for name in names}}))
            generation = root / "datasets/generation_00001"
            task = root / "tasks/task_1"
            for target, identity, name in ((generation, "generation_1", "generated.rds"),
                                           (task, "task_1", "completed.rds")):
                target.mkdir(parents=True); (target / name).write_bytes(b"completed bytes")
                (target / "terminal.json").write_text(json.dumps({"schema": 1, "terminal": True,
                    "id": identity, "plan_sha256": "a" * 64,
                    "checkpoint_sha256": DISPATCH.sha256(target / name)}))
            command = [sys.executable, str(MODULE_PATH), "--plan", str(plan), "--repo", str(root),
                       "--bootstrap", str(bootstrap), "--dry-run", "--once"]
            first = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            state = json.loads((root / "dispatch/state.json").read_text())
            self.assertEqual(state["complete"], {"generation": [1], "methods": [1]})
            (generation / "terminal.json").unlink()
            second = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(second.returncode, 0, second.stderr)
            state = json.loads((root / "dispatch/state.json").read_text())
            self.assertEqual(state["complete"], {"generation": [], "methods": [1]})
            self.assertIn('"would_submit": "generation"', second.stdout)
            (task / "completed.rds").write_bytes(b"altered but cached complete")
            third = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(third.returncode, 0)
            self.assertIn("Terminal checkpoint changed", third.stderr)


if __name__ == "__main__":
    unittest.main()
