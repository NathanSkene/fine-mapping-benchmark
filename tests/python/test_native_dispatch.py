"""Native orchestration must preserve panel dependence and scheduler evidence."""
import contextlib
import csv
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("native_dispatch", REPO / "scripts/server_setup/native_dispatch.py")
NATIVE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(NATIVE)
HELPERS = NATIVE.load_helpers(REPO)


class NativeDispatchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.plan = self.root / "native"
        self.preflight = self.root / "preflight"
        self.plan.mkdir(); self.preflight.mkdir()
        for root in (self.plan, self.preflight):
            (root / "plan.rds").write_bytes(b"plan")
            (root / "plan.sha256").write_text("a" * 64)
        for name in ("tasks.csv", "datasets.csv", "generation.csv"):
            (self.preflight / name).write_text("frozen preflight projection\n")
        (self.preflight / "dispatch-manifest.json").write_text(json.dumps(dict(
            schema=1, plan_sha256="a" * 64,
            plan_rds_sha256=HELPERS.sha256(self.preflight / "plan.rds"),
            files={name: HELPERS.sha256(self.preflight / name) for name in
                   ("tasks.csv", "datasets.csv", "generation.csv")})))
        self.bootstrap = self.root / "environment.sh"
        self.bootstrap.write_text("#!/bin/bash\n")
        self.write_table("draws", [dict(draw_id="draw_001")])
        self.write_table("inputs", [dict(input_id="input_00001", draw_id="draw_001")])
        self.write_table("chains", [dict(chain_id=f"chain_{i:06d}", draw_id="draw_001") for i in (1, 2)])
        self.write_table("regional", [dict(regional_id=f"regional_{i:06d}", input_id="input_00001", method=m)
                                      for i, m in enumerate(("susie", "abf"), 1)])
        (self.plan / "dispatch-manifest.json").write_text(json.dumps(dict(
            schema=1, plan_sha256="a" * 64, plan_rds_sha256=HELPERS.sha256(self.plan / "plan.rds"),
            files={name: HELPERS.sha256(self.plan / name) for name in
                   ("draws.csv", "inputs.csv", "chains.csv", "regional.csv")})))
        self.argv = ["native_dispatch", "--plan", str(self.plan / "plan.rds"), "--repo", str(REPO),
                     "--bootstrap", str(self.bootstrap), "--preflight-plan", str(self.preflight / "plan.rds"),
                     "--once"]

    def write_table(self, name, rows):
        with (self.plan / f"{name}.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)

    def receipt(self, root, folder, identity):
        directory = root / folder / identity
        directory.mkdir(parents=True, exist_ok=True)
        checkpoint = directory / "completed.rds"
        checkpoint.write_bytes(b"completed fitted artifact")
        (directory / "terminal.json").write_text(json.dumps(dict(
            schema=1, terminal=True, id=identity, plan_sha256="a" * 64,
            checkpoint_sha256=HELPERS.sha256(checkpoint))))

    def dry_run(self):
        output = io.StringIO()
        with patch.object(sys, "argv", self.argv + ["--dry-run"]), contextlib.redirect_stdout(output):
            NATIVE.main()
        return [json.loads(line) for line in output.getvalue().splitlines() if line.startswith("{")]

    def test_preparation_gates_chains_and_regional_panels(self):
        self.assertEqual([x["would_submit"] for x in self.dry_run() if "would_submit" in x], ["prepare"])
        self.receipt(self.plan, "prepared", "draw_001")
        self.assertEqual([x["would_submit"] for x in self.dry_run() if "would_submit" in x], ["chains", "chains"])
        for method in ("susie", "abf"):
            self.receipt(self.preflight, "preflight", method)
        submissions = [x for x in self.dry_run() if "would_submit" in x]
        self.assertEqual([x["would_submit"] for x in submissions], ["chains", "chains", "regional"])
        self.assertIn("--input-id input_00001", Path(submissions[-1]["script"]).read_text())

    def test_generated_jobs_use_default_or_explicit_queue(self):
        base = self.argv.copy()
        for arguments, queue in (([], None), (["--queue", ""], None),
                                 (["--queue", "approved_queue"], "approved_queue")):
            with self.subTest(queue=queue, arguments=arguments):
                self.argv = base + arguments
                submissions = [x for x in self.dry_run() if "would_submit" in x]
                self.assertEqual(len(submissions), 1)
                script = Path(submissions[0]["script"]).read_text()
                lines = [line for line in script.splitlines() if line.startswith("#PBS -q")]
                self.assertEqual(lines, [] if queue is None else [f"#PBS -q {queue}"])

    def test_new_controller_can_target_an_unchanged_worker_snapshot(self):
        worker_repo = self.root / "unchanged-workers"
        helper = worker_repo / "scripts/campaign/dispatch.py"
        helper.parent.mkdir(parents=True)
        helper.write_text((REPO / "scripts/campaign/dispatch.py").read_text())
        worker = worker_repo / "scripts/native/campaign/prepare.R"
        worker.parent.mkdir(parents=True); worker.write_text("# frozen scientific worker\n")
        self.argv[self.argv.index("--repo") + 1] = str(worker_repo)
        submissions = [x for x in self.dry_run() if "would_submit" in x]
        self.assertEqual(len(submissions), 1)
        script = Path(submissions[0]["script"]).read_text()
        self.assertIn(str(worker.resolve()), script)
        self.assertIn(f"cd {worker_repo.resolve()}", script)
        state = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(state["fingerprint"]["controller_sha256"], HELPERS.sha256(Path(NATIVE.__file__)))
        self.assertEqual(state["fingerprint"]["observation_helpers_sha256"], HELPERS.sha256(helper))
        self.assertEqual(worker.read_text(), "# frozen scientific worker\n")

    def test_partial_regional_panel_does_not_become_complete(self):
        unit = NATIVE.make_units(self.plan)["regional"]["input_00001"]
        self.receipt(self.plan, "regional", "regional_000001")
        self.assertFalse(NATIVE.unit_complete(self.plan, "regional", unit, "a" * 64, HELPERS))
        self.receipt(self.plan, "regional", "regional_000002")
        self.assertTrue(NATIVE.unit_complete(self.plan, "regional", unit, "a" * 64, HELPERS))

    def test_corrupted_csv_projection_is_rejected_before_dispatch(self):
        with (self.plan / "chains.csv").open("a") as handle:
            handle.write("foreign_chain,draw_001\n")
        with self.assertRaisesRegex(RuntimeError, "projection changed: chains.csv"):
            self.dry_run()

    def test_resume_revalidates_previously_completed_checkpoint_bytes(self):
        self.receipt(self.plan, "prepared", "draw_001")
        self.dry_run()
        (self.plan / "prepared/draw_001/completed.rds").write_bytes(b"altered valid-looking checkpoint")
        with self.assertRaisesRegex(RuntimeError, "checkpoint changed"):
            self.dry_run()

    def test_ambiguous_submission_is_durable_and_never_double_submitted(self):
        with patch.object(sys, "argv", self.argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", side_effect=subprocess.TimeoutExpired("qsub", 60)) as submit:
            with self.assertRaisesRegex(RuntimeError, "outcome unknown"):
                NATIVE.main()
            with self.assertRaisesRegex(RuntimeError, "Unresolved submission"):
                NATIVE.main()
            self.assertEqual(submit.call_count, 1)
        state = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(next(iter(state["jobs"].values()))["status"], "submission_unknown")

    def test_failed_scheduler_observation_retains_assignment_and_lock(self):
        with patch.object(sys, "argv", self.argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "123.pbs\n", "")):
            NATIVE.main()
        worker_lock = self.plan / "prepared/draw_001.lock"
        worker_lock.mkdir(parents=True)
        with patch.object(sys, "argv", self.argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=dict(state="observation_error")), \
             patch.object(NATIVE.subprocess, "run") as submit:
            NATIVE.main()
            submit.assert_not_called()
        self.assertTrue(worker_lock.exists())
        state = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(next(iter(state["jobs"].values()))["status"], "live")

    def test_terminal_scheduler_failure_allows_bounded_retry_after_lock_cleanup(self):
        with patch.object(sys, "argv", self.argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "123.pbs\n", "")):
            NATIVE.main()
        worker_lock = self.plan / "prepared/draw_001.lock"
        worker_lock.mkdir(parents=True)
        with patch.object(sys, "argv", self.argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=dict(state="terminal", exit_status=137)), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "124.pbs\n", "")) as submit:
            NATIVE.main()
            self.assertEqual(submit.call_count, 1)
        self.assertFalse(worker_lock.exists())
        state = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(state["attempts"]["prepare:draw_001"], 2)


if __name__ == "__main__":
    unittest.main()
