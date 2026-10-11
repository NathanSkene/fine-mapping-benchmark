"""Native orchestration must preserve panel dependence and scheduler evidence."""
import contextlib
import csv
import importlib.util
import io
import json
import os
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

    def test_shared_directory_claim_precedes_controller_body_and_releases(self):
        control = self.plan / "dispatch"
        owner_path = control / "controller.owner" / "owner.json"
        captured = []

        def inspect_claim(*args):
            owner = json.loads(owner_path.read_text())
            captured.append(owner)
            self.assertEqual(owner["schema"], 1)
            self.assertEqual(owner["pid"], os.getpid())
            self.assertTrue(owner["hostname"])
            self.assertTrue(owner["nonce"])
            self.assertGreater(owner["started_ns"], 0)
            self.assertFalse((control / "state.json").exists())
            with self.assertRaises(SystemExit):
                NATIVE.main()
            self.assertEqual(json.loads(owner_path.read_text()), owner)

        with patch.object(sys, "argv", self.argv), patch.object(NATIVE, "run_controller", side_effect=inspect_claim) as body:
            NATIVE.main()
            self.assertEqual(body.call_count, 1)
        self.assertEqual(len(captured), 1)
        self.assertFalse(owner_path.parent.exists())

    def test_existing_or_incomplete_controller_claim_is_never_cleared(self):
        control = self.plan / "dispatch"
        control.mkdir()
        claim = control / "controller.owner"
        for metadata in (None, dict(schema=1, hostname="other-node", pid=987654,
                                    nonce="foreign-owner", started_ns=1)):
            with self.subTest(metadata=metadata):
                claim.mkdir()
                owner_path = claim / "owner.json"
                if metadata is not None:
                    owner_path.write_text(json.dumps(metadata))
                before = {path.name: path.read_bytes() for path in claim.iterdir()}
                with patch.object(sys, "argv", self.argv), patch.object(NATIVE, "run_controller") as body, \
                     patch.object(NATIVE.subprocess, "run") as submit:
                    with self.assertRaises(SystemExit):
                        NATIVE.main()
                    body.assert_not_called(); submit.assert_not_called()
                self.assertEqual({path.name: path.read_bytes() for path in claim.iterdir()}, before)
                self.assertFalse((control / "state.json").exists())
                if owner_path.exists(): owner_path.unlink()
                claim.rmdir()

    def test_exceptional_controller_exit_releases_only_its_own_claim(self):
        claim = self.plan / "dispatch/controller.owner"
        with patch.object(sys, "argv", self.argv), \
             patch.object(NATIVE, "run_controller", side_effect=RuntimeError("test controller failure")):
            with self.assertRaisesRegex(RuntimeError, "test controller failure"):
                NATIVE.main()
        self.assertFalse(claim.exists())

    def test_mutated_controller_claim_is_preserved_on_exit(self):
        owner_path = self.plan / "dispatch/controller.owner/owner.json"

        def replace_owner(*args):
            owner = json.loads(owner_path.read_text())
            owner["nonce"] = "replacement-owner"
            owner_path.write_text(json.dumps(owner))

        with patch.object(sys, "argv", self.argv), \
             patch.object(NATIVE, "run_controller", side_effect=replace_owner):
            with self.assertRaisesRegex(RuntimeError, "ownership claim changed; preserved"):
                NATIVE.main()
        self.assertEqual(json.loads(owner_path.read_text())["nonce"], "replacement-owner")
        with patch.object(sys, "argv", self.argv), patch.object(NATIVE, "run_controller") as body:
            with self.assertRaises(SystemExit):
                NATIVE.main()
            body.assert_not_called()

    def test_owner_metadata_sync_failure_retains_claim_before_controller_body(self):
        claim = self.plan / "dispatch/controller.owner"
        with patch.object(sys, "argv", self.argv), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS.os, "fsync", side_effect=OSError("owner sync failed")), \
             patch.object(NATIVE, "run_controller") as body:
            with self.assertRaisesRegex(OSError, "owner sync failed"):
                NATIVE.main()
            body.assert_not_called()
        self.assertTrue(claim.is_dir())
        self.assertFalse((self.plan / "dispatch/state.json").exists())
        with patch.object(sys, "argv", self.argv), patch.object(NATIVE, "run_controller") as body:
            with self.assertRaises(SystemExit):
                NATIVE.main()
            body.assert_not_called()
        self.assertTrue(claim.is_dir())

    def test_worker_helpers_without_directory_locking_are_refused(self):
        with patch.object(sys, "argv", self.argv), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "exclusive_controller", None, create=True), \
             patch.object(NATIVE, "run_controller") as body:
            with self.assertRaisesRegex(SystemExit, "lack shared atomic-directory controller locking"):
                NATIVE.main()
            body.assert_not_called()
        self.assertFalse((self.plan / "dispatch/state.json").exists())

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
        self.assertEqual(state["attempts"], {})

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

    def test_interrupted_terminal_lock_cleanup_replays_without_duplicate_retry(self):
        argv = self.argv + ["--max-attempts", "2"]
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "123.pbs\n", "")):
            NATIVE.main()
        worker_lock = self.plan / "prepared/draw_001.lock"
        worker_lock.mkdir(parents=True)
        observation = dict(state="terminal", exit_status=137)
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=observation), \
             patch.object(NATIVE.shutil, "rmtree", side_effect=OSError("interrupted lock cleanup")), \
             patch.object(NATIVE.subprocess, "run") as submit:
            with self.assertRaisesRegex(OSError, "interrupted lock cleanup"):
                NATIVE.main()
            submit.assert_not_called()
        interrupted = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(interrupted["attempts"], {"prepare:draw_001": 1})
        self.assertEqual(next(iter(interrupted["jobs"].values()))["status"], "live")
        self.assertTrue(worker_lock.exists())
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=observation), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "124.pbs\n", "")) as submit:
            NATIVE.main()
            self.assertEqual(submit.call_count, 1)
        resumed = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(resumed["attempts"], {"prepare:draw_001": 2})
        self.assertEqual({job["job_id"]: job["status"] for job in resumed["jobs"].values()},
                         {"123.pbs": "terminal", "124.pbs": "live"})
        self.assertFalse(worker_lock.exists())

    def test_many_exhausted_units_bound_console_examples_and_preserve_full_report(self):
        self.write_table("draws", [dict(draw_id=f"draw_{i:03d}") for i in range(1, 26)])
        manifest_path = self.plan / "dispatch-manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["files"]["draws.csv"] = HELPERS.sha256(self.plan / "draws.csv")
        manifest_path.write_text(json.dumps(manifest))
        argv = self.argv + ["--max-attempts", "1", "--max-workers", "25", "--max-prepare", "25"]
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", side_effect=[
                 subprocess.CompletedProcess([], 0, f"{i}.pbs\n", "") for i in range(501, 526)]) as submit:
            NATIVE.main()
            self.assertEqual(submit.call_count, 25)
        output = io.StringIO()
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(output), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=dict(state="terminal", exit_status=1)), \
             patch.object(NATIVE.subprocess, "run") as submit:
            with self.assertRaisesRegex(RuntimeError, "retry budget exhausted for 25 unresolved units"):
                NATIVE.main()
            submit.assert_not_called()
        report = [json.loads(line) for line in output.getvalue().splitlines() if line.startswith("{")][-1]
        self.assertEqual(report["unresolved_count"], 25)
        self.assertEqual(len(report["unresolved_examples"]), 20)
        self.assertEqual(report["unresolved_examples"][0], "prepare:draw_001")
        self.assertEqual(report["unresolved_examples"][-1], "prepare:draw_020")
        state = json.loads(Path(report["full_report"]).read_text())
        self.assertEqual(len(state["exhausted_work"]), 25)
        self.assertEqual(state["exhausted_work"][-1]["id"], "draw_025")
        self.assertEqual(len(state["jobs"]), 25)
        self.assertTrue(all(job["status"] == "terminal" for job in state["jobs"].values()))

    def test_preworker_failures_exhaust_confirmed_submissions_and_stop_once(self):
        argv = self.argv + ["--max-attempts", "2"]
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "101.pbs\n", "")):
            NATIVE.main()
        # A bootstrap error can terminate PBS before any worker or receipt starts.
        observation = dict(state="terminal", exit_status=1, resources_used={"walltime": "00:00:01"})
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=observation), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "102.pbs\n", "")) as submit:
            NATIVE.main()
            self.assertEqual(submit.call_count, 1)
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=observation), \
             patch.object(NATIVE.subprocess, "run") as submit:
            with self.assertRaisesRegex(RuntimeError, "retry budget exhausted.*prepare:draw_001"):
                NATIVE.main()
            submit.assert_not_called()
        state = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(state["attempts"], {"prepare:draw_001": 2})
        self.assertEqual(len(state["jobs"]), 2)
        self.assertTrue(all(job["status"] == "terminal" for job in state["jobs"].values()))
        self.assertTrue(all(job["last_observation"] == observation for job in state["jobs"].values()))
        self.assertEqual(state["campaign_state"], "retry_budget_exhausted")
        self.assertEqual(state["exhausted_work"], [dict(kind="prepare", id="draw_001",
            confirmed_submissions=2, max_attempts=2, missing_receipts=["draw_001"])])
        prior_jobs = state["jobs"]
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job") as observe, \
             patch.object(NATIVE.subprocess, "run") as submit:
            with self.assertRaisesRegex(RuntimeError, "retry budget exhausted"):
                NATIVE.main()
            observe.assert_not_called(); submit.assert_not_called()
        self.assertEqual(json.loads((self.plan / "dispatch/state.json").read_text())["jobs"], prior_jobs)

    def test_exhausted_chain_stops_even_while_another_assignment_is_live(self):
        self.receipt(self.plan, "prepared", "draw_001")
        argv = [x for x in self.argv if x != "--once"] + ["--once", "--max-attempts", "1"]
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", side_effect=[
                 subprocess.CompletedProcess([], 0, f"{i}.pbs\n", "") for i in (201, 202)]):
            NATIVE.main()
        locks = [self.plan / f"chains/chain_{i:06d}.lock" for i in (1, 2)]
        for lock in locks: lock.mkdir(parents=True)
        output = io.StringIO()
        with patch.object(sys, "argv", [x for x in argv if x != "--once"]), contextlib.redirect_stdout(output), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", side_effect=lambda job: dict(
                 state="terminal", exit_status=137) if job == "201.pbs" else dict(state="live", job_state="R")), \
             patch.object(NATIVE.subprocess, "run") as submit, patch.object(NATIVE.time, "sleep") as sleep:
            with self.assertRaisesRegex(RuntimeError, "retry budget exhausted.*chains:chain_000001"):
                NATIVE.main()
            submit.assert_not_called(); sleep.assert_not_called()
        state = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(len(state["jobs"]), 2)
        self.assertEqual({job["job_id"]: job["status"] for job in state["jobs"].values()},
                         {"201.pbs": "terminal", "202.pbs": "live"})
        self.assertFalse(locks[0].exists()); self.assertTrue(locks[1].exists())
        report = [json.loads(line) for line in output.getvalue().splitlines() if line.startswith("{")][-1]
        self.assertEqual(report["active"]["chains"], 1)
        self.assertEqual(report["unresolved_examples"], ["chains:chain_000001"])

    def test_receipt_completion_overrides_a_failed_terminal_handle(self):
        argv = self.argv + ["--max-attempts", "1"]
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "301.pbs\n", "")):
            NATIVE.main()
        self.receipt(self.plan, "prepared", "draw_001")
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=dict(state="terminal", exit_status=1)), \
             patch.object(NATIVE.subprocess, "run", side_effect=[
                 subprocess.CompletedProcess([], 0, f"{i}.pbs\n", "") for i in (302, 303)]) as submit:
            NATIVE.main()
            self.assertEqual(submit.call_count, 2)
        state = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(state["complete"]["prepare"], ["draw_001"])
        self.assertEqual(state["exhausted_work"], [])

    def test_partial_regional_exhaustion_lists_only_missing_receipts(self):
        self.receipt(self.plan, "prepared", "draw_001")
        for i in (1, 2): self.receipt(self.plan, "chains", f"chain_{i:06d}")
        for method in ("susie", "abf"): self.receipt(self.preflight, "preflight", method)
        argv = self.argv + ["--max-attempts", "1"]
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE.subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "401.pbs\n", "")):
            NATIVE.main()
        self.receipt(self.plan, "regional", "regional_000001")
        directory = self.plan / "regional/regional_000001"
        before = {p.name: p.read_bytes() for p in directory.iterdir()}
        with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
             patch.object(NATIVE, "load_helpers", return_value=HELPERS), \
             patch.object(HELPERS, "observe_job", return_value=dict(state="terminal", exit_status=137)), \
             patch.object(NATIVE.subprocess, "run") as submit:
            with self.assertRaisesRegex(RuntimeError, "retry budget exhausted.*regional:input_00001"):
                NATIVE.main()
            submit.assert_not_called()
        state = json.loads((self.plan / "dispatch/state.json").read_text())
        self.assertEqual(state["exhausted_work"][0]["missing_receipts"], ["regional_000002"])
        self.assertEqual(state["complete"]["regional"], [])
        self.assertEqual({p.name: p.read_bytes() for p in directory.iterdir()}, before)


if __name__ == "__main__":
    unittest.main()
