"""Scheduler observation is evidence, never an optimistic restart hint."""
import importlib.util
import contextlib
import io
import json
from pathlib import Path
import subprocess
import sqlite3
import sys
import tempfile
import time
import unittest
from unittest import mock

MODULE_PATH = Path(__file__).resolve().parents[2] / "scripts/campaign/dispatch.py"
SPEC = importlib.util.spec_from_file_location("campaign_dispatch", MODULE_PATH)
DISPATCH = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DISPATCH)


@contextlib.contextmanager
def sqlite_connection(path):
    connection = sqlite3.connect(path)
    try:
        with connection:
            yield connection
    finally:
        connection.close()


class DispatchTests(unittest.TestCase):
    def priority_fixture(self, root):
        plan = root / "plan.rds"; plan.write_bytes(b"immutable full plan")
        (root / "plan.sha256").write_text("a" * 64)
        bootstrap = root / "bootstrap.sh"; bootstrap.write_text("# frozen\n")
        (root / "generation.csv").write_text("generation_id\n1\n2\n3\n4\n")
        (root / "tasks.csv").write_text("index,task_id,stage,method,dataset_ids\n")
        (root / "datasets.csv").write_text("dataset_id,sim_sha256,sim_file\n")
        names = ("tasks.csv", "datasets.csv", "generation.csv")
        manifest = root / "dispatch-manifest.json"
        manifest.write_text(json.dumps({"schema": 1, "plan_sha256": "a" * 64,
            "plan_rds_sha256": DISPATCH.sha256(plan),
            "files": {name: DISPATCH.sha256(root / name) for name in names}}))
        priority = root / "priority.json"
        priority.write_text(json.dumps({"schema": 1, "plan_sha256": "a" * 64,
            "plan_rds_sha256": DISPATCH.sha256(plan),
            "dispatch_manifest_sha256": DISPATCH.sha256(manifest),
            "method_task_indices": [], "generation_indices": [3, 1]}))
        command = [sys.executable, str(MODULE_PATH), "--plan", str(plan), "--repo", str(root),
                   "--bootstrap", str(bootstrap), "--once", "--max-workers", "3",
                   "--max-generation", "3", "--priority-manifest", str(priority)]
        return plan, manifest, priority, command

    def test_priority_manifest_rejects_foreign_invalid_and_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); plan, manifest, path, _ = self.priority_fixture(root)
            good = json.loads(path.read_text())
            validate = lambda: DISPATCH.validate_priority_manifest(
                path, "a" * 64, {1, 2}, {1, 2, 3, 4}, plan, DISPATCH.sha256(manifest))
            self.assertEqual(validate()[:2], ([], [3, 1]))
            for name, value in (("schema", True), ("plan_sha256", "b" * 64),
                                ("method_task_indices", [1, 1]),
                                ("method_task_indices", [3]), ("generation_indices", [0]),
                                ("generation_indices", [True]), ("generation_indices", [1.0]),
                                ("generation_indices", ["1"]),
                                ("plan_rds_sha256", "b" * 64),
                                ("dispatch_manifest_sha256", "b" * 64)):
                with self.subTest(name=name, value=value):
                    path.write_text(json.dumps(dict(good, **{name: value})))
                    with self.assertRaises(RuntimeError): validate()

    def test_priority_is_a_prefix_with_unchanged_complete_remainder(self):
        self.assertEqual(DISPATCH.scheduling_order([1, 2, 3, 4], [3, 1]), [3, 1, 2, 4])
        self.assertEqual(DISPATCH.scheduling_order([1, 2, 3, 4], []), [1, 2, 3, 4])
        tasks = {i: {"method": "carma" if i < 4 else "susie", "stage": "challenge"}
                 for i in range(1, 7)}
        self.assertEqual(DISPATCH.priority_batches([3, 1, 2, 6, 4], tasks,
            lambda kind, index: True, 4, 3), [[3, 1], [2], [6, 4]])
        # A blocked prerequisite remains blocked and does not prevent ready work.
        self.assertEqual(DISPATCH.priority_batches([3, 1, 2, 6, 4], tasks,
            lambda kind, index: index not in {3, 2}, 4, 2), [[1], [6, 4]])
        self.assertEqual(DISPATCH.priority_batches([1], tasks,
            lambda kind, index: True, 4, 0), [])

    def test_priority_generation_order_hash_resume_and_full_remainder(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); _, _, path, command = self.priority_fixture(root)
            result = subprocess.run(command + ["--dry-run"], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            jobs = [json.loads(x) for x in result.stdout.splitlines() if '"would_submit"' in x]
            self.assertEqual([x["indices"] for x in jobs], [[3], [1], [2]])
            state = json.loads((root / "dispatch/state.json").read_text())
            self.assertEqual(state["priority_manifest_sha256"], DISPATCH.sha256(path))
            value = json.loads(path.read_text()); value["generation_indices"] = [1, 3]
            path.write_text(json.dumps(value))
            changed = subprocess.run(command + ["--dry-run"], text=True, capture_output=True)
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn("priority changed", changed.stderr)
            # Removing an already committed schedule is also incompatible.
            removed = command[:-2]
            changed = subprocess.run(removed + ["--dry-run"], text=True, capture_output=True)
            self.assertNotEqual(changed.returncode, 0)

    def test_live_priority_tamper_stops_before_second_submission(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); _, _, path, command = self.priority_fixture(root)
            calls = []
            def submit(arguments, **kwargs):
                calls.append(arguments)
                path.write_text(path.read_text() + "\n")
                return subprocess.CompletedProcess(arguments, 0, "123.pbs\n", "")
            opened = []; original_open = Path.open
            def tracked_open(path, *arguments, **kwargs):
                handle = original_open(path, *arguments, **kwargs)
                if path.name == "controller.lock": opened.append(handle)
                return handle
            with mock.patch.object(sys, "argv", command[1:]), \
                    mock.patch.object(Path, "open", new=tracked_open), \
                    mock.patch.object(DISPATCH.subprocess, "run", side_effect=submit), \
                    contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(RuntimeError, "changed during dispatch"):
                    DISPATCH.main()
            for handle in opened: handle.close()
            self.assertEqual(len(calls), 1)
            state = json.loads((root / "dispatch/state.json").read_text())
            self.assertEqual(next(iter(state["jobs"].values()))["job_id"], "123.pbs")

    def test_ready_priority_methods_precede_remaining_tasks_without_bypassing_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); plan, manifest, path, command = self.priority_fixture(root)
            (root / "tasks.csv").write_text("index,task_id,stage,method,dataset_ids\n"
                "1,t1,challenge,exact,d1\n2,t2,challenge,forced_one,d2\n"
                "3,t3,original_refit,exact,d3\n4,t4,corrected_refit,exact,d3\n"
                "5,t5,original_refit,forced_one,d5\n6,t6,corrected_refit,forced_one,d5\n")
            (root / "datasets.csv").write_text("dataset_id,sim_sha256,sim_file\n"
                "d1,frozen,path\nd2,frozen,path\nd3,frozen,path\nd5,frozen,path\n")
            value = json.loads(manifest.read_text())
            for name in value["files"]: value["files"][name] = DISPATCH.sha256(root / name)
            manifest.write_text(json.dumps(value))
            value = json.loads(path.read_text()); value.update(method_task_indices=[4, 5, 2],
                generation_indices=[], dispatch_manifest_sha256=DISPATCH.sha256(manifest))
            path.write_text(json.dumps(value))
            for i in (1, 2, 3, 4):
                target = root / f"datasets/generation_{i:05d}"; target.mkdir(parents=True)
                checkpoint = target / "generated.rds"; checkpoint.write_bytes(b"sealed generation")
                (target / "terminal.json").write_text(json.dumps({"schema": 1, "terminal": True,
                    "id": f"generation_{i}", "plan_sha256": "a" * 64,
                    "checkpoint_sha256": DISPATCH.sha256(checkpoint)}))
            command[command.index("--max-workers") + 1] = "4"
            result = subprocess.run(command + ["--dry-run"], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            jobs = [json.loads(x) for x in result.stdout.splitlines() if '"would_submit"' in x]
            self.assertEqual([x["indices"] for x in jobs], [[5], [2], [1], [3]])

    def test_generated_jobs_use_default_or_explicit_queue_and_frozen_workers(self):
        for arguments, queue in (([], None), (["--queue", ""], None),
                                 (["--queue", "approved_queue"], "approved_queue")):
            with self.subTest(queue=queue, arguments=arguments), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); campaign = root / "campaign"; campaign.mkdir()
                worker_repo = root / "unchanged-workers"; worker_repo.mkdir()
                worker = worker_repo / "scripts/campaign/generate.R"
                worker.parent.mkdir(parents=True); worker.write_text("# frozen scientific worker\n")
                plan = campaign / "plan.rds"; plan.write_bytes(b"immutable scientific plan")
                (campaign / "plan.sha256").write_text("a" * 64)
                bootstrap = worker_repo / "bootstrap.sh"; bootstrap.write_text("# frozen environment\n")
                (campaign / "generation.csv").write_text("generation_id\n1\n")
                (campaign / "tasks.csv").write_text("index,task_id,stage,method,dataset_ids\n")
                (campaign / "datasets.csv").write_text("dataset_id,sim_sha256,sim_file\n")
                names = ("tasks.csv", "datasets.csv", "generation.csv")
                (campaign / "dispatch-manifest.json").write_text(json.dumps({
                    "schema": 1, "plan_sha256": "a" * 64,
                    "plan_rds_sha256": DISPATCH.sha256(plan),
                    "files": {name: DISPATCH.sha256(campaign / name) for name in names}}))
                command = [sys.executable, str(MODULE_PATH), "--plan", str(plan),
                           "--repo", str(worker_repo), "--bootstrap", str(bootstrap),
                           "--dry-run", "--once", *arguments]
                result = subprocess.run(command, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                submissions = [json.loads(line) for line in result.stdout.splitlines()
                               if line.startswith("{") and "would_submit" in json.loads(line)]
                self.assertEqual(len(submissions), 1)
                script = Path(submissions[0]["script"]).read_text()
                queue_lines = [line for line in script.splitlines() if line.startswith("#PBS -q")]
                self.assertEqual(queue_lines, [] if queue is None else [f"#PBS -q {queue}"])
                self.assertIn(str(worker.resolve()), script)
                self.assertIn(f"cd {worker_repo.resolve()}", script)
                self.assertIn(f"source {bootstrap.resolve()}", script)
                state = json.loads((campaign / "dispatch/state.json").read_text())
                self.assertEqual(state["controller_sha256"], DISPATCH.sha256(MODULE_PATH))
                self.assertEqual(worker.read_text(), "# frozen scientific worker\n")

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
            self.assertEqual(state["complete_counts"], {"generation": 1, "methods": 1})
            (generation / "terminal.json").unlink()
            second = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(second.returncode, 0, second.stderr)
            state = json.loads((root / "dispatch/state.json").read_text())
            self.assertEqual(state["complete_counts"], {"generation": 0, "methods": 1})
            self.assertIn('"would_submit": "generation"', second.stdout)
            (task / "completed.rds").write_bytes(b"altered but cached complete")
            third = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(third.returncode, 0)
            self.assertIn("Terminal checkpoint changed", third.stderr)

    def run_main(self, command, submit, observation=None):
        with mock.patch.object(sys,"argv",command[1:]), \
                mock.patch.object(DISPATCH.subprocess,"run",side_effect=submit), \
                contextlib.redirect_stdout(io.StringIO()):
            if observation is None:
                DISPATCH.main()
            else:
                with mock.patch.object(DISPATCH,"observe_job",return_value=observation):
                    DISPATCH.main()

    def test_submission_intent_is_transactional_and_unknown_result_stays_assigned(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,_,_,command=self.priority_fixture(root)
            calls=[]
            def uncertain(arguments,**kwargs):
                calls.append(arguments)
                # A second SQLite reader sees committed intent BEFORE qsub.
                with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                    token,status=db.execute("SELECT token,status FROM jobs").fetchone()
                    self.assertEqual(status,"submission_pending")
                    self.assertEqual(db.execute("SELECT active_job FROM units WHERE kind='generation' AND unit_index=3").fetchone()[0],token)
                    self.assertEqual(db.execute("SELECT event FROM events WHERE token=?",(token,)).fetchone()[0],"submission_intent")
                raise subprocess.TimeoutExpired(arguments,60)
            with self.assertRaisesRegex(RuntimeError,"outcome unknown"):
                self.run_main(command,uncertain)
            self.assertEqual(len(calls),1)
            state=json.loads((root/"dispatch/state.json").read_text())
            self.assertEqual(next(iter(state["jobs"].values()))["status"],"submission_unknown")
            with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                self.assertEqual(db.execute("SELECT attempts FROM units WHERE kind='generation' AND unit_index=3").fetchone()[0],0)
                self.assertEqual(db.execute("SELECT COUNT(*) FROM jobs").fetchone()[0],1)
            restart=subprocess.run(command,text=True,capture_output=True)
            self.assertNotEqual(restart.returncode,0)
            self.assertIn("Unresolved submission",restart.stderr)
            with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                self.assertEqual(db.execute("SELECT COUNT(*) FROM jobs").fetchone()[0],1)

    def test_observation_errors_and_new_missing_ids_never_release_assignments(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,_,_,command=self.priority_fixture(root)
            command[command.index("--max-workers")+1]="1"
            good=lambda a,**k:subprocess.CompletedProcess(a,0,"123.pbs\n","")
            self.run_main(command,good)
            lock=root/"datasets/generation_00003.lock"; lock.mkdir(parents=True)
            for observation in (dict(state="observation_error",detail="network error"),dict(state="missing")):
                with self.subTest(observation=observation):
                    self.run_main(command,lambda *a,**k:self.fail("Duplicate qsub"),observation)
                    self.assertTrue(lock.exists())
                    with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                        self.assertEqual(db.execute("SELECT status FROM jobs").fetchone()[0],"live")
                        self.assertEqual(db.execute("SELECT attempts FROM units WHERE kind='generation' AND unit_index=3").fetchone()[0],1)
            with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                last=json.loads(db.execute("SELECT payload FROM events WHERE event='scheduler_observation' ORDER BY sequence DESC LIMIT 1").fetchone()[0])
                self.assertEqual(last["raw"]["state"],"missing")
                self.assertEqual(last["effective"]["state"],"observation_error")

    def test_authoritative_terminal_retry_keeps_history_and_cleans_lock_first(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,_,_,command=self.priority_fixture(root)
            command[command.index("--max-workers")+1]="1"
            self.run_main(command,lambda a,**k:subprocess.CompletedProcess(a,0,"123.pbs\n",""))
            lock=root/"datasets/generation_00003.lock"; lock.mkdir(parents=True)
            before=json.loads((root/"dispatch/state.json").read_text())
            old=next(iter(before["jobs"].values()))
            def retry(arguments,**kwargs):
                self.assertFalse(lock.exists())
                with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                    self.assertEqual(db.execute("SELECT status FROM jobs WHERE token=?",(old["token"],)).fetchone()[0],"terminal")
                return subprocess.CompletedProcess(arguments,0,"124.pbs\n","")
            self.run_main(command,retry,dict(state="terminal",exit_status=137))
            state=json.loads((root/"dispatch/state.json").read_text())
            self.assertEqual(len(state["jobs"]),1)
            self.assertEqual(next(iter(state["jobs"].values()))["job_id"],"124.pbs")
            self.assertTrue(Path(old["script"]).exists())
            with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                self.assertEqual(db.execute("SELECT COUNT(*) FROM jobs").fetchone()[0],2)
                self.assertEqual(db.execute("SELECT COUNT(*) FROM job_units").fetchone()[0],2)
                self.assertEqual(db.execute("SELECT attempts FROM units WHERE kind='generation' AND unit_index=3").fetchone()[0],2)
                self.assertEqual(db.execute("SELECT completed FROM cleanup").fetchone()[0],1)
                event=json.loads(db.execute("SELECT payload FROM events WHERE event='job_authoritatively_finished'").fetchone()[0])
                self.assertEqual(event["exit_status"],137)

    def test_restart_finishes_committed_terminal_cleanup_before_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,_,_,command=self.priority_fixture(root)
            command[command.index("--max-workers")+1]="1"
            self.run_main(command,lambda a,**k:subprocess.CompletedProcess(a,0,"123.pbs\n",""))
            lock=root/"datasets/generation_00003.lock"; lock.mkdir(parents=True)
            with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                token=db.execute("SELECT token FROM jobs").fetchone()[0]
                # Simulate death after terminal commit, before filesystem cleanup.
                db.execute("UPDATE jobs SET status='terminal' WHERE token=?",(token,))
                db.execute("UPDATE units SET active_job=NULL WHERE active_job=?",(token,))
                db.execute("INSERT INTO cleanup(token,kind,unit_index) VALUES (?,'generation',3)",(token,))
            def retry(a,**k):
                self.assertFalse(lock.exists())
                return subprocess.CompletedProcess(a,0,"124.pbs\n","")
            self.run_main(command,retry)
            with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                self.assertEqual(db.execute("SELECT completed FROM cleanup").fetchone()[0],1)

    def test_legacy_state_or_missing_ledger_never_silently_discards_history(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,_,_,command=self.priority_fixture(root)
            control=root/"dispatch"; control.mkdir()
            original=json.dumps(dict(schema=1,jobs={"historical":{"status":"live"}}))
            state=control/"state.json"; state.write_text(original)
            result=subprocess.run(command+["--dry-run"],text=True,capture_output=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn("not automatically migrated",result.stderr)
            self.assertEqual(state.read_text(),original)
            self.assertFalse((control/"ledger.sqlite3").exists())
            state.write_text(json.dumps(dict(schema=2,jobs={"committed":{"status":"submission_unknown"}})))
            result=subprocess.run(command+["--dry-run"],text=True,capture_output=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn("refusing to discard",result.stderr)

    def test_resume_rejects_changed_immutable_ledger_projection(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,_,_,command=self.priority_fixture(root)
            result=subprocess.run(command+["--dry-run"],text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr)
            path=root/"dispatch/ledger.sqlite3"
            with sqlite_connection(path) as db:
                self.assertEqual(db.execute("PRAGMA journal_mode").fetchone()[0],"delete")
                db.execute("UPDATE units SET payload='{}' WHERE kind='generation' AND unit_index=3")
            changed=subprocess.run(command+["--dry-run"],text=True,capture_output=True)
            self.assertNotEqual(changed.returncode,0)
            self.assertIn("Frozen scheduler ledger projection changed",changed.stderr)
            with sqlite_connection(path) as db:
                self.assertEqual(db.execute("SELECT COUNT(*) FROM jobs").fetchone()[0],0)

    def test_completion_status_counts_are_rebuilt_from_revalidated_receipts(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,_,_,command=self.priority_fixture(root)
            first=subprocess.run(command+["--dry-run"],text=True,capture_output=True)
            self.assertEqual(first.returncode,0,first.stderr)
            with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                db.execute("UPDATE stats SET completed=total")
            resumed=subprocess.run(command+["--dry-run"],text=True,capture_output=True)
            self.assertEqual(resumed.returncode,0,resumed.stderr)
            state=json.loads((root/"dispatch/state.json").read_text())
            self.assertEqual(state["complete_counts"]["generation"],0)
            self.assertIn('"would_submit": "generation"',resumed.stdout)

    def test_resume_rejects_lost_live_assignment_before_any_submission(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,_,_,command=self.priority_fixture(root)
            command[command.index("--max-workers")+1]="1"
            self.run_main(command,lambda a,**k:subprocess.CompletedProcess(a,0,"123.pbs\n",""))
            with sqlite_connection(root/"dispatch/ledger.sqlite3") as db:
                db.execute("UPDATE units SET active_job=NULL WHERE active_job IS NOT NULL")
            with self.assertRaisesRegex(RuntimeError,"assignment ledger is inconsistent"):
                self.run_main(command,lambda *a,**k:self.fail("Duplicate qsub"),dict(state="live"))

    def test_large_history_stays_on_disk_and_status_has_only_counts_and_live_jobs(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); _,manifest,path,command=self.priority_fixture(root)
            n=10000
            with (root/"tasks.csv").open("w") as output:
                output.write("index,task_id,stage,method,dataset_ids\n")
                for i in range(1,n+1): output.write(f"{i},t{i},challenge,exact,d1\n")
            (root/"datasets.csv").write_text("dataset_id,sim_sha256,sim_file\nd1,frozen,path\n")
            value=json.loads(manifest.read_text())
            for name in value["files"]: value["files"][name]=DISPATCH.sha256(root/name)
            manifest.write_text(json.dumps(value))
            value=json.loads(path.read_text()); value["dispatch_manifest_sha256"]=DISPATCH.sha256(manifest)
            path.write_text(json.dumps(value))
            result=subprocess.run(command+["--dry-run"],text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stderr)
            state=json.loads((root/"dispatch/state.json").read_text())
            fingerprint={key:state[key] for key in ("plan_sha256","bootstrap_sha256","dispatch_manifest_sha256","priority_manifest_sha256","controller_sha256")}
            ledger=DISPATCH.ControllerLedger(root/"dispatch/ledger.sqlite3",fingerprint)
            try:
                with ledger.transaction():
                    ledger.db.executemany("INSERT INTO jobs VALUES (?,?,?,NULL,?,'terminal',?,?,NULL,NULL)",
                        ((f"history{i}",f"old{i}","methods",time.time(),"archived-script","a"*64) for i in range(n)))
                    ledger.db.executemany("INSERT INTO events(epoch,event,payload) VALUES (?,'historical','{}')",
                                         ((time.time(),) for _ in range(n)))
                    ledger.db.execute("UPDATE units SET completed=1,attempts=1 WHERE kind='methods'")
                snapshot=ledger.snapshot(fingerprint)
                self.assertEqual(snapshot["jobs"],{})
                self.assertEqual(snapshot["complete_counts"]["methods"],n)
                self.assertLess(len(json.dumps(snapshot)),3000)
                self.assertNotIn("complete",snapshot)
                self.assertNotIn("attempts",snapshot)
                plan=ledger.db.execute("EXPLAIN QUERY PLAN SELECT * FROM jobs WHERE status IN ('submission_pending','submission_unknown','live')").fetchall()
                self.assertTrue(any("jobs_active" in x[3] for x in plan))
                with self.assertRaisesRegex(RuntimeError,"rollback"):
                    with ledger.transaction():
                        ledger.db.execute("UPDATE units SET completed=0 WHERE kind='methods' AND unit_index=1")
                        ledger.event("must_rollback",{})
                        raise RuntimeError("rollback")
                self.assertEqual(ledger.counts()["methods"]["completed"],n)
                self.assertEqual(ledger.db.execute("SELECT COUNT(*) FROM events WHERE event='must_rollback'").fetchone()[0],0)
            finally: ledger.close()
            # The 10k cached completions are only hints: no receipt exists.
            again=subprocess.run(command+["--dry-run"],text=True,capture_output=True)
            self.assertEqual(again.returncode,0,again.stderr)
            snapshot=json.loads((root/"dispatch/state.json").read_text())
            self.assertEqual(snapshot["complete_counts"]["methods"],0)
            self.assertLess(len(json.dumps(snapshot)),3000)


if __name__ == "__main__":
    unittest.main()
