#!/usr/bin/env python3
"""Bounded PBS dispatch with durable assignments and verified terminal receipts.

This process is lightweight scheduler orchestration. R inference runs only in
PBS jobs. A failed observation never authorizes restarting an existing job.
"""
from __future__ import annotations
import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import time
import uuid


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    with temporary.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def terminal_receipt(directory: Path, identity: str, plan_sha: str,
                     checkpoint_name: str = "completed.rds") -> bool:
    receipt = directory / "terminal.json"
    if not receipt.exists():
        return False
    value = json.loads(receipt.read_text())
    if value.get("schema") != 1 or value.get("terminal") is not True \
            or value.get("id") != identity or value.get("plan_sha256") != plan_sha:
        raise RuntimeError(f"Foreign or malformed terminal receipt: {receipt}")
    checkpoint = directory / checkpoint_name
    if not checkpoint.exists() or sha256(checkpoint) != value.get("checkpoint_sha256"):
        raise RuntimeError(f"Terminal checkpoint changed: {checkpoint}")
    return True


def validate_dispatch_manifest(root: Path, plan: Path, plan_sha: str) -> str:
    manifest = root / "dispatch-manifest.json"
    value = json.loads(manifest.read_text())
    if value.get("schema") != 1 or value.get("plan_sha256") != plan_sha \
            or value.get("plan_rds_sha256") != sha256(plan):
        raise RuntimeError("Frozen dispatch plan commitment changed")
    files = value.get("files", {})
    if set(files) != {"tasks.csv", "datasets.csv", "generation.csv"}:
        raise RuntimeError("Malformed frozen dispatch projections")
    for name, expected in files.items():
        if sha256(root / name) != expected:
            raise RuntimeError(f"Frozen dispatch projection changed: {name}")
    return sha256(manifest)


def observe_job(job_id: str, runner=subprocess.run) -> dict:
    try:
        result = runner(["qstat", "-xf", "-F", "json", job_id],
                        text=True, capture_output=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as error:
        return {"state": "observation_error", "detail": str(error)}
    if result.returncode:
        text = result.stdout + result.stderr
        if re.search(r"Unknown Job Id|Unknown Job ID|Job has finished", text, re.I):
            return {"state": "missing", "detail": text[-2000:]}
        return {"state": "observation_error", "detail": text[-2000:]}
    try:
        data = json.loads(result.stdout)
        jobs = data.get("Jobs", {})
        if job_id not in jobs:
            # A syntactically valid empty response is still weaker than an
            # explicit server-side Unknown Job Id; retain the assignment.
            return {"state": "observation_error", "detail": "Requested ID absent from qstat JSON"}
        value = jobs[job_id]
        return {"state": "terminal" if value.get("job_state") in {"F", "C", "X"} else "live",
                "pbs_state": value.get("job_state"), "exit_status": value.get("Exit_status"),
                "resources_used": value.get("resources_used", {})}
    except (ValueError, TypeError) as error:
        return {"state": "observation_error", "detail": str(error)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--bootstrap", type=Path, required=True,
                        help="Frozen shell environment script that loads R/Python/library paths")
    parser.add_argument("--queue", default="v1_small72a")
    parser.add_argument("--select", default="1:ncpus=1:mem=32gb")
    parser.add_argument("--walltime", default="72:00:00")
    parser.add_argument("--max-workers", type=int, default=64)
    parser.add_argument("--max-generation", type=int, default=16)
    parser.add_argument("--batch-tasks", type=int, default=4)
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if min(args.max_workers, args.max_generation, args.batch_tasks, args.max_attempts) < 1:
        raise SystemExit("Worker, generation, batch and attempt limits must be positive")
    root = args.plan.resolve().parent
    repo = args.repo.resolve()
    bootstrap = args.bootstrap.resolve(strict=True)
    plan_sha = (root / "plan.sha256").read_text().strip()
    if not re.fullmatch(r"[a-f0-9]{64}", plan_sha):
        raise SystemExit("Missing valid immutable plan commitment")
    dispatch_manifest_sha = validate_dispatch_manifest(root, args.plan.resolve(), plan_sha)
    with (root / "tasks.csv").open(newline="") as handle:
        tasks = list(csv.DictReader(handle))
    with (root / "generation.csv").open(newline="") as handle:
        generation = list(csv.DictReader(handle))
    with (root / "datasets.csv").open(newline="") as handle:
        datasets = {row["dataset_id"]: row for row in csv.DictReader(handle)}
    control = root / "dispatch"
    control.mkdir(exist_ok=True)
    lock = (control / "controller.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit("Another controller owns this campaign")
    state_path = control / "state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {
        "schema": 1, "plan_sha256": plan_sha, "bootstrap_sha256": sha256(bootstrap),
        "dispatch_manifest_sha256": dispatch_manifest_sha,
        "controller_sha256": sha256(Path(__file__)),
        "jobs": {}, "complete": {"generation": [], "methods": []}, "attempts": {}}
    if state["plan_sha256"] != plan_sha or state["bootstrap_sha256"] != sha256(bootstrap) \
            or state.get("dispatch_manifest_sha256") != dispatch_manifest_sha \
            or state.get("controller_sha256") != sha256(Path(__file__)):
        raise SystemExit("Plan or bootstrap changed; refusing incompatible scheduler resume")
    # Persisted IDs are scheduling hints. Revalidate every receipt once on
    # startup, including IDs marked complete in a previous controller run.
    done = {kind: set() for kind in ("generation", "methods")}
    by_index = {"generation": {int(g["generation_id"]): g for g in generation},
                "methods": {int(t["index"]): t for t in tasks}}
    original_units = {(t["method"], t["dataset_ids"]): int(t["index"])
                      for t in tasks if t["stage"] == "original_refit"}

    def completed(kind: str, index: int) -> bool:
        if index in done[kind]:
            return True
        row = by_index[kind][index]
        if kind == "generation":
            directory = root / "datasets" / f"generation_{index:05d}"
            value = terminal_receipt(directory, f"generation_{index}", plan_sha, "generated.rds")
        else:
            directory = root / "tasks" / row["task_id"]
            value = terminal_receipt(directory, row["task_id"], plan_sha)
        if value:
            done[kind].add(index)
        return value

    # Resume untracked terminal work once. Subsequent polls inspect only live
    # assignments and candidate units, not every file in a 500k-task campaign.
    for kind, rows in by_index.items():
        for index in rows:
            completed(kind, index)

    def persist() -> None:
        state["complete"] = {kind: sorted(values) for kind, values in done.items()}
        state["updated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        atomic_json(state_path, state)

    def submit(kind: str, indices: list[int]) -> None:
        token = uuid.uuid4().hex[:10]
        name = "fmc" + token
        script = control / f"job-{token}.sh"
        commands = ["#!/bin/bash", f"#PBS -N {name}", f"#PBS -q {args.queue}",
                    f"#PBS -l select={args.select}", f"#PBS -l walltime={args.walltime}",
                    f"#PBS -o {control}/", f"#PBS -e {control}/", "set -euo pipefail",
                    f"source {shlex.quote(str(bootstrap))}", f"cd {shlex.quote(str(repo))}"]
        worker = "generate.R" if kind == "generation" else "run_task.R"
        for index in indices:
            commands.append(" ".join(shlex.quote(x) for x in [args.rscript, "--vanilla",
                str(repo / "scripts/campaign" / worker), "--plan", str(args.plan.resolve()),
                "--index", str(index)]))
        script.write_text("\n".join(commands) + "\n")
        if args.dry_run:
            print(json.dumps({"would_submit": kind, "indices": indices, "script": str(script)}), flush=True)
            return
        # Persist pending assignment before dispatch. An ambiguous qsub outcome
        # remains assigned and needs explicit recovery; never double-submit it.
        job = {"name": name, "kind": kind, "indices": indices, "job_id": None,
               "submitted_epoch": time.time(), "status": "submission_pending", "script": str(script)}
        state["jobs"][token] = job
        persist()

        try:
            result = subprocess.run(["qsub", str(script)], text=True, capture_output=True, timeout=60)
        except (OSError, subprocess.TimeoutExpired) as error:
            job["status"] = "submission_unknown"; job["detail"] = str(error); persist()
            raise RuntimeError("qsub outcome unknown; inspect qselect -N " + name) from error
        if result.returncode or not re.fullmatch(r"[0-9]+(?:\[[^]]*\])?\.[A-Za-z0-9_.-]+", result.stdout.strip()):
            job["status"] = "submission_unknown"; job["detail"] = result.stdout + result.stderr; persist()
            raise RuntimeError("qsub did not return a confirmed handle; inspect qselect -N " + name)
        job["job_id"] = result.stdout.strip(); job["status"] = "live"
        for index in indices:
            key = f"{kind}:{index}"
            state["attempts"][key] = state["attempts"].get(key, 0) + 1
        persist()
        print(json.dumps({"submitted": job["job_id"], "kind": kind, "indices": indices}), flush=True)

    while True:
        assigned = {"generation": set(), "methods": set()}
        live_counts = {"generation": 0, "methods": 0}
        for token, job in state["jobs"].items():
            if job["status"] in {"terminal", "missing"}:
                continue
            if not job.get("job_id"):
                raise RuntimeError("Unresolved submission: inspect qselect -N " + job["name"])
            observation = observe_job(job["job_id"])
            job["last_observation"] = observation
            if observation["state"] == "missing" and time.time() - job["submitted_epoch"] < 120:
                observation = {"state": "observation_error", "detail": "Unknown new ID inside propagation grace period"}
            if observation["state"] in {"terminal", "missing"}:
                job["status"] = observation["state"]
                for index in job["indices"]:
                    if not completed(job["kind"], index):
                        # Scheduler has authoritatively stopped the job. Retain
                        # its failure before removing an abandoned worker lock.
                        row = by_index[job["kind"]][index]
                        directory = root / ("datasets" if job["kind"] == "generation" else "tasks") / (
                            f"generation_{index:05d}" if job["kind"] == "generation" else row["task_id"])
                        lock_dir = Path(str(directory) + ".lock")
                        if lock_dir.exists():
                            import shutil
                            shutil.rmtree(lock_dir)
                print(json.dumps({"finished_handle": job["job_id"], "observation": observation}), flush=True)
            else:
                assigned[job["kind"]].update(job["indices"])
                live_counts[job["kind"]] += 1
        persist()

        ready_methods = {t["method"] for t in tasks if t["method"] in {"exact", "forced_one"}}
        for method in {t["method"] for t in tasks} - ready_methods:
            if terminal_receipt(root / "preflight" / method, method, plan_sha):
                ready_methods.add(method)

        def candidate(kind: str, index: int) -> bool:
            if index in assigned[kind] or completed(kind, index):
                return False
            if state["attempts"].get(f"{kind}:{index}", 0) >= args.max_attempts:
                return False
            if kind == "methods":
                task = by_index[kind][index]
                method = task["method"]
                if method not in ready_methods:
                    return False
                if task["stage"] == "corrected_refit":
                    original_index = original_units.get((method, task["dataset_ids"]))
                    if original_index is not None and original_index not in done["methods"]:
                        # Identical authored code/inputs can reuse the first
                        # realized fit; serialize this pair to avoid wasted
                        # concurrent duplicate computation.
                        return False
                for identity in task["dataset_ids"].split("|"):
                    dataset = datasets[identity]
                    if dataset["sim_sha256"] == "generated":
                        match = re.search(r"generation_([0-9]+)", dataset["sim_file"])
                        if not match or int(match.group(1)) not in done["generation"]:
                            return False
            return True

        candidates = {kind: [] for kind in by_index}
        slots = {"generation": max(0, min(args.max_generation - live_counts["generation"],
                        args.max_workers - sum(live_counts.values()))),
                 "methods": 0}
        for index in by_index["generation"]:
            if candidate("generation", index):
                candidates["generation"].append(index)
                if len(candidates["generation"]) >= slots["generation"]:
                    break
        if slots["generation"] == 0:
            candidates["generation"] = []
        for index in candidates["generation"]:
            submit("generation", [index]); live_counts["generation"] += 1
        slots["methods"] = max(0, args.max_workers - sum(live_counts.values()))
        # Keep expensive methods in separate batches so CARMA cannot hold a
        # otherwise fast job beyond walltime. A batch is bounded, never an
        # enormous scheduler array of all planned task units.
        buckets: dict[tuple[str, str], list[int]] = {}
        for index, task in by_index["methods"].items():
            if slots["methods"] == 0:
                break
            key = (task["method"], task["stage"])
            bucket = buckets.setdefault(key, [])
            per_batch = min(args.batch_tasks, 2) if task["method"] == "carma" else args.batch_tasks
            if len(bucket) >= per_batch:
                continue
            if not candidate("methods", index):
                continue
            bucket.append(index)
        submissions = 0
        # Dispatch one bounded batch per method/stage before returning to the
        # first method, so optional/native implementations make progress early.
        for bucket in buckets.values():
            if bucket and submissions < slots["methods"]:
                submit("methods", bucket); submissions += 1
        persist()
        print(json.dumps({"completed_generation": len(done["generation"]),
                          "planned_generation": len(generation),
                          "completed_methods": len(done["methods"]), "planned_methods": len(tasks),
                          "active_before_submission": live_counts}), flush=True)
        if len(done["generation"]) == len(generation) and len(done["methods"]) == len(tasks):
            print("Campaign terminal receipts complete; run the R collector and requirement audit.", flush=True)
            return
        if args.once or args.dry_run:
            return
        time.sleep(args.poll_seconds)


if __name__ == "__main__":
    main()
