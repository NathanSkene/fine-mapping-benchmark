#!/usr/bin/env python3
"""Bounded native PBS campaign controller with verified per-worker receipts.

Preparation, native chains and regional panels have separate resource limits.
The shared legacy controller supplies conservative scheduler observations;
an unavailable scheduler never permits restarting an assigned job.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import importlib.util
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time
import uuid


def load_helpers(repo: Path):
    path = repo / "scripts/campaign/dispatch.py"
    spec = importlib.util.spec_from_file_location("fmb_dispatch", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def validate_projection_manifest(root: Path, plan: Path, plan_sha: str, helpers) -> str:
    path = root / "dispatch-manifest.json"
    value = json.loads(path.read_text())
    if value.get("schema") != 1 or value.get("plan_sha256") != plan_sha or \
            value.get("plan_rds_sha256") != helpers.sha256(plan):
        raise RuntimeError("Frozen native dispatch plan changed")
    files = value.get("files", {})
    if set(files) != {"draws.csv", "inputs.csv", "chains.csv", "regional.csv"}:
        raise RuntimeError("Malformed frozen native projections")
    for name, expected in files.items():
        if helpers.sha256(root / name) != expected:
            raise RuntimeError(f"Frozen native projection changed: {name}")
    return helpers.sha256(path)


def make_units(root: Path) -> dict[str, dict[str, dict]]:
    inputs = {row["input_id"]: row for row in read_rows(root / "inputs.csv")}
    units = {"prepare": {}, "chains": {}, "regional": {}}
    for i, row in enumerate(read_rows(root / "draws.csv"), 1):
        units["prepare"][row["draw_id"]] = {
            "index": i, "draw_id": row["draw_id"], "receipts": [row["draw_id"]]}
    for i, row in enumerate(read_rows(root / "chains.csv"), 1):
        units["chains"][row["chain_id"]] = {
            "index": i, "draw_id": row["draw_id"], "receipts": [row["chain_id"]]}
    for row in read_rows(root / "regional.csv"):
        identity = row["input_id"]
        unit = units["regional"].setdefault(identity, {
            "input_id": identity, "draw_id": inputs[identity]["draw_id"],
            "receipts": [], "methods": []})
        unit["receipts"].append(row["regional_id"])
        unit["methods"].append(row["method"])
    return units


def unit_directories(root: Path, kind: str, unit: dict) -> list[Path]:
    folder = "prepared" if kind == "prepare" else kind
    return [root / folder / identity for identity in unit["receipts"]]


def unit_complete(root: Path, kind: str, unit: dict, plan_sha: str, helpers) -> bool:
    return all(helpers.terminal_receipt(path, identity, plan_sha)
               for path, identity in zip(unit_directories(root, kind, unit), unit["receipts"]))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--repo", type=Path, required=True)
    p.add_argument("--bootstrap", type=Path, required=True)
    p.add_argument("--preflight-plan", type=Path, required=True,
                   help="Regional campaign plan with successful actual-method preflight receipts")
    p.add_argument("--queue", default="v1_small72a")
    p.add_argument("--max-workers", type=int, default=12)
    p.add_argument("--max-prepare", type=int, default=2)
    p.add_argument("--max-chains", type=int, default=6)
    p.add_argument("--max-regional", type=int, default=4)
    p.add_argument("--prepare-select", default="1:ncpus=1:mem=32gb")
    p.add_argument("--chains-select", default="1:ncpus=1:mem=8gb")
    p.add_argument("--regional-select", default="1:ncpus=1:mem=32gb")
    p.add_argument("--walltime", default="72:00:00")
    p.add_argument("--poll-seconds", type=int, default=60)
    p.add_argument("--max-attempts", type=int, default=3)
    p.add_argument("--rscript", default="Rscript")
    p.add_argument("--once", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if min(args.max_workers, args.max_prepare, args.max_chains, args.max_regional,
           args.max_attempts, args.poll_seconds) < 1:
        raise SystemExit("Resource and retry limits must be positive")
    root = args.plan.resolve(strict=True).parent
    repo = args.repo.resolve(strict=True)
    bootstrap = args.bootstrap.resolve(strict=True)
    helpers = load_helpers(repo)
    plan_sha = (root / "plan.sha256").read_text().strip()
    preflight_root = args.preflight_plan.resolve(strict=True).parent
    preflight_sha = (preflight_root / "plan.sha256").read_text().strip()
    if any(not re.fullmatch(r"[a-f0-9]{64}", x) for x in (plan_sha, preflight_sha)):
        raise SystemExit("Invalid immutable plan commitment")
    preflight_manifest_sha = helpers.validate_dispatch_manifest(
        preflight_root, args.preflight_plan.resolve(), preflight_sha)
    projection_sha = validate_projection_manifest(root, args.plan.resolve(), plan_sha, helpers)
    units = make_units(root)
    control = root / "dispatch"
    control.mkdir(exist_ok=True)
    lock = (control / "controller.lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise SystemExit("Another controller owns this native campaign")
    try:
        run_controller(args, root, repo, bootstrap, helpers, plan_sha,
                       preflight_root, preflight_sha, preflight_manifest_sha,
                       projection_sha, units, control)
    finally:
        lock.close()


def run_controller(args, root, repo, bootstrap, helpers, plan_sha,
                   preflight_root, preflight_sha, preflight_manifest_sha,
                   projection_sha, units, control) -> None:
    path = control / "state.json"
    fingerprint = {
        "plan_sha256": plan_sha, "preflight_sha256": preflight_sha,
        "preflight_manifest_sha256": preflight_manifest_sha,
        "projection_manifest_sha256": projection_sha,
        "bootstrap_sha256": helpers.sha256(bootstrap),
        "controller_sha256": helpers.sha256(Path(__file__)),
        "observation_helpers_sha256": helpers.sha256(repo / "scripts/campaign/dispatch.py")}
    state = json.loads(path.read_text()) if path.exists() else dict(
        schema=1, fingerprint=fingerprint, jobs={}, complete={k: [] for k in units}, attempts={})
    if state.get("fingerprint") != fingerprint:
        raise SystemExit("Plan, bootstrap or controller changed; refusing incompatible resume")
    # Persisted completion is a hint. Revalidate every receipt/checkpoint once
    # on startup; then retain verified identities in memory during this run.
    done = {kind: set() for kind in units}
    caps = {"prepare": args.max_prepare, "chains": args.max_chains,
            "regional": args.max_regional}
    resources = {"prepare": args.prepare_select, "chains": args.chains_select,
                 "regional": args.regional_select}

    def completed(kind: str, identity: str) -> bool:
        if identity in done[kind]:
            return True
        if unit_complete(root, kind, units[kind][identity], plan_sha, helpers):
            done[kind].add(identity)
            return True
        return False

    for kind in units:
        for identity in units[kind]:
            completed(kind, identity)

    def persist() -> None:
        state["complete"] = {kind: sorted(values) for kind, values in done.items()}
        state["updated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        helpers.atomic_json(path, state)

    def submit(kind: str, identity: str) -> None:
        token = uuid.uuid4().hex[:10]
        name = "fmn" + token
        script = control / f"job-{token}.sh"
        unit = units[kind][identity]
        worker = {"prepare": "prepare.R", "chains": "run_chain.R",
                  "regional": "run_regional.R"}[kind]
        command = [args.rscript, "--vanilla", str(repo / "scripts/native/campaign" / worker),
                   "--plan", str(args.plan.resolve())]
        command += ["--input-id", identity] if kind == "regional" else ["--index", str(unit["index"])]
        script.write_text("\n".join([
            "#!/bin/bash", f"#PBS -N {name}", f"#PBS -q {args.queue}",
            f"#PBS -l select={resources[kind]}", f"#PBS -l walltime={args.walltime}",
            f"#PBS -o {control}/", f"#PBS -e {control}/", "set -euo pipefail",
            f"source {shlex.quote(str(bootstrap))}", f"cd {shlex.quote(str(repo))}",
            " ".join(shlex.quote(x) for x in command)]) + "\n")
        if args.dry_run:
            print(json.dumps({"would_submit": kind, "id": identity, "script": str(script)}), flush=True)
            return
        job = {"name": name, "kind": kind, "id": identity, "job_id": None,
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
            raise RuntimeError("Unconfirmed qsub; inspect qselect -N " + name)
        job["job_id"] = result.stdout.strip(); job["status"] = "live"
        key = f"{kind}:{identity}"
        state["attempts"][key] = state["attempts"].get(key, 0) + 1
        persist()
        print(json.dumps({"submitted": job["job_id"], "kind": kind, "id": identity}), flush=True)

    while True:
        assigned = {kind: set() for kind in units}
        live = {kind: 0 for kind in units}
        for job in state["jobs"].values():
            if job["status"] in {"terminal", "missing"}:
                continue
            if not job.get("job_id"):
                raise RuntimeError("Unresolved submission: inspect qselect -N " + job["name"])
            observation = helpers.observe_job(job["job_id"])
            if observation["state"] == "missing" and time.time() - job["submitted_epoch"] < 120:
                observation = {"state": "observation_error", "detail": "Unknown new ID inside propagation grace period"}
            job["last_observation"] = observation
            if observation["state"] in {"terminal", "missing"}:
                job["status"] = observation["state"]
                kind, identity = job["kind"], job["id"]
                if not completed(kind, identity):
                    for directory in unit_directories(root, kind, units[kind][identity]):
                        worker_lock = Path(str(directory) + ".lock")
                        if worker_lock.exists():
                            shutil.rmtree(worker_lock)
                print(json.dumps({"finished_handle": job["job_id"], "observation": observation}), flush=True)
            else:
                assigned[job["kind"]].add(job["id"])
                live[job["kind"]] += 1
        persist()
        ready_methods = {method for unit in units["regional"].values() for method in unit["methods"]
                         if helpers.terminal_receipt(preflight_root / "preflight" / method, method, preflight_sha)}
        submissions = 0
        for kind in units:
            for identity, unit in units[kind].items():
                if live[kind] >= caps[kind] or sum(live.values()) >= args.max_workers:
                    break
                if identity in assigned[kind] or completed(kind, identity):
                    continue
                if state["attempts"].get(f"{kind}:{identity}", 0) >= args.max_attempts:
                    continue
                if kind != "prepare" and unit["draw_id"] not in done["prepare"]:
                    continue
                if kind == "regional" and not set(unit["methods"]).issubset(ready_methods):
                    continue
                submit(kind, identity); live[kind] += 1; submissions += 1
        persist()
        print(json.dumps({"complete": {k: len(done[k]) for k in units},
                          "planned": {k: len(units[k]) for k in units}, "active": live}), flush=True)
        if all(len(done[k]) == len(units[k]) for k in units):
            print("Native terminal receipts complete; run collector and diagnostic audit.", flush=True)
            return
        if args.once or args.dry_run:
            return
        if not sum(live.values()) and not submissions:
            raise RuntimeError("Native campaign stalled: inspect exhausted retries and regional preflight gates")
        time.sleep(args.poll_seconds)


if __name__ == "__main__":
    main()
