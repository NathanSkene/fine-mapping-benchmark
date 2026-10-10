import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import socket
import sqlite3
import subprocess
import sys
import time

EXPECTED_HELPER = "cc2da80c47856fba3aeff5b0817ec63ad13990c168d34e052a35f8caec957510"
root = Path(sys.argv[2]).resolve()
helper_path = Path(sys.argv[3]).resolve()
helper_hash = hashlib.sha256(helper_path.read_bytes()).hexdigest()
if helper_hash != EXPECTED_HELPER:
    raise RuntimeError("Unexpected frozen production ownership helper")
spec = importlib.util.spec_from_file_location("ownership_proof_helper", helper_path)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def write_json(path, value):
    temporary = path.with_name(path.name + ".pending")
    with temporary.open("x") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.link(temporary, path)
    temporary.unlink()


def wait_file(path, seconds):
    deadline = time.monotonic() + seconds
    while not path.exists():
        if time.monotonic() >= deadline:
            raise TimeoutError(f"No proof file arrived: {path}")
        time.sleep(2)
    return json.loads(path.read_text())


def database():
    db = sqlite3.connect(root / "ledger.sqlite3", isolation_level=None)
    assert db.execute("PRAGMA journal_mode=DELETE").fetchone()[0] == "delete"
    db.execute("PRAGMA synchronous=FULL")
    db.execute("PRAGMA foreign_keys=ON")
    assert db.execute("PRAGMA synchronous").fetchone()[0] == 2
    assert db.execute("PRAGMA foreign_keys").fetchone()[0] == 1
    return db


if sys.argv[1] == "child":
    try:
        with helper.exclusive_controller(root):
            raise RuntimeError("Compute node acquired the login controller claim")
    except SystemExit as error:
        if "controller ownership claim exists" not in str(error):
            raise
    write_json(root / "child-blocked.json", dict(host=socket.gethostname(),
        pbs_job_id=os.environ.get("PBS_JOBID"), production_claim_refused=True,
        helper_sha256=helper_hash, script_sha256=script_hash))
    wait_file(root / "parent-released.json", 300)
    with helper.exclusive_controller(root) as owner:
        assert owner["hostname"] == socket.gethostname() and owner["pid"] == os.getpid()
        assert json.loads((root / "controller.owner/owner.json").read_text()) == owner
        db = database()
        assert db.execute("SELECT value FROM units WHERE id=1").fetchone()[0] == "login-committed"
        db.execute("BEGIN IMMEDIATE")
        db.execute("UPDATE units SET value='compute-committed' WHERE id=1")
        db.execute("COMMIT")
        assert db.execute("PRAGMA quick_check").fetchone()[0] == "ok"
        db.close()
    assert not (root / "controller.owner").exists()
    assert helper.sha256(helper_path) == helper_hash
    write_json(root / "child-complete.json", dict(host=socket.gethostname(),
        production_claim_acquired_and_released=True, synced_owner_metadata=True,
        parent_commit_read=True, child_commit=True, sqlite_version=sqlite3.sqlite_version,
        helper_sha256=helper_hash, script_sha256=script_hash))
elif sys.argv[1] == "parent":
    root.mkdir(parents=True, exist_ok=False)
    parent_host = socket.gethostname()
    with helper.exclusive_controller(root) as owner:
        assert owner["hostname"] == parent_host and owner["pid"] == os.getpid()
        db = database()
        db.execute("CREATE TABLE units(id INTEGER PRIMARY KEY,value TEXT NOT NULL)")
        db.execute("BEGIN IMMEDIATE")
        db.execute("INSERT INTO units VALUES(1,'login-committed')")
        db.execute("COMMIT")
        db.close()
        job_script = root / "probe.pbs"
        job_script.write_text("#!/usr/bin/env bash\n#PBS -N fmownproof\n"
            "#PBS -l select=1:ncpus=1:mem=2gb\n#PBS -l walltime=00:10:00\n"
            "set -euo pipefail\nexec " + shlex.join([sys.executable,
                str(Path(__file__).resolve()), "child", str(root), str(helper_path)]) + "\n")
        write_json(root / "submission-intent.json", dict(parent_host=parent_host,
            helper_sha256=helper_hash, script_sha256=script_hash,
            job_script_sha256=helper.sha256(job_script), role="Actual production ownership helper; no inference"))
        submitted = subprocess.run(["/opt/pbs/bin/qsub", "-o", str(root / "pbs.stdout"),
            "-e", str(root / "pbs.stderr"), str(job_script)], capture_output=True, text=True)
        write_json(root / "submission-response.json", dict(returncode=submitted.returncode,
            stdout=submitted.stdout, stderr=submitted.stderr))
        if submitted.returncode or not re.fullmatch(r"[0-9]+\.[A-Za-z0-9_.-]+", submitted.stdout.strip()):
            raise RuntimeError("Inspect submission outcome; no automatic retry")
        job_id = submitted.stdout.strip()
        blocked = wait_file(root / "child-blocked.json", 900)
        assert blocked["script_sha256"] == script_hash and blocked["helper_sha256"] == helper_hash
        assert blocked["host"] != parent_host
    assert not (root / "controller.owner").exists()
    write_json(root / "parent-released.json", dict(released=True))
    completed = wait_file(root / "child-complete.json", 300)
    assert completed["script_sha256"] == script_hash and completed["helper_sha256"] == helper_hash
    with helper.exclusive_controller(root):
        db = database()
        assert db.execute("SELECT value FROM units WHERE id=1").fetchone()[0] == "compute-committed"
        assert db.execute("PRAGMA quick_check").fetchone()[0] == "ok"
        db.close()
    assert not (root / "controller.owner").exists() and helper.sha256(helper_path) == helper_hash
    proof = dict(status="PASS", pbs_job_id=job_id, parent_host=parent_host,
        child_host=completed["host"], sqlite_version=completed["sqlite_version"],
        production_claim_refused_cross_node=True, production_claim_acquired_after_release=True,
        production_claim_reacquired_on_login=True, graceful_claim_release=True,
        synced_owner_metadata=True, commits_visible_in_both_directions=True,
        journal_mode="delete", synchronous=2, foreign_keys=1,
        helper_sha256=helper_hash, script_sha256=script_hash,
        role="Actual frozen production ownership helper on login and PBS compute hosts",
        limitations="No storage-outage or killed-owner simulation; no simultaneous database writers are permitted")
    write_json(root / "proof.json", proof)
    print(json.dumps(proof, indent=2), flush=True)
else:
    raise ValueError("Expected parent or child mode")
