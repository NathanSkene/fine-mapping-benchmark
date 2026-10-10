#!/usr/bin/env python3
"""Bounded PBS dispatch with durable assignments and verified terminal receipts.

This process is lightweight scheduler orchestration. R inference runs only in
PBS jobs. A failed observation never authorizes restarting an existing job.
"""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import csv
import fcntl
import hashlib
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


@contextmanager
def exclusive_controller(control: Path):
    """Shared-storage ownership: atomic mkdir first, local flock second.

    No existing claim, including an incomplete claim left by a killed process,
    is automatically cleared. The owner can remove only the unchanged directory
    and exact metadata it created. Workers never enter this context.
    """
    claim = control / "controller.owner"
    try:
        claim.mkdir()
    except FileExistsError as error:
        raise SystemExit("Another or stale controller ownership claim exists: "+str(claim)+
                         "; inspect its owner and scheduler state before manual recovery") from error
    identity = claim.stat()
    owner = dict(schema=1,hostname=socket.gethostname(),pid=os.getpid(),
                 nonce=uuid.uuid4().hex,started_ns=time.time_ns())
    raw = (json.dumps(owner,sort_keys=True)+"\n").encode()
    # A failure before metadata has been fully synced leaves an incomplete
    # claim deliberately. It cannot be safely auto-reconciled by another node.
    owner_path = claim / "owner.json"
    with owner_path.open("xb") as handle:
        handle.write(raw); handle.flush(); os.fsync(handle.fileno())
    for path in (claim,control):
        descriptor = os.open(path,os.O_RDONLY)
        try: os.fsync(descriptor)
        finally: os.close(descriptor)
    try:
        with (control/"controller.lock").open("a") as lock:
            try:
                fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise SystemExit("Another local controller owns this campaign") from error
            yield owner
    finally:
        # Never recursively delete a claim, and never remove foreign, missing,
        # edited or replaced metadata even after our own process encounters an
        # exception. Such a claim needs explicit operator reconciliation.
        try:
            current = claim.lstat()
            unchanged = (current.st_dev,current.st_ino)==(identity.st_dev,identity.st_ino) \
                and not claim.is_symlink() and owner_path.read_bytes()==raw \
                and set(x.name for x in claim.iterdir())=={"owner.json"}
        except OSError:
            unchanged = False
        if not unchanged:
            raise RuntimeError("Controller ownership claim changed; preserved for manual recovery: "+str(claim))
        owner_path.unlink(); claim.rmdir()
        descriptor = os.open(control,os.O_RDONLY)
        try: os.fsync(descriptor)
        finally: os.close(descriptor)


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


def validate_priority_manifest(path: Path, plan_sha: str, method_ids: set[int],
                               generation_ids: set[int], plan: Path,
                               dispatch_manifest_sha: str) -> tuple[list[int], list[int], str]:
    """Read a design-only scheduling prefix, committed to this exact full plan."""
    raw = path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict) or type(value.get("schema")) is not int \
            or value["schema"] != 1 or value.get("plan_sha256") != plan_sha:
        raise RuntimeError("Foreign or malformed priority manifest")
    if "plan_rds_sha256" in value and value["plan_rds_sha256"] != sha256(plan):
        raise RuntimeError("Priority manifest plan bytes changed")
    if "dispatch_manifest_sha256" in value \
            and value["dispatch_manifest_sha256"] != dispatch_manifest_sha:
        raise RuntimeError("Priority manifest dispatch commitment changed")
    orders = []
    for name, allowed in (("method_task_indices", method_ids),
                          ("generation_indices", generation_ids)):
        order = value.get(name)
        # JSON booleans compare equal to integers in Python; reject them (and
        # fractional/string coercions) rather than silently changing task IDs.
        if not isinstance(order, list) or any(type(x) is not int for x in order) \
                or len(set(order)) != len(order) or any(x not in allowed for x in order):
            raise RuntimeError(f"Invalid, duplicate or out-of-range priority {name}")
        orders.append(order)
    return orders[0], orders[1], hashlib.sha256(raw).hexdigest()


def scheduling_order(indices, prefix: list[int]) -> list[int]:
    """Retain the complete original remainder after the explicit prefix."""
    selected = set(prefix)
    return prefix + [index for index in indices if index not in selected]


def priority_batches(order: list[int], tasks: dict[int, dict], candidate,
                     batch_tasks: int, slots: int) -> list[list[int]]:
    """Batch adjacent ready prefix tasks without moving a later method ahead."""
    batches: list[list[int]] = []
    previous_key = None
    for index in order:
        if not candidate("methods", index):
            continue
        task = tasks[index]
        key = (task["method"], task["stage"])
        limit = min(batch_tasks, 2) if task["method"] == "carma" else batch_tasks
        if batches and key == previous_key and len(batches[-1]) < limit:
            batches[-1].append(index)
        elif len(batches) < slots:
            batches.append([index])
        else:
            break
        previous_key = key
    return batches


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


def start_ledger_path(script: str | Path) -> Path:
    return Path(script).with_suffix(".starts.jsonl")


def read_start_ledger(path: Path, token: str, kind: str, plan_sha: str,
                      controller_sha: str, indices: list[int] | None = None,
                      max_records: int | None = None) -> list[tuple[dict, str]]:
    """A complete prefix of this job's ordered start intents, never a hint."""
    import stat
    if not path.exists() or path.is_symlink() or not stat.S_ISREG(path.stat().st_mode):
        raise RuntimeError(f"Missing or invalid worker start ledger: {path}")
    limit = len(indices) if indices is not None else max_records
    if limit is None or limit < 0:
        raise RuntimeError("Missing worker start-ledger bound")
    fields = {"schema", "token", "kind", "index", "ordinal", "plan_sha256",
              "controller_sha256", "started_ns"}
    def object_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result: raise ValueError("Duplicate start-ledger field")
            result[key] = value
        return result
    records = []
    with path.open("rb") as handle:
        while True:
            line = handle.readline(2049)
            if not line: break
            if len(line) > 2048 or not line.endswith(b"\n") or len(records) >= limit:
                raise RuntimeError(f"Truncated, oversized or overlong worker start ledger: {path}")
            try:
                value = json.loads(line, object_pairs_hook=object_pairs)
            except (ValueError, UnicodeDecodeError) as error:
                raise RuntimeError(f"Malformed worker start ledger: {path}") from error
            ordinal = len(records) + 1
            if not isinstance(value, dict) or set(value) != fields \
                    or type(value["schema"]) is not int or value["schema"] != 1 \
                    or value["token"] != token or value["kind"] != kind \
                    or value["plan_sha256"] != plan_sha or value["controller_sha256"] != controller_sha \
                    or type(value["ordinal"]) is not int or value["ordinal"] != ordinal \
                    or type(value["index"]) is not int or value["index"] < 1 \
                    or (indices is not None and value["index"] != indices[ordinal-1]) \
                    or type(value["started_ns"]) is not int or not 0 < value["started_ns"] < 2**63:
                raise RuntimeError(f"Foreign or out-of-order worker start ledger: {path}")
            records.append((value, hashlib.sha256(line).hexdigest()))
    return records


def record_task_start(path: Path, token: str, kind: str, index: int, ordinal: int,
                      plan_sha: str, controller_sha: str) -> None:
    """Worker-only mode: fsync intent before shell is permitted to invoke R."""
    if not re.fullmatch(r"[a-f0-9]{32}", token) or kind not in {"generation", "methods"} \
            or type(index) is not int or index < 1 or type(ordinal) is not int or ordinal < 1 \
            or not re.fullmatch(r"[a-f0-9]{64}", plan_sha) \
            or sha256(Path(__file__)) != controller_sha:
        raise RuntimeError("Invalid worker task-start contract")
    prior = read_start_ledger(path, token, kind, plan_sha, controller_sha,
                              max_records=ordinal-1)
    if len(prior) != ordinal-1:
        raise RuntimeError("Worker start ledger is not the expected preceding prefix")
    value = dict(schema=1, token=token, kind=kind, index=index, ordinal=ordinal,
                 plan_sha256=plan_sha, controller_sha256=controller_sha,
                 started_ns=time.time_ns())
    line = (json.dumps(value, sort_keys=True, separators=(",", ":"))+"\n").encode()
    descriptor = os.open(path, os.O_WRONLY | os.O_APPEND | getattr(os, "O_NOFOLLOW", 0))
    try:
        if os.write(descriptor, line) != len(line):
            raise RuntimeError("Incomplete task-start append; refuse to invoke R")
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


class ControllerLedger:
    """Disk-backed plan, assignments and history; completion rows are hints.

    A single flock owner writes this database. A FULL-synchronous rollback
    journal preserves committed submission intent before qsub; JSON is only
    a bounded status projection. WAL is avoided on shared cluster filesystems.
    """
    SCHEMA = 3

    def __init__(self, path: Path, fingerprint: dict):
        self.path = path
        self.db = sqlite3.connect(path, isolation_level=None)
        try:
            self.db.row_factory = sqlite3.Row
            self.db.execute("PRAGMA foreign_keys=ON")
            self.db.execute("PRAGMA busy_timeout=30000")
            self.db.execute("PRAGMA synchronous=FULL")
            if self.db.execute("PRAGMA journal_mode=DELETE").fetchone()[0] != "delete":
                raise RuntimeError("Scheduler requires SQLite rollback journaling")
            value = self.db.execute("SELECT value FROM metadata WHERE key='fingerprint'").fetchone()
            schema = self.db.execute("SELECT value FROM metadata WHERE key='schema'").fetchone()
            if schema is None or json.loads(schema[0]) != self.SCHEMA or value is None \
                    or json.loads(value[0]) != fingerprint:
                raise RuntimeError("Plan, bootstrap, controller or priority changed; refusing incompatible scheduler resume")
            expected = self.db.execute("SELECT value FROM metadata WHERE key='projection_sha256'").fetchone()
            if expected is None or self.projection_digest(self.db) != expected[0]:
                raise RuntimeError("Frozen scheduler ledger projection changed")
        except Exception:
            self.db.close()
            raise

    @contextmanager
    def transaction(self):
        self.db.execute("BEGIN IMMEDIATE")
        try:
            yield
            self.db.execute("COMMIT")
        except BaseException:
            self.db.execute("ROLLBACK")
            raise

    @classmethod
    def projection_digest(cls, db) -> str:
        """Canonical immutable rows; live assignments and completion hints vary."""
        digest = hashlib.sha256()
        for label,query in (
                ("units","SELECT kind,unit_index,position,payload,method,stage,dataset_ids,priority_rank,original_index FROM units ORDER BY kind,unit_index"),
                ("datasets","SELECT dataset_id,payload,generation_index FROM datasets ORDER BY dataset_id"),
                ("dependencies","SELECT method_index,generation_index FROM dependencies ORDER BY method_index,generation_index"),
                ("method_groups","SELECT method,stage,first_position FROM method_groups ORDER BY method,stage")):
            digest.update(label.encode()+b"\n")
            for row in db.execute(query):
                digest.update(json.dumps(list(row),separators=(",",":"),ensure_ascii=True).encode()+b"\n")
        return digest.hexdigest()

    @classmethod
    def create(cls, path: Path, root: Path, fingerprint: dict,
               priority_methods: list[int], priority_generation: list[int]) -> None:
        # No jobs can be submitted until this complete projection is installed.
        # A crash may leave a staging file; it never becomes resume authority.
        temporary = path.with_name(path.name + f".initializing.{os.getpid()}")
        if temporary.exists():
            raise RuntimeError(f"Preserve and inspect interrupted ledger initialization: {temporary}")
        db = sqlite3.connect(temporary)
        db.execute("PRAGMA foreign_keys=ON")
        db.execute("PRAGMA synchronous=FULL")
        try:
            db.executescript("""
                CREATE TABLE metadata(key TEXT PRIMARY KEY, value TEXT NOT NULL);
                CREATE TABLE stats(kind TEXT PRIMARY KEY, total INTEGER NOT NULL DEFAULT 0,
                    completed INTEGER NOT NULL DEFAULT 0, attempts INTEGER NOT NULL DEFAULT 0,
                    prestart_failures INTEGER NOT NULL DEFAULT 0);
                INSERT INTO stats(kind) VALUES ('generation'),('methods');
                CREATE TABLE jobs(token TEXT PRIMARY KEY, name TEXT NOT NULL UNIQUE,
                    kind TEXT NOT NULL, job_id TEXT, submitted_epoch REAL NOT NULL,
                    status TEXT NOT NULL CHECK(status IN ('submission_pending','submission_unknown','live','terminal','missing')),
                    script TEXT NOT NULL, script_sha256 TEXT NOT NULL,
                    last_observation TEXT, detail TEXT);
                CREATE TABLE units(kind TEXT NOT NULL, unit_index INTEGER NOT NULL,
                    position INTEGER NOT NULL, payload TEXT NOT NULL, method TEXT, stage TEXT,
                    dataset_ids TEXT, priority_rank INTEGER, original_index INTEGER,
                    completed INTEGER NOT NULL DEFAULT 0 CHECK(completed IN (0,1)),
                    attempts INTEGER NOT NULL DEFAULT 0, active_job TEXT REFERENCES jobs(token),
                    prestart_failures INTEGER NOT NULL DEFAULT 0,
                    PRIMARY KEY(kind,unit_index));
                CREATE INDEX units_ready ON units(kind,completed,active_job,position);
                CREATE INDEX units_group ON units(kind,method,stage,completed,active_job,position);
                CREATE INDEX units_original ON units(method,stage,dataset_ids,position);
                CREATE INDEX units_priority ON units(kind,completed,active_job,priority_rank,position);
                CREATE INDEX units_attempt_limit ON units(completed,active_job,attempts);
                CREATE INDEX units_prestart_limit ON units(completed,active_job,prestart_failures);
                CREATE TABLE method_groups(method TEXT NOT NULL, stage TEXT NOT NULL,
                    first_position INTEGER NOT NULL, PRIMARY KEY(method,stage));
                CREATE TABLE job_units(token TEXT NOT NULL REFERENCES jobs(token),
                    kind TEXT NOT NULL, unit_index INTEGER NOT NULL,
                    PRIMARY KEY(token,kind,unit_index),
                    FOREIGN KEY(kind,unit_index) REFERENCES units(kind,unit_index));
                CREATE TABLE datasets(dataset_id TEXT PRIMARY KEY, payload TEXT NOT NULL,
                    generation_index INTEGER);
                CREATE TABLE dependencies(method_index INTEGER NOT NULL, generation_index INTEGER NOT NULL,
                    PRIMARY KEY(method_index,generation_index));
                CREATE INDEX dependencies_method ON dependencies(method_index);
                CREATE TABLE events(sequence INTEGER PRIMARY KEY AUTOINCREMENT, epoch REAL NOT NULL,
                    event TEXT NOT NULL, token TEXT, kind TEXT, unit_index INTEGER, payload TEXT NOT NULL);
                CREATE TABLE cleanup(token TEXT NOT NULL REFERENCES jobs(token), kind TEXT NOT NULL,
                    unit_index INTEGER NOT NULL, completed INTEGER NOT NULL DEFAULT 0,
                    PRIMARY KEY(token,kind,unit_index),
                    FOREIGN KEY(kind,unit_index) REFERENCES units(kind,unit_index));
                CREATE INDEX cleanup_pending ON cleanup(completed);
                CREATE INDEX jobs_active ON jobs(status);
                CREATE TABLE task_starts(token TEXT NOT NULL REFERENCES jobs(token),
                    kind TEXT NOT NULL, unit_index INTEGER NOT NULL, ordinal INTEGER NOT NULL,
                    record_sha256 TEXT NOT NULL, started_ns INTEGER NOT NULL,
                    PRIMARY KEY(token,kind,unit_index), UNIQUE(token,ordinal),
                    FOREIGN KEY(kind,unit_index) REFERENCES units(kind,unit_index));
                CREATE INDEX task_starts_unit ON task_starts(kind,unit_index);
                CREATE TRIGGER unit_insert AFTER INSERT ON units BEGIN
                    UPDATE stats SET total=total+1 WHERE kind=NEW.kind;
                END;
                CREATE TRIGGER unit_completion AFTER UPDATE OF completed ON units BEGIN
                    UPDATE stats SET completed=completed+NEW.completed-OLD.completed WHERE kind=NEW.kind;
                END;
                CREATE TRIGGER unit_attempt AFTER UPDATE OF attempts ON units BEGIN
                    UPDATE stats SET attempts=attempts+NEW.attempts-OLD.attempts WHERE kind=NEW.kind;
                END;
                CREATE TRIGGER unit_prestart AFTER UPDATE OF prestart_failures ON units BEGIN
                    UPDATE stats SET prestart_failures=prestart_failures+NEW.prestart_failures-OLD.prestart_failures WHERE kind=NEW.kind;
                END;
            """)
            db.execute("INSERT INTO metadata VALUES ('schema',?)", (json.dumps(cls.SCHEMA),))
            db.execute("INSERT INTO metadata VALUES ('fingerprint',?)", (json.dumps(fingerprint, sort_keys=True),))
            for filename, kind, key, priority in (
                    ("generation.csv", "generation", "generation_id", priority_generation),
                    ("tasks.csv", "methods", "index", priority_methods)):
                ranks = {index: rank for rank, index in enumerate(priority)}
                with (root / filename).open(newline="") as handle:
                    for position, row in enumerate(csv.DictReader(handle)):
                        index = int(row[key])
                        db.execute("INSERT INTO units(kind,unit_index,position,payload,method,stage,dataset_ids,priority_rank) "
                                   "VALUES (?,?,?,?,?,?,?,?)", (kind,index,position,json.dumps(row),
                                   row.get("method"),row.get("stage"),row.get("dataset_ids"),ranks.get(index)))
            with (root / "datasets.csv").open(newline="") as handle:
                for row in csv.DictReader(handle):
                    match = re.search(r"generation_([0-9]+)",row["sim_file"]) if row["sim_sha256"] == "generated" else None
                    if row["sim_sha256"] == "generated" and not match:
                        raise RuntimeError("Generated dataset has no generation identity")
                    db.execute("INSERT INTO datasets VALUES (?,?,?)", (row["dataset_id"],json.dumps(row),
                               int(match.group(1)) if match else None))
            # SQL joins retain pairing/dependency constraints without keeping a
            # million Python dictionaries or scanning old job history per poll.
            db.execute("UPDATE units AS u SET original_index=(SELECT o.unit_index FROM units o "
                       "WHERE o.kind='methods' AND o.stage='original_refit' AND o.method=u.method "
                       "AND o.dataset_ids=u.dataset_ids ORDER BY o.position DESC LIMIT 1) "
                       "WHERE u.kind='methods' AND u.stage='corrected_refit'")
            db.execute("INSERT INTO method_groups SELECT method,stage,MIN(position) FROM units "
                       "WHERE kind='methods' GROUP BY method,stage")
            for index,payload in db.execute("SELECT unit_index,payload FROM units WHERE kind='methods'"):
                for identity in json.loads(payload)["dataset_ids"].split("|"):
                    row = db.execute("SELECT generation_index FROM datasets WHERE dataset_id=?",(identity,)).fetchone()
                    if row is None:
                        raise RuntimeError(f"Method task references missing dataset: {identity}")
                    if row[0] is not None:
                        if db.execute("SELECT 1 FROM units WHERE kind='generation' AND unit_index=?",(row[0],)).fetchone() is None:
                            raise RuntimeError("Dataset references missing generation unit")
                        db.execute("INSERT OR IGNORE INTO dependencies VALUES (?,?)",(index,row[0]))
            db.execute("INSERT INTO metadata VALUES ('projection_sha256',?)",(cls.projection_digest(db),))
            db.commit()
            db.close()
            os.replace(temporary,path)
        except BaseException:
            db.close()
            raise

    def close(self):
        self.db.close()

    def event(self, name: str, payload: dict, token=None, kind=None, index=None):
        self.db.execute("INSERT INTO events(epoch,event,token,kind,unit_index,payload) VALUES (?,?,?,?,?,?)",
                        (time.time(),name,token,kind,index,json.dumps(payload,sort_keys=True)))

    def unit(self, kind: str, index: int) -> dict:
        row = self.db.execute("SELECT payload FROM units WHERE kind=? AND unit_index=?",(kind,index)).fetchone()
        if row is None:
            raise RuntimeError("Unknown planned unit")
        return json.loads(row[0])

    def counts(self) -> dict:
        return {row["kind"]:dict(row) for row in self.db.execute("SELECT * FROM stats")}

    def rebuild_counts(self):
        # Counters are a fast status projection, not completion authority.
        self.db.execute("UPDATE stats SET total=(SELECT COUNT(*) FROM units u WHERE u.kind=stats.kind), "
                        "completed=(SELECT COALESCE(SUM(completed),0) FROM units u WHERE u.kind=stats.kind), "
                        "attempts=(SELECT COALESCE(SUM(attempts),0) FROM units u WHERE u.kind=stats.kind), "
                        "prestart_failures=(SELECT COALESCE(SUM(prestart_failures),0) FROM units u WHERE u.kind=stats.kind)")

    def validate_assignments(self):
        active = "('submission_pending','submission_unknown','live')"
        mismatch = self.db.execute(
            "SELECT j.token FROM jobs j JOIN job_units a ON a.token=j.token "
            "JOIN units u ON u.kind=a.kind AND u.unit_index=a.unit_index "
            "WHERE j.status IN "+active+" AND u.active_job IS NOT j.token LIMIT 1").fetchone()
        orphan = self.db.execute(
            "SELECT u.active_job FROM units u JOIN jobs j ON j.token=u.active_job "
            "WHERE u.active_job IS NOT NULL AND (j.status NOT IN "+active+
            " OR NOT EXISTS (SELECT 1 FROM job_units a WHERE a.token=u.active_job "
            "AND a.kind=u.kind AND a.unit_index=u.unit_index)) LIMIT 1").fetchone()
        empty = self.db.execute(
            "SELECT j.token FROM jobs j WHERE j.status IN "+active+
            " AND NOT EXISTS (SELECT 1 FROM job_units a WHERE a.token=j.token) LIMIT 1").fetchone()
        if mismatch or orphan or empty:
            raise RuntimeError("Scheduler assignment ledger is inconsistent; refusing possible duplicate submission")

    def active_jobs(self) -> list[dict]:
        jobs = []
        for row in self.db.execute("SELECT * FROM jobs WHERE status IN ('submission_pending','submission_unknown','live') ORDER BY submitted_epoch,token"):
            job = dict(row)
            job["indices"] = [x[0] for x in self.db.execute("SELECT unit_index FROM job_units WHERE token=? ORDER BY rowid",(job["token"],))]
            job["last_observation"] = json.loads(job["last_observation"]) if job["last_observation"] else None
            jobs.append(job)
        return jobs

    def ingest_starts(self, job: dict, fingerprint: dict) -> int:
        records = read_start_ledger(start_ledger_path(job["script"]), job["token"],
            job["kind"], fingerprint["plan_sha256"], fingerprint["controller_sha256"], job["indices"])
        known = list(self.db.execute("SELECT ordinal,record_sha256 FROM task_starts WHERE token=? ORDER BY ordinal", (job["token"],)))
        if len(known) > len(records) or any(row["ordinal"] != i+1 or row["record_sha256"] != records[i][1]
                                            for i,row in enumerate(known)):
            raise RuntimeError("Previously counted worker starts were truncated or rewritten")
        with self.transaction():
            for value, digest in records[len(known):]:
                self.db.execute("INSERT INTO task_starts VALUES (?,?,?,?,?,?)",
                    (job["token"],job["kind"],value["index"],value["ordinal"],digest,value["started_ns"]))
                self.db.execute("UPDATE units SET attempts=attempts+1 WHERE kind=? AND unit_index=?", (job["kind"],value["index"]))
                self.event("task_start_counted", dict(ordinal=value["ordinal"],record_sha256=digest),job["token"],job["kind"],value["index"])
        return len(records)

    def snapshot(self, fingerprint: dict) -> dict:
        counts = self.counts()
        live = self.active_jobs()
        for job in live:
            if job["detail"]: job["detail"] = job["detail"][-2000:]
            if job["last_observation"] and job["last_observation"].get("detail"):
                job["last_observation"]["detail"] = job["last_observation"]["detail"][-2000:]
        return dict(fingerprint, schema=self.SCHEMA, ledger=self.path.name,
                    ledger_projection_sha256=self.db.execute("SELECT value FROM metadata WHERE key='projection_sha256'").fetchone()[0],
                    complete_counts={kind:x["completed"] for kind,x in counts.items()},
                    planned_counts={kind:x["total"] for kind,x in counts.items()},
                    task_start_intents={kind:x["attempts"] for kind,x in counts.items()},
                    prestart_infrastructure_failures={kind:x["prestart_failures"] for kind,x in counts.items()},
                    retry_blockers=json.loads((self.db.execute("SELECT value FROM metadata WHERE key='retry_blockers'").fetchone() or ["{}"]) [0]),
                    jobs={job["token"]:job for job in live},
                    event_sequence=self.db.execute("SELECT COALESCE(MAX(sequence),0) FROM events").fetchone()[0],
                    updated_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()))

    def ready(self, kind: str, max_attempts: int, ready_methods: set[str] | None = None,
              priority: bool | None = None, group: tuple[str,str] | None = None,
              limit: int = 1, max_prestart_failures: int = 3) -> list[int]:
        where = ["u.kind=?","u.completed=0","u.active_job IS NULL","u.attempts<?","u.prestart_failures<?",
                 "NOT EXISTS (SELECT 1 FROM cleanup c WHERE c.kind=u.kind AND c.unit_index=u.unit_index AND c.completed=0)"]
        values: list = [kind,max_attempts,max_prestart_failures]
        if kind == "methods":
            if not ready_methods:
                return []
            where += ["u.method IN ("+",".join("?" for _ in ready_methods)+")",
                      "(u.original_index IS NULL OR EXISTS (SELECT 1 FROM units o WHERE o.kind='methods' AND o.unit_index=u.original_index AND o.completed=1))",
                      "NOT EXISTS (SELECT 1 FROM dependencies d JOIN units g ON g.kind='generation' AND g.unit_index=d.generation_index WHERE d.method_index=u.unit_index AND g.completed=0)"]
            values.extend(sorted(ready_methods))
        if priority is not None:
            where.append("u.priority_rank IS "+("NOT NULL" if priority else "NULL"))
        if group is not None:
            where.extend(["u.method=?","u.stage=?"]); values.extend(group)
        order = "u.priority_rank,u.position" if priority else \
                "(u.priority_rank IS NULL),u.priority_rank,u.position" if kind == "generation" else "u.position"
        values.append(limit)
        return [row[0] for row in self.db.execute("SELECT u.unit_index FROM units u WHERE "+
                   " AND ".join(where)+" ORDER BY "+order+" LIMIT ?",values)]

    def report_retry_blockers(self, path: Path, max_attempts: int, max_prestart_failures: int) -> dict:
        # Most polls have no exhausted units. Avoid walking every dependency
        # in a full campaign to establish that fact; these two lookups use the
        # dedicated retry-limit indexes. A blocked prerequisite must itself be
        # an incomplete, unassigned exhausted unit.
        exhausted = self.db.execute(
            "SELECT 1 FROM units WHERE completed=0 AND active_job IS NULL AND attempts>=? "
            "UNION ALL SELECT 1 FROM units WHERE completed=0 AND active_job IS NULL "
            "AND prestart_failures>=? LIMIT 1", (max_attempts,max_prestart_failures)).fetchone()
        failure = "({alias}.attempts>=? OR {alias}.prestart_failures>=?)"
        own = failure.format(alias="u")
        original = "(u.kind='methods' AND EXISTS (SELECT 1 FROM units o WHERE o.kind='methods' AND o.unit_index=u.original_index AND o.completed=0 AND o.active_job IS NULL AND "+failure.format(alias="o")+"))"
        generation = "(u.kind='methods' AND EXISTS (SELECT 1 FROM dependencies d JOIN units g ON g.kind='generation' AND g.unit_index=d.generation_index WHERE d.method_index=u.unit_index AND g.completed=0 AND g.active_job IS NULL AND "+failure.format(alias="g")+"))"
        sql = "SELECT u.kind,u.unit_index,u.method,u.stage,u.attempts,u.prestart_failures, CASE WHEN u.attempts>=? THEN 'task_start_limit' WHEN u.prestart_failures>=? THEN 'prestart_infrastructure_limit' WHEN "+generation+" THEN 'generation_prerequisite_exhausted' WHEN "+original+" THEN 'original_prerequisite_exhausted' END AS reason FROM units u WHERE u.completed=0 AND u.active_job IS NULL AND ("+own+" OR "+generation+" OR "+original+") ORDER BY u.kind,u.unit_index"
        bounds = [max_attempts,max_prestart_failures]*6
        temporary = path.with_name(path.name+f".tmp.{os.getpid()}")
        counts = {}
        with temporary.open("w",newline="") as handle:
            writer = csv.DictWriter(handle,fieldnames=["kind","unit_index","method","stage","attempts","prestart_failures","reason"])
            writer.writeheader()
            for row in self.db.execute(sql,bounds) if exhausted else ():
                writer.writerow(dict(row)); counts[row["reason"]] = counts.get(row["reason"],0)+1
            handle.flush(); os.fsync(handle.fileno())
        os.replace(temporary,path)
        result = dict(counts=counts,total=sum(counts.values()),report=path.name,
                      report_sha256=sha256(path),max_task_start_intents=max_attempts,
                      max_prestart_infrastructure_failures=max_prestart_failures)
        with self.transaction():
            self.db.execute("INSERT INTO metadata VALUES ('retry_blockers',?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",(json.dumps(result,sort_keys=True),))
            self.event("retry_blockers_reported",result)
        return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--bootstrap", type=Path, required=True,
                        help="Frozen shell environment script that loads R/Python/library paths")
    parser.add_argument("--priority-manifest", type=Path,
                        help="Immutable design-only scheduling prefix; every remaining task is retained")
    parser.add_argument("--queue", default=None,
                        help="Explicit PBS queue; omit to use the site's permitted default")
    parser.add_argument("--select", default="1:ncpus=1:mem=32gb")
    parser.add_argument("--walltime", default="72:00:00")
    parser.add_argument("--max-workers", type=int, default=64)
    parser.add_argument("--max-generation", type=int, default=16)
    parser.add_argument("--batch-tasks", type=int, default=4)
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--max-prestart-failures", type=int, default=3,
                        help="Separate cap on confirmed terminal worker failures before its first task start")
    parser.add_argument("--start-python", type=Path, default=Path(sys.executable),
                        help="Controller-compatible Python executable available on compute nodes for durable start markers")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if min(args.max_workers,args.max_generation,args.batch_tasks,args.max_attempts,args.max_prestart_failures,args.poll_seconds) < 1:
        raise SystemExit("Worker, generation, batch, poll and attempt limits must be positive")
    root,repo = args.plan.resolve().parent,args.repo.resolve()
    bootstrap = args.bootstrap.resolve(strict=True)
    start_python = args.start_python.resolve(strict=True)
    plan_sha = (root/"plan.sha256").read_text().strip()
    if not re.fullmatch(r"[a-f0-9]{64}",plan_sha):
        raise SystemExit("Missing valid immutable plan commitment")
    dispatch_sha = validate_dispatch_manifest(root,args.plan.resolve(),plan_sha)
    priority_methods,priority_generation,priority_sha = [],[],None
    priority_path = args.priority_manifest.resolve(strict=True) if args.priority_manifest else None
    if priority_path is not None:
        # Only retain IDs requested by the prefix while streaming the complete
        # frozen CSV. This keeps validation memory independent of total tasks.
        requested = json.loads(priority_path.read_bytes())
        allowed = []
        for filename,key,field in (("tasks.csv","index","method_task_indices"),
                                   ("generation.csv","generation_id","generation_indices")):
            order = requested.get(field,[]) if isinstance(requested,dict) else []
            wanted = set(order) if isinstance(order,list) and all(type(x) is int for x in order) else set()
            with (root/filename).open(newline="") as handle:
                allowed.append({int(row[key]) for row in csv.DictReader(handle) if int(row[key]) in wanted})
        priority_methods,priority_generation,priority_sha = validate_priority_manifest(
            priority_path,plan_sha,allowed[0],allowed[1],args.plan.resolve(),dispatch_sha)
    fingerprint = dict(plan_sha256=plan_sha,bootstrap_sha256=sha256(bootstrap),
                       dispatch_manifest_sha256=dispatch_sha,priority_manifest_sha256=priority_sha,
                       controller_sha256=sha256(Path(__file__)),start_python=str(start_python),
                       start_python_sha256=sha256(start_python))
    control = root/"dispatch"; control.mkdir(exist_ok=True)
    with exclusive_controller(control):
        state_path,ledger_path = control/"state.json",control/"ledger.sqlite3"
        if state_path.exists() and json.loads(state_path.read_text()).get("schema") != ControllerLedger.SCHEMA:
            raise SystemExit("Legacy JSON scheduler state is not automatically migrated; preserve/archive it explicitly before starting a new ledger")
        if not ledger_path.exists():
            if state_path.exists():
                raise RuntimeError("Scheduler ledger is missing; refusing to discard committed history")
            ControllerLedger.create(ledger_path,root,fingerprint,priority_methods,priority_generation)
            if validate_dispatch_manifest(root,args.plan.resolve(),plan_sha) != dispatch_sha:
                raise RuntimeError("Dispatch projection changed during ledger initialization")
        ledger = ControllerLedger(ledger_path,fingerprint)
        try:
            _dispatch(args,root,repo,bootstrap,control,state_path,ledger,fingerprint,priority_path)
        finally:
            ledger.close()


def _dispatch(args,root,repo,bootstrap,control,state_path,ledger,fingerprint,priority_path):
    plan_sha = fingerprint["plan_sha256"]
    ledger.validate_assignments()
    def persist():
        atomic_json(state_path,ledger.snapshot(fingerprint))

    def directory(kind,index,row=None):
        row = row or ledger.unit(kind,index)
        return root/("datasets" if kind=="generation" else "tasks")/(
            f"generation_{index:05d}" if kind=="generation" else row["task_id"])

    def receipt(kind,index,row=None):
        row = row or ledger.unit(kind,index)
        return terminal_receipt(directory(kind,index,row),
                                f"generation_{index}" if kind=="generation" else row["task_id"],
                                plan_sha,"generated.rds" if kind=="generation" else "completed.rds")

    def mark_completed(kind,index,complete,phase):
        ledger.db.execute("UPDATE units SET completed=? WHERE kind=? AND unit_index=?",(int(complete),kind,index))
        ledger.event("completion_hint_changed",dict(complete=complete,phase=phase),kind=kind,index=index)

    # A ledger completion is never receipt authority. Verify every planned unit
    # once on startup, streaming rows and writing changed hints in small batches.
    changes,total_verified = [],0
    for row in ledger.db.execute("SELECT kind,unit_index,payload,completed FROM units ORDER BY kind,unit_index"):
        value = receipt(row["kind"],row["unit_index"],json.loads(row["payload"]))
        total_verified += 1
        if value != bool(row["completed"]):
            changes.append((row["kind"],row["unit_index"],value))
        if len(changes)>=1024:
            with ledger.transaction():
                for kind,index,value in changes: mark_completed(kind,index,value,"startup")
            changes.clear()
    with ledger.transaction():
        for kind,index,value in changes: mark_completed(kind,index,value,"startup")
        ledger.rebuild_counts()
        ledger.event("startup_receipts_verified",dict(n_units=total_verified))
    persist()

    def completed(kind,index):
        hint = ledger.db.execute("SELECT completed FROM units WHERE kind=? AND unit_index=?",(kind,index)).fetchone()[0]
        if hint: return True
        if not receipt(kind,index): return False
        with ledger.transaction(): mark_completed(kind,index,True,"poll")
        return True

    def guard_priority():
        if priority_path is not None and sha256(priority_path)!=fingerprint["priority_manifest_sha256"]:
            raise RuntimeError("Priority manifest changed during dispatch; refusing additional submissions")
        if sha256(bootstrap)!=fingerprint["bootstrap_sha256"] or sha256(Path(__file__))!=fingerprint["controller_sha256"]:
            raise RuntimeError("Bootstrap or controller changed during dispatch; refusing additional submissions")
        if sha256(Path(fingerprint["start_python"])) != fingerprint["start_python_sha256"]:
            raise RuntimeError("Task-start Python executable changed during dispatch")

    def cleanup_terminal_locks():
        import shutil
        for item in list(ledger.db.execute("SELECT c.* FROM cleanup c JOIN jobs j ON j.token=c.token "
                                           "WHERE c.completed=0 AND j.status IN ('terminal','missing')")):
            kind,index,token = item["kind"],item["unit_index"],item["token"]
            active = ledger.db.execute("SELECT active_job FROM units WHERE kind=? AND unit_index=?",(kind,index)).fetchone()[0]
            if active is not None:
                raise RuntimeError("Pending terminal cleanup overlaps a current assignment")
            is_complete = completed(kind,index)
            lock_dir = Path(str(directory(kind,index))+".lock")
            removed = not is_complete and lock_dir.exists()
            if removed: shutil.rmtree(lock_dir)
            with ledger.transaction():
                ledger.db.execute("UPDATE cleanup SET completed=1 WHERE token=? AND kind=? AND unit_index=?",(token,kind,index))
                ledger.event("authoritative_terminal_lock_cleanup",dict(removed=removed,receipt_complete=is_complete),token,kind,index)

    def submit(kind,indices):
        guard_priority()
        token = uuid.uuid4().hex
        name = "fmc"+token[:10]
        script = control/f"job-{token}.sh"
        commands = ["#!/bin/bash",f"#PBS -N {name}",
                    *([f"#PBS -q {args.queue}"] if args.queue else []),
                    f"#PBS -l select={args.select}",f"#PBS -l walltime={args.walltime}",
                    f"#PBS -o {control}/",f"#PBS -e {control}/","set -euo pipefail",
                    f"source {shlex.quote(str(bootstrap))}",f"cd {shlex.quote(str(repo))}"]
        worker = "generate.R" if kind=="generation" else "run_task.R"
        for ordinal,index in enumerate(indices,1):
            commands.append(" ".join(shlex.quote(x) for x in [fingerprint["start_python"],str(Path(__file__).resolve()),
                "--record-start",str(start_ledger_path(script)),token,kind,str(index),str(ordinal),plan_sha,fingerprint["controller_sha256"]]))
            commands.append(" ".join(shlex.quote(x) for x in [args.rscript,"--vanilla",
                str(repo/"scripts/campaign"/worker),"--plan",str(args.plan.resolve()),"--index",str(index)]))
        with script.open("w") as handle:
            handle.write("\n".join(commands)+"\n"); handle.flush(); os.fsync(handle.fileno())
        if args.dry_run:
            print(json.dumps(dict(would_submit=kind,indices=indices,script=str(script))),flush=True)
            return
        with start_ledger_path(script).open("xb") as handle:
            handle.flush(); os.fsync(handle.fileno())
        descriptor = os.open(control, os.O_RDONLY)
        try: os.fsync(descriptor)
        finally: os.close(descriptor)
        # Reserve all indices and commit intent atomically BEFORE qsub. A crash
        # or ambiguous result cannot make them eligible for another submission.
        with ledger.transaction():
            ledger.db.execute("INSERT INTO jobs(token,name,kind,submitted_epoch,status,script,script_sha256) VALUES (?,?,?,?,?,?,?)",
                              (token,name,kind,time.time(),"submission_pending",str(script),sha256(script)))
            for index in indices:
                cursor = ledger.db.execute("UPDATE units SET active_job=? WHERE kind=? AND unit_index=? AND completed=0 AND active_job IS NULL AND attempts<? AND prestart_failures<?",
                                           (token,kind,index,args.max_attempts,args.max_prestart_failures))
                if cursor.rowcount!=1: raise RuntimeError("Attempted duplicate or ineligible assignment")
                ledger.db.execute("INSERT INTO job_units VALUES (?,?,?)",(token,kind,index))
            ledger.event("submission_intent",dict(name=name,indices=indices),token)
        persist()
        try:
            result = subprocess.run(["qsub",str(script)],text=True,capture_output=True,timeout=60)
        except (OSError,subprocess.TimeoutExpired) as error:
            with ledger.transaction():
                ledger.db.execute("UPDATE jobs SET status='submission_unknown',detail=? WHERE token=?",(str(error),token))
                ledger.event("submission_unknown",dict(error=str(error)),token)
            persist()
            raise RuntimeError("qsub outcome unknown; inspect qselect -N "+name) from error
        if result.returncode or not re.fullmatch(r"[0-9]+(?:\[[^]]*\])?\.[A-Za-z0-9_.-]+",result.stdout.strip()):
            with ledger.transaction():
                ledger.db.execute("UPDATE jobs SET status='submission_unknown',detail=? WHERE token=?",(result.stdout+result.stderr,token))
                ledger.event("submission_unknown",dict(returncode=result.returncode,stdout=result.stdout,stderr=result.stderr),token)
            persist()
            raise RuntimeError("qsub did not return a confirmed handle; inspect qselect -N "+name)
        with ledger.transaction():
            ledger.db.execute("UPDATE jobs SET job_id=?,status='live' WHERE token=?",(result.stdout.strip(),token))
            ledger.event("submission_confirmed",dict(job_id=result.stdout.strip(),stdout=result.stdout,stderr=result.stderr),token)
        persist()
        print(json.dumps(dict(submitted=result.stdout.strip(),kind=kind,indices=indices)),flush=True)

    def verified_ready(kind,limit,ready_methods=None,priority=None,group=None):
        selected = []
        while len(selected)<limit:
            ids = ledger.ready(kind,args.max_attempts,ready_methods,priority,group,limit+len(selected),args.max_prestart_failures)
            changed = False
            for index in ids:
                if index in selected: continue
                if completed(kind,index): changed=True
                else: selected.append(index)
                if len(selected)==limit: break
            if not changed: break
        return selected

    cleanup_terminal_locks()
    def ingest_starts(job):
        try:
            return ledger.ingest_starts(job,fingerprint)
        except RuntimeError as error:
            with ledger.transaction():
                ledger.event("worker_start_ledger_error",dict(error=str(error),assignments_retained=True),job["token"])
            persist()
            raise RuntimeError("Unverifiable worker start ledger; assignments retained: "+str(error)) from error

    while True:
        guard_priority()
        live_counts = dict(generation=0,methods=0)
        for job in ledger.active_jobs():
            starts = ingest_starts(job)
            if not job.get("job_id"):
                raise RuntimeError("Unresolved submission: inspect qselect -N "+job["name"])
            raw = observe_job(job["job_id"])
            observation = raw
            if raw["state"]=="missing" and time.time()-job["submitted_epoch"]<120:
                observation = dict(state="observation_error",detail="Unknown new ID inside propagation grace period")
            with ledger.transaction():
                ledger.db.execute("UPDATE jobs SET last_observation=? WHERE token=?",(json.dumps(observation),job["token"]))
                ledger.event("scheduler_observation",dict(raw=raw,effective=observation),job["token"])
            if observation["state"] in {"terminal","missing"}:
                # A worker can append while qstat is in flight. The terminal
                # observation now establishes that its final prefix is stable;
                # account it again before releasing any assignment.
                starts = ingest_starts(job)
                outcomes = [(index,completed(job["kind"],index)) for index in job["indices"]]
                with ledger.transaction():
                    ledger.db.execute("UPDATE jobs SET status=? WHERE token=?",(observation["state"],job["token"]))
                    for index,value in outcomes:
                        cursor = ledger.db.execute("UPDATE units SET active_job=NULL WHERE kind=? AND unit_index=? AND active_job=?",
                                                   (job["kind"],index,job["token"]))
                        if cursor.rowcount!=1: raise RuntimeError("Terminal job lost its committed assignment")
                        if not value:
                            ledger.db.execute("INSERT INTO cleanup(token,kind,unit_index) VALUES (?,?,?)",(job["token"],job["kind"],index))
                    if starts==0 and not outcomes[0][1]:
                        ledger.db.execute("UPDATE units SET prestart_failures=prestart_failures+1 WHERE kind=? AND unit_index=?",(job["kind"],job["indices"][0]))
                        ledger.event("prestart_infrastructure_failure",dict(first_index=job["indices"][0],observation=observation),job["token"])
                    ledger.event("job_authoritatively_finished",dict(observation,started_task_count=starts,
                        start_ledger_sha256=sha256(start_ledger_path(job["script"]))),job["token"])
                cleanup_terminal_locks()
                print(json.dumps(dict(finished_handle=job["job_id"],observation=observation)),flush=True)
            else:
                live_counts[job["kind"]]+=1
        persist()
        methods = {x[0] for x in ledger.db.execute("SELECT DISTINCT method FROM method_groups")}
        ready_methods = methods & {"exact","forced_one"}
        for method in methods-ready_methods:
            if terminal_receipt(root/"preflight"/method,method,plan_sha): ready_methods.add(method)
        generation_slots = max(0,min(args.max_generation-live_counts["generation"],args.max_workers-sum(live_counts.values())))
        for index in verified_ready("generation",generation_slots):
            submit("generation",[index]); live_counts["generation"]+=1
        method_slots = max(0,args.max_workers-sum(live_counts.values()))
        submissions = 0
        if method_slots:
            ids = verified_ready("methods",method_slots*args.batch_tasks,ready_methods,True)
            task_rows = {index:ledger.unit("methods",index) for index in ids}
            for batch in priority_batches(ids,task_rows,lambda kind,index:True,args.batch_tasks,method_slots):
                submit("methods",batch); submissions+=1
            groups = ledger.db.execute("SELECT method,stage FROM method_groups ORDER BY first_position")
            for group in groups:
                if submissions>=method_slots: break
                size = min(args.batch_tasks,2) if group[0]=="carma" else args.batch_tasks
                batch = verified_ready("methods",size,ready_methods,False,(group[0],group[1]))
                if batch: submit("methods",batch); submissions+=1
        persist()
        counts = ledger.counts()
        blockers = ledger.report_retry_blockers(control/"retry-blockers.csv",args.max_attempts,args.max_prestart_failures)
        persist()
        if blockers["total"]:
            print(json.dumps(dict(incomplete_retry_blockers=blockers)),flush=True)
        print(json.dumps(dict(completed_generation=counts["generation"]["completed"],
              planned_generation=counts["generation"]["total"],completed_methods=counts["methods"]["completed"],
              planned_methods=counts["methods"]["total"],active_before_submission=live_counts)),flush=True)
        if all(x["completed"]==x["total"] for x in counts.values()):
            print("Campaign terminal receipts complete; run the R collector and requirement audit.",flush=True)
            return
        if blockers["total"] and not ledger.active_jobs() and not args.dry_run:
            raise RuntimeError("Campaign incomplete: retry-exhausted or dependency-blocked work remains; inspect "+str(control/"retry-blockers.csv"))
        if args.once or args.dry_run: return
        time.sleep(args.poll_seconds)

if __name__ == "__main__":
    if len(sys.argv)>1 and sys.argv[1]=="--record-start":
        if len(sys.argv)!=9: raise SystemExit("Invalid task-start arguments")
        record_task_start(Path(sys.argv[2]),sys.argv[3],sys.argv[4],int(sys.argv[5]),int(sys.argv[6]),sys.argv[7],sys.argv[8])
    else:
        main()
