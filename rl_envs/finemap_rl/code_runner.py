"""Opt-in Docker execution: agent code never shares a Python process with truth.

Docker is a containment layer, not a verified boundary against kernel exploits.
Use dedicated workers and independently reviewed infrastructure for production.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import secrets
import subprocess
import tempfile
import threading
import time

from .tasks import canonical, public_task
from .environments import MethodDesignEnv
from .scoring import brier_loss

STARTER_CODE = '''import json, sys
import numpy as np
tasks = json.load(sys.stdin)
outputs = []
for x in tasks:
    # A deliberately limited forced-one approximation: improve the null model,
    # effect prior, multiple-cause treatment and LD dependence.
    w = np.exp(np.minimum(np.asarray(x["z"])**2 / 2, 500))
    w = w / w.sum()
    outputs.append({"task_id": x["task_id"], "variant_id": x["variant_id"], "pip": w.tolist(), "null_probability": 0.0})
json.dump(outputs, sys.stdout, allow_nan=False)
'''

ABF_CODE = '''import json, sys
import numpy as np
tasks = json.load(sys.stdin)
outputs = []
for x in tasks:
    p = len(x["z"])
    prior = .7 / p
    v = x["n"] * .05**2
    log_bf = -.5 * np.log1p(v) + .5 * np.asarray(x["z"])**2 * v / (1 + v)
    log_odds = np.log(prior / (1 - prior)) + log_bf
    pip = 1 / (1 + np.exp(-np.clip(log_odds, -700, 700)))
    outputs.append({"task_id": x["task_id"], "variant_id": x["variant_id"],
                    "pip": pip.tolist(), "null_probability": float(np.prod(1-pip))})
json.dump(outputs, sys.stdout, allow_nan=False)
'''


class DockerPredictorRunner:
    def __init__(self, image="finemap-rl-predictor:0.1", timeout=20, output_limit=1_048_576):
        if not 1 <= timeout <= 300 or not 1024 <= output_limit <= 4_194_304:
            raise ValueError("Invalid sandbox limits")
        self.timeout, self.output_limit = timeout, output_limit
        result = subprocess.run(["docker", "image", "inspect", image, "--format", "{{.Id}}"],
                                capture_output=True, text=True, timeout=15, check=True)
        self.image_id = result.stdout.strip()
        if not self.image_id.startswith("sha256:") or len(self.image_id) != 71:
            raise RuntimeError("Cannot resolve immutable predictor image")

    def run(self, source, inputs):
        if not isinstance(source, str) or not 1 <= len(source.encode()) <= 65_536:
            raise ValueError("Predictor source exceeds contract")
        payload = canonical([public_task(x) for x in inputs])
        name = "fmb-rl-" + secrets.token_hex(8)
        with tempfile.TemporaryDirectory(prefix="fmb-agent-") as directory:
            root = Path(directory); root.chmod(0o755)
            program = root / "predictor.py"; program.write_text(source); program.chmod(0o444)
            command = ["docker", "run", "--rm", "--name", name, "--log-driver", "none", "--network", "none", "--read-only",
                       "--user", "65534:65534", "--memory", "256m", "--memory-swap", "256m",
                       "--cpus", "1", "--pids-limit", "64", "--cap-drop", "ALL",
                       "--security-opt", "no-new-privileges", "--tmpfs", "/tmp:rw,noexec,nosuid,size=16m",
                       "--mount", f"type=bind,src={directory},dst=/agent,readonly",
                       "-i", self.image_id, "python", "-I", "/agent/predictor.py"]
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            buffers = {"stdout": bytearray(), "stderr": bytearray()}; overflow = threading.Event()

            def drain(stream, key):
                while True:
                    chunk = stream.read(4096)
                    if not chunk: break
                    if len(buffers[key]) + len(chunk) > self.output_limit:
                        # Keep draining without buffering until kill completes.
                        # Stopping reads can backpressure Docker's attach stream
                        # and interfere with timely daemon-side cleanup.
                        overflow.set()
                    elif not overflow.is_set():
                        buffers[key].extend(chunk)

            readers = [threading.Thread(target=drain, args=(stream, key), daemon=True)
                       for stream, key in ((process.stdout, "stdout"), (process.stderr, "stderr"))]
            for reader in readers: reader.start()
            try:
                # Input volume is bounded by panel_size and task schema; writing
                # occurs in a thread so a non-reading malicious script can time out.
                def feed():
                    try: process.stdin.write(payload); process.stdin.close()
                    except (BrokenPipeError, OSError, ValueError): pass
                feeder = threading.Thread(target=feed, daemon=True); feeder.start()
                start = time.monotonic()
                while process.poll() is None:
                    if overflow.is_set() or time.monotonic() - start > self.timeout:
                        try:
                            subprocess.run(["docker", "rm", "--force", name], stdout=subprocess.DEVNULL,
                                           stderr=subprocess.DEVNULL, timeout=10)
                        except subprocess.TimeoutExpired:
                            process.kill(); process.wait(timeout=5)
                            raise RuntimeError("Predictor cleanup timed out; worker requires inspection")
                        if process.poll() is None: process.kill()
                        process.wait(timeout=5)
                        raise RuntimeError("Predictor exceeded output or time limit")
                    time.sleep(.02)
                for reader in readers: reader.join(timeout=2)
                if overflow.is_set() or process.returncode:
                    raise RuntimeError("Predictor failed or exceeded output limit")
                outputs = json.loads(buffers["stdout"])
                if not isinstance(outputs, list) or len(outputs) != len(inputs):
                    raise ValueError("Predictor must return one prediction per attempted task")
                return outputs
            finally:
                if process.poll() is None: process.kill(); process.wait(timeout=5)
                try:
                    subprocess.run(["docker", "rm", "--force", name], stdout=subprocess.DEVNULL,
                                   stderr=subprocess.DEVNULL, timeout=10)
                except subprocess.TimeoutExpired:
                    raise RuntimeError("Predictor cleanup timed out; worker requires inspection")
                for stream in (process.stdin, process.stdout, process.stderr):
                    if stream is not None: stream.close()


class CodeImprovementEnv(MethodDesignEnv):
    """Actual Python predictor revision with budgeted development and frozen holdout."""

    def __init__(self, runner, **kwargs):
        super().__init__(**kwargs); self.runner = runner

    def _observation(self):
        obs = super()._observation()
        del obs["candidate_schema"]
        obs.update({"starter_code": STARTER_CODE, "source_limit_bytes": 65_536,
                    "predictor_contract": "JSON list of public tasks in; aligned prediction list out",
                    "sandbox_image": self.runner.image_id})
        return obs

    def _score_source(self, source, tasks):
        try:
            outputs = self.runner.run(source, [t.inputs for t in tasks])
        except (Exception,):
            return {"loss": 1.0, "failures": len(tasks), "attempted": len(tasks)}
        losses, failures = [], 0
        for prediction, task in zip(outputs, tasks):
            try: losses.append(brier_loss(prediction, task))
            except (ValueError, TypeError, KeyError, OverflowError): losses.append(1.0); failures += 1
        return {"loss": sum(losses) / len(losses), "failures": failures, "attempted": len(tasks)}

    def step(self, action):
        if self._done: raise RuntimeError("Reset required before stepping")
        if (not isinstance(action, dict) or set(action) != {"type", "source"}
                or action["type"] not in ("evaluate", "submit") or not isinstance(action["source"], str)
                or not 1 <= len(action["source"].encode()) <= 65_536
                or (action["type"] == "evaluate" and len(self._history) >= self.query_budget)):
            self._done = True
            return self._observation(), -1.0, True, False, {"status": "invalid_action"}
        report = self._score_source(action["source"], self._dev if action["type"] == "evaluate" else self._final)
        if action["type"] == "evaluate":
            self._history.append({"source": action["source"], "metrics": report})
            return self._observation(), -self.query_cost, False, False, {"status": "evaluated"}
        self._done = True
        return self._observation(), 1 - report["loss"], True, False, {"status": "submitted", "attempted": report["attempted"], "failures": report["failures"]}
