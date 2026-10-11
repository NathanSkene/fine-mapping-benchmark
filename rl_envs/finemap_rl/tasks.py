"""Procedural score experiments with a narrow, explicit public schema.

This is an independent NumPy implementation of the repository's Gaussian
score-statistic experiment, not a replacement for its haplotype generator.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import hmac
import json
import math
import os
from pathlib import Path
from typing import Any

import numpy as np

SPECIFICATION = "finemap-rl-gaussian-v1"
SCENARIOS = ("matched_prior", "null", "multi_sparse", "many_sparse",
             "opposing_effects", "sparse_inf", "ld_reference_shift",
             "perfect_ld", "annotation_informative", "annotation_reversal",
             "allele_error")
PUBLIC_KEYS = {"task_id", "variant_id", "z", "LD", "annotations", "n", "statistic"}


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def public_task(value: dict[str, Any]) -> dict[str, Any]:
    """Validate, copy and reject extra fields rather than passing latent metadata."""
    if not isinstance(value, dict) or set(value) != PUBLIC_KEYS:
        raise ValueError("Public task must contain exactly the documented input fields")
    ids = value["variant_id"]
    if not isinstance(ids, list) or not 2 <= len(ids) <= 512 or len(set(ids)) != len(ids):
        raise ValueError("Invalid variant identities or task size")
    if any(not isinstance(i, str) or not i or len(i) > 160 for i in ids):
        raise ValueError("Invalid variant identity")
    if not isinstance(value["task_id"], str) or not value["task_id"]:
        raise ValueError("Invalid task identity")
    p = len(ids)
    z = np.asarray(value["z"], dtype=float)
    ld = np.asarray(value["LD"], dtype=float)
    anno = np.asarray(value["annotations"], dtype=float)
    if z.shape != (p,) or ld.shape != (p, p) or anno.shape != (p, 2):
        raise ValueError("Input dimensions do not align")
    if not all(np.isfinite(a).all() for a in (z, ld, anno)):
        raise ValueError("Non-finite public input")
    if not np.allclose(ld, ld.T, atol=1e-8) or not np.allclose(np.diag(ld), 1, atol=1e-8):
        raise ValueError("LD must be a symmetric correlation matrix")
    if np.linalg.eigvalsh(ld).min() < -1e-7:
        raise ValueError("LD must be positive semidefinite")
    n = value["n"]
    if isinstance(n, bool) or not isinstance(n, (int, float)) or not math.isfinite(n) or n < 3:
        raise ValueError("Invalid n")
    if value["statistic"] != "known-noise-score":
        raise ValueError("Unsupported statistic convention")
    return json.loads(canonical(value))


@dataclass(frozen=True)
class PrivateTask:
    inputs: dict[str, Any]
    labels: np.ndarray
    null: bool
    scenario: str
    group_id: str

    def __post_init__(self):
        public_task(self.inputs)
        labels = np.asarray(self.labels, dtype=float)
        if labels.shape != (len(self.inputs["z"]),) or not np.isin(labels, (0, 1)).all():
            raise ValueError("Invalid sparse-membership labels")
        if self.null and labels.any():
            raise ValueError("A null task cannot have an observed sparse target")
        labels = labels.copy()
        labels.flags.writeable = False
        object.__setattr__(self, "labels", labels)


class TaskFactory:
    """Keep the key assessor-side; split namespaces yield independent loci.

    A public episode number never suffices to reconstruct a private RNG seed.
    Fixed keys support development reproducibility. Final keys must be fresh,
    unreleased and inaccessible to agents; a Python object is not a sandbox.
    """

    def __init__(self, key: bytes | None = None, p: int = 16, n: int = 10000,
                 scenarios: tuple[str, ...] = SCENARIOS):
        if not isinstance(p, int) or not 5 <= p <= 128 or not isinstance(n, int) or n < 3:
            raise ValueError("Unsupported generator dimensions")
        if not scenarios or any(s not in SCENARIOS for s in scenarios):
            raise ValueError("Unsupported scenarios")
        if key is not None and (not isinstance(key, bytes) or len(key) < 16):
            raise ValueError("Assessor key needs at least 128 bits")
        self._key = key if key is not None else os.urandom(32)
        self.p, self.n, self.scenarios = p, n, tuple(scenarios)

    def rng(self, domain: str, episode: int) -> np.random.Generator:
        if domain not in ("train", "validation", "test", "assay", "design-dev", "design-final"):
            raise ValueError("Invalid split/domain")
        if isinstance(episode, bool) or not isinstance(episode, int) or not 0 <= episode < 2**63:
            raise ValueError("Episode must be a nonnegative integer")
        digest = hmac.new(self._key, f"{SPECIFICATION}:{domain}:{episode}".encode(), hashlib.sha256).digest()
        return np.random.default_rng(int.from_bytes(digest[:16], "big"))

    def make(self, split: str, episode: int) -> PrivateTask:
        if split not in ("train", "validation", "test", "design-dev", "design-final"):
            raise ValueError("Invalid task split")
        rng = self.rng(split, episode)
        p, n = self.p, self.n
        scenario = self.scenarios[int(rng.integers(len(self.scenarios)))]
        rho = .97 if scenario == "opposing_effects" else rng.uniform(.35, .85)
        ld = rho ** np.abs(np.subtract.outer(np.arange(p), np.arange(p)))
        if scenario == "perfect_ld":
            ld[:, 1] = ld[:, 0]; ld[1, :] = ld[0, :]; ld[1, 1] = 1
        annotations = rng.normal(size=(p, 2))
        k = int(rng.choice(3, p=[.5, .3, .2])) if scenario == "matched_prior" else 1
        if scenario == "null": k = 0
        if scenario in ("multi_sparse", "sparse_inf", "opposing_effects"): k = 2
        if scenario == "many_sparse": k = 5
        weights = np.ones(p)
        if scenario in ("annotation_informative", "annotation_reversal"):
            eta = annotations[:, 0] * (1 if scenario == "annotation_informative" else -1)
            weights = np.exp(eta - eta.max())
        causal = rng.choice(p, size=k, replace=False, p=weights / weights.sum())
        if scenario == "perfect_ld": causal = np.array([int(rng.integers(2))])
        if scenario == "opposing_effects": causal = np.array([0, 1])
        beta = np.zeros(p)
        if scenario == "matched_prior": beta[causal] = rng.normal(0, .05, len(causal))
        else: beta[causal] = rng.choice([-1, 1], len(causal)) * 4.5 / np.sqrt(n)
        if scenario == "opposing_effects": beta[causal] = np.array([1, -1]) * 4.5 / np.sqrt(n)
        background = rng.normal(0, 1.25 / np.sqrt(n), p) if scenario == "sparse_inf" else np.zeros(p)
        eig, vectors = np.linalg.eigh(ld)
        noise = vectors @ (np.sqrt(np.maximum(eig, 0)) * rng.normal(size=p))
        z = np.sqrt(n) * ld @ (beta + background) + noise
        fit_ld = .2 ** np.abs(np.subtract.outer(np.arange(p), np.arange(p))) if scenario == "ld_reference_shift" else ld
        if scenario == "allele_error": z[causal[0]] *= -1
        # Random ordering prevents positional priors for opposing/perfect-LD targets.
        order = rng.permutation(p)
        token = hmac.new(self._key, f"identity:{split}:{episode}".encode(), hashlib.sha256).hexdigest()[:24]
        labels = np.zeros(p); labels[causal] = 1
        inputs = dict(task_id=token, variant_id=[f"{token}:v{i}" for i in range(p)],
                      z=z[order].tolist(), LD=fit_ld[np.ix_(order, order)].tolist(),
                      annotations=annotations[order].tolist(), n=n, statistic="known-noise-score")
        return PrivateTask(public_task(inputs), labels[order], not len(causal), scenario, token)


def write_release(tasks: list[PrivateTask], public_dir: str | Path, private_dir: str | Path) -> dict:
    """Freeze JSON-only inputs and separate assessor truth, refusing overlap/overwrite."""
    pub, priv = Path(public_dir).resolve(), Path(private_dir).resolve()
    if pub == priv or pub in priv.parents or priv in pub.parents:
        raise ValueError("Public and private directories overlap")
    if pub.exists() or priv.exists():
        raise ValueError("Refusing to overwrite release")
    if not tasks or len({t.inputs["task_id"] for t in tasks}) != len(tasks):
        raise ValueError("Empty or duplicate task release")
    pub.mkdir(parents=True); priv.mkdir(parents=True, mode=0o700)
    inputs = [public_task(t.inputs) for t in tasks]
    truth = [{"task_id": t.inputs["task_id"], "labels": t.labels.tolist(), "null": t.null,
              "scenario": t.scenario, "group_id": t.group_id} for t in tasks]
    manifests = {}
    for directory, data, name in ((pub, inputs, "inputs.json"), (priv, truth, "truth.json")):
        payload = canonical(data); (directory / name).write_bytes(payload)
        manifest = {"specification": SPECIFICATION, "n_tasks": len(tasks),
                    "files": {name: hashlib.sha256(payload).hexdigest()}}
        (directory / "manifest.json").write_bytes(canonical(manifest))
        manifests[str(directory)] = manifest
    return manifests


def read_release(public_dir: str | Path, private_dir: str | Path | None = None):
    def read(directory, filename):
        directory = Path(directory)
        manifest = json.loads((directory / "manifest.json").read_text())
        if set(manifest) != {"specification", "n_tasks", "files"} or set(manifest["files"]) != {filename}:
            raise ValueError("Invalid manifest")
        if manifest["specification"] != SPECIFICATION:
            raise ValueError("Unsupported release specification")
        data = (directory / filename).read_bytes()
        if hashlib.sha256(data).hexdigest() != manifest["files"][filename]:
            raise ValueError("Release hash mismatch")
        rows = json.loads(data)
        if not isinstance(rows, list) or len(rows) != manifest["n_tasks"]:
            raise ValueError("Invalid release row count")
        return rows
    inputs = [public_task(t) for t in read(public_dir, "inputs.json")]
    if len({t["task_id"] for t in inputs}) != len(inputs):
        raise ValueError("Duplicate public task identity")
    if private_dir is None:
        return inputs
    truth = read(private_dir, "truth.json")
    if [t.get("task_id") for t in truth] != [t["task_id"] for t in inputs]:
        raise ValueError("Public/private identities or order do not align")
    if any(set(t) != {"task_id", "labels", "null", "scenario", "group_id"} for t in truth):
        raise ValueError("Invalid private schema")
    if any(not isinstance(t["null"], bool) for t in truth):
        raise ValueError("Invalid null labels")
    return [PrivateTask(x, t["labels"], t["null"], t["scenario"], t["group_id"]) for x, t in zip(inputs, truth)]
