"""Gym-style core API without a Gymnasium dependency or execution privileges."""

from __future__ import annotations

import copy
import hashlib
import hmac
import secrets
from typing import Any

import numpy as np

from .tasks import TaskFactory, PrivateTask, public_task, canonical
from .scoring import brier_loss, panel_loss
from .baselines import DEFAULT, candidate, predict


class InferenceEnv:
    """One prediction, one terminal reward; assessor objects stay in this process."""

    def __init__(self, factory=None, split="train"):
        if split not in ("train", "validation", "test"):
            raise ValueError("Invalid split")
        self.factory = factory if factory is not None else TaskFactory()
        self.split = split
        self._episode = -1
        self._task = None
        self._done = True

    def reset(self, *, seed=None, options=None):
        if options:
            raise ValueError("No reset options accepted")
        self._episode = self._episode + 1 if seed is None else seed
        self._task = self.factory.make(self.split, self._episode)
        self._done = False
        return self._observation(), {"contract": "finemap-rl-inference-v1"}

    def reset_task(self, task: PrivateTask):
        """Assessor-only frozen challenge hook, never exposed as a client operation."""
        if not isinstance(task, PrivateTask):
            raise ValueError("Expected an assessor task")
        self._task, self._done = task, False
        return self._observation(), {"contract": "finemap-rl-inference-v1"}

    def _observation(self):
        return {"inputs": public_task(self._task.inputs)}

    def _require_active(self):
        if self._task is None or self._done:
            raise RuntimeError("Reset required before stepping")

    def _invalid(self):
        self._done = True
        return self._observation(), -1.0, True, False, {"status": "invalid_action"}

    def step(self, action):
        self._require_active()
        try:
            if not isinstance(action, dict) or set(action) != {"type", "prediction"} or action["type"] != "submit":
                raise ValueError("Expected submit action")
            loss = brier_loss(action["prediction"], self._task)
        except (ValueError, TypeError, KeyError, OverflowError):
            return self._invalid()
        self._done = True
        # Do not expose label-dependent submetrics even in terminal info.
        return self._observation(), 1 - loss, True, False, {"status": "submitted"}

    def state(self):
        """Public lifecycle metadata only, never a serialised object/__dict__."""
        return {"active": not self._done, "task_id": self._task.inputs["task_id"] if self._task else None}

    def close(self):
        self._done = True


class ExperimentEnv(InferenceEnv):
    """Allocate noisy sparse-membership assays, then report calibrated probabilities.

    This is a stylised confirmatory assay channel, not a biological MPRA oracle.
    Costs are utility units, not empirical money or wall-clock measurements.
    """

    def __init__(self, factory=None, split="train", budget=6, max_repeats=2,
                 sensitivity=.85, specificity=.9, cost=.003):
        super().__init__(factory, split)
        if not isinstance(budget, int) or not 0 <= budget <= 64 or not isinstance(max_repeats, int) or not 1 <= max_repeats <= 16:
            raise ValueError("Invalid assay budget or repetition bound")
        if not .5 < sensitivity < 1 or not .5 < specificity < 1 or not 0 <= cost <= .1:
            raise ValueError("Invalid assay model")
        self.budget, self.max_repeats = budget, max_repeats
        self.sensitivity, self.specificity, self.cost = sensitivity, specificity, cost

    def reset(self, *, seed=None, options=None):
        # Initialisation must precede the base observation call.
        self._attempts = np.zeros(self.factory.p, dtype=int)
        self._successes = np.zeros(self.factory.p, dtype=int)
        self._remaining = self.budget
        obs, info = super().reset(seed=seed, options=options)
        digest = hmac.new(self.factory._key, f"assay:{self.split}:{self._episode}".encode(), hashlib.sha256).digest()
        self._assay_rng = np.random.default_rng(int.from_bytes(digest[:16], "big"))
        return obs, {"contract": "finemap-rl-experiment-v1"}

    def reset_task(self, task):
        self._attempts = np.zeros(len(task.inputs["z"]), dtype=int)
        self._successes = self._attempts.copy(); self._remaining = self.budget
        self._assay_rng = np.random.default_rng(secrets.randbits(128))
        return super().reset_task(task)

    def _observation(self):
        return {"inputs": public_task(self._task.inputs), "budget_remaining": self._remaining,
                "assay_attempts": self._attempts.tolist(), "assay_successes": self._successes.tolist(),
                "legal_assays": ((self._attempts < self.max_repeats) & (self._remaining > 0)).tolist(),
                "assay_model": {"sensitivity": self.sensitivity, "specificity": self.specificity,
                                "cost_per_assay": self.cost, "target": "sparse-membership"}}

    def step(self, action):
        self._require_active()
        if isinstance(action, dict) and action.get("type") == "submit":
            return super().step(action)
        try:
            if not isinstance(action, dict) or set(action) != {"type", "variant_id"} or action["type"] != "assay":
                raise ValueError("Expected assay or submit")
            j = self._task.inputs["variant_id"].index(action["variant_id"])
            if self._remaining < 1 or self._attempts[j] >= self.max_repeats:
                raise ValueError("Assay unavailable")
        except (ValueError, TypeError, KeyError):
            return self._invalid()
        probability = self.sensitivity if self._task.labels[j] else 1 - self.specificity
        self._successes[j] += int(self._assay_rng.random() < probability)
        self._attempts[j] += 1; self._remaining -= 1
        return self._observation(), -self.cost, False, False, {"status": "assayed"}


class MethodDesignEnv:
    """Budgeted optimisation of executable declarative inference programs.

    Development losses can guide improvement. Final loci are separate and never
    appear in observations. Every final task contributes, including failures.
    Free-form code is handled only by the separate Docker runner, not eval/exec.
    """

    def __init__(self, factory=None, panel_size=12, query_budget=4, query_cost=.002):
        self.factory = factory if factory is not None else TaskFactory()
        if not isinstance(panel_size, int) or not 2 <= panel_size <= 128 or not isinstance(query_budget, int) or not 1 <= query_budget <= 16:
            raise ValueError("Invalid design budget or panel size")
        if not 0 <= query_cost <= .1:
            raise ValueError("Invalid validation query cost")
        self.panel_size, self.query_budget, self.query_cost = panel_size, query_budget, query_cost
        self._episode, self._done = -1, True

    def reset(self, *, seed=None, options=None):
        if options: raise ValueError("No reset options accepted")
        self._episode = self._episode + 1 if seed is None else seed
        if isinstance(self._episode, bool) or not isinstance(self._episode, int) or not 0 <= self._episode < 2**40:
            raise ValueError("Invalid design episode")
        indexes = [self._episode * self.panel_size + i for i in range(self.panel_size)]
        self._dev = [self.factory.make("design-dev", i) for i in indexes]
        self._final = [self.factory.make("design-final", i) for i in indexes]
        self._history, self._done = [], False
        return self._observation(), {"contract": "finemap-rl-method-design-v1"}

    def _observation(self):
        return {"development_inputs": [public_task(t.inputs) for t in self._dev],
                "query_budget_remaining": self.query_budget - len(self._history),
                "validation_history": copy.deepcopy(self._history),
                "candidate_schema": copy.deepcopy(DEFAULT), "query_cost": self.query_cost,
                "final_panel_size": self.panel_size}

    def step(self, action):
        if self._done: raise RuntimeError("Reset required before stepping")
        try:
            if not isinstance(action, dict) or set(action) != {"type", "candidate"} or action["type"] not in ("evaluate", "submit"):
                raise ValueError("Invalid design action")
            cfg = candidate(action["candidate"])
            if action["type"] == "evaluate" and len(self._history) >= self.query_budget:
                raise ValueError("Validation query budget exhausted")
        except (ValueError, TypeError, KeyError, OverflowError):
            self._done = True
            return self._observation(), -1.0, True, False, {"status": "invalid_action"}
        predictor = lambda inputs: predict(inputs, cfg)
        if action["type"] == "evaluate":
            report = panel_loss(predictor, self._dev)
            self._history.append({"candidate": cfg, "metrics": report})
            return self._observation(), -self.query_cost, False, False, {"status": "evaluated"}
        report = panel_loss(predictor, self._final)
        self._done = True
        return self._observation(), 1 - report["loss"], True, False, {"status": "submitted", "attempted": report["attempted"], "failures": report["failures"]}

    def state(self):
        return {"active": not self._done, "queries_used": len(getattr(self, "_history", []))}

    def close(self): self._done = True
