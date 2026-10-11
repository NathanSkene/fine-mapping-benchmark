"""Gymnasium 1.3 wrappers with numeric spaces and deterministic seeded resets."""

from itertools import product
import copy

import gymnasium as gym
from gymnasium import spaces
import numpy as np

from .tasks import TaskFactory
from .environments import InferenceEnv, ExperimentEnv, MethodDesignEnv
from .baselines import DEFAULT


def task_vector(inputs):
    return np.concatenate([inputs["z"], np.asarray(inputs["LD"]).ravel(),
                           np.asarray(inputs["annotations"]).ravel()]).astype(np.float64)


class FineMappingGym(gym.Env):
    metadata = {"render_modes": []}

    def __init__(self, p=16, experiment=False, factory=None, budget=6):
        self.p, self.experiment = p, experiment
        factory = factory if factory is not None else TaskFactory(p=p)
        if factory.p != p: raise ValueError("Factory dimension mismatch")
        self.core = ExperimentEnv(factory, budget=budget) if experiment else InferenceEnv(factory)
        length = p*p + 3*p + (2*p + 1 if experiment else 0)
        self.observation_space = spaces.Box(-np.inf, np.inf, shape=(length,), dtype=np.float64)
        self.action_space = spaces.Dict({"operation": spaces.Discrete(2 if experiment else 1),
                                         "variant": spaces.Discrete(p),
                                         "probabilities": spaces.Box(0, 1, shape=(p + 1,), dtype=np.float64)})

    def _vector(self, observation):
        x = task_vector(observation["inputs"])
        if self.experiment:
            x = np.concatenate([x, observation["assay_attempts"], observation["assay_successes"],
                                [observation["budget_remaining"]]])
        return x

    def reset(self, *, seed=None, options=None):
        super().reset(seed=seed)
        obs, info = self.core.reset(seed=seed, options=options)
        return self._vector(obs), {**info, "public_observation": obs}

    def step(self, action):
        if not self.action_space.contains(action):
            raise ValueError("Action outside declared Gymnasium space")
        if action["operation"] == 0:
            probabilities = action["probabilities"]
            command = {"type": "submit", "prediction": {
                "task_id": self.core._task.inputs["task_id"], "variant_id": list(self.core._task.inputs["variant_id"]),
                "pip": probabilities[:-1].tolist(), "null_probability": float(probabilities[-1])}}
        else:
            command = {"type": "assay", "variant_id": self.core._task.inputs["variant_id"][int(action["variant"])]}
        obs, reward, terminated, truncated, info = self.core.step(command)
        return self._vector(obs), reward, terminated, truncated, {**info, "public_observation": obs}

    def close(self): self.core.close()


def candidate_catalog():
    result = []
    for algorithm, sd, counts, weight in product(
            ("exact", "independent_abf"), (.02, .05, .1),
            ([.5, .3, .2], [.8, .15, .05], [.1, .45, .45]), (0., 1., -1.)):
        result.append({**DEFAULT, "algorithm": algorithm, "effect_sd": sd,
                       "count_prior": list(counts), "annotation_weight": weight})
    return result


class MethodDesignGym(gym.Env):
    metadata = {"render_modes": []}

    def __init__(self, p=8, panel_size=6, query_budget=4, factory=None):
        self.catalog = candidate_catalog()
        self.core = MethodDesignEnv(factory if factory is not None else TaskFactory(p=p), panel_size, query_budget)
        self.action_space = spaces.Discrete(2 * len(self.catalog))
        length = panel_size * (p*p + 3*p) + query_budget * 4 + 1
        self.observation_space = spaces.Box(-np.inf, np.inf, shape=(length,), dtype=np.float64)

    def _vector(self, obs):
        history = np.zeros((self.core.query_budget, 4))
        for j, item in enumerate(obs["validation_history"]):
            cfg = item["candidate"]
            index = self.catalog.index(cfg)
            metrics = item["metrics"]
            history[j] = [1, index / len(self.catalog), metrics["loss"], metrics["failures"] / metrics["attempted"]]
        return np.concatenate([*(task_vector(x) for x in obs["development_inputs"]),
                               history.ravel(), [obs["query_budget_remaining"]]])

    def reset(self, *, seed=None, options=None):
        super().reset(seed=seed)
        obs, info = self.core.reset(seed=seed, options=options)
        return self._vector(obs), {**info, "public_observation": obs}

    def step(self, action):
        if not self.action_space.contains(action): raise ValueError("Invalid discrete action")
        index = int(action) % len(self.catalog)
        operation = "evaluate" if action < len(self.catalog) else "submit"
        obs, reward, terminated, truncated, info = self.core.step(
            {"type": operation, "candidate": copy.deepcopy(self.catalog[index])})
        return self._vector(obs), reward, terminated, truncated, {**info, "public_observation": obs}

    def close(self): self.core.close()


def register_environments():
    for env_id, entry, kwargs in (
        ("FineMappingInference-v0", "finemap_rl.gym_adapter:FineMappingGym", {}),
        ("FineMappingExperiment-v0", "finemap_rl.gym_adapter:FineMappingGym", {"experiment": True}),
        ("FineMappingMethodDesign-v0", "finemap_rl.gym_adapter:MethodDesignGym", {})):
        if env_id not in gym.envs.registry:
            gym.register(env_id, entry_point=entry, kwargs=kwargs)
