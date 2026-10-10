"""Optional OpenEnv 0.3 adapter using its installed typed server interfaces."""

from typing import Any
from pydantic import ConfigDict
from openenv.core.env_server import Environment, create_app
from openenv.core.env_server.types import Action, Observation, State

from .tasks import TaskFactory, canonical
from .environments import InferenceEnv, ExperimentEnv, MethodDesignEnv


class FineMappingAction(Action):
    model_config = ConfigDict(extra="forbid")
    command: dict[str, Any]


class FineMappingObservation(Observation):
    model_config = ConfigDict(extra="forbid")
    payload: dict[str, Any]


class FineMappingState(State):
    model_config = ConfigDict(extra="forbid")
    active: bool = False


class FineMappingOpenEnv(Environment[FineMappingAction, FineMappingObservation, FineMappingState]):
    SUPPORTS_CONCURRENT_SESSIONS = True
    def __init__(self, environment="inference", p=16, split="train", max_episodes=100):
        super().__init__()
        if environment not in ("inference", "experiment", "design") or not 1 <= max_episodes <= 100000:
            raise ValueError("Invalid OpenEnv configuration")
        # A fresh session key prevents a restart/new session from recreating old
        # final labels. A frozen multi-agent competition needs coordinator custody.
        factory = TaskFactory(p=p)
        self.core = (InferenceEnv(factory, split) if environment == "inference" else
                     ExperimentEnv(factory, split) if environment == "experiment" else MethodDesignEnv(factory))
        self.split, self.max_episodes = split, max_episodes
        self._issued, self._steps, self._task_id = 0, 0, None

    def reset(self, seed=None, episode_id=None, **kwargs):
        if (kwargs or (self.split == "test" and seed is not None) or self._issued >= self.max_episodes
                or (self.split == "test" and self.core.state()["active"])):
            raise ValueError("Invalid reset or episode cap")
        observation, info = self.core.reset(seed=seed)
        self._issued += 1; self._steps = 0
        self._task_id = observation.get("inputs", {}).get("task_id")
        if self._task_id is None:
            self._task_id = observation["development_inputs"][0]["task_id"]
        return FineMappingObservation(payload=observation, reward=0., done=False, metadata=info)

    def step(self, action, timeout_s=None, **kwargs):
        if kwargs or action.metadata or len(canonical(action.model_dump())) > 1_048_576:
            raise ValueError("Unsupported action metadata or request size")
        observation, reward, terminated, truncated, info = self.core.step(action.command)
        self._steps += 1
        return FineMappingObservation(payload=observation, reward=reward, done=terminated or truncated,
                                      metadata={**info, "terminated": terminated, "truncated": truncated})

    @property
    def state(self):
        return FineMappingState(episode_id=self._task_id, step_count=self._steps,
                                active=self.core.state()["active"])

    def close(self): self.core.close()


def create_finemap_app(environment="inference", p=16, split="train", max_episodes=100):
    return create_app(lambda: FineMappingOpenEnv(environment, p, split, max_episodes),
                      FineMappingAction, FineMappingObservation, env_name="finemap_rl", max_concurrent_envs=8)


app = create_finemap_app()
