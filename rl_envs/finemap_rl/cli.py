"""Bounded JSON-lines assessor protocol and reproducible development smoke runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

from .tasks import TaskFactory, SPECIFICATION, canonical, SCENARIOS
from .environments import InferenceEnv, ExperimentEnv, MethodDesignEnv
from .baselines import DEFAULT, predict, uncertainty_assay, value_information_assay


def serve(env, input_stream, output_stream, split="train", max_episodes=100):
    """Run in an assessor process; remote agents receive this protocol only.

    No arbitrary attribute access, filesystem commands, truth retrieval or code
    execution is exposed. Terminal score feedback inherently conveys information:
    never reuse a final episode for tuning and cap issued final episodes.
    """
    issued = 0
    while True:
        line = input_stream.readline(1_048_577)
        if not line: break
        try:
            if len(line.encode()) > 1_048_576 or not line.endswith("\n"):
                raise ValueError("Request too large")
            request = json.loads(line, parse_constant=lambda value: (_ for _ in ()).throw(ValueError("Nonfinite JSON")))
            if not isinstance(request, dict): raise ValueError("Invalid request")
            if request.get("op") == "reset" and set(request) <= {"op", "seed"}:
                if issued >= max_episodes: raise ValueError("Episode cap reached")
                if split == "test" and "seed" in request: raise ValueError("Final seeds are assessor-controlled")
                if split == "test" and env.state()["active"]:
                    raise ValueError("A final task must be submitted before issuing another")
                obs, info = env.reset(seed=request.get("seed"))
                issued += 1
                result = dict(observation=obs, reward=0., terminated=False, truncated=False, info=info)
            elif request.get("op") == "step" and set(request) == {"op", "action"}:
                obs, reward, terminated, truncated, info = env.step(request["action"])
                result = dict(observation=obs, reward=reward, terminated=terminated, truncated=truncated, info=info)
            elif request == {"op": "state"}:
                result = {"state": env.state()}
            else: raise ValueError("Unknown request")
        except (ValueError, TypeError, KeyError, RuntimeError, OverflowError):
            result = {"error": "invalid_request"}
        output_stream.write(canonical(result).decode() + "\n"); output_stream.flush()
        # Oversized unterminated input closes the protocol, preventing the
        # remainder from being parsed as an independent request.
        if not line.endswith("\n"): break
    env.close()


def features(inputs):
    z = np.abs(np.asarray(inputs["z"]))
    ld = np.asarray(inputs["LD"])
    return np.asarray([1., min(z.max(), 15) / 15, min(z.mean(), 10) / 10,
                       np.mean(np.abs(ld - np.eye(len(z))))])


class ContextualBandit:
    """A real one-step RL policy update; no claim of LLM post-training."""
    def __init__(self, catalog, seed=7):
        self.catalog, self.rng = catalog, np.random.default_rng(seed)
        self.a = np.tile(np.eye(4), (len(catalog), 1, 1))
        self.b = np.zeros((len(catalog), 4)); self.counts = np.zeros(len(catalog), dtype=int)

    def choose(self, inputs, explore=False):
        phi = features(inputs)
        weights = np.stack([np.linalg.solve(a, b) for a, b in zip(self.a, self.b)])
        values = weights @ phi
        if explore:
            uncertainty = np.asarray([np.sqrt(phi @ np.linalg.solve(a, phi)) for a in self.a])
            values += .15 * uncertainty
            if self.rng.random() < .12: return int(self.rng.integers(len(self.catalog)))
        return int(values.argmax())

    def update(self, inputs, action, reward):
        phi = features(inputs); self.a[action] += np.outer(phi, phi)
        self.b[action] += phi * reward; self.counts[action] += 1


def smoke(output_dir, episodes=250, test_episodes=100, p=12, include_docker=False):
    directory = Path(output_dir)
    if directory.exists(): raise ValueError("Refusing to overwrite smoke artifact directory")
    directory.mkdir(parents=True)
    # Deliberately public development key. This smoke is not a blind final panel.
    factory = TaskFactory(b"public-development-smoke-key-v1", p=p)
    catalog = [dict(DEFAULT), {**DEFAULT, "algorithm": "independent_abf"},
               {**DEFAULT, "count_prior": [.8, .15, .05]},
               {**DEFAULT, "count_prior": [.1, .45, .45]},
               {**DEFAULT, "annotation_weight": 1.}, {**DEFAULT, "annotation_weight": -1.},
               {**DEFAULT, "effect_sd": .02}, {**DEFAULT, "effect_sd": .1}]
    learner = ContextualBandit(catalog)
    training = []; env = InferenceEnv(factory)
    start = time.monotonic()
    for episode in range(episodes):
        obs, _ = env.reset(seed=episode)
        index = learner.choose(obs["inputs"], explore=True)
        _, reward, done, _, _ = env.step({"type": "submit", "prediction": predict(obs["inputs"], catalog[index])})
        learner.update(obs["inputs"], index, reward)
        training.append({"episode": episode, "action": index, "reward": reward})
        if not done: raise AssertionError("Inference episode did not terminate")
    rows, trajectory = [], []
    for episode in range(test_episodes):
        for method in ("working_exact", "independent_abf", "prior", "learned_bandit", "active_uncertainty", "active_model_voi"):
            core = ExperimentEnv(factory, "validation") if method.startswith("active_") else InferenceEnv(factory, "validation")
            obs, _ = core.reset(seed=episode)
            total, steps = 0., 0
            while True:
                if method == "active_uncertainty": action = uncertainty_assay(obs)
                elif method == "active_model_voi": action = value_information_assay(obs)
                else:
                    cfg = (catalog[learner.choose(obs["inputs"])] if method == "learned_bandit" else
                           {**DEFAULT, "algorithm": "independent_abf"} if method == "independent_abf" else
                           {**DEFAULT, "algorithm": "prior"} if method == "prior" else DEFAULT)
                    action = {"type": "submit", "prediction": predict(obs["inputs"], cfg)}
                next_obs, reward, done, truncated, info = core.step(action)
                total += reward; steps += 1
                if episode == 0:
                    trajectory.append({"method": method, "observation": obs, "action": action,
                                       "reward": reward, "terminated": done, "truncated": truncated, "info": info})
                obs = next_obs
                if done or truncated: break
            rows.append({"episode": episode, "method": method, "utility": total, "steps": steps,
                         "terminal_brier": 1 - reward, "query_cost": reward - total})
    design = MethodDesignEnv(factory, panel_size=6, query_budget=3)
    obs, _ = design.reset(seed=0)
    design_trajectory = []; budgeted = catalog[:3]
    for cfg in budgeted:
        action = {"type": "evaluate", "candidate": cfg}
        nxt, reward, done, truncated, info = design.step(action)
        design_trajectory.append(dict(observation=obs, action=action, reward=reward, terminated=done, info=info)); obs = nxt
    best = min(obs["validation_history"], key=lambda x: x["metrics"]["loss"])["candidate"]
    action = {"type": "submit", "candidate": best}
    nxt, reward, done, _, info = design.step(action)
    design_trajectory.append(dict(observation=obs, action=action, reward=reward, terminated=done, info=info))
    docker_report = None
    if include_docker:
        from .code_runner import CodeImprovementEnv, DockerPredictorRunner, STARTER_CODE, ABF_CODE
        runner = DockerPredictorRunner()
        code_env = CodeImprovementEnv(runner, factory=factory, panel_size=6, query_budget=2)
        obs, _ = code_env.reset(seed=0)
        validation_rewards = 0.
        for source in (STARTER_CODE, ABF_CODE):
            nxt, validation_reward, _, _, _ = code_env.step({"type": "evaluate", "source": source})
            validation_rewards += validation_reward
        chosen = min(nxt["validation_history"], key=lambda h: h["metrics"]["loss"])["source"]
        _, reward, done, _, info = code_env.step({"type": "submit", "source": chosen})
        docker_report = {"image_id": runner.image_id,
                         "development": [h["metrics"] for h in nxt["validation_history"]],
                         "selected": "independent_abf" if chosen == ABF_CODE else "forced_one_starter",
                         "heldout_utility": reward + validation_rewards, "terminated": done, "info": info}
    for filename, data in (("training.csv", training), ("evaluation.csv", rows)):
        with (directory / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0])); writer.writeheader(); writer.writerows(data)
    rng = np.random.default_rng(19); summaries = []
    for method in sorted({r["method"] for r in rows}):
        values = np.asarray([r["utility"] for r in rows if r["method"] == method])
        boot = np.asarray([values[rng.integers(len(values), size=len(values))].mean() for _ in range(1000)])
        summaries.append(dict(method=method, n=len(values), mean_utility=float(values.mean()),
                              bootstrap_low=float(np.quantile(boot, .025)), bootstrap_high=float(np.quantile(boot, .975))))
    (directory / "summary.json").write_bytes(canonical(summaries))
    reference = np.asarray([r["utility"] for r in rows if r["method"] == "working_exact"])
    comparisons = []
    for method in sorted({r["method"] for r in rows} - {"working_exact"}):
        difference = np.asarray([r["utility"] for r in rows if r["method"] == method]) - reference
        draws = np.asarray([difference[rng.integers(len(difference), size=len(difference))].mean() for _ in range(1000)])
        comparisons.append({"method": method, "reference": "working_exact", "paired_loci": test_episodes,
                            "mean_utility_difference": float(difference.mean()),
                            "bootstrap_low": float(np.quantile(draws,.025)), "bootstrap_high": float(np.quantile(draws,.975))})
    (directory / "paired-comparisons.json").write_bytes(canonical(comparisons))
    (directory / "example-trajectories.json").write_bytes(canonical(trajectory + design_trajectory))
    provenance = {"specification": SPECIFICATION, "development_key_public": True, "blind": False,
                  "episodes": episodes, "test_episodes": test_episodes, "p": p, "scenarios": SCENARIOS,
                  "numpy_version": np.__version__, "python_version": sys.version,
                  "policy": "LinUCB-style contextual bandit; numerical RL, not LLM training",
                  "action_counts": learner.counts.tolist(), "parameter_norm": float(np.linalg.norm(learner.b)),
                  "elapsed_seconds": time.monotonic() - start, "docker": docker_report,
                  "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                    for path in sorted(Path(__file__).parent.glob("*.py"))}}
    (directory / "provenance.json").write_bytes(canonical(provenance))
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("smoke"); s.add_argument("--output", required=True)
    s.add_argument("--episodes", type=int, default=250); s.add_argument("--test-episodes", type=int, default=100)
    s.add_argument("--p", type=int, default=12); s.add_argument("--docker", action="store_true")
    v = sub.add_parser("serve"); v.add_argument("--environment", choices=("inference", "experiment", "design", "code"), default="inference")
    v.add_argument("--split", choices=("train", "validation", "test"), default="train")
    v.add_argument("--key-file"); v.add_argument("--p", type=int, default=16)
    v.add_argument("--max-episodes", type=int, default=100)
    args = parser.parse_args()
    if args.command == "smoke":
        if args.episodes < 1 or args.test_episodes < 2: parser.error("Positive training and at least two evaluation episodes required")
        print(json.dumps(smoke(args.output, args.episodes, args.test_episodes, args.p, args.docker), indent=2))
    else:
        if not 1 <= args.max_episodes <= 100000: parser.error("Invalid episode cap")
        key = Path(args.key_file).read_bytes() if args.key_file else None
        factory = TaskFactory(key, p=args.p)
        if args.environment == "inference": env = InferenceEnv(factory, args.split)
        elif args.environment == "experiment": env = ExperimentEnv(factory, args.split)
        elif args.environment == "design": env = MethodDesignEnv(factory)
        else:
            from .code_runner import DockerPredictorRunner, CodeImprovementEnv
            env = CodeImprovementEnv(DockerPredictorRunner(), factory=factory)
        serve(env, sys.stdin, sys.stdout, args.split, args.max_episodes)


if __name__ == "__main__": main()
