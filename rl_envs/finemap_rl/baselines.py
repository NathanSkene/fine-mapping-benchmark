"""Inspectable CPU baselines and a restricted, safe method-design language."""

from itertools import combinations
import math
import numpy as np

from .tasks import public_task

DEFAULT = {"algorithm": "exact", "effect_sd": .05, "count_prior": [.5, .3, .2],
           "annotation_weight": 0.0, "ld_ridge": 0.0}


def candidate(value):
    if not isinstance(value, dict) or set(value) != set(DEFAULT):
        raise ValueError("Candidate must match the documented declarative program schema")
    if value["algorithm"] not in ("exact", "independent_abf", "prior"):
        raise ValueError("Unknown algorithm")
    for name, lo, hi in (("effect_sd", .001, .5), ("annotation_weight", -3, 3), ("ld_ridge", 0, .5)):
        x = value[name]
        if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or not lo <= x <= hi:
            raise ValueError("Invalid candidate parameter")
    q = np.asarray(value["count_prior"], dtype=float)
    if q.shape != (3,) or not np.isfinite(q).all() or (q < 0).any() or not np.isclose(q.sum(), 1, atol=1e-9):
        raise ValueError("Count prior requires three nonnegative probabilities summing to one")
    return {**value, "count_prior": q.tolist()}


def configuration_posterior(inputs, config=None, successes=None, attempts=None,
                            sensitivity=.85, specificity=.9):
    """Enumerate a working 0/1/2-effect posterior with optional noisy assays.

    Exact for the matched-prior score experiment and these assay parameters;
    a misspecified working model in stress tracks. Never a universal ceiling.
    """
    x = public_task(inputs); cfg = candidate(DEFAULT if config is None else config)
    p = len(x["z"])
    if p > 64:
        raise ValueError("Exact baseline limited to 64 variants for bounded evaluation cost")
    z, ld, n = np.asarray(x["z"]), np.asarray(x["LD"]), x["n"]
    ld = (1 - cfg["ld_ridge"]) * ld + cfg["ld_ridge"] * np.eye(p)
    eta = cfg["annotation_weight"] * np.asarray(x["annotations"])[:, 0]
    configurations = [()] + [(i,) for i in range(p)] + list(combinations(range(p), 2))
    sizes = np.asarray([len(s) for s in configurations])
    logits = np.full(len(configurations), -np.inf)
    w = cfg["effect_sd"] ** 2
    for k in range(3):
        indexes = np.flatnonzero(sizes == k)
        if cfg["count_prior"][k] == 0:
            continue
        log_weights = np.asarray([sum(eta[list(configurations[i])]) for i in indexes])
        m = log_weights.max(); denominator = m + np.log(np.exp(log_weights - m).sum())
        for i, prior_weight in zip(indexes, log_weights):
            s = configurations[i]
            log_bf = 0.0
            if s and cfg["algorithm"] != "prior":
                h = np.sqrt(n) * z[list(s)]
                q = n * ld[np.ix_(s, s)]
                a = np.eye(k) + w * q
                sign, logdet = np.linalg.slogdet(a)
                if sign <= 0:
                    raise ValueError("Invalid working information matrix")
                log_bf = -.5 * logdet + .5 * w * float(h @ np.linalg.solve(a, h))
            logits[i] = np.log(cfg["count_prior"][k]) + prior_weight - denominator + log_bf
    indicators = np.zeros((len(configurations), p), dtype=float)
    for i, s in enumerate(configurations): indicators[i, list(s)] = 1
    if successes is not None or attempts is not None:
        succ, att = np.asarray(successes), np.asarray(attempts)
        if succ.shape != (p,) or att.shape != (p,) or (succ < 0).any() or (att < succ).any():
            raise ValueError("Invalid assay history")
        prob = indicators * sensitivity + (1 - indicators) * (1 - specificity)
        logits += (np.log(prob) * succ + np.log1p(-prob) * (att - succ)).sum(axis=1)
    finite = np.isfinite(logits)
    if not finite.any(): raise ValueError("No supported configurations")
    m = logits[finite].max(); probs = np.exp(logits - m); probs /= probs.sum()
    prediction = {"task_id": x["task_id"], "variant_id": x["variant_id"], "pip": (probs @ indicators).tolist(),
                  "null_probability": float(probs[0])}
    return prediction, indicators, probs


def predict(inputs, config=None, **assays):
    cfg = candidate(DEFAULT if config is None else config)
    if cfg["algorithm"] != "independent_abf":
        return configuration_posterior(inputs, cfg, **assays)[0]
    if assays:
        raise ValueError("Independent ABF baseline does not implement assay integration")
    x = public_task(inputs); p = len(x["z"])
    prior = (cfg["count_prior"][1] + 2 * cfg["count_prior"][2]) / p
    if prior == 0: pi = np.zeros(p)
    else:
        eta = np.log(prior / (1 - prior)) + cfg["annotation_weight"] * np.asarray(x["annotations"])[:, 0]
        v = x["n"] * cfg["effect_sd"] ** 2
        log_bf = -.5 * np.log1p(v) + .5 * np.asarray(x["z"])**2 * v / (1 + v)
        t = np.clip(eta + log_bf, -700, 700)
        pi = 1 / (1 + np.exp(-t))
    return {"task_id": x["task_id"], "variant_id": x["variant_id"], "pip": pi.tolist(), "null_probability": float(np.prod(1 - pi))}


def uncertainty_assay(observation, config=None):
    """Acquire the highest marginal posterior-variance affordable unexhausted assay."""
    prediction = predict(observation["inputs"], config,
                         successes=observation["assay_successes"], attempts=observation["assay_attempts"],
                         sensitivity=observation["assay_model"]["sensitivity"],
                         specificity=observation["assay_model"]["specificity"])
    pip = np.asarray(prediction["pip"])
    legal = np.asarray(observation["legal_assays"], dtype=bool)
    if not legal.any() or observation["budget_remaining"] < 1:
        return {"type": "submit", "prediction": prediction}
    utility = np.where(legal, pip * (1 - pip), -1)
    return {"type": "assay", "variant_id": observation["inputs"]["variant_id"][int(utility.argmax())]}


def value_information_assay(observation, config=None):
    """Myopic Bayes value of an assay under the declared working model.

    Compute the expected reduction in weighted posterior Brier risk, including
    the null event and correlated labels. Stop when no legal measurement's
    expected gain exceeds its declared cost. This is not an omniscient optimum
    under a shifted generator or for the full remaining acquisition horizon.
    """
    model = observation["assay_model"]
    prediction, indicators, probabilities = configuration_posterior(
        observation["inputs"], config, successes=observation["assay_successes"],
        attempts=observation["assay_attempts"], sensitivity=model["sensitivity"],
        specificity=model["specificity"])
    pip = np.asarray(prediction["pip"]); null = prediction["null_probability"]
    current_risk = .5*np.mean(pip*(1-pip)) + .5*null*(1-null)
    gains = np.full(len(pip), -np.inf)
    for j, legal in enumerate(observation["legal_assays"]):
        if not legal: continue
        likelihood = (indicators[:,j]*model["sensitivity"]
                      + (1-indicators[:,j])*(1-model["specificity"]))
        probability_yes = float(probabilities @ likelihood)
        future_risk = 0.
        for channel, chance in ((likelihood,probability_yes),(1-likelihood,1-probability_yes)):
            posterior = probabilities*channel/chance
            marginals = posterior @ indicators; null_probability = posterior[0]
            risk = .5*np.mean(marginals*(1-marginals)) + .5*null_probability*(1-null_probability)
            future_risk += chance*risk
        gains[j] = current_risk - future_risk
    best = int(gains.argmax())
    if observation["budget_remaining"] < 1 or gains[best] <= model["cost_per_assay"]:
        return {"type":"submit","prediction":prediction}
    return {"type":"assay","variant_id":observation["inputs"]["variant_id"][best]}
