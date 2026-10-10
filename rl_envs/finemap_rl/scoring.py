"""Deterministic verification: no model judge, ranking proxy or hidden score bins."""

import math
import numpy as np


def validate_prediction(prediction, inputs):
    if not isinstance(prediction, dict) or set(prediction) != {"task_id", "variant_id", "pip", "null_probability"}:
        raise ValueError("Prediction requires task_id, variant_id, pip and null_probability")
    if prediction["task_id"] != inputs["task_id"]:
        raise ValueError("Prediction task identity does not align")
    ids = prediction["variant_id"]
    if not isinstance(ids, list) or len(ids) != len(inputs["variant_id"]) or len(set(ids)) != len(ids):
        raise ValueError("Duplicate, missing or invalid variant identities")
    if set(ids) != set(inputs["variant_id"]):
        raise ValueError("Prediction identities differ from supplied candidates")
    if any(isinstance(x, bool) for x in prediction["pip"]):
        raise ValueError("PIPs must be real probabilities, not booleans")
    pip = np.asarray(prediction["pip"], dtype=float)
    null = prediction["null_probability"]
    if isinstance(null, bool) or not isinstance(null, (int, float)) or not math.isfinite(null) or not 0 <= null <= 1:
        raise ValueError("Invalid null probability")
    if pip.shape != (len(ids),) or not np.isfinite(pip).all() or ((pip < 0) | (pip > 1)).any():
        raise ValueError("Invalid PIPs")
    # Necessary event coherence, not a proof of a unique joint distribution.
    if pip.max() > 1 - null + 1e-7:
        raise ValueError("PIPs and null probability violate event bounds")
    lookup = dict(zip(ids, pip.tolist()))
    return np.asarray([lookup[i] for i in inputs["variant_id"]]), float(null)


def brier_loss(prediction, task):
    """Positive weighted proper score for observed sparse membership and whole-locus null.

    Equal weight makes the null event consequential despite sparse candidate labels.
    Neither term asserts a complete joint posterior or credible-set probability.
    """
    pip, null = validate_prediction(prediction, task.inputs)
    return float(.5 * np.mean((pip - task.labels) ** 2) + .5 * (null - task.null) ** 2)


def panel_loss(predictor, tasks):
    """Score every attempted task; exceptions get maximum loss one."""
    losses, failures = [], 0
    for task in tasks:
        try:
            losses.append(brier_loss(predictor(task.inputs), task))
        except (Exception,):
            losses.append(1.0); failures += 1
    if not losses:
        raise ValueError("Cannot score an empty panel")
    return {"loss": float(np.mean(losses)), "failures": failures, "attempted": len(tasks)}
