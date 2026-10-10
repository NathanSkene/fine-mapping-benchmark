"""Versioned BEATRICE configuration likelihoods and final-prior probabilities.

The search cache contains log Bayes factors only. Learned priors may change
while configurations are discovered; posterior probabilities are computed under
one final independent-Bernoulli prior, including the null configuration.

This module uses only the standard library so the probability algebra can be
tested independently of neural-network and plotting dependencies.
"""
import math
import numbers

MEMO_FORMAT = "beatrice-log-bf-final-prior-v1"


def _configuration(key, p):
    if not isinstance(key, tuple) or any(not isinstance(j, numbers.Integral) or isinstance(j, bool) for j in key):
        raise ValueError("Configuration keys must be tuples of integer indices")
    if tuple(sorted(set(key))) != key or any(j < 0 or j >= p for j in key):
        raise ValueError("Configuration indices must be unique, sorted and in range")
    return tuple(int(j) for j in key)


def make_context(log_bf, final_prior):
    """Freeze a pure likelihood cache and final prior; null BF is exactly one."""
    prior = [float(x) for x in final_prior]
    if not prior or any(not math.isfinite(x) or x < 0 or x > 1 for x in prior):
        raise ValueError("The final prior must contain finite probabilities in [0,1]")
    likelihoods = {}
    for key, value in log_bf.items():
        key = _configuration(key, len(prior))
        value = float(value)
        if not math.isfinite(value):
            raise ValueError("Cached log Bayes factors must be finite")
        likelihoods[key] = value
    if () in likelihoods and likelihoods[()] != 0:
        raise ValueError("The null configuration must have log Bayes factor zero")
    likelihoods[()] = 0.0
    return {"format": MEMO_FORMAT, "log_bf": likelihoods, "final_prior": prior}


def validate_context(context):
    if (not isinstance(context, dict) or context.get("format") != MEMO_FORMAT or
            not {"log_bf", "final_prior"}.issubset(context)):
        raise ValueError("Legacy or unknown posterior memo format; regenerate this fit")
    return make_context(context["log_bf"], context["final_prior"])


def posterior_weights(context):
    """Normalize configuration probabilities by log-sum-exp, without clipping."""
    context = validate_context(context)
    prior = context["final_prior"]
    mandatory = {j for j, x in enumerate(prior) if x == 1}
    forbidden = {j for j, x in enumerate(prior) if x == 0}
    interior = [(j, x) for j, x in enumerate(prior) if 0 < x < 1]
    base = math.fsum(math.log1p(-x) for _, x in interior)
    odds = {j: math.log(x) - math.log1p(-x) for j, x in interior}
    log_joint = {}
    for key, log_bf in context["log_bf"].items():
        members = set(key)
        if not mandatory.issubset(members) or members.intersection(forbidden):
            log_joint[key] = -math.inf
        else:
            log_joint[key] = log_bf + base + math.fsum(odds.get(j, 0.0) for j in key)
    maximum = max(log_joint.values())
    if not math.isfinite(maximum):
        raise ValueError("No explored configuration has positive probability under the final prior")
    weights = {key: math.exp(value - maximum) for key, value in log_joint.items()}
    total = math.fsum(weights.values())
    return {key: value / total for key, value in weights.items()}


def freeze_memo(memo, final_prior):
    """Replace a training log-BF map with probabilities and return its context."""
    context = make_context(memo, final_prior)
    weights = posterior_weights(context)
    memo.clear()
    memo.update(weights)
    return context


def add_configuration(memo, context, key, log_bf):
    """Extend the search under the same final prior and renormalize every model."""
    checked = validate_context(context)
    key = _configuration(tuple(key), len(checked["final_prior"]))
    log_bf = float(log_bf)
    if not math.isfinite(log_bf):
        raise ValueError("New log Bayes factor must be finite")
    if key == () and log_bf != 0:
        raise ValueError("The null configuration must have log Bayes factor zero")
    checked["log_bf"][key] = log_bf
    weights = posterior_weights(checked)
    context.clear()
    context.update(checked)
    memo.clear()
    memo.update(weights)
    return memo[key]


class _ConfigurationBatch:
    """One complete search round with no intervening probability queries.

    Working likelihoods are private until a successful exit. Existing context
    validation/copying is done once per round; each new key and likelihood is
    checked immediately. Commit normalizes the complete support under the same
    final prior. Aborted rounds do not publish a partially explored posterior.
    """
    def __init__(self, memo, context):
        self._memo = memo
        self._context = context
        self._checked = validate_context(context)
        self._prior = tuple(self._checked["final_prior"])
        self._source_likelihoods = dict(context["log_bf"])
        self._source_memo = dict(memo)
        self._closed = False
        self._entered = False
        self._dirty = False

    def __enter__(self):
        if self._closed or self._entered:
            raise RuntimeError("A configuration batch cannot be reused")
        self._entered = True
        return self

    def _key(self, key):
        if self._closed or not self._entered:
            raise RuntimeError("Configuration additions require an active batch")
        return _configuration(tuple(key), len(self._prior))

    def contains(self, key):
        """Check explored likelihood support; this is not a probability query."""
        return self._key(key) in self._checked["log_bf"]

    def add(self, key, log_bf):
        """Buffer a validated pure log BF; no posterior probability is returned."""
        key = self._key(key)
        log_bf = float(log_bf)
        if not math.isfinite(log_bf):
            raise ValueError("New log Bayes factor must be finite")
        if key == () and log_bf != 0:
            raise ValueError("The null configuration must have log Bayes factor zero")
        self._checked["log_bf"][key] = log_bf
        self._dirty = True

    def __exit__(self, exception_type, exception, traceback):
        self._closed = True
        if exception_type is not None:
            return False
        # Reject concurrent/reentrant changes instead of silently overwriting a
        # different prior, likelihood cache or publicly observable posterior.
        if (self._context.get("format") != MEMO_FORMAT or
                tuple(self._context.get("final_prior", ())) != self._prior or
                self._context.get("log_bf") != self._source_likelihoods or
                self._memo != self._source_memo):
            raise ValueError("Frozen posterior context changed during a configuration batch")
        if self._dirty:
            weights = posterior_weights(self._checked)
            self._context.clear()
            self._context.update(self._checked)
            self._memo.clear()
            self._memo.update(weights)
        return False


def configuration_batch(memo, context):
    """Explicitly batch a complete candidate round, then publish probabilities.

    Only use when no selection/credible-set decision reads ``memo`` inside the
    round. The existing ``add_configuration`` API retains eager normalization.
    """
    return _ConfigurationBatch(memo, context)


def load_context(result):
    """Reject old mixed-prior archives: they cannot be scientifically converted."""
    if not isinstance(result, dict) or result.get("memo_format") != MEMO_FORMAT:
        raise ValueError("Legacy posterior archive lacks pure likelihoods; regenerate this fit")
    return validate_context(result.get("posterior_context"))
