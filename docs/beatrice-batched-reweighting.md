# BEATRICE posterior expansion and conditional sets

This fork fixes two scientific errors in credible-set generation. Conditional
membership probabilities now accumulate in floating point: the original integer
arrays truncated fractional model probabilities to zero, suppressing additional
signals. Credible-set stopping now uses the actual conditional union-event
probability, rather than a sum of membership marginals that can count the same
model several times. These changes intentionally alter inference from the
original implementation. Training, likelihood formulas, random-number generation
and the persisted annotation-prior contract are unchanged.

Separately, candidate expansion now normalizes the posterior once after each
complete add/delete/change round. Batching retains candidate traversal, pure log
Bayes factors and configuration insertion order. Its equivalence reference is
the **corrected eager implementation**, with both scientific fixes applied. A
round publishes probabilities before every subsequent conditional decision or
final output; the batching change itself leaves successful inference unchanged.

`configuration_batch(memo, context)` validates and copies the existing likelihood
context once on entry. Every added configuration and log Bayes factor is checked
immediately. Successful exit normalizes the complete explored support using the
existing `posterior_weights` function; failed rounds publish no pending additions.
Changes to the frozen prior, existing likelihood cache or observable memo during a
round are rejected. The context is for candidate enumeration without probability
queries; reading the memo inside it would read the previous committed posterior.
The public single-add `add_configuration` API still normalizes eagerly and returns
the added model's current probability.

This is valid because all candidates in these three helpers are generated without
using intermediate probabilities. At the end of a round, both implementations
evaluate the same expression for every explored configuration C:

    P(C | data, explored support) ∝ exp(log_BF[C])
                                  × ∏(j in C) prior[j]
                                  × ∏(j not in C) (1 − prior[j]).

The empty configuration remains in that support. This computational equivalence
does not establish calibration outside the explored support or address omissions
by the search procedure.

For a candidate set A and the other selected variants G, the default
`credible_set_rule="conditional_union"` stops when

    Σ(C: G ⊆ C and C ∩ A ≠ ∅) P(C | data, explored support)
    -----------------------------------------------------  ≥ target.
          Σ(C: G ⊆ C) P(C | data, explored support)

This conditions on G being causal and normalizes within the final explored
configuration support. It is not an unconditional 95% guarantee or a claim about
unexplored models. A zero-mass conditioning event is undefined and raises an
error. The no-duplicate path masks previously used members, without renormalizing
the remaining member marginals into a choice distribution. The outer stepwise
loop controls the selected signal count; the retained inner `n_caus` argument
does not impose a cardinality filter.

For example, models `{0,1}`, `{0}`, `{1}`, `{2}` with probabilities
0.4, 0.3, 0.2, 0.1 give member probabilities 0.7, 0.6, 0.1. Set `{0,1}` has
membership sum 1.3 but union-event coverage 0.9. At target 0.95 the corrected rule
adds variant 2 where selection and purity constraints permit. The explicit
`credible_set_rule="legacy_membership_sum"` option preserves sum stopping as an
ablation with floating-point accumulators; it does not restore the integer bug
or represent a replay of the whole original algorithm. Its sum is a diagnostic,
never reported as coverage.

`credible_set_diagnostics.json` records actual conditional union coverage,
conditioning and set indices, requested target and `target_met` for each exported
set. If constraints prevent reaching the target, the set is retained with its
actual coverage and `target_met=false`. All search expansions finish before these
events are evaluated against the same frozen support persisted in `res` and used
by exported PIPs. Deduplication preserves the first set's matching conditioning
record. The versioned `scripts/credible_set_contract.json` marker declares the
default rule. Native variant names must exactly match the ordered Z-file IDs;
JSON uses zero-based indices and rejects nonfinite numbers. The reported null
probability is also normalized within explored support under the final prior.

`tests/python/test_posterior_batch.py` compares the production batch path with
the retained eager API. It checks exact likelihood support and order, weights,
every traced conditional selection and probability vector, PIPs, credible sets,
null posterior and frozen prior. Cases include seeded 12- and 300-variant real
Torch expansions with genuine multi-signal selection, all three round helpers,
500 candidate additions, duplicate
updates, zero evidence, extreme finite log Bayes factors and priors at zero and
one. NaN/infinite log Bayes factors, malformed keys/context, changed prior/cache
and aborted or reused batches are rejected. Independent probability-table tests
cover the conditional denominator, masked marginals, overlapping models, an
unattainable target and deduplication. A production `main` roundtrip checks the
sidecar against the final persisted posterior, output sets and ordered IDs.

Local measurements on 10 October 2026 used Python 3.14.7 and Torch 2.14.1 with one
Torch thread. Timings include test instrumentation and are illustrative; the
operation counts and exact-output assertions are the regression contracts.

| Scope | Eager normalization calls | Batch calls | Eager time | Batch time |
| --- | ---: | ---: | ---: | ---: |
| 500-variant memo round, 500 additions, 501 models | 500 | 1 | 0.642 s | 0.00175 s |
| Corrected 300-variant Torch expansion, 898 models | 895 | 3 | 3.661 s | 0.115 s |

The 500-addition test also counts complete context validations: 1,000 eagerly
versus two in a batch. These are expansion measurements, excluding model training,
data loading and the full server campaign. No full-size campaign runtime or peak
memory claim follows from them. New source requires fresh source fingerprints and
scheduled preflights.

Reproduce the 14 focused tests in an environment with the BEATRICE dependencies:

```sh
FMB_NEURAL_SMOKE=1 python -m unittest discover \
    -s tests/python -p test_posterior_batch.py -v
```

Set `FMB_BATCH_METRICS_PATH` to a JSON output path to retain the operation counts
and timings. The broader dependency-free suite remains runnable without Torch;
the seven Torch-dependent equivalence, event and sidecar tests require the
explicit opt-in above.
