# Executable challenge protocol

Specification: `conditional-gaussian-v1`. This is a synthetic development
suite with an interface suitable for independent assessment. It is not already
a blind contest. Read [challenge-science.md](challenge-science.md) for the
likelihood, target and limits, and [REVIEW.md](REVIEW.md) for original defects.

## Mechanisms and assessed questions

| Mechanism | Deliberate change | What a good result establishes |
|---|---|---|
| `matched_prior` | K ~ (.5,.3,.2) for 0/1/2; random normal effects match reference prior | Bayesian calibration positive control under a declared experiment |
| `null` | No sparse or background effects | False discoveries and null uncertainty |
| `single_sparse`, `multi_sparse` | One/two sparse effects of declared score magnitude | Conditional localisation and recovery |
| `many_sparse` | Five causes beyond exact reference's model support | Sensitivity to architecture/model-size misspecification |
| `opposing_effects` | Adjacent effects of opposite signs in strong LD | Behaviour when marginal signals cancel |
| `sparse_inf` | Sparse effects plus independent effect draws on every variant | Sparse-target recovery with a nonzero background; store both truths |
| `ld_reference_small` | Independent 30-observation Gaussian reference by default | Reference sampling error, including singular LD when p >= reference_n |
| `ld_reference_shift` | Reference correlation differs from generating population | Likelihood robustness under declared LD mismatch |
| `annotation_informative` | Causal sampling favours activity | Benefit of prespecified useful annotations |
| `annotation_irrelevant` | Replicate one-signal control with supplied nuisance annotations independent of truth | Resistance to fitting-only nuisance annotations; not a new architecture |
| `annotation_reversal` | Causal sampling favours low activity | Negative transfer from an assumed positive annotation relation |
| `missing_causal` | Remove a generating cause from observed candidates | Recovery ceiling, confidence in proxy variants, and outside-model truth |
| `perfect_ld` | Two variants generate the same score variable | Inability to resolve identical likelihoods; class ranking diagnostic |
| `allele_flip` | Change one score sign only | Sensitivity to inconsistent allele coding |
| `allele_recoding` | Change z, coefficients and matching LD rows/columns | Coding invariance positive control |
| `selected_null` | Null scores retained only when max(abs(z)) >= 3.5 | Repeated-sampling performance conditional on discovery selection |
| `selected_signal` | Weak causal mean magnitude 3; retain max(abs(z)) >= 5 | Selected weak-signal robustness; no universal selection-correction claim |
| `meta_missingness` | Alternating variants have one or two independent contributing studies | Variable n with insufficient overlap metadata; adapter should expose unsupported likelihood |
| `meta_known_overlap` | Same experiment with known study-overlap correlation | Covariance-aware inference under the declared information set |
| `multi_ancestry_shared` | Independent cohorts, different Gaussian LD, shared coefficients | Joint shared-effect localisation with complementary LD |
| `multi_ancestry_heterogeneous` | Reverse causal effects in second cohort | Shared-effect model misspecification |

The ancestry-labelled mechanisms vary covariance, not real human genotype
frequencies, demographics, admixture or imputation. Do not claim population
transportability from them. Both discovery-selection arms use rejection
sampling with an explicit cap and assessor-side selection counts. Stress
signals usually scale as `signal_z/sqrt(n)`, so varying n without changing
that policy does not assess power at fixed effects.

Tracks include invariant and placebo controls, not 22 independent architecture
perturbations. Nuisance annotations are also supplied in the baseline tracks;
`annotation_irrelevant` replicates the `single_sparse` distribution under fresh
draws. Its label identifies an assessment question, and does not confer extra
problem diversity. The suite's per-track reports avoid silently doubling that
architecture's weight in an overall score.

## Predictor contract

`generate_challenge()` returns `inputs`, `truth` and `config`. A predictor is
called with one `input` only. Mechanisms are shuffled; task IDs contain no
scenario name. Assessor seeds and truth are excluded from the public tree.

Public input contains:

- `task_id`, ordered unique `variant_id`, and `build="synthetic-v1"`;
- `z`, `n`, score-scale `beta_hat`, `se`, positive-semidefinite `LD`;
- an `annotations` matrix with activity and decoy columns;
- optional independent `cohorts` with z/LD/n and declared cohort IDs;
- optional `n_per_variant` and, in the known-overlap track, `sample_overlap`.

Synthetic IDs are deliberately labelled as synthetic. A future real-data
release needs canonical chromosome, position, REF, ALT and genome-build
identities, provenance, allele reconciliation and variant-wise effective n.
The summary challenge's score convention is not the marginal-OLS convention
used by every existing method; this is part of comparison interpretation.

Output: `task_id`, all ordered or matchable `variant_id`, finite [0,1] `pip`,
and `status="ok"`. Optional native `credible_sets` contain variant-ID sets.
A `credible_configurations` collection must declare `configuration_coverage`;
its elements are complete configurations, including the empty null set when
applicable. Optional `null_probability` must obey the PIP union bounds.
Warnings, nonconvergence and elapsed time remain diagnostic metadata.

Example wrapper:

```r
predictor <- function(input) {
  # Use input only; do not close over development/final truth or tune here.
  result <- my_method(z = input$z, LD = input$LD, n = input$n)
  list(task_id = input$task_id, variant_id = input$variant_id,
       pip = result$pip, status = "ok", credible_sets = list())
}
```

`challenge_region_adapter()` connects summary-based legacy region wrappers.
It rejects multi-cohort and variable-n inputs because discarding these fields
would silently change the declared method. Existing oracle wrappers depend on
simulator truth and do not belong in this participant API. The source-only
runner is an interface convention, not a security sandbox: trusted predictor
code in the same R process can access globals or files. Real untrusted
submissions require isolated execution and access-controlled assessor storage.

## Frozen release and prediction sequence

1. Fix the scientific specification, permitted information, method configuration,
   splits, selection procedure, discovery thresholds and primary contrasts.
2. Generate development data; publish its truth for learning and debugging.
3. Independent assessors generate final inputs and store truth outside the
   participant-accessible tree. Do not publish final seeds or generator state.
4. `write_challenge_release()` writes public input/metadata and separately held
   truth/config with SHA-256 file manifests; overwriting is refused.
5. `read_challenge_release()` verifies hashes and contracts before execution.
6. Seal participant predictions and record code/environment/source hashes,
   hardware, resource budget, retries and completion status before unblinding.
7. Score the frozen predictions; publish failures and availability alongside
   per-mechanism results, prespecified uncertainty and assessor commentary.

SHA-256 manifests detect accidental artifact changes if the manifest is trusted.
They do not authenticate an independently replaceable manifest. Sign or publish
manifest commitments through an independent channel for a live competition.
RDS is convenient for trusted R development; use an interoperable schema and
safe deserialisation policy for untrusted community submissions.

`split="final"` uses disjoint task IDs and an independent seed offset. It does
not create new empirical loci, prevent training-data recognition, or establish
that a final set was kept unseen. Actual exposure history and custody matter.

## Scores and uncertainty

`evaluate_challenge()` produces one task row per attempted method/task, a
mechanism-specific summary and variant-level reliability records. Unknown,
duplicate or inconsistent IDs are rejected. Invalid or missing outputs remain
failures; they are not converted into purported all-zero posterior forecasts.

- AP uses exact distinct-score thresholds and admits ties as a group. It is
  undefined in an observed candidate set with no positive targets.
- Raw Brier and clipped log loss are successful-fit probability scores for
  observed sparse-membership truth. Log clipping is at 1e-15 and labelled.
- Total-target discovery recall and distinct-set recall retain missing causes
  and give failed attempts zero operational recovery credit.
- FDR is the pooled fraction of non-target calls among calls, undefined with
  no calls. Null false-call rate is conditional on completed null tasks.
- Native set hit rate is an existential event. Maximum matching assesses
  distinct sparse-target capture; it does not establish native signal identity.
- Configuration hit rate assesses the full sparse configuration. Nominal
  coverage is shown only when a common declared level exists, with collection
  size alongside hit rate. Mixed declared levels are not pooled as a calibrated
  nominal event; hit rates remain descriptive.
- Perfect-LD class top-hit is a tie-averaged ranking diagnostic, not a posterior
  probability for any member of the class.
- Availability, warning counts, nonconvergence, successful runtime and total
  measured attempted runtime accompany quality. CPU/GPU time, peak memory,
  installation and retries still require an isolated resource runner.

`challenge_compare()` resamples independent locus clusters while preserving
method pairing. It reports A minus B, a percentile interval, complete-pair
availability and attempted pairs for one chosen mechanism/metric. Select
contrasts before looking at results; exploratory multiple comparisons are
not family-wise discoveries. Few clusters yield imprecise intervals; a single
cluster has no reported interval.

The demo intentionally does not optimise a single global score or discard hard
mechanisms. A future competitive track must separately define same-input
regional comparisons, genome-wide prior-learning/strategy comparisons and
real-data robustness/assay predictions. The exact reference's two-effect model
is an explicit limitation on the five-effect track.
