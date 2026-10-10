# Evaluation methodology

`evaluate_methods()` scores the simulation design against method outputs. Version
2 adds exact ranking metrics, explicit failure denominators, null-locus metrics,
proper probability scores, validation and duplicate-aware signal recovery. The
legacy threshold-grid `auprc` remains available for comparison with previous
runs. Use `ap_macro` with `ap_macro_attempted`, success rate and calibration for
new locus-ranking comparisons; use pooled threshold counts for decisions made
at one PIP threshold across loci. These answer different questions.

## Inputs and the population being scored

The two required inputs are outputs of `run_simulation()` and `run_methods()`.
A fit identifies one `scenario_id` and `region_id`. Every scenario × region in
the supplied simulation is an intended attempt. Missing outputs are synthesized
as failures; duplicate identities or identities outside the simulation are
rejected. Evaluate an intentional subset using a correspondingly subsetted,
consistently indexed simulation, rather than silently dropping its difficult
fits.

The simulated genotype matrix/LD matrix defines the variant universe. Submitted
PIP length cannot redefine it. Scenario metadata defines the stratification and
replicate labels. Each successful output must have finite numeric PIPs in
`[0,1]`, one per simulated variant, and credible sets must be nonempty vectors
of in-range integer variant indices. Malformed output becomes a failure, counted
in `n_invalid` and listed in `validation_errors`. Scores are never rescued by
clipping probabilities or dropping nonfinite variants. A method's explicit
error is counted in `n_failed`; missing results also contribute `n_missing`.

Truth is the simulation's `truth$causal_indices`. In a sparse-plus-infinitesimal
architecture these are the designated sparse effects, not every nonzero
infinitesimal effect. The benchmark therefore evaluates recovery of the sparse
component. A different scientific target requires different truth labels.

Optional arguments retain their original defaults:

| Argument | Default | Meaning |
|---|---|---|
| `pip_thresholds` | `seq(0, 1, by = 0.005)` | Finite probability thresholds for operating-point curves; does not affect exact AP |
| `n_pip_cal_bins` | `10L` | Positive integer number of equal-width calibration bins |
| `save` | `FALSE` | Write the complete RDS and flat global CSV |
| `output_dir` | `"results"` | Destination when saving |
| `verbose` | `TRUE` | Progress messages |

## Failure denominators and null regions

Historical fields `fdr_power_curve`, `cs_power`, `ap`, `ap_macro` and proper
probability scores condition on valid outputs. That conditioning must be visible
when comparing methods: a method may preferentially fail on difficult loci.
Version 2 additionally reports:

- `success_rate`, `n_fits`, `n_failed`, `n_invalid`, `n_missing`;
- `fdr_power_curve_attempted`: the same successful predictions, but with every
  simulated causal in its recall denominator; failures make no calls and leave
  their causal variants as false negatives;
- `cs_power_attempted`: every failed locus has an empty collection of sets;
- `ap_macro_attempted`: failed signal-containing loci contribute zero AP, and
  all intended signal-containing loci contribute to the denominator.

A failure has no probability forecast. Consequently it cannot receive an honest
Brier or logarithmic score without a declared fallback forecasting rule. These
scores remain conditional and are accompanied by success rate. Treating failure
as no calls for recovery is a declared operational convention, not a posterior
probability of zero.

For a null locus, every call is false: FP and FDR remain observable, including
at threshold zero. Recall and AP are undefined when there are no true positives,
and are reported as `NA`. They are not filled with reassuring zeros. The
successful-null denominators `n_null_fits`, `null_region_false_call_rate` (at
the predeclared PIP threshold 0.5) and `null_region_cs_rate` expose methods that
invent signals. The complete threshold curve allows analysis at other
predeclared cutoffs. Report failure rate on null loci as well.

If every fit fails, the successful-fit curve and calibration table are `NULL`;
conditional ranking/probability scores are `NA`. Attempted recall and attempted
CS power are zero if the simulation contains causal variants. No sets were
reported, so `n_cs_reported = 0`.

## Operating points: pooled precision, false discovery proportion and recall

Within a stratum, successful PIPs and their truth labels are concatenated. At
threshold `t` every PIP `>= t` is selected:

```text
TP(t) = selected causal variants
FP(t) = selected noncausal variants
FN(t) = target causal variants - TP(t)
precision(t) = TP(t) / (TP(t) + FP(t))
fdr(t) = FP(t) / (TP(t) + FP(t))
recall(t) = power(t) = TP(t) / target causal variants
```

The stored `fdr` is the realized pooled false discovery proportion (FDP). Its
name is retained for compatibility; formal FDR is the expectation of FDP over
repeated data. Methods are compared on these data at a predeclared threshold,
not by optimizing the threshold using hidden test truth.

With no selections, precision is 1 and FDP is 0 by convention; this is
abstention and must be assessed alongside selection volume and recall. At
`t=1` variants with PIP exactly 1 are selected. At `t=0` every successful-fit
variant is selected. Sorting and cumulative counts handle tied scores together;
thresholds may be unsorted or repeated without changing their corresponding
operating points.

## Exact average precision and legacy AUPRC

Exact AP groups variants at each distinct observed PIP value. If `P_b` and
`R_b` are precision and recall after admitting the entire tied block `b`, then

```text
AP = sum_b (R_b - R_(b-1)) * P_b, with R_0 = 0.
```

There is no arbitrary PIP grid, interpolation or genomic-position tie break.
A causal at rank 61 has AP `1/61` in a single-causal locus even if its PIP is
0.0032. A constant predictor gives AP equal to the causal prevalence; a perfect
ranking gives 1, including when its causal PIP is exactly 1. This is the
non-interpolated estimator in the
[scikit-learn AP documentation](https://scikit-learn.org/1.8/modules/generated/sklearn.metrics.average_precision_score.html).

`ap` is exact **pooled** AP. It evaluates a global ordering across loci and
therefore depends on between-locus score scales and the variant/causal mixture.
`ap_macro` is the unweighted mean of exact AP across successful
signal-containing loci; it asks about a typical locus's local ordering and is
invariant to strictly monotone score transformations applied separately within
each locus. `ap_n_fits` gives its non-null successful denominator.
`ap_macro_attempted` adds failed signal-containing loci with score zero. Null
loci do not have a ranking-recovery task and are scored separately.

Pooled and macro AP can legitimately rank methods differently. Neither is
universally correct for every intended deployment. Their names state the
estimand. Large loci and different causal counts also change the task; report
strata and baseline prevalence rather than treating one global scalar as a
scientific conclusion.

`auprc` continues to use trapezoidal integration of the supplied threshold grid
for historical compatibility. It can miss observed score transitions and the
recall-zero endpoint, and its linear interpolation can bias the summary. It is
not exact AP. Previous numerical claims of a fixed percentage bias or a
universally downward bias are unsupported and should not be used.

## Probability calibration and proper scores

Calibration compares observed causal frequency to predicted probability in
`n_pip_cal_bins` equal-width bins. The final bin includes PIP 1. Empty bins
have undefined means. The returned table contains `n`, `n_causal`, `mean_pip`,
`frac_causal` and `frac_causal_se`, plus `sum_pip`, `sum_pip_sq` and
`sum_causal_pip` so the raw quadratic score can be reconstructed exactly.

Three scalar diagnostics are computed directly from the original predictions:

```text
brier_score = mean((PIP - truth)^2)
log_loss = -mean(truth*log(PIP) + (1-truth)*log(1-PIP))
ece = sum_b n_b * abs(mean_PIP_b - causal_fraction_b) / N
```

The logarithmic implementation avoids `0 * log(0)` and preserves infinite loss
for certainty assigned to an event contradicted by truth. It does not clip
probabilities to disguise impossible-event forecasts. Brier and logarithmic
loss reward accurate probability forecasts; ECE is a bin-dependent descriptive
calibration diagnostic, not a proper scoring rule. Calibration alone can reward
an uninformative constant forecast.

`brier_skill_score` compares raw Brier to the empirical constant-prevalence
reference `prevalence*(1-prevalence)`. It is undefined for an all-null or
all-causal population. The reference uses test truth as a descriptive baseline,
not a probability estimate available to the method. Raw Brier remains the
proper scoring target; skill scores and mixtures across changing prevalences
need careful interpretation.

The optional `.murphy_decomposition()` helper reports reliability, resolution
and uncertainty following
[Murphy's probability-score decomposition](https://journals.ametsoc.org/view/journals/apme/12/4/1520-0450_1973_012_0595_anvpot_2_0_co_2.xml).
Its historical `bs = rel-res+unc` and `bss` score **bin-mean predictions**.
They need not equal raw Brier because binning loses within-bin prediction
variation and association with truth. With squared-PIP and causal-PIP moments,
`bs_raw` and `bss_raw` recover the original score. Without both moments, raw
fields are `NA`; no exact raw score is claimed from insufficient counts.

For example PIPs `(0.1,0.4)` with truth `(1,0)` in one bin give raw Brier
0.485 but bin-mean Brier 0.3125. The discrepancy is scientific information
lost by aggregation, not roundoff.

## Native credible sets and signal-recovery diagnostics

Native sets are kept in their submitted collection. Distinct variant indices
within a set define its size. The original summaries are descriptive:

| Field | Definition |
|---|---|
| `cs_coverage` | Fraction of reported sets containing at least one true causal |
| `cs_power` | Fraction of successful-fit causals appearing in any reported set |
| `cs_power_attempted` | The same recovery over every intended causal |
| `cs_size_median`, `cs_size_mean` | Number of distinct variants in each reported set |
| `n_cs_reported` | Number of reported sets, including repeated reports |
| `cs_duplicate_sets` | Number of exact repeated set memberships within loci |

A method reporting two copies of `{causal_1}` has hit fraction 1 but has not
recovered two distinct signals. `cs_matched_precision` and `cs_matched_power`
use **maximum bipartite matching** between reported sets and causal variants:
each set can receive one causal match and each causal can be credited once.
Their numerators are the number of matches and denominators are reported sets
and successful-fit causals, respectively. Overlapping sets require an
augmenting-path algorithm; greedily claiming the first causal can undercount
recovery. Repeated sets remain in the denominator and cannot earn additional
credit for the same causal.

The matching convention is useful for methods claiming one set per signal.
It is not a valid penalty for a method claiming one *global* set for all
signals. Similarly the native hit fraction is not automatically nominal
95% coverage: a SuSiE signal-specific set, a posterior model/configuration set,
a global set and a greedy prefix of marginal PIPs represent different events.
Compare nominal coverage only after establishing that the methods claim the
same event and conditioning regime. Repeated sets, selection of reportable
signals and uncertain signal count must also be accounted for. Native
coverage/power fields alone cannot establish calibrated signal inference.

The analysis helper `.rank_set_metrics()` retains historical `set_*` names
for compatibility. Its prefixes reach an absolute **sum of marginal PIPs**,
which is expected causal count, not the probability of recovering a causal
signal. These are diagnostic candidate lists, not credible sets with nominal
coverage. All members of a boundary tie are included. Fixed-budget
`prec_at_1`, `prec_at_S` and `prec_at_2S` instead report expected precision under
a uniform tie break at that budget. An unreachable mass target is reported
with `set_reached_* = FALSE` and whole-locus size.

When no sets are reported, hit fraction and set sizes are undefined; causal
recovery is zero if causal variants exist. Runtime summaries use nonnegative,
finite `runtime_seconds` from valid fits, and omit missing timing values.

## Stratification and aggregation

Output includes `global`, `by_S`, `by_phi`, `by_p_causal`, `by_region_size`,
`by_causal_maf` and `by_true_annotation_type` when those variables exist. MAF
bins use minimum causal MAF: rare `<=0.01`, low `<=0.05`, common `>0.05`.
Null loci have no causal MAF and are excluded from that stratification.

Nominal region length is a design axis; realized post-filter length still
defines each fit's actual variant universe. A global summary pools design
cells and is a descriptive mixture. It is not a replacement for difficulty
strata or a predeclared challenge weighting scheme. One-variable strata also
mix the other variables; use a complete scenario-cell analysis to investigate
interactions or fairness across region sizes.

## Uncertainty targets the reported aggregate

Every scalar `_se` field and the pointwise curve/calibration SEs use a
delete-one-`iter` **cluster jackknife**. Each deletion recomputes the actual
pooled or macro statistic on all remaining replicate clusters:

```text
SE = sqrt((K-1)/K * sum_r (theta_(-r) - mean(theta_(-r)))^2).
```

This replaces `sd(per_replicate_metric)/sqrt(K)`, which generally estimates
uncertainty of a mean of replicate statistics rather than of a pooled nonlinear
ratio, AP or median. Variants in LD and multiple fits inside the same replicate
are not treated as independent observations. `n_replicate_clusters`,
`se_method` and `se_conditioning` declare the uncertainty unit.

If fewer than two clusters exist, a cluster label is missing, or a required
leave-one-out statistic is undefined, its SE is `NA`. Undefined deletions are
not removed to manufacture a smaller denominator. With very few clusters,
non-smooth rankings/medians or sparse calibration bins, jackknife SEs can be
unstable; they are descriptive uncertainty estimates, not automatic confidence
intervals or evidence of a winner.

Replicates share the simulated loci and genotype/LD backgrounds in the original
simulation design. These SEs are conditional on those backgrounds and the fixed
scenario grid. They do not quantify generalization to new genomic loci,
ancestries or study designs. An independent locus holdout with paired
locus-level resampling is needed for that claim. Global pooling of fixed S/phi
values does not by itself introduce random parameter uncertainty or guarantee
that the global SE exceeds every stratified SE.

## Saving and existing plots

`save=TRUE` writes `evaluation.rds` with all fields and
`evaluation_summary.csv` with core global scores, failures, proper scores,
matching and null metrics. Metadata includes `evaluation_version = "2.0"`,
`simulation_params` and `pip_thresholds_used`. Historical columns remain.

```r
ev <- evaluate_methods(simulation, results, save = TRUE,
                       output_dir = "results/my_run")
ev$susie$global[c("ap_macro", "ap_macro_attempted", "success_rate")]
ev$susie$by_region_size[["500"]]$pip_calibration
```

`plot_results()` retains the historical PR/calibration and AUPRC summary plots.
These legacy plots do not automatically display all version-2 validity metrics.
For new comparisons, inspect exact AP, intended-attempt recovery, failures,
null-region calls and probability scores together. Selecting nothing is not
proof of useful error control, and calibrated constant probabilities are not
proof of useful fine-mapping.
