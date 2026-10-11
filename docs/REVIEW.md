# Scientific and benchmark validity review

Reviewed on **10 October 2026**. The original source is upstream commit
[`f053fdf318a3abe2c44eeb0ac3f3c99a03979b81`](https://github.com/lucassiu14/fine-mapping-benchmark/tree/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81).
Source line links below refer to that immutable revision, not to the changing
fork. Repair descriptions refer to code inspected in this fork during the
review; they do not certify that every external method has been installed and
rerun. Regression-test results belong in the validation record for the released
revision.

The original framework is useful infrastructure, but its existing experiments
cannot establish a general ranking of fine-mapping methods. Several scoring
and truth-generation defects change the meaning of the results. The new
challenge is a second, explicit simulation experiment with better separation of
prediction and assessment. It is not an experimental causal-variant gold
standard. See the [scientific contract](challenge-science.md) and
[challenge protocol](challenge-protocol.md).

## Review scope and severity

- **P1:** can invalidate a performance or calibration conclusion without an
  obvious runtime failure.
- **P2:** materially limits interpretation, comparative fairness or
  generalisation; may be acceptable if explicitly scoped.
- **P3:** explanation or reproducibility debt that makes an otherwise usable
  result harder to assess.

This review covers simulation, evaluation, comparator meaning, experiment
design and explanatory claims. It does not amount to a line-by-line audit of
all bundled Python training code or every external binary.

## Findings and repair state

| Finding | Severity | State in this fork |
|---|---|---|
| AP evaluated on a coarse PIP grid | P1 | Exact tied-rank AP used in the repaired evaluator |
| Successful fits form method-specific denominators | P1 | Missing/invalid fits retained; conditional and attempted metrics separated |
| Sparse/background covariance omitted from PVE calibration | P1 | Total genetic variance used; variance ledger and background coefficients retained |
| Sum of PIPs presented as credible-set probability | P1 | Historical helper retained as a candidate-list diagnostic; new reference reports configuration collections |
| Binned Murphy score called raw Brier score | P1 | Raw moments and direct raw scoring added |
| Oracle prior described as a universal ceiling | P1 | Comparator explicitly labelled information-advantaged and not a ceiling |
| LD-score intercept interpreted as genetic baseline; restricted solver | P1 | Baseline LD score, null intercept and KKT-checked solver separated |
| Duplicate sets repeatedly credited for the same signal | P2 | Distinct signal matching and duplicate diagnostics added |
| Genome-wide probabilities lose requested prevalence after clipping | P2 | Global probability calibration and actual probabilities retained |
| Mixed-population haplotypes imply unjustified population claims | P2 | Still an open empirical-population gap |
| Focal-method tuning confused with competitive method evaluation | P2 | New protocol separates comparisons; original competitive run not rerun |
| Phenotype-replicate SE interpreted as locus generalisation | P2 | Conditional uncertainty labelled; new challenge supports paired locus clusters |
| Singular LD and weak signal removed from the frontier | P2 | Restored as challenge questions, with explicit interpretation limits |
| Equal-effect sign policy differs between generators | P2 | Random signs aligned; opposing-effect challenge added |
| Null loci receive zero false positives by construction of the evaluator | P1 | Null discoveries counted; conditional false-call rates exposed |
| Marginal regression omits centering/intercept | P1 | Centered regression and independently checked OLS standard errors |
| Learned-prior Python cache mixes different training steps | P1 | Pure likelihood cache; one final prior applied to all explored configurations |
| Conditional BEATRICE probabilities accumulated as integers | P1 | Floating-point accumulation and analytical conditional-probability tests |
| BEATRICE set stopping double-counts configurations | P1 | Conditional union-event stopping, attained-probability metadata and an explicit legacy diagnostic |
| Tuning objective, sampling design and study reuse change the experiment | P1 | Balanced panel, attempted macro AP, failure checks and fingerprinted seeded studies |
| Historical README and iteration claims outlive their scope | P3 | Revised explanations distinguish historical runs from repaired code |

### 1. P1 — coarse-grid AP can reverse low-signal comparisons

**Original source:** [`R/evaluate.R:344–381`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/evaluate.R#L344-L381).

The evaluator pools PIPs and computes `ap` from a fixed threshold curve, even
though an exact helper already exists. Distinct predictions in the same grid
interval enter the selected set simultaneously. A causal variant at rank 61
with PIP 0.0032 can receive precision 1/1000 at the zero threshold rather than
1/61 at its actual rank. The error grows when scores compress towards zero;
it is not a harmless constant offset. Trapezoidal PR area also differs from
step-integral average precision.

**Consequence:** weak-signal methods and strata can be misranked for a reason
unrelated to fine-mapping ability.

**Repair:** exact distinct-score AP, explicit tie handling and per-fit macro AP.
The historical trapezoidal quantity can remain for compatibility, but it must
not inherit the meaning of exact AP. Micro AP measures a global ranking across
loci; macro AP measures within-locus localisation. Neither substitutes for the
other.

### 2. P1 — missing and failed hard loci disappear from the competition

**Original source:** [`R/evaluate.R:346–405`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/evaluate.R#L346-L405).

Only fits without an `error` field enter scores and runtime summaries. A
method that fails on difficult loci can look better than a method that returns
usable, uncertain predictions on every locus. Missing result objects are not
necessarily counted as attempted tasks. The variant universe is inferred from
the submitted PIP length, which also makes malformed short outputs dangerous.

**Repair:** complete the expected task grid, validate predictions against the
simulation universe and record missing/invalid attempts explicitly. Report
conditional successful-fit probability scores alongside availability and
recall over all attempts. A failed fit is not an all-zero posterior. Zero
discovery credit after a failure describes an operational policy, not a
probabilistic prediction.

**Remaining limit:** successful-fit runtime does not describe total deployment
cost. Failed work, retries, installation, CPU/GPU use and memory need separate
accounting. An algorithmically unsupported track must also be distinguished
from an attempted fit that failed.

### 3. P1 — sparse and infinitesimal variances do not generally add

**Original sources:** [`R/simulate_phenotypes.R:742–767`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/simulate_phenotypes.R#L742-L767)
and [`R/simulate_gwfm_data.R:979–1027`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/simulate_gwfm_data.R#L979-L1027).

Independent effect draws do not make the realised genetic signals orthogonal
under LD. If `g = g_sparse + g_inf`, then

```text
Var(g) = Var(g_sparse) + Var(g_inf) + 2 Cov(g_sparse, g_inf).
```

The original per-locus path sets component variances to `q * phi` and
`(1-q) * phi`, then residual variance to `1-phi`, omitting covariance. Its
reported realised PVE numerator omits covariance too. The sparse-only path
labels a ratio using the noise distribution's variance as realised PVE, though
that ratio is the target conditional expectation.

**Repair:** calibrate residual variance from `Var(g)`; retain sparse and
background coefficients and an explicit variance ledger. Expected PVE and the
finite-sample ratio `Var(g)/Var(y)` are different quantities. The latter is not
sample regression R-squared and can exceed one when realised signal/noise
covariance is negative.

**Remaining interpretation:** `p_causal` is a target ratio of component
variances before their covariance contribution. It is not a uniquely defined
fraction of total genetic variance in every draw. Background variants have
nonzero effects; calling them biologically noncausal confuses sparse-component
membership with all-effect causality.

### 4. P1 — cumulative marginal PIPs are not a multi-causal credible set

**Original source:** [`R/evaluate_extras.R:79–113`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/evaluate_extras.R#L79-L113).

The helper accumulates marginal PIPs until their sum reaches a requested
`alpha`, then evaluates whether the candidate list contains any cause. However,
the sum is an expected count. It is not the probability that the list contains
at least one cause, or all causes. Two independent variants with PIPs 0.6 and
0.6 have sum 1.2 but probability 0.84 of containing any cause. A list selected
at mass 0.95 can therefore lack 95% event coverage even if its PIPs are exact.
Normalising PIPs to sum to one does not fix the problem.

**Repair:** retain this historical calculation only as a cumulative-mass
candidate-list diagnostic. Use native signal-specific credible sets, or a
credible collection of complete configurations whose event is explicitly
defined. Report list size and hit rate without attaching a nominal guarantee
the method did not calculate.

### 5. P1 — the binned Murphy calculation is not the raw Brier score

**Original source:** [`R/evaluate_extras.R:117–153`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/evaluate_extras.R#L117-L153).

`REL - RES + UNC`, using one mean PIP per broad bin, scores the bin-collapsed
forecasts. Raw binary Brier score requires

```text
BS = [sum(p_j^2) - 2 sum(gamma_j p_j) + sum(gamma_j)] / J.
```

The original bin counts and PIP sums cannot reconstruct the squared-PIP or
label/PIP cross-product terms. Rankings can change when predictions vary
within bins. ECE and mass ratios are also insufficient: a diffuse predictor
can have little calibration error while providing little useful resolution.

**Repair:** direct raw Brier and log loss; sufficient raw moments where only
aggregates are stored; binned diagnostics labelled separately. Probability
scores reward a different property from ranking or candidate-set size. Proper
scoring rules motivate assessing the stated probabilities themselves.
[Gneiting and Raftery, 2007](https://doi.org/10.1198/016214506000001437).

### 6. P1 — the annotation oracle is neither an exact prior nor a ceiling

**Original sources:** [`R/wrapper_polyfun_oracle.R:15–29`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/wrapper_polyfun_oracle.R#L15-L29)
and [`R/simulate_phenotypes.R:600–609`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/simulate_phenotypes.R#L600-L609).

For more than one cause, sequential weighted sampling without replacement
does not give marginal inclusion probabilities equal to the normalised
sampling weights. Passing those weights to SuSiE also does not reconstruct
the simulator's joint fixed-count prior. In the genome-wide path, Bernoulli
probability saturation creates a further distinction.

Even a correct posterior minimises its corresponding expected Bayes loss;
it need not win every finite replicate or every ranking metric. This wrapper
also uses an approximate inference backbone and can face likelihood
misspecification. Better performance by another method is not proof of
cheating.

**Repair:** retain the method as an enrichment-informed, information-advantaged
diagnostic. Label its truth access and exclude it from ordinary participant
rankings. The new small-locus exact reference is exact only for its declared
working model, not a universal biological optimum.

### 7. P1 — LD-aware prior estimation still needs a genetic baseline

**Original sources:** [`R/wrapper_polyfun_est.R:9–19`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/wrapper_polyfun_est.R#L9-L19)
and [`R/wrapper_polyfun_ldsc.R:307–384`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/wrapper_polyfun_ldsc.R#L307-L384).

In the independent, zero-mean random-effect score model,

```text
E[z_j^2] = 1 + n sum_k R_jk^2 E[beta_k^2].
```

Regressing `z_j^2` on variant j's own annotations attributes tagged effects to
the tag's annotation. Simulation does not remove LD. The LD-score wrapper
improves this but originally converts its regression intercept into a genetic
baseline prior. The null sampling intercept of one and the all-variant
genetic LD score are different terms. The original sign-dropping solver also
cannot establish general NNLS optimality because removed coordinates never
re-enter and KKT conditions are not checked.

**Repair:** separate fixed null intercept, sample-size scaling and all-variant
LD score; use a solver with a convergence criterion; hold out chromosomes when
metadata supports it and otherwise label region holding out explicitly.

**Remaining limits:** this restricted comparator is not upstream PolyFun. A
nonnegative-coefficient constraint is a modelling restriction; negative
conditional annotation coefficients can represent depletion while total
per-variant variance remains positive. Canonical PolyFun's cross-fitting and
prior-estimation procedure should not be replaced by a local approximation
under its name. [Weissbrod et al., 2020](https://www.nature.com/articles/s41588-020-00735-5).

### 8. P2 — credible-set hits can count one signal repeatedly

**Original source:** [`R/evaluate.R:646–673`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/evaluate.R#L646-L673).

Every reported set is credited if it contains any cause. Repeating one correct
set, or reporting several overlapping sets around one cause, can produce
excellent set hit coverage without discovering independent signals. A large
set can hit multiple causes without resolving any of them.

**Repair:** retain the historical hit fraction for compatibility, but add
duplicate diagnostics and maximum bipartite matching between reported sets
and distinct true variants. Include signal recall, set size and reporting
frequency. Matching is a resolution diagnostic, not a new nominal posterior
coverage theorem.

### 9. P2 — probability clipping changes genome-wide causal prevalence

**Original source:** [`R/simulate_gwfm_data.R:881–900`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/simulate_gwfm_data.R#L881-L900).

Normalising weights to mean `pi` within each region, then clipping values
above one, reduces the expected causal count under concentrated annotations.
Different regions are separately renormalised, so between-region annotation
composition does not influence prevalence as a common global prior would.

**Repair:** globally calibrate the saturation multiplier so mean selection
probability matches the requested `pi`; record actual probabilities. Saturated
weights no longer preserve literal enrichment ratios. Realised counts still
vary under Bernoulli sampling and must not be called exactly fixed.

### 10. P2 — the reference population is not defined by a window download

**Original sources:** [`inst/scripts/prepare_vcfs.R:7–13`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/inst/scripts/prepare_vcfs.R#L7-L13)
and [`R/simulate_genotypes.R:332–340`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/simulate_genotypes.R#L332-L340).

The downloaded VCF contains all 2,504 Phase 3 individuals; the reader has no
population-subsetting step. Recombining a mixed source does not establish a
European reference, a realistic admixed cohort or a multi-ancestry analysis.
Drawing main and reference samples from the same source tests finite-panel
sampling error, not all forms of population mismatch. One reused example
region is useful for smoke tests, not evidence of generalisation across loci.

**State:** open empirical gap. Add explicit sample-panel metadata, ancestry
selection, allele/build alignment and distinct locus panels. Synthetic
cohorts in the new challenge test differing LD domains and effect assumptions;
they do not replace real ancestry-specific allele-frequency and LD evidence.
Multi-ancestry methods make materially different effect-sharing assumptions.
[SuSiEx, 2024](https://www.nature.com/articles/s41588-024-01870-z),
[MESuSiE, 2024](https://www.nature.com/articles/s41588-023-01604-7).

### 11. P2 — a controlled ablation is not a competitive tuning comparison

**Original sources:** [`docs/autoresearch/iteration-004.md:139–165`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-004.md#L139-L165)
and [`scripts/tuning/make_tuning_sim.R:90–106`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/scripts/tuning/make_tuning_sim.R#L90-L106).

The focal method was selected from 962 trials on a small fixed tuning panel.
Applying its selected hyperparameters to BEATRICE can isolate a prior change
at that setting. It does not establish which method performs best after each
has had a comparable opportunity to tune. A shared setting selected for one
method can disadvantage the other. Searching many times on a tiny validation
panel also overfits that panel without direct test-label leakage.

**Required future comparison:** report author defaults, independently tuned
competitive configurations under equal declared budgets, and controlled
ablations as different tracks. Separate development, validation and final
loci or chromosomes, not only noise seeds. Annotation fitting on test summary
statistics is a transductive analysis choice that must be declared; oracle
truth must remain unavailable.

### 12. P2 — uncertainty is conditional on reused genotype panels

**Original source:** [`R/evaluate.R:423–468`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/evaluate.R#L423-L468).

The SE is computed over phenotype iteration summaries after pooling their
loci and parameter cells. For nonlinear pooled metrics, the SE of iteration
metrics is not automatically the SE of the reported pooled point estimate.
Repeated phenotypes on the same LD panels also do not create new independent
genomic loci.

**Repair/limit:** the legacy evaluator labels conditional uncertainty and uses
a delete-one-iteration estimate. The new challenge compares methods on paired
locus clusters. Resampling must preserve shared genotype/LD ancestry,
repeated-phenotype and study-overlap structure. A ten-locus panel still has ten
locus-level units for generalisation, however many phenotype draws it has.

### 13. P2 — excluding hard regimes narrows the question rather than solves it

**Original source:** [`docs/autoresearch/iteration-004.md:35–55`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-004.md#L35-L55).

Raising local PVE and removing reference LD makes an in-sample experiment
easier to interpret, but excludes the uncertainty regimes most relevant to
many summary-statistic applications. Rank-deficient LD is a genuine
fine-mapping condition, especially when variants outnumber reference samples.
Low numbers of confident calls call for more replicates and an honest
uncertainty statement; they do not justify general claims from stronger-signal
cells alone.

**Repair:** distinct matched-input and stress tracks, including singular
reference LD, incomplete causal variant observation, nulls, selection and
annotation reversal. Methods designed for in-sample LD should have their
intended-use ranking separated from deliberate robustness tests. The
fine-mapping-inf authors explicitly scope their methods to single-cohort
in-sample LD. [Author software](https://github.com/FinucaneLab/fine-mapping-inf).

### 14. P2/P3 — effect-scale and sign choices need precise explanations

**Original sources:** [`R/simulate_phenotypes.R:304–313`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/simulate_phenotypes.R#L304-L313)
and [`R/simulate_phenotypes.R:628–659`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/simulate_phenotypes.R#L628-L659).

The per-locus equal-effect path assigns positive unit effects, while the
genome-wide path uses random signs. This changes cancellation among correlated
causes. With fixed target PVE, scaling all normal effects by a constant and
scaling residual noise accordingly leaves marginal t statistics unchanged;
`effect_variance` is then an effect/phenotype-scale parameter rather than an
independent signal-strength axis. Fixed-effect priors can still react to the
scale change if supplied beta/SE units are not harmonised.

**Repair:** both phenotype generators use random signs for equal-magnitude
effects; the new challenge includes opposing effects in highly correlated
variants. Do not interpret duplicated scale settings as new independent
evidence of difficulty coverage.

### 15. P1 — null loci cannot have their false calls removed

**Original source:** [`R/evaluate.R:509–514`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/evaluate.R#L509-L514).

The curve builder returns zero true and false positives whenever a locus has
no labelled causes. High-PIP predictions on null loci consequently disappear
from the discovery denominator. This can reward a method that confidently
invents a signal, precisely the behaviour a null benchmark should detect.

**Repair:** count predicted null variants as false positives, and report null
false-call rates separately. Selected null loci form an additional conditional
challenge: association-based selection can make an unselected-null control
look misleadingly easy.

### 16. P1 — the marginal OLS calculation needs an intercept

**Original source:** [`R/simulate_phenotypes.R:792–798`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/R/simulate_phenotypes.R#L792-L798).

Using uncentred genotype cross-products with an intercept-style residual
degrees of freedom does not calculate the usual marginal GWAS regression.
Adding a constant to genotype coding or phenotype values can change the
result. Correct standard errors also depend on the centred residual sum of
squares, rather than only on the total phenotype variance.

**Repair:** centre both variables, fit the marginal slope with an intercept,
and calculate its residual standard error. Regression tests compare against
independent `lm()` fits and shifts of genotype/phenotype origins. This repair
requires fresh simulation summaries; it cannot be inferred from historical
PIPs alone.

### 17. P1 — cached neural posterior weights use incompatible learned priors

**Original sources:** [`trainer_annot.py:565–580`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/BEATRICE_annot_sparse/scripts/trainer_annot.py#L565-L580),
[`joint_trainer.py:112–123`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/BEATRICE_annot_sparse/scripts/joint_trainer.py#L112-L123),
and [`generate_credible_sets.py:326–331`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/BEATRICE_annot_sparse/scripts/generate_credible_sets.py#L326-L331).

Cached configuration weights include the learned prior at the time each model
was visited. Finalisation updates the null model but leaves historical priors
inside other weights. Credible-set expansion can then add models under a
uniform prior. The resulting normalised numbers are not one posterior under
any declared final prior; they can depend on search order and training history.

**Repair:** cache pure log Bayes factors, persist the final Bernoulli inclusion
prior, apply it once to every explored configuration including the null, and
normalise with log-sum-exp. Expanded archives and CSV PIPs now agree. Old
mixed-prior archives are rejected because their likelihoods cannot be reliably
reconstructed. Independent Gaussian likelihood and Bernoulli calculations,
including actual single-locus and joint neural optimisation, verify the repair.

**Remaining limit:** normalisation is conditional on the explored model
support; it does not certify negligible mass in unvisited configurations.
The selected-index KL in the training objective is an approximation, rather
than a full-dimensional exact ELBO. Correct cache mathematics does not prove
global optimisation or posterior coverage.

### 18. P1 — a tuning study can silently change its estimand and dataset

**Original sources:** [`scripts/tuning/fb_objective.R`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/scripts/tuning/fb_objective.R)
and [`scripts/tuning/optuna_fb.py`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/scripts/tuning/optuna_fb.py).

The R objective claims macro AP while reading pooled AP. Rotating one locus
per parameter cell confounds locus size with particular architecture settings,
and leaves per-region metadata out of alignment. The Python caller discards
failure counts, while unseeded sampling and unverified study reuse can combine
different datasets or objective implementations in one apparent campaign.
An AP-only pruning rule can also discard nondominated candidates from a
multi-objective probability/calibration search.

**Repair:** balanced fixed-locus subsets across the complete scenario grid,
aligned metadata, attempted macro AP and explicit partial-failure rejection;
seeded sampler settings and dataset/code fingerprints guard study reuse.
Multi-objective campaigns must retain their declared tradeoffs. This is a
repair to the study design, not a rerun of the historical 962-trial campaign.

### 19. P1 — integer conditional accumulators suppress additional BEATRICE signals

**Original source:** [`generate_credible_sets.py:41–85`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/BEATRICE_annot_sparse/scripts/generate_credible_sets.py#L41-L85).

Both constrained conditional-probability helpers initialise their membership
accumulator with `np.array([0]*bp)`. Scalar additions of fractional model
probabilities truncate to integers. With posterior probabilities 0.1 on the
empty model, 0.4 on `{0}`, 0.2 on `{1}` and 0.3 on `{0,1}`, the conditional
probability of variant 1 given variant 0 is 0.3/0.7. The original helpers
instead return zero. Stepwise signal selection then rejects an additional
signal even when its actual conditional probability clears the threshold.

**Repair:** floating-point accumulation, explicit failure for a zero-mass
conditioning event, and independently calculated conditional examples. A
real two-signal expansion tests selection across multiple candidate rounds.
Actual vanilla and functional CLI runs pass the strict R importer with two
sets and fractional conditional probabilities. These are execution checks,
not evidence of population-wide accuracy or a historical ranking reversal.

### 20. P1 — BEATRICE's summed conditional marginals do not measure set coverage

**Original source:** [`generate_credible_sets.py:125–175`](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/BEATRICE_annot_sparse/scripts/generate_credible_sets.py#L125-L175).

After repairing the accumulator, the set builder still stops using a sum of
conditional membership marginals. A configuration containing two members
contributes twice. For model probabilities 0.4 on `{0,1}`, 0.3 on `{0}`, 0.2 on
`{1}` and 0.1 on `{2}`, set `{0,1}` has marginal sum 1.3 but union probability
0.9. Stopping at a requested 0.95 on the sum can therefore stop below the
intended event probability. This is the same mathematical distinction as
finding 4, inside the method's own conditional set calculation.

**Repair:** stop on `P(C intersects A | G is a subset of C, explored support)`
using the final frozen posterior, with unchanged candidate order and purity
constraints. Record the actual conditional probability, conditioning indices
and target-attainment flag, including when constraints leave the target
unmet. The R importer checks identity and necessary probability bounds without
claiming those checks prove the complete joint posterior. The old sum rule is
available as an explicit diagnostic ablation with corrected accumulators;
the original source replay retains the full original algorithm.

**Remaining limit:** conditioning on other selected variants and restricting
to explored configurations does not guarantee unconditional nominal coverage.
Missing causes and unexplored configurations can still matter. The reported
null is the empty candidate configuration in this working posterior; comparing
it to a generating sparse-null label is an operational stress assessment when
those events differ. See the [mathematics and batch validation](beatrice-batched-reweighting.md).

## Conclusions the repaired suite can support

1. An evaluator defect can be established with an analytic counterexample and
   a regression test without rerunning every external method.
2. Previously saved PIPs can often be rescored, but new PVE/probability
   calibration changes require new simulations for a corrected comparison.
3. Exact matched-model inference provides a positive control for the
   simulation and score definitions. Stress performance measures robustness
   to explicit model or input changes.
4. A correct Bayesian posterior need not be calibrated conditional on the
   realised, hidden true causal count. Such truth strata are diagnostics,
   not universal calibration obligations; the scientific contract explains
   the distinction.
5. No result here certifies biological causal truth, all real ancestries or
   a universal winning method. A prospectively blinded, independently
   assessed experimental round remains future work.

## Important work still open

The new synthetic missingness and multiple-cohort cases are mathematical
mechanism checks. Real cohort variant availability, allele frequencies,
dosage/imputation quality, covariate-adjusted LD and sample overlap still need
validated empirical fixtures. The five-effect track tests misspecification
beyond the exact reference's two-effect limit; it has no exact-reference
coverage guarantee. Other frontier questions include non-additive haplotypes, correlated
traits, rare/structural variants and chromosome-boundary effects. Each needs
an explicit estimand and an input track with enough information to make the
question answerable.

An MPRA-positive allele is evidence of activity in the measured assay context,
not a complete disease-causal label; an untested or underpowered allele is not
a negative. The prospective protocol must retain those distinctions. Public
assay resources support an external evidence track, but past selection by
fine-mapping PIP must be recorded rather than repackaged as independent truth.

