# Scientific contract for the fine-mapping challenge

This document defines the mathematical claims behind the challenge. The
[protocol](challenge-protocol.md) describes its executable interfaces and
mechanisms; the [review](REVIEW.md) records original-source defects. The local
challenge is a reproducible synthetic experiment. File separation and held-out
seeds help prevent accidental leakage, but do not themselves constitute a
blind community assessment.

## 1. Define the scientific target before scoring

For a sparse model, the target indicator can be
`gamma_j = 1{beta_j != 0}`. In a sparse-plus-background model, distinguish

```text
beta_total = beta_sparse + beta_background
gamma_sparse,j = 1{j belongs to the sparse component}
gamma_nonzero,j = 1{beta_total,j != 0}.
```

The current challenge scores **sparse-component membership**. With continuous
background effects, almost every variant can have a nonzero total effect. The
term "false discovery" means a discovery outside the declared sparse target,
not proof that the allele has no biological effect. Store both coefficient
vectors. Sparse membership is also not identical to a biologically meaningful
effect-size threshold; background effects can occasionally exceed sparse
effects. Changing to a threshold target requires changing labels, priors and
claims together. Modelling broad effects is an established alternative to
strict sparsity. [Cui et al., 2024](https://www.nature.com/articles/s41588-023-01597-3).

Each score must say whether it concerns the observed variant universe or all
generating targets. Hidden causes remain in total-target recall denominators.
Observed-variant AP answers the ranking question among supplied candidates;
it is undefined when none of those candidates are targets. Report the
missing-target fraction and attainable recall ceiling separately.

## 2. The score-statistic experiment and its units

The synthetic likelihood is

```text
z | beta, R, n ~ Normal(sqrt(n) R beta, R),
R symmetric positive semidefinite, diag(R) = 1.
```

One exact construction is a centred design `X` with column squared norms `n`,
`R = X'X/n`, phenotype `y = X beta + epsilon` and independent unit-variance
Gaussian residuals. Define the score statistic `z = X'y/sqrt(n)`. Its mean is
`sqrt(n) R beta` and its covariance is `R`. Thus beta is in units of residual
standard deviations per genotype standard deviation. The public
`beta_hat = z/sqrt(n)` and `se = 1/sqrt(n)` are summaries of this **known-noise
score experiment**.

These are not exact finite-sample marginal OLS t statistics with estimated
residual standard errors. In particular, strong effects alter marginal
residual estimates. Methods that convert z to sufficient statistics using a
standardised-phenotype OLS convention can therefore be misspecified on this
input even with matched LD. An adapter must state its assumed statistic and
variance convention. The haplotype/phenotype generator is a different family
and supplies marginal OLS statistics. Agreement across both families is more
persuasive than dominance under either one alone.

The summary-data likelihood and its connection to sufficient statistics are
developed in [Zou et al., 2022](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010299).
The challenge's scale declaration is part of the experiment, not an assertion
that every published wrapper uses that convention.

`signal_z` holds noncentral signal strength approximately fixed. Because the
stress effects scale as `signal_z/sqrt(n)`, changing n alone does not create a
power-versus-sample-size experiment. For that question, hold effect sizes
fixed and vary n instead.

## 3. Exact reference means a declared small working model

The reference enumerates configurations with zero, one or two causes and
normal effect priors. Its count prior is declared, includes a null model and
must be independent of hidden final labels. For configuration S and prior
effect variance W, define

```text
h_S = sqrt(n) z_S
Q_S = n R_SS
log BF_S = -0.5 log det(I + W Q_S)
           + 0.5 W h_S' (I + W Q_S)^(-1) h_S.
```

For independent cohorts with shared effects, sum cohort score vectors `h_S`
and information matrices `Q_S` before evaluating the expression. Multiplying
independently integrated cohort Bayes factors instead corresponds to different
effect-sharing assumptions.

The `matched_prior` generator draws causal counts and effects from the
reference's declared prior, with compatible score covariance. It is the
positive control for probability calibration and configuration coverage.
Fixed-count, fixed-effect and shifted-input arms deliberately change this
distribution. The exact reference's two-effect cap is also a model limit; it
is not evidence that larger causal counts are unimportant.

There is an additional support condition for singular covariance. Under
`z ~ Normal(sqrt(n) R beta, R)`, scores lie in the range of R. If independently
sampled reference LD is singular while scores come from a full-rank study
matrix, the observed z will generally lie outside that range. No
positive-density observation then exists under the literal reference-score
Gaussian model. The same formula can still enumerate a **working RSS
likelihood/pseudo-posterior**, but its result is not exact inference from that
misspecified Gaussian generative distribution. Adding regularisation is a
method choice that changes the likelihood and must be recorded.

No reference in this suite is a universal performance ceiling. Even a true
Bayesian posterior has an expected-loss optimality claim tied to its prior,
likelihood and score, not a guarantee of winning each finite realisation or
every utility metric.

## 4. Calibration: what is and is not promised

For a correct marginal posterior `p_j = E[gamma_j | D]`, the tower property
gives `E[gamma_j | p_j] = p_j` under the joint generating model. The same
reasoning applies when conditioning on an observed input property included in
the posterior information.

It does **not** imply

```text
E[gamma_j | p_j, realised true S] = p_j
```

when S is hidden and the method marginalises over it. Revealing true S gives
additional information. Stratifying on true effect size, sparse/background
membership or realised causal count can similarly change a calibration plot.
Those strata diagnose prior shift and practical failure modes, but should not
be presented as universal obligations of a correctly marginalised Bayesian
posterior. Conversely, averaging incompatible true-S strata can conceal
practically important over- and under-confidence; retain both the marginal
positive control and the diagnostic strata.

Use raw marginal Brier score and log loss alongside reliability tables,
precision/recall and set size. These probability scores assess honest
marginal predictions even when target indicators are correlated. They do not
score the entire joint configuration distribution. Broad-bin ECE and Murphy
components are descriptive summaries, not substitutes for raw scoring.
Report the counts behind high-PIP bins and the sampling uncertainty; a bin
with a handful of calls cannot establish calibration.

The challenge currently uses explicitly clipped log loss for finite reporting;
the legacy evaluator also exposes exact log loss. A wrong endpoint prediction
has infinite exact loss. A clipping constant must be declared and held fixed
across methods. It changes extreme-event penalties; clipped reporting should
not be represented as an unqualified strictly proper score over the full
closed probability interval. Proper-score theory supports evaluating
probabilities without arbitrary rewards for confidence alone.
[Gneiting and Raftery, 2007](https://doi.org/10.1198/016214506000001437).

## 5. Credible sets, configurations and distinct signals

Marginal PIPs obey `sum_j p_j = E[|gamma| | D]`. They do not form a probability
distribution over mutually exclusive variant identities when multiple causes
are possible. Consequently,

```text
P(any cause in C | D) != sum_{j in C} p_j, in general.
P(all causes in C | D) is a different event again.
```

A credible collection of complete configurations has a clean event:
the true configuration belongs to the collection. Its declared probability
is the sum of posterior configuration masses. A signal-specific credible set
can instead concern one single-effect component. These outputs must not be
silently converted into one another.

For every nominal coverage claim, record the level, event, variant universe
and postprocessing rule. Mixing 50%, 95% and 99% collections into a single hit
fraction does not measure calibration at a common level. Native outputs from
different methods may warrant different tracks. A collection union can be a
useful follow-up list, but its size and coverage answer another question.

Also count independent discoveries. Maximum bipartite matching between
reported sets and distinct true variants prevents multiple overlapping sets
from claiming the same cause repeatedly. Report duplicate sets, reporting
frequency, distinct-signal recall and set size alongside raw set hit rate.
Neither matching nor high within-set LD proves a nominal probability claim.

## 6. Identifiability and unavailable targets

With perfectly identical genotype columns, `R[,j] = R[,k]`, moving one effect
between j and k produces the same score distribution. If the priors and
annotations are exchangeable, the data cannot identify which twin is causal.
Signed equivalent columns create the same issue when effect signs can change.
An exact variant metric should still reward correct uncertainty: for an
exchangeable one-cause pair, PIPs of one half have a positive but unavoidable
loss. It must not treat an arbitrary certain guess as an attainable scientific
ideal.

The LD-class top-rank diagnostic in this suite assesses localisation to an
observationally equivalent class. It is not a class posterior probability.
Summing class-member PIPs gives an event probability only under additional
mutual-exclusivity assumptions; otherwise a joint posterior is needed.
Near-perfect LD is not exact equivalence and should be reported separately.
Resolution depends on sample size, effect strength and alternative
configurations as well as pairwise LD.

Generate the full signal before hiding causal variants. The observed mean must
remain `sqrt(n) R_observed,all beta_all`, including the omitted contribution.
Deleting the effect before generating z would create a null locus rather than
a missing-cause problem. Without extra information, an omitted causal allele
in exact LD with an observed allele may be intrinsically undetectable.
Evaluate false certainty, observed ranking, total recall and attainable recall
separately; do not demand an impossible diagnosis in every missing-cause case.

## 7. LD mismatch, orientation and heterogeneous samples

Finite reference sampling, population LD shift and inconsistent allele
orientation are different mechanisms. Matched coding changes have the form

```text
z_new = D z
R_new = D R D
beta_new = D beta,    D diagonal with entries +1 or -1.
```

PIPs should be invariant to this simultaneous recoding. Changing z alone is
an error mechanism. Variant permutations must preserve predictions after IDs
are aligned; tie breaking must not depend on genomic index. Truth covariance
and coefficients must use the same declared orientation too.

With study-specific variant inclusion, a common scalar sample size and nominal
LD matrix are insufficient. For independent study scores with standardised
meta-analysis weights `w_ja`, the noise covariance is

```text
Cov(z_j, z_k) = sum_a w_ja w_ka R_a,jk.
```

In the two equal-size, same-LD study construction, a variant included in both
studies and one included only in the first have covariance `R_jk/sqrt(2)`.
Each variant still has variance one, and its mean scales with the square root
of its own effective sample size. Effect/SE summaries must therefore use that
variant's sample size. An exact heterogeneous model needs the appropriate
mean map and covariance, not merely insertion of `n_j` into a scalar-n formula.
Exposing `n_j` does not expose all study-overlap information.

This mechanism has real relevance: heterogeneous cohort processing and
variant availability affect meta-analysis fine-mapping, while QC can diagnose
some inconsistencies. A QC method can also reject genuine multiple-signal
loci; detection performance and retained-locus fine-mapping performance need
separate denominators. [Kanai et al., 2022](https://pubmed.ncbi.nlm.nih.gov/36643910/),
[DENTIST, 2021](https://www.nature.com/articles/s41467-021-27438-7).

## 8. Discovery selection is a conditional experiment

The selection arms retain scores conditional on a declared maximum-|z| rule.
They test selected nulls and weak signals that cross a discovery threshold.
The target pool, attempted draws and selection rule must remain in the
assessor record. Null selection is intentionally a tail experiment; it is
not the population frequency of false GWAS loci without the discovery
denominator.

Selection does not automatically invalidate a correct Bayesian posterior. If
selection is a deterministic function of the complete observed data D, it
adds no information after conditioning on D. However, fixed-null or
fixed-effect selected stress distributions differ from a mixed prior, and
selection on discovery data omitted from the predictor input can add
information. Distinguish those cases from universally requiring a selection
correction. Winner's curse in effect estimates is a related but separate
question from causal membership calibration.

A future paired discovery/replication track should specify overlapping versus
independent participants, the association selection policy and which study
summaries are provided. Replication of a high PIP is an external consistency
endpoint, not complete causal truth.

## 9. Annotation learning and fair comparisons

Annotation enrichment, neutrality, correlation, nonlinear interaction,
depletion and transfer failure are separate hypotheses. Positive annotations
alone create an easy home-model advantage. An informative-versus-reversed
comparison assesses a prescribed prior's failure under transfer; demonstrating
that an annotation-learning method actually transfers requires training on
one set and predicting on another.

Avoid using the target locus to learn its own prior in a track advertised as
held-out. Hold out chromosomes or independent locus families where possible.
Pooling test summary statistics is allowed in an explicitly transductive
track, with equal data exposure for competitors. A causal annotation model
and an effect-variance annotation model are not generally interchangeable.

Report author defaults, individually tuned competitive models and controlled
ablations separately. Give comparable declared tuning budgets, freeze
hyperparameters before final assessment and report all attempted methods.
Oracle enrichment can illuminate an information bottleneck, but is not a
normal participant. A baseline that forces one cause is deliberately useful
for exposing null overconfidence, not a recommended general model.

## 10. Uncertainty, aggregation and resource fairness

Pair methods on the same targets. Bootstrap independent locus clusters,
keeping repeated phenotypes and shared LD or study structures together.
When all generated loci are independent, one task per cluster is appropriate.
When the generator reuses one LD panel, merely assigning a new task ID does
not create an independent cluster. Conditional noise uncertainty and
generalisation to new genomic loci are different claims.

Macro localisation gives loci equal weight. Micro discovery precision weights
calls, and pooled probability scores weight variants. These are different
estimands. Mechanism weights in an overall score are scientific choices;
there is no unique universal aggregate winner across an arbitrary stress mix.
Report trade-offs between power, calibration, resolution, availability and
resources. Do not select a PIP threshold on final truth, or infer control at
all thresholds from the best observed operating point.

Conditional successful-fit Brier/AP/log loss require the success column.
Attempted-task recall can implement a failure-means-no-discovery policy.
Neither should silently replace the other. The current runner records elapsed
time, but it does not enforce a container resource limit, measure peak memory
or fully charge failure/retry/setup costs. Those remain requirements for a
fair public competition, including separate CPU and accelerator tracks.

## 11. A credible prospective 2027 assessment

CASP's transferable contribution is its institutional experiment: targets
whose reference measurements are unavailable during prediction, fixed
submission deadlines and independent assessment. The organisers also seek
targets where current methods struggle. An open synthetic generator and a
`final` filename reproduce neither confidential truth custody nor independent
biological measurement. The target embargo includes preprints and public
images, not just journal publication. [CASP17 target instructions](https://predictioncenter.org/casp17/targets_submission.cgi),
[CASP assessor policy](https://www.predictioncenter.org/index.cgi?page=assessors_list).

A proposed fine-mapping round should proceed as follows:

1. **Pre-register endpoints and select target loci.** Select loci using a
   declared association or sampling policy, without choosing variants by the
   methods being scored. Define trait, tissue, condition, ancestry, variant
   universe, label uncertainty and desired resolution before predictions.
2. **Recruit independent experimental providers.** Obtain unreleased,
   context-matched measurements with positive controls, biological replicates,
   editing verification and assay-power assessment. Keep target labels and
   manuscripts under an agreed embargo through the submission deadline.
3. **Release frozen predictor inputs.** Supply aligned variant identities,
   alleles, build, cohort-level statistics, LD provenance and allowed
   annotations. Distinguish summary-only, annotation-informed, multiple-cohort
   and genotype/haplotype-input tracks. Do not expose truth-bearing target
   mechanism names or investigator annotations derived from the withheld assay.
4. **Measure candidates without method-specific verification bias.** Where
   feasible, test the candidate universe within a locus. Otherwise use a
   preregistered probability sample with known inclusion probabilities and a
   fixed follow-up budget. Untested variants remain unknown. Partial-label
   scoring needs a justified sampling/measurement model rather than treating
   all omitted variants as negatives.
5. **Freeze submissions and assess independently.** Archive the container,
   software, training sources, timestamp and resource use. Assessors should
   not tune competing entries on the final labels. Publish every task,
   withdrawal and failure, uncertainty intervals and method trade-offs.
6. **Release evidence and retain uncertainty.** Publish measurement records,
   selection provenance and adjudication alongside scores after embargo.
   Maintain an ongoing time-split evidence track so methods cannot repeatedly
   optimise against one static final panel.

The first plausible experimental endpoint is a **context-specific molecular
effect**, paired with an association dataset in the same system and endogenous
single-variant editing. This can give much stronger variant evidence than a
generic reporter assay. It still does not automatically establish that a
molecular perturbation mediates a human disease phenotype. Editing efficiency,
bystander edits and phenotype measurement need their own uncertainty model.
[Ryu et al., 2024](https://www.nature.com/articles/s41588-024-01726-6).

An existing large external resource is the 2026 Siraj et al. MPRA study of
221,412 fine-mapped variants and additional controls in five cell types. It is
a valuable assay-evidence track, but its candidate selection used PIPs and
credible sets, and high/low PIPs also served as causal-status proxies.
Reusing those proxies as independent fine-mapping labels would be circular.
MPRA-negative does not mean disease-noncausal; reporter context and assay power
matter. Its measured effects can be scored under an explicit assay endpoint
and recorded ascertainment. [Siraj et al., 2026](https://www.nature.com/articles/s41586-026-10121-6).

## 12. Frontier work remains a research programme

The present challenge offers explicit mechanism tests, not complete coverage
of fine-mapping. The following extensions require new validated inputs and
new experiment contracts:

| Frontier | Required scientific distinction |
|---|---|
| Real multiple-ancestry cohorts | Shared per-allele effects versus shared standardised effects; ancestry-specific MAF and LD; shared versus distinct causal support |
| Real meta-analysis | Variant-wise sample masks, imputation/dosage quality, covariates, participant overlap and allele/build QC |
| Multiple traits | Any-trait versus trait-specific causal events; residual trait correlation and effect sharing |
| Haplotype and non-additive effects | Marginal additive summaries may contain insufficient information; use a genotype/haplotype input track |
| Rare, indel and structural causes | Variant representation, ascertainment and LD accuracy rather than only common SNP toy matrices |
| Boundary and long-range LD | Missing outside-window effects and alternative region definitions |
| Very many effects | Model-cap failure versus inference failure; scalable references instead of exhaustive two-effect enumeration |
| Prospective experimental truth | Assay-defined positives, adequately powered negatives and unknowns; independent target selection and embargo |

For real multiple-ancestry sharing, a useful score mean is
`sqrt(n_a) R_a D_a beta_per_allele / sigma_a`, where `D_a` contains genotype
standard deviations. Sharing standardised beta directly is a different
assumption. Multiple-trait models likewise need residual trait covariance;
simply repeating a univariate locus does not benchmark joint inference.
[mvSuSiE, 2026](https://www.nature.com/articles/s41588-025-02486-7).

Progress should be measured by which uncertainty or information bottleneck a
new target resolves, not by the number of scenarios in a manifest. The next
release can add empirical fixtures and participant adapters while retaining
this synthetic suite as a fast, falsifiable set of scientific controls.
