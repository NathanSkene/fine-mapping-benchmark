---
title: "What makes a fine-mapping benchmark trustworthy? A comparative reassessment and a prospective challenge framework"
subtitle: "Working comparative manuscript — campaign results not yet complete"
date: "10 October 2026"
lang: en-GB
fontsize: 10pt
geometry: margin=0.9in
colorlinks: true
---

**Manuscript status.** This is a research draft for Nathan Skene's review.
It credits Lucas Siu's original software and methodological contributions,
but does not assign a publication author list or claim that collaborators have
approved its contents. Nathan's comparison target is the corrected suite
developed here, not a separate pre-existing fine-mapping method. The corrected
full server campaign is in preparation.
The results section therefore records verified development evidence and the
historical claims to be tested, not fabricated final rankings. The completion
contract is [results-requirements.md](results-requirements.md).

# Abstract

Fine-mapping methods turn association evidence and linkage disequilibrium into
probabilities about genetic variants. Comparing those probabilities requires
more than a leaderboard: the benchmark must define the causal event, retain
failed and difficult targets, preserve each method's intended information set,
and distinguish localisation from calibrated confidence. We reassess Lucas
Siu's open fine-mapping benchmarking framework at a frozen source revision.
The framework provides a useful common interface, haplotype-based simulation,
multiple method adapters, cross-region functional-prior experiments and an
explicit record of iterative criticism. Our review identifies issues
in scoring, simulation, posterior calculation, comparator identity and
experimental design. We repair these issues and add a 22-track synthetic
challenge with explicit nulls, missing causes, LD ambiguity and mismatch,
annotation transfer, discovery selection and multiple-cohort inputs. We
separate historical rescoring, implementation repair and new-data comparisons
to make changes in method ranking attributable. Completed development runs
and native-method integration checks establish that the infrastructure works;
they do not yet establish a full method ranking. The final comparative campaign
will report localisation, probability quality, uncertainty coverage, attempted
recovery and resource use by setting. We also propose a prospective assessment
based on independently withheld molecular evidence and reusable environments
for training agents to make statistically defensible fine-mapping decisions.
The central contribution is an assessable comparison contract, rather than an
assumption that one method should win every problem.

# Introduction

Statistical fine-mapping asks which variants explain an association, given
that neighbouring variants often carry similar information. A successful
method should identify a causal variant when the data distinguish it, and
acknowledge ambiguity when they do not. Ranking variants and quantifying that
uncertainty are related tasks, but they are not the same task. A method can
place the correct variant near the top while assigning much too much
probability to the remaining variants. It can also be appropriately uncertain
at a locus where several variants have identical likelihoods. A benchmark that
rewards confidence or a sharp-looking credible set without checking its event
probability can favour the wrong behaviour.

This problem is particularly clear in high LD. SuSiE was designed in part to
report uncertainty over correlated candidate variables instead of arbitrarily
selecting one. Its signal-level credible sets and variant PIPs summarise
different aspects of a fitted posterior. [Wang et al., 2020](https://doi.org/10.1111/rssb.12388).
Summary-statistic inference adds another layer: score, marginal effect and LD
inputs must correspond to a coherent sampling experiment. [Zou et al., 2022](https://doi.org/10.1371/journal.pgen.1010299).
Comparisons can fail before the fine-mapping algorithm is even reached if the
input scales, covariances or candidate identities do not agree.

Siu's repository is a useful starting point because it brings methods into one
pipeline, supplies haplotype-derived genotype and phenotype simulations, and
already examines annotation learning and reference-panel LD. The original
reports identify substantial limitations, including poorly calibrated
functional-prior inference, incorrect PAINTOR pooling in early runs, and
changes in experiment design that forbid naive comparisons between
iterations. [Original software and reports](https://github.com/lucassiu14/fine-mapping-benchmark/tree/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81).
Our work continues that critical process. We aim to make the scientific
question and its evidence clear enough that an unexpected result can teach us
something about inference, rather than merely reflect an unnoticed change in
the evaluator.

The reassessment has three purposes. First, it establishes which original
quantitative conclusions can be checked from archived evidence or reproduced
inputs. Second, it tests the corrected methods across mechanisms that the
original grids did not include. Third, it turns the benchmark into a renewable
assessment and agent-training interface, with explicit reward and information
boundaries. These purposes share infrastructure but need separate evidence.
A new synthetic stress suite is not a historical rescoring, and a runnable
agent environment is not proof of commercial demand or scientific discovery.

# Questions and comparison scope

The paper evaluates five scientific questions.

1. How much do evaluator corrections change measured performance when
   predictions and truth are held fixed?
2. How much do repaired inference calculations change predictions on the same
   inputs, particularly when a learned prior changes during optimisation?
3. When does sharing functional-prior information across regions improve
   localisation, and when does it produce overconfidence or negative transfer?
4. Which methods preserve useful uncertainty under nulls, LD ambiguity,
   missing candidates, nonsparse effects, selection and input mismatch?
5. What evidence is required to move from development simulation to a blind
   community assessment, or to a training environment with meaningful rewards?

"Original" refers to the frozen upstream source
`f053fdf318a3abe2c44eeb0ac3f3c99a03979b81`. "Corrected" refers to the repaired
fork and its recorded source hashes. These are benchmark and implementation
versions; they are not rival ownership claims. Functional BEATRICE and its
shared-prior variants existed in the original repository. Nathan's supervision
records describe a collaborative programme, and the current benchmark repairs
should not be described as inventing those methods independently.

The comparison includes deployable methods, annotation-free and
annotation-informed variants, local surrogate implementations and explicit
controls. We identify each of these roles. A simulator-informed oracle can
show the value of information, but it is not an ordinary participant. A local
mixture sampler called SBayesRC is not interchangeable with native GCTB, and
the standard `--sbayes RC` command is not automatically the newer GWFM
implementation described in the 2026 fine-mapping study. Configuration and
version are part of the method identity.

The full campaign, including optional and native methods, is required for the
final comparative conclusions. Methods that cannot accept a track's declared
inputs remain visible as unsupported attempts. For example, a single-cohort
method does not become a multiple-cohort method by silently discarding the
second cohort. Its performance can be described within its supported domain,
and its availability can be reported across the operational task set, without
confusing the two.

# Historical evidence and the original experiments

## What the original reports establish

The original Git tree contains analysis code and narrative reports; it does
not contain the historical per-locus PIP archives or collected tables. The
notebook has no stored output. We therefore distinguish documented historical
claims from independently recomputed findings. The complete claim ledger and
immutable links are in [the original-results audit](../docs/original-results-audit.md).
Access to an original archive would permit rescoring. Recreating the experiment
from frozen source permits a reproduction, but not a claim that the exact
historical outputs have been recovered.

The Iteration 003 report describes improvements from sharing a functional
prior across loci. It reports AP 0.539 for `fb_xregion` on continuous
annotations with in-sample LD, and describes it as the best deployable method
in that setting. It also reports that the joint LassoNet head outperformed
the linear head on AP in 62–97% of cells across the reported strata. With a
reference panel of 200, reported mass ratios fell from 11.45 to 5.04 for
continuous annotations and from 9.02 to 5.41 for binary annotations. Yet the
same report says `fb_xregion` beat FINEMAP on AP in only 2–10% of cells outside
in-sample LD. These observations motivate hypotheses about shared-prior
learning and LD sensitivity; they do not yet settle them under the repaired
scoring and inference definitions. [Iteration 003 findings](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L246-L265).

We also retain the original report's self-corrections. Early PAINTOR results
are explicitly withdrawn because its enrichment estimation was incorrectly
run one locus at a time. CARMA was registered but absent from early executed
comparisons. Neither an invalid configuration nor a missing installation is
evidence of poor performance by the published method. [Original method audit](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L276-L374).

## Experiments that must not be pooled

Iterations 002/003 use 135 grid rows, 250 scenarios per row and ten regions
per scenario. The genotype sample size is 1,000; region sizes run from 100
to 1,000; causal counts are 1, 2, 3, 5 and 10; and LD is in-sample or estimated
from panels of 500 or 200. The sparse-plus-infinitesimal arm varies the sparse
proportion of genetic variance. Annotation type and enrichment are additional
axes. The joint methods reuse the original simulations, allowing paired data
exposure against earlier baselines.

Iteration 004 changes sample size to 5,000, region sizes to 500–2,000, the
heritability and sparse-variance grids, and removes reference-panel LD. It
contains 11,250 scenarios and nineteen method IDs. Its document is marked
running and contains early execution QC, rather than a completed final
ranking. It explicitly states that its results are not comparable to
Iterations 001–003. [Iteration 004 specification](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-004.md).

The new challenge changes the data family again. A difference between an old
haplotype/OLS result and a new Gaussian-score result can reflect the input
convention, architecture, sample information or method implementation. It
must not be called a ranking reversal caused by an evaluator repair unless
the paired evidence isolates that repair.

The planned expanded campaign contains the complete 45-row Iteration 004 grid
and 135 Iteration 003 control rows, with 250 scenarios per row (45,000 scenarios
in total). It will evaluate original-source refits, corrected refits of those
same inputs, and corrected phenotype resimulations on the same frozen genotype
panels. Both evaluators will score every archived prediction. These are new-seed
upstream-code replays, not recovered historical outputs. The planned challenge
adds 13,200 tasks: 22 tracks, 100 draws, three variant counts (40, 200 and 1,000),
and two sample sizes (1,000 and 10,000). These counts describe the requested run
design until manifests and completed output archives verify execution.

# Methods

## Four attribution contrasts

We organise comparison evidence around the following contrasts.

| Contrast | What is fixed | What the result can establish |
|---|---|---|
| Historical rescoring | Predictions, truth, candidate IDs and locus identities | Change due to the evaluator |
| Inference repair | Inputs, method settings and declared stochastic budget | Change due to implementation |
| Simulation repair | Source genotype family and matched generation policy | Change due to the generating experiment |
| New challenge | Frozen new inputs and participant information contract | Behaviour under newly tested mechanisms |

Only the first contrast is a direct update to an old ranking with predictions
held fixed. The second can test a posterior bug without changing the
phenotype. The third addresses whether the originally labelled model was
actually simulated. The fourth extends the questions the benchmark can ask.
When historical artifacts are unavailable, the paper must state that it uses
a reproduction at the original source, and record any dependency or source
panel differences.

## Causal target and probability event

For a sparse model, the target may be written as

$$\gamma_j=\mathbf{1}\{\beta_j\ne0\}.$$

With broad background effects, we instead record both components,

$$\beta_{\mathrm{total}}=\beta_{\mathrm{sparse}}+
\beta_{\mathrm{background}},\qquad
\gamma_{\mathrm{sparse},j}=\mathbf{1}\{j\text{ belongs to the sparse component}\}.$$

The present challenge scores sparse-component membership. Almost all variants
can have nonzero total effects when a continuous background is present. A
"false discovery" outside the sparse target therefore does not prove that
the variant has no biological effect. Nor is sparse membership an effect-size
threshold: a background draw can occasionally exceed a sparse effect. A
different scientific target requires a different truth and scoring contract.
The value of modelling nonsparse effects is established in the SuSiE-inf and
FINEMAP-inf work; our benchmark should assess that target clearly rather than
assume that all small effects are literally noncausal. [Cui et al., 2024](https://doi.org/10.1038/s41588-023-01597-3).

Mixture-model PIPs require particular care here. Both standard native SBayesRC
and the repository's local SBayesRC surrogate
reports the probability of membership in any nonzero effect component; that
is not necessarily the probability of membership in the simulator's planted
sparse component. With continuous background effects, the full nonzero-effect
label can be positive for nearly every variant. The campaign therefore
reports native-PIP ranking against the sparse target as an operational
localisation diagnostic, and separately scores the declared nonzero-effect
event. It does not call native PIPs miscalibrated sparse probabilities merely
because they assign mass to small background effects. Sparse-target FDR and
native-PIP calibration must remain distinct in conclusions.

Candidate-universe and generating-universe scores also differ. AP ranks
observed candidates. Total-target recall retains causes omitted from those
candidates, so the attainable recall ceiling is explicit. A hidden cause is
not erased from the denominator just because no method can predict its ID.
The benchmark records that information bottleneck instead of blaming an
algorithm for failing to identify an unprovided variant.

Null forecasts require the same distinction. The Gaussian assessor's null
label means that the entire generating universe has no sparse targets.
BEATRICE's reported null probability instead concerns the empty modeled
configuration among observed candidates, conditional on explored support and
the frozen learned prior. Comparing that number with the generating-null
label is an operational forecast assessment under model misspecification.
It cannot establish calibration for the same event when an omitted cause or
background architecture separates the model's event from the generating
target. The sparse-target null also differs from absence of every nonzero
total effect under a continuous background.

The current R challenge accepts a candidate-configuration null forecast and
checks both union bounds, $\max_jp_j\le1-q\le\sum_jp_j$. Its optional null
Brier score uses the generating sparse-target label; that interpretation is
stated explicitly. The RL interface below instead targets the whole
generating locus and requires only the upper bound on individual observed
PIPs. Future submissions should declare the null event directly, rather than
using one untyped field for these distinct claims.
[Challenge event validation and scoring](../R/challenge_evaluation.R),
[Generating truth](../R/challenge_simulation.R).

## Corrected phenotype construction

The repaired genotype/phenotype family centres the marginal regression and
checks effect estimates against an intercept-containing ordinary least-squares
fit. The original uncentred calculation could distort summaries when finite
sample means were nonzero. The variance ledger retains the cross-covariance
between sparse and background genetic components:

$$\operatorname{Var}(g_s+g_b)=\operatorname{Var}(g_s)+
\operatorname{Var}(g_b)+2\operatorname{Cov}(g_s,g_b).$$

Separately rescaling the two components does not make their realised
cross-covariance vanish. We therefore distinguish intended component scales
from realised total PVE, retain background coefficients, and verify the total
phenotype construction. Causal-selection probabilities likewise require a
joint prior with the stated prevalence. Clipping individual probabilities
after normalisation changes that prevalence and can change the experiment.

These repairs do not imply that old conclusions necessarily reverse. They
establish that the new model label corresponds to a checked construction, so
a comparison is interpretable.

## Gaussian score challenge and exact reference

The synthetic summary-data family specifies

$$z\mid\beta,R,n\sim N(\sqrt{n}R\beta,R),\qquad
R=R^\mathsf{T}\succeq0,\quad\operatorname{diag}(R)=1.$$

One construction uses centred columns of a design matrix with squared norms
equal to n, unit residual variance and $z=X^\mathsf{T}y/\sqrt n$. The effect
unit is a residual standard deviation per genotype standard deviation. The
public summaries $\hat\beta=z/\sqrt n$ and $\mathrm{se}=1/\sqrt n$ belong to
this known-noise score experiment. They are not exact finite-sample marginal
OLS t-statistics with estimated residual variance. A wrapper that assumes a
different convention can be misspecified even when its LD matches.

The corrected SuSiE score adapter supplies $X^\mathsf{T}X=R$ and
$X^\mathsf{T}y=z$, treating its coefficient as $\theta=\sqrt n\beta$.
It fixes the known residual variance and scales the coefficient prior variance
by n, then converts posterior means back to the effect scale. Passing these
scores through the ordinary `susie_rss(z,R,n)` interface would instead invoke
its OLS-t adjustment and alter the experiment. The current SuSiE-inf adapter
explicitly rejects this score input: disabling residual-variance estimation
also disables its infinitesimal-variance update, whereas enabling it would
estimate variance from fictitious phenotype sufficient statistics. Neither
route implements the declared fixed-noise background model. Individual-data
and compatible OLS fits remain supported. Earlier Gaussian SuSiE/SuSiE-inf
prototype results are superseded by this input-contract correction and do not
support the final comparison. [Typed score adapter](../R/gaussian_score_adapter.R),
[SuSiE-inf adapter](../R/wrapper_susie_inf.R).

The finite-sample OLS track has a different reconstruction. For marginal
t-statistic t, the genotype–phenotype correlation is
$r=t/\sqrt{n-2+t^2}$. The FINEMAP-inf adapter now supplies
$z_{\mathrm{native}}=\sqrt{n\,\mathrm{meansq}}\,r$, so the author implementation
constructs $X^\mathsf{T}X=nR$, $X^\mathsf{T}y=\sqrt n\,z_{\mathrm{native}}$
and $y^\mathsf{T}y=n\,\mathrm{meansq}$ consistently. The default
`meansq=1` declares population-SD normalisation of genotype and phenotype;
it does not assert that the archived raw phenotype has variance one. This
likelihood repair is part of the same-phenotype implementation contrast.
[FINEMAP-inf adapter](../R/wrapper_finemap_inf.R),
[Author implementation](https://github.com/FinucaneLab/fine-mapping-inf).

The exact reference enumerates zero-, one- and two-effect configurations.
For configuration C and normal effect variance W, its Bayes factor is

$$\log\mathrm{BF}_C=-\tfrac12\log\det(I+WQ_C)
+\tfrac12W h_C^\mathsf{T}(I+WQ_C)^{-1}h_C,$$

where $h_C=\sqrt n z_C$ and $Q_C=nR_{CC}$. A declared causal-count prior
includes the empty model. For independent cohorts with shared effects, scores
and information matrices are summed before integration. Multiplying separately
integrated cohort Bayes factors would imply a different effect-sharing model.

Enumeration is exact only for its declared working model and compatible
sampling inputs. It is not an upper bound for every comparator. Five-effect
data violate its support; unknown overlap violates its input contract; and
incompatible singular reference LD can invalidate the presumed Gaussian
support. Results in those settings are working-model or pseudo-posterior
diagnostics, not calibrated ground truth.

The matched-prior track draws counts and effects from the exact reference
prior and is a positive control for calibration. Fixed-count and fixed-effect
stress tracks assess conditional behaviour. A calibrated Bayesian posterior
under its mixed prior need not have exactly nominal calibration after
conditioning on hidden S or a particular effect size. We report that change
of distribution explicitly.

## Challenge mechanisms

The suite contains 22 named tracks, including invariant and placebo controls.
It does not contain 22 independent genetic architectures. In particular,
the irrelevant-annotation track repeats the one-signal control distribution
with fresh draws; nuisance annotations are also present in baseline tracks.
No aggregate score silently gives this architecture extra weight.

| Family | Mechanisms | Scientific question |
|---|---|---|
| Positive and null controls | Matched prior, null, one/two signals | Does probability quality agree with the declared model and null event? |
| Architecture stress | Five signals, opposing high-LD effects, broad background | Do model-size limits and masked marginal signals change recovery? |
| LD and representation | Small/shifted reference, perfect LD, allele error/recoding | Does the method handle mismatch and preserve coding invariance? |
| Annotation evidence | Informative, irrelevant, reversed | Does useful evidence help without unsupported confidence or transfer? |
| Candidate and discovery selection | Missing cause, selected null, selected weak signal | What remains recoverable, and which denominator defines performance? |
| Cohort information | Variable n, known overlap, shared/heterogeneous effects | Are covariance and effect-sharing assumptions respected? |

An ancestry-labelled track changes Gaussian covariance between cohorts. It
does not represent real population allele frequencies, admixture, imputation
or demographic history. It is a covariance-domain experiment. Population
transportability needs empirically grounded inputs and per-allele effect
conventions, not an ancestry label on a correlation matrix.

Selected-null and selected-signal inputs condition on a declared discovery
threshold. Deterministic selection using all provided data adds no new
information after conditioning on those data under a correct Bayesian model.
The stress distribution can nevertheless differ from the prior predictive
distribution. If the discovery data are omitted or participants overlap across
studies, additional modelling may be needed. The benchmark does not assume
that every selection requires the same correction.

## Corrected scoring

We compute AP at distinct observed score thresholds, admitting exact ties as
a group. This avoids the original fixed-grid compression of low PIPs and
method-order tie breaking. AP remains undefined where no observed target is
positive. It is a localisation score; it does not validate PIP magnitudes.
Any auxiliary failure-zero AP uses the same observed-positive locus universe
for successful and failed attempts. Null loci remain in completion and
false-call assessments instead of entering an AP denominator only on failure.

For observed sparse-membership labels, probability quality uses the raw Brier
score and endpoint-preserving log loss,

$$B=\frac1p\sum_j(p_j-\gamma_j)^2,\qquad
L=-\frac1p\sum_j\{\gamma_j\log p_j+(1-\gamma_j)\log(1-p_j)\}.$$

An incorrect probability of zero or one has infinite log loss; infinite
outcomes and their frequency are retained. A separately labelled numerical
diagnostic clips probabilities at $10^{-15}$. The raw Brier score is not renamed
from a coarse binned Murphy decomposition. Binning is useful for reliability
plots, but approximation and within-bin variation have to be retained when
interpreting a decomposition. Proper scores incentivise honest probabilities
under their target distribution. [Gneiting and Raftery, 2007](https://doi.org/10.1198/016214506000001437).
With sparse positives, an all-zero predictor can have a small average Brier
score; discovery, recall and positive-event diagnostics therefore accompany
the proper score instead of treating it as a complete leaderboard.

At prespecified PIP thresholds, FDR is the number of non-target calls divided
by all calls. We pool counts when estimating discovery precision rather than
average locus rates with different denominators. No calls means the FDR is
undefined, not zero. Null loci remain in false-call assessment. Failed attempts
remain in operational recovery denominators and are not silently discarded.
They are also not fabricated as valid all-zero posteriors; successful-fit
proper scores and attempted-task completion/recovery are reported separately.

For a candidate set A, the sum of PIPs has the meaning

$$\sum_{j\in A}p_j=E\!\left[\sum_{j\in A}\gamma_j\mid D\right],$$

not generally $P(\exists j\in A:\gamma_j=1\mid D)$. A sum can exceed one.
Native credible sets are therefore evaluated under their declared event.
Set hit rate concerns at least one sparse target; maximum matching between
sets and targets assesses distinct recovery without double credit for
duplicates. A credible collection of full configurations is a separate
posterior event and is assessed by whether the complete truth configuration
appears in the collection. Its size and declared nominal level accompany
coverage. These distinctions prevent one summary from being judged as
another. [SuSiE credible-set definition](https://doi.org/10.1111/rssb.12388).

CAVIAR returns a variant superset intended to contain all causal variants.
We assess complete containment, causal-variant recall and support size; a
single correct superset is not penalised as recovering only one signal.
Its distinct-per-signal set metric is consequently inapplicable.
[CAVIAR](https://pmc.ncbi.nlm.nih.gov/articles/PMC4196608/).
CARMA returns several LD-grouped sets. The original wrapper merged these into
one set, losing the grouping needed to assess distinct recovery. We preserve
the native groups and describe their construction: within candidate LD
neighbourhoods, members are added in descending PIP order until their summed
PIP exceeds the chosen level. The source selects the LD threshold producing
the most groups. We assess empirical set hit rates, sizes and distinct-target
recovery, as appropriate to the authors' intended at-least-one-target event.
The PIP-mass rule alone does not establish an exact joint coverage probability.
For example, two variants that are jointly causal with probability .6 each
have PIP .6; their sum exceeds .99, while the probability of any causal variant
is still .6. This matters when several causes occupy one LD group.
[CARMA set construction](https://github.com/ZikunY/CARMA/blob/51f0fcf6a74ee74c75fd2ce70912fdd6cdd214fe/R/CARMA.R#L146-L200),
[CARMA assessment definitions](https://www.nature.com/articles/s41588-023-01392-0).
FiniMOM's signal sets instead target a causal variant, with their stated
coverage conditioned on the selected causal-count model. Target recovery and
nominal conditional coverage are separate questions.
[FiniMOM methods](https://oulurepo.oulu.fi/bitstream/handle/10024/45769/nbnfi-fe20230908121884.pdf?isAllowed=y&sequence=1).

## Learned priors and posterior consistency

A learned annotation prior can change during training. The corrected BEATRICE
fork memoises the data likelihood or log Bayes factor for a configuration,
then applies a single final learned prior consistently across every retained
configuration, including the null. The original mixed cache could combine
weights evaluated under different historical priors. That collection cannot
represent one coherent final posterior.

The correction does not make the search exhaustive. The posterior remains
conditional on the retained configuration support, and the selected-index KL
term remains an approximation. Numerical guards prevent nonfinite or singular
samples from being silently counted as valid training evidence. The method's
original deep variational design is retained and identified; the repair
addresses posterior coherence rather than claims to invent a replacement
inference family. [BEATRICE, 2024](https://doi.org/10.1093/bioinformatics/btae590).

An independent audit found an additional scientific bug in conditional signal
selection: integer-valued accumulator arrays truncated fractional posterior
weights to zero. For example, configuration probabilities .1 for the null,
.4 for {0}, .2 for {1}, and .3 for {0,1} imply
$P(1\in C\mid0\in C)=.3/.7$, while the original eligible-variant accumulator
returned zero. Floating-point accumulation restores this conditional
probability and can change later signal selection. This is an intentional
inference correction. The separate batched reweighting optimisation is tested
against the corrected eager calculation, not against the faulty original
selection. [Correction and analytical tests](../docs/beatrice-batched-reweighting.md).

BEATRICE's sequential sets also need an explicit probability event. Let
$\mathcal T$ denote the explored configuration support, $G$ the selected MAP
variants conditioned on at that step, and $A$ a candidate set. The relevant
conditional union mass is

$$q(A;G,\mathcal T)=
\frac{\sum_{C\in\mathcal T:G\subseteq C,\ C\cap A\ne\varnothing}w_C}
{\sum_{C\in\mathcal T:G\subseteq C}w_C}.$$

Summing conditional PIPs instead gives the expected number of variants in
$A$, which can exceed this union probability. The corrected selector targets
the conditional union event and retains the original sum rule as a labelled
diagnostic ablation. Set-level metadata must state the conditioning variants,
explored support and whether the target was reached. All sets are assessed
against the final memo subsequently persisted with the PIPs; zero supported
conditioning mass is rejected. A .95 conditional target under searched
support is not an unconditional .95 guarantee under the full model. Search
omissions and uncertainty about the conditioned MAP variants remain
scientifically relevant.

Local end-to-end checks ran the repaired vanilla and functional methods for
500 training iterations on eight-variant inputs. Each emitted two singleton
sets; the R wrappers imported their ordered variant identities, conditioning
variants, final-support null probability and union-target attainment flags.
This verifies the execution and event-import contract on that small example,
not empirical coverage or a method ranking.
[Recorded two-method check](../docs/validation-artifacts/campaign-hardening-20261010/beatrice-events-proof.json).

## Functional-prior and native-method comparisons

Functional annotations can improve localisation if they carry relevant,
transferable information. A useful-only generator creates a favourable model
match but cannot establish resistance to irrelevant or reversed annotations.
We separately test learning, fixed-prior transfer, and cross-region sharing.
Held-out chromosome/locus-family training is distinct from transductive prior
learning that sees every test-region summary. Both can be legitimate tracks,
provided competitors receive the same allowed information.

PolyFun estimates priors using genome-wide annotation information; its
published strategy is not equivalent to every local nonnegative annotation
regression followed by SuSiE. [Weissbrod et al., 2020](https://doi.org/10.1038/s41588-020-00735-5).
The repaired local surrogate uses a coherent LDSC baseline, constrained fitting
checks and held-out chromosome calculations, but retains a surrogate label.
Native PolyFun comparisons require its actual implementation and appropriate
genome-wide inputs.

Native SBayesRC likewise needs enough genome-wide evidence to estimate a
shared architecture and functional prior. The 2026 GWFM paper evaluates a
genome-wide mixture model and includes improved mixing and multichain
convergence assessment. Its advantage cannot be established from a small
identity-LD integration fixture. [Wu et al., 2026](https://doi.org/10.1038/s41588-026-02549-3).
Native binary, command mode, LD eigen representation, variant alignment,
effect scale, iteration/burn-in settings and diagnostics are archived.
Successful exit is one necessary condition for a valid comparison, not proof
of convergence.

The expanded native comparison plans five independent genotype draws per
architecture, each with 128 EUR source blocks, 200 variants per block and 5,000
individuals. One phenotype is shared across the blocks. Sparse and
sparse-plus-infinitesimal architectures are crossed with polygenicity .001/.01,
heritability .1/.3/.5 and three phenotype replicates. Causal selection uses
three informative binary annotations with fivefold enrichment. For each
phenotype, the same observed statistics are supplied with no categories,
informative categories or jointly permuted annotation rows within each block.
The design therefore contains 540 inputs and 2,160 native chains. Regional
competitors receive matching summary statistics and annotation exposures.

Each input starts with four independent chains of 25,000 iterations, 5,000
burn-in and thinning 10. A real release-binary check established that its raw
`CoreParameters` file writes a cached row on every iteration, without iteration
labels, including burn-in and repeated values between retained draws. The
parser is bound to that executable's SHA-256 and observed header. It requires
exactly the planned raw row count and retains only one-based iterations greater
than burn-in and divisible by thinning. This yields 2,000 retained draws per
full-design chain, with explicit iteration numbers and sealed extraction
metadata. Burn-in and stale repeats do not contribute to effective sample size.
The retained core-parameter draws are assessed
using the official `posterior` implementation of rank-normalised/folded split
R-hat and bulk/tail effective sample size. The prespecified thresholds are
R-hat at most 1.01 and both effective sample sizes at least 400. The declared
RC/default-gamma structural constant `Vg1 = 0` must be exactly zero in every
retained draw, is reported separately with missing R-hat/ESS, and is excluded
from sampled-parameter pass criteria. The source code and actual release
log/draws support this one prespecified constant. Other constant or nonfinite
sampled-parameter traces cannot establish convergence. Core-parameter
assessability and threshold failures are distinguished: the source can clamp
residual variance to trait variance, so a constant ancillary trace need not
mean the posterior predictions are wrong, and longer chains may not resolve
that assessability limitation. Raw imported probability-score outcomes remain
available. Ordinary native accuracy summaries and paired contrasts require
every input in the draw/setting to have four coherent, successful chains passing
the strict core criterion; no favourable subset is compared against a regional
method's full setting. Core-parameter
diagnostics do not establish mixing of every unrecorded SNP indicator;
truth-free across-chain PIP variation is reported separately through median,
95th-percentile and maximum SNP ranges, disagreements on 0.95 calls, and ranges
of block-level PIP sums. Equal block sums can conceal chains selecting different
LD-equivalent variants. These descriptions carry no invented convergence or
precision threshold. No convergence claim is made
for short-chain integration or resource canaries. A prespecified optional
extension retains the initial outcomes and applies equally to every input
failing these diagnostics, independently of accuracy.

The audited source commit is a reference for interpreting the writer, not proof
of the compiled release commit. A source-level robust-restart branch raises a
possible mismatch between replacement-run traces and final summaries; this
has not been demonstrated in the release binary. Observed restart/disruption
markers therefore withhold ordinary trace diagnostics and flag coherence as
unverified while preserving the execution record. The native executable is
unmodified. Its seven full-variant text arrays are redirected to `/dev/null`;
raw core, small annotation/block traces and final SNP summaries remain. This
policy saves storage without claiming per-variant trajectory diagnostics or
avoiding native formatting cost. Matched release-binary runs with and without
text output gave byte-identical final SNP summaries on the short validation
input. The all/core/no-text retention check also matched parameter,
parameter-set and enrichment summaries, and all/core retained core draws were
identical. These are output-contract checks, not scientific performance results.
[Observed trace contract](../scripts/native/campaign/trace-contract.json),
[GCTB trace writer](https://github.com/jianzeng/GCTB/blob/cc7fa7d765c83a89c6375946cf77fe50ba1a317e/scr/mcmc.cpp),
[conditional residual-variance update](https://github.com/jianzeng/GCTB/blob/cc7fa7d765c83a89c6375946cf77fe50ba1a317e/scr/model.cpp#L2141-L2155),
[source restart branch](https://github.com/jianzeng/GCTB/blob/cc7fa7d765c83a89c6375946cf77fe50ba1a317e/scr/main.cpp#L255-L276).

Native effects and standard errors are converted to original per-A1-allele
units using empirical raw-dosage standard deviations, with verified ALT
orientation and stable coordinate/allele IDs. The within-block empirical LD is
archived as a checked float32 eigendecomposition. The joint block approximation
omits finite-sample cross-block LD; its 128-region spread is a representative
subset, not dense genome-wide coverage. Intervals resample independent genotype
draws to retain shared-phenotype and learned-prior dependence. Five draws give
limited generalisation precision, so draw-level effects accompany intervals.
[Native campaign design](../scripts/native/campaign/README.md),
[posterior R-hat](https://mc-stan.org/posterior/reference/rhat.html),
[posterior bulk ESS](https://mc-stan.org/posterior/reference/ess_bulk.html).

FINEMAP, CARMA, FiniMOM, PAINTOR, Funmap and SparsePro are compared using
identified upstream software. CARMA's discrepancy modelling is particularly
relevant to external/meta-analysis LD, while SuSiE-inf and FINEMAP-inf are
explicitly developed for single-cohort in-sample LD. Stress results outside
that scope are labelled accordingly. [CARMA](https://doi.org/10.1038/s41588-023-01392-0),
[FiniMOM](https://doi.org/10.1093/bioinformatics/btad396),
[fine-mapping-inf input scope](https://github.com/FinucaneLab/fine-mapping-inf).
Funmap and SparsePro are functional-prior comparators with their own
information and fitting conventions. [Funmap](https://doi.org/10.1093/bioinformatics/btaf017),
[SparsePro software](https://github.com/zhwm/SparsePro).

Native probability support is also part of the comparison. The audited PAINTOR
enumerator advances from the all-zero configuration before evaluating it, so
its PIPs condition on at least one modelled cause. At a two-cause limit its
support contains $p + {p\choose2}$ nonempty configurations. Its distinct MCMC
path includes null in its search, but neither output supplies a locus-null
posterior. Forced-nonnull assumptions remain assessable through operational
Brier/log loss on null tasks; those outcomes show their consequences under the
benchmark's population. They do not establish calibration of a conditional
probability as an unconditional one. ABF's forced-single-cause baseline and
FINEMAP's observed default zero prior on the null need the same explicit
support disclosure. [PAINTOR source enumeration](https://github.com/gkichaev/PAINTOR_V3.0/blob/ffce82a7114eeb98220db39ba5d6b69e5a419a9e/Functions_model.cpp#L124-L151),
[native support audit](../docs/paintor-native-contract.md).

## Tuning, pairing and uncertainty

We distinguish author defaults, competitive tuning and controlled ablations.
Giving an annotation-enabled model a tuned configuration while leaving its
annotation-free counterpart at unrelated defaults does not isolate the prior.
Conversely, sharing settings selected for only one model can disadvantage the
other. Both comparisons can be useful if their interpretation is explicit.
Hyperparameters, training data, study fingerprints and tuning budgets are
frozen before final assessment. Trial failures and selective pruning are
reported, not converted into favourable missing values.

Methods are paired on the same inputs. Uncertainty resamples independent
locus/panel families and preserves shared structures. In the legacy generator,
genotypes are simulated once per job row and reused across all causal-count,
heritability and replicate scenarios. Repeated phenotypes do not create new
empirical loci. Joint prior learning also links regions within a training
panel. Intervals that condition on those LD/training families must be described
as conditional uncertainty; generalisation to new genomic loci requires
additional independent families.

Macro localisation assigns equal weight to loci. Micro discovery precision
assigns weight to calls; pooled probability scores assign weight to variants.
These are different estimands. The final report presents them together with
paired effect differences, intervals, completed-pair counts and attempted
availability. Exploratory multiple contrasts are not presented as family-wise
confirmatory discoveries. No universal winner is computed from arbitrary
equal weighting of the named stress tracks.

# Results: verified evidence and pending comparative assessment

**The full optional/native-method campaign is not yet complete.** This section
must be replaced with its frozen tables, plots and setting-specific comparison
before the manuscript is described as a full results paper. The following
evidence establishes implementation and integration, not the requested final
ranking.

The earlier corrected release completed 964 R test expectations and 21 Python tests;
`R CMD check` returned OK. A 220-task local challenge evaluated six predictors,
and an Imperial development run evaluated four predictors on 1,100 tasks.
The latter spans 50 draws per each of the 22 tracks with 40 observed variants
before track-specific changes. These are small-model development runs.
Unsupported multiple-cohort/variant-wise-sample inputs remain explicit rather
than being converted into standard single-cohort analyses. The validated
implementation revision is `3ac28638a8e4249c2dec2b76eb131c6ff6bf27aa`;
those checks describe that revision. Later campaign likelihood adapters and
execution controls have their own source commitments and validation records.
[Validation record](../docs/validation.md).

A native FINEMAP 1.4.2 integration run found exact PIP agreement between the
wrapper and direct binary output on its synthetic fixture, including a
consistent variant permutation. A later actual-binary parser check completed
two loci whose reported causal-count probabilities summed to 0.9996 and
1.0011. These values are compatible with FINEMAP's three-significant-digit
output precision; the parser checks the corresponding rounding intervals and
complete causal-count coverage without renormalising the imported values.
This verifies an output contract, not posterior calibration.
Native GCTB 2.5.5 standard SBayesRC executed
on a 256-variant four-block synthetic fixture. Its short chain had PIP sum
61.86 against sixteen planted sparse effects. That difference demonstrates
why successful execution is not a calibration claim. The run was not designed
to establish convergence or genome-wide accuracy. Both raw evidence and
provenance are retained. [Native integration evidence](../docs/validation-artifacts/README.md).

The final native-wrapper retention check also compared actual Linux R API
imports with the seven documented text sinks enabled. Imported PIPs, posterior
effect means and effect standard deviations were identical in this short
256-variant run. The underlying same-seed all/core/no-text checks described
above additionally matched native final output bytes. These checks support
the recording policy; they do not establish full-sized runtime or mixing.

The RL extension completed 30 local tests, including actual Gymnasium checks,
Docker source execution and attack fixtures, and persistent OpenEnv WebSocket
sessions. A public-key synthetic development demonstration trained a contextual
bandit for 250 updates and evaluated six policies on 100 distinct validation
loci. The learned policy achieved utility 0.928781, compared with 0.898580 for
the working 0/1/2-effect reference: a paired difference of 0.030201, with a
percentile bootstrap interval [0.015384, 0.048274]. The active marginal-uncertainty
policy achieved 0.926501, and a myopic working-model value-of-information policy
achieved 0.907386. Their intervals against the working reference included zero.
These utilities include acquisition costs. They demonstrate the value of
testing decisions under a declared stress distribution; they do not establish
a blinded competition win, biological assay performance or LLM post-training.

An actual Docker code-improvement example selected an independently written
ABF predictor using six development loci, reducing development loss from
0.149947 to 0.112943 relative to a forced-one starter. A distinct six-locus
held-out panel gave utility 0.881110 after two query costs, with no failed tasks.
This is evidence of executable development feedback and held-out scoring,
with labels outside the predictor container. It is a small worked example,
not an estimate of trained coding-agent capability.
[RL validation and frozen artifacts](../rl_envs/docs/validation.md).

The final results need to answer the historical claims individually. For
example, if corrected AP changes the lead on the original continuous,
in-sample inputs, the paper should show the paired change and whether it
arises from the evaluator, inference or both. If a new annotation-reversal
track penalises a functional prior, that demonstrates transfer sensitivity,
not a direct contradiction of a useful-annotation result. If a method
performs well in matching LD but degrades with reference mismatch, the
comparison should specify what information was lost and whether the model's
intended scope includes that input.

Setting-specific guidance will be written from completed evidence in this
form: which method or strategy leads on localisation; which controls
high-confidence false discoveries at the chosen threshold; whether credible
uncertainty is useful and assessed under the correct event; whether methods
are distinguishable given independent clusters; and what compute or
availability cost accompanies the result. No winner is assigned here before
those results exist.

# From development benchmark to prospective assessment

A renewable benchmark needs new evidence that participants could not already
have used. CASP's useful lesson is independent target provision and blind
assessment, rather than merely a large number of targets. Experimental
structures remain confidential through the modelling deadline, and independent
assessors evaluate frozen predictions. [CASP17 target policy](https://predictioncenter.org/casp17/callForTargets.cgi).

A fine-mapping counterpart should begin with a tractable, assay-defined
endpoint. For example, context-matched endogenous single-variant editing can
assess an allele's molecular effect alongside appropriate association data.
That endpoint is stronger than a generic reporter readout, but still does
not automatically prove mediation of a human disease phenotype. Trait,
tissue, condition, ancestry, candidate universe, experimental power and
measurement uncertainty must be recorded before predictions.

Providers should select targets under a declared policy that does not depend
on the competing fine-mappers. Where an entire candidate universe cannot be
tested, a preregistered probability sample with known inclusion probabilities
permits a justified partial-label analysis. Untested variants are unknown;
an assay-negative call can be an adequately powered negative only for the
specified endpoint. Reporter-negative is not equivalent to disease-noncausal.

The synthetic implementation separates public predictor inputs from private
assessor truth, shuffles task mechanisms before opaque IDs are assigned,
validates output IDs and probabilities, and records file hashes. These features
reduce accidental leakage and support reproducibility. They are not a
security sandbox or independent custody. An actual community round requires
access-controlled assessor storage, isolated untrusted submissions, independent
commitment of manifests, a frozen deadline and post-assessment evidence release.

The first prospective round should retain several input tracks: same-input
regional inference, annotation-informed strategies, multiple-cohort analyses,
and genotype/haplotype methods where marginal summaries are insufficient.
Genome-wide prior learning should be assessed as a strategy with its declared
information advantage, not hidden inside a nominally same-input locus contest.
The round should publish every failure and withdrawal, complete artifacts,
uncertainty intervals, and assessor commentary on methods' useful operating
domains.

# Agent-training environments and frontier-lab use

Fine-mapping supplies several complementary decision problems for agents.
The accompanying implementation provides four environments: `InferenceEnv`
for probability forecasts, `ExperimentEnv` for budgeted evidence acquisition,
`MethodDesignEnv` for choosing and validating an inference strategy, and
`CodeImprovementEnv` for evaluating submitted implementation changes. Their
interfaces have passed local execution and isolation tests, with the measured
development demonstration reported above. Independently administered final
evaluation remains a separate requirement.
A further diagnostic task could ask an agent to identify allele or LD
inconsistency and request missing covariance information. These are different
action spaces and should not share an ill-defined reward.

Rewards must respect the probability event and information boundary. A proper
score rewards sparse-membership forecasts and the whole-locus null event:

$$r=1-\left[\frac{1}{2p}\sum_j(p_j-\gamma_j)^2+
\frac12(p_{\mathrm{null}}-\mathbf1\{C=\varnothing\})^2\right].$$

The forecast contract requires exact variant alignment and the necessary
event bound $\max_jp_j\le1-p_{\mathrm{null}}$. Here the null event concerns
the whole generating locus, while PIPs concern observed candidates. If a
causal candidate is missing, all observed PIPs can be zero with null
probability zero. The lower union bound based on the sum of observed PIPs
therefore cannot be imposed without a public completeness guarantee.
Malformed actions receive an explicit failure penalty rather than
a counterfeit all-zero forecast. The experimental environment supplies a
stylised measurement channel with sensitivity 0.85, specificity 0.90, a
six-query budget, at most two repeated measurements per variant, and a cost
of 0.003 per acquisition. Its measurements concern the synthetic sparse label,
not an assertion that a real MPRA or editing assay has those operating
characteristics. The method-design environment has four development
validation queries, each costing 0.002. Explicit operational costs account
for measurement requests and validation.
AP alone is invariant to many probability distortions and is unsuitable as
the sole reward for calibrated inference. A hidden assessor label must not
appear in an observation, tool file, seed, mechanism-bearing task ID or reusable
final-data directory. Failed and malformed actions remain failures instead of
receiving invented posterior probabilities.

Development episodes can expose labels after completion for learning.
Evaluation episodes require separately held truth and an exposure record;
different seeds alone do not guarantee secrecy. Holdouts should include
new locus families and mechanisms, with enough baseline agents and tested
reward behaviour to reveal simple reward hacking. A credible package must
specify environment APIs, reset/step semantics, observation and action schemas,
termination rules, reference policies, resource limits, scoring tests,
provenance and licensing.

The implemented assessor uses private keyed development/final locus streams,
opaque IDs and randomised observation order, with bounded issued evaluation
episodes. The local Python/Gym interface is for trusted development. The
separate JSON-lines assessment interface and code-submission Docker path are
the intended boundary for untrusted policies: public JSON is supplied to a
nonroot, read-only, network-disabled container, while truth stays in the host
assessor. This distinction is part of the package contract, and isolation
claims require actual container execution and adversarial tests before
delivery. A future independently administered assessment remains distinct
from a locally held private key.

Such environments could provide frontier laboratories with tasks that demand
numerical reasoning, tool use, coherent uncertainty and experimental decisions.
That is a proposed use, not an observed sales channel. Readiness to offer a
package requires clean rights to its code and input data, reproducible held-out
evaluation, documented safety/isolation boundaries, and an honest account of
what the reward measures. Any future biological validation or commercial
agreement would be additional evidence, not inferred from passing an API test.

# Discussion

The principal lesson is that fine-mapping performance is a relationship
between a method, a target event and an information set. A ranking can change
because a better algorithm uses the same information, because a strategy
learns a prior from more data, because the data construction changes, or
because the evaluator stops rewarding an artefact. These explanations should
remain visible rather than be collapsed into one leaderboard.

The original framework contributes useful machinery and a thoughtful line of
functional-prior experiments. Its reports already take several limitations
seriously. The repaired suite preserves that contribution while making the
remaining assumptions and failure denominators testable. A fair comparison
can revise a historical conclusion without diminishing the value of the
infrastructure that made the question inspectable.

Several benchmark design choices are necessarily scientific choices. Sparse
membership is useful for comparing sparse-plus-background models, but it is
not a complete biological definition of causality. Fixed-effect stress tests
probe conditional behaviour, but do not provide the same calibration claim as
a matched prior-predictive experiment. Perfect-LD labels test honest ambiguity,
but successful localisation within an indistinguishable class does not prove
that a particular member is causal. These distinctions should guide the
interpretation of final method-specific results.

There are also clear limits. The Gaussian challenge is small and mechanistic;
it is not a substitute for genome-wide genotype realism, variant-wise imputation
and sample masks, real ancestry transfer, multiple traits, nonadditive effects,
long-range LD or experimental truth. Native MCMC methods require sufficient
iterations, independent chains and suitable architecture-learning inputs.
Tuning can consume or leak assessment information if training families and
budgets are not recorded. The final campaign must quantify these limits where
it has evidence and leave unsupported claims unasserted.

The useful endpoint is practical guidance: when does a method produce
recoverable targets with trustworthy probabilities at a tolerable cost, and
what should an analyst do when its assumptions fail? A future blinded round
and agent environments can extend that guidance with renewable evidence and
decision-based assessment. They should inherit the same discipline about
events, inputs, provenance and uncertainty.

# Data, code and reproducibility

A source-bound local orchestration-only check constructed the schema-3
649,220-method-task and 240-generation-unit scheduler projection with about
32 MiB peak process memory. Initialisation took 202 seconds with Python
allocation tracing enabled; resume verification of its frozen rows took
22 seconds. Adding 100,000 synthetic terminal job records and events left
the live-state JSON at 1,280 bytes. With no exhausted units, the indexed
retry-blocker check and empty-report publication took 2.4 milliseconds.
These measurements test disk-backed scheduling feasibility on that local
filesystem. They do not measure inference throughput, time a full exhausted-
dependency report or establish the shared server filesystem's locking and sync
behaviour. The record preserves the earlier schema-2 source-bound probe;
its timings are not a controlled before/after speed comparison.
[Recorded scheduler scale check](../docs/validation-artifacts/campaign-hardening-20261010/sqlite-controller-scale.json).

Interrupted sequential batches require separate submission and task-start records.
The repaired controller syncs a per-job ordered start-intent ledger before each
R command and counts each intent once transactionally; untouched suffix tasks
retain their retry budgets. Exhausted units and dependent comparisons remain
explicitly incomplete. Cross-node ownership also needs a storage-specific
contract: a deployment probe allowed two nodes to acquire the same `flock`.
The controller therefore claims an atomic ownership directory before opening
scheduler state and preserves unknown stale claims for manual reconciliation.
These are execution-integrity repairs, not changes to a fine-mapping model or
its scientific score. [Controller execution contract](../scripts/campaign/README.md).

Original software: [Lucas Siu's repository at the audited revision](https://github.com/lucassiu14/fine-mapping-benchmark/tree/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81).
Corrected implementation: [Nathan Skene's fork](https://github.com/NathanSkene/fine-mapping-benchmark),
branch `codex/challenge-benchmark-suite`, [draft comparison PR](https://github.com/NathanSkene/fine-mapping-benchmark/pull/1).
The [critical review](../docs/REVIEW.md), [scientific contract](../docs/challenge-science.md),
[execution protocol](../docs/challenge-protocol.md), [original-results audit](../docs/original-results-audit.md)
and [validation record](../docs/validation.md) define the present evidence.

The final paper will cite frozen campaign manifests, prediction and evaluation
archives, method/environment hashes, PBS exit evidence, resource and convergence
records, analysis scripts, tables and rendered figures. Historical artifact
availability and reproduction differences will be disclosed. The current
draft does not claim that these full-campaign artifacts already exist.
The runtime guard records package versions and hashes configured external
source, executable and interpreter files. It does not hash every installed R
or Python distribution file; same-version package changes require an immutable
installation environment and its retained build records. Derived analysis
manifests seal the actual output tables and disclose incomplete collections.

# Contribution and disclosure statement

Lucas Siu is credited for the original benchmarking framework and the
Functional BEATRICE/shared-prior code present at the audited revision.
The reassessment implements validity repairs, expanded challenge mechanisms,
provenance and assessment interfaces in a fork for Nathan Skene. The vault's
programme records indicate collaborative supervision and method development.
Specific publication authorship and contribution roles require human
agreement before submission. Automated coding and writing assistance was used
in this review and build; claims require independent scientific review before
publication.
