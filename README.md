# Fine-Mapping Benchmarking Framework

A framework for testing **which fine-mapping conclusions remain reliable when
LD, genetic architecture, annotations or the available variants change**.
Originally developed by Lucas Siu; this fork adds a scientific validity review,
repairs to simulation and scoring, and an executable challenge suite.

Start with the [critical review](docs/REVIEW.md). Historical results produced
before the repairs require reevaluation; passing new tests does not validate
old method rankings. The [scientific contract](docs/challenge-science.md)
explains exactly what the new experiments measure.

## Two complementary simulation families

| Family | Inputs and purpose | Boundary |
|---|---|---|
| Haplotype/phenotype | sim1000G genotypes, individual phenotypes, marginal OLS summaries; existing regional and shared-phenotype workflows | Population composition and reusable reference loci need explicit provenance |
| Conditional Gaussian challenge | Declared score likelihood, separated predictor inputs and truth, exact small-model reference, adversarial scoring tests | Synthetic experiments; neither empirical genotype realism nor prospective causal truth |

The new generator supports 22 mechanisms: matched-prior calibration, nulls,
one/two/five sparse effects, correlated opposing effects, sparse plus background,
finite and shifted LD references, informative/irrelevant/reversed annotations,
missing causes, identical-LD twins, inconsistent allele flips and consistent
recoding, discovery-selected nulls and weak signals, variant-wise sample sizes
with missing or known overlap, and shared or heterogeneous effects across two
independent cohorts. See the [scenario and input protocol](docs/challenge-protocol.md).

## Run the challenge without genotype downloads

Requirements: R >= 4.1 and `openssl` on PATH for SHA-256 release integrity.
The generator, exact reference and ABF baseline use base/recommended R packages.

```bash
git clone https://github.com/NathanSkene/fine-mapping-benchmark.git
cd fine-mapping-benchmark
git switch codex/challenge-benchmark-suite
Rscript --vanilla scripts/challenge/run_demo.R \
  --replicates 10 --variants 24 --out results/development-v1
```

The script writes separate public/assessor directories, validated predictions,
per-task and per-mechanism CSVs, paired bootstrap contrasts, a PDF report and run
provenance. If susieR is installed it also runs SuSiE and SuSiE-inf wrappers.
Missing support for an input type remains visible in the availability column.
A larger locus count is needed for a precise performance claim; this example
is an executed development diagnostic.

There is **no single overall rank**. Read availability, the target definition,
proper scores and operational recovery together. Methods that return sharp
PIPs on the loci they complete can still fail the tasks that matter.

```r
# Source checkout: does not require loading the full package's genotype imports.
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
panel <- generate_challenge(n_replicates = 10, p = 24)
predictions <- run_challenge(panel$inputs, list(
  exact = challenge_exact,
  forced_one = challenge_forced_one,
  abf = challenge_region_adapter(run_abf_region)
))
report <- evaluate_challenge(panel, predictions)
report$summary
challenge_compare(report, "exact", "forced_one", "brier", "null")
```

`challenge_exact()` enumerates a declared zero/one/two-effect Gaussian model.
It is an exact probabilistic reference in the compatible matched-prior track,
and exact enumeration of a working model in misspecified tracks. It is not a
universal performance ceiling. An entirely missing causal variant remains in
total-target recall; no method can recover an ID it was never supplied.

## Use the existing genotype-based pipeline

Restore the pinned R environment with `renv::restore()`, then load the source
package with `pkgload::load_all()` or install it normally. `sim1000G` and
`hapsim` are archived dependencies. Some optional methods require external
binaries, Python environments or GitHub packages; their requirements are in
[methods.md](docs/methods.md) and [environment.yml](environment.yml).

```r
pkgload::load_all()
sim <- run_simulation(n_regions = 2, n = 500, p = 100,
                      n_iter = 5, S = c(1, 2), phi = c(.05, .2), seed = 42)
results <- run_methods(sim, methods = c("susie", "abf"))
evaluation <- evaluate_methods(sim, results)
evaluation$susie$global$ap         # exact tied-rank average precision
plot_results(evaluation, output_dir = "results/haplotype-example")
```

For meaningful locus generalisation, supply multiple documented reference
VCFs rather than repeatedly using the small bundled example. The genotype
[getting-started vignette](vignettes/getting-started.Rmd),
[genome-wide simulation specification](docs/gw_simulation_documentation.md) and
[HPC instructions](scripts/hpc/README.md) cover those workflows.

`simulate_gwfm_data()` generates one shared phenotype across regions. This
alone does not make every comparator a faithful genome-wide method. The legacy
`"sbayesrc"` key is explicitly an **SBayesRC-style surrogate**, and the local
PolyFun-style variants are not the upstream PolyFun implementation. Compare
upstream software at its intended information scale before making a claim
about those named methods. Oracle-prior methods belong in a labelled diagnostic
track; they must not receive the same competitive rank as deployable methods.
The separate [official GCTB/SBayesRC interface](docs/native-sbayesrc.md) jointly
fits a prepared full input panel and records its declared genomic scale,
allele alignment, upstream version and file commitments. Its small native
smoke test establishes executable integration; a genome-wide accuracy claim
requires a genome-wide experiment.

The [tuning contract](docs/tuning-validity.md) requires balanced locus subsets,
complete scenario grids, attempted macro AP, explicit fit failures and
fingerprinted studies. Historical unfingerprinted studies cannot be reused
as if they came from the repaired experiment.

The bundled Functional BEATRICE Python code now caches pure log Bayes factors
and applies one final learned prior. Historical mixed-prior `res` archives
must be regenerated. Its posterior remains conditional on explored models,
and its selected-index training KL is approximate. Dependency-free tests run
with `python -m unittest discover -s tests/python`; opt-in actual neural tests
use `FMB_NEURAL_SMOKE=1` in an environment with the documented Python packages.

## What the scores mean

- Exact AP measures ranking; macro and pooled AP answer different questions.
  The old threshold-grid `auprc` is retained as a labelled legacy quantity.
- Raw Brier and clipped logarithmic scores evaluate marginal sparse-membership
  forecasts. Binned reliability diagrams are descriptive, not proof of calibration.
- Native credible sets, configuration collections and cumulative-PIP candidate
  lists represent different events. A sum of marginal PIPs is an expected count.
- Null false calls, missing outputs, invalid predictions, convergence diagnostics
  and all-attempt recovery are reported explicitly.
- Shared genotype regions require clustered uncertainty. The new challenge
  supplies independent locus clusters and paired mechanism-specific bootstrap
  comparisons; the legacy phenotype-replicate jackknife is conditional on its
  reference loci.

See [evaluation.md](docs/evaluation.md) for formulas, tie conventions,
denominators, uncertainty and limitations. Sparse-component membership is the
primary target in sparse-plus-background simulations; a non-target variant can
still have a nonzero background effect.

## Test the evaluator and generator

```r
install.packages("testthat", repos = "https://cloud.r-project.org")
```

```bash
Rscript --vanilla scripts/challenge/test_core.R
```

These scientific invariant tests include hand-computed posterior odds and
proper scores, independent prior-predictive calibration, singular LD,
recoding/permutation invariance, missing-target limits, failure denominators,
input/truth separation, artifact tampering and adversarial submissions.
The dedicated CI workflow runs them on Linux and macOS.

With the full pinned environment restored, run the existing package suite:

```r
pkgload::load_all()
Sys.setenv(NOT_CRAN = "true")
testthat::test_dir("tests/testthat")
```

External method tests can skip when their dependencies are absent. Report
those skips; source parsing or mocked wrappers do not establish that a real
binary ran. Validation results and executed examples are recorded in
[the validation record](docs/validation.md).

## Towards a community challenge

The software supports frozen public-input/private-truth releases and distinct
specification, dataset and run identities. A real blind round additionally
needs independent custody, concealed fresh targets, prespecified metrics and
resource budgets, separate development data, and independent assessors.

A plausible 2027 round would combine a simulation calibration division,
real-genotype robustness tasks and a prospective endogenous variant-effect
assay division. Experimental assay effects must keep their own endpoint;
reporter activity, replication and high PIPs are not complete GWAS causal truth.
The [scientific contract](docs/challenge-science.md) specifies this proposal
and the evidence still needed before claiming a CASP-equivalent assessment.

## Reuse and citation

MIT licence: [LICENSE](LICENSE). Preserve Lucas Siu's authorship and cite the
original framework using [CITATION.cff](CITATION.cff). Source code and benchmark
specifications can be improved independently; a new code revision does not
retroactively change the meaning of a frozen dataset or historical run.
