# Fine-Mapping Methods

This document describes the registered fine-mapping wrappers, their scientific
interpretation, and the arguments available for each. Registration is not evidence
that an upstream dependency was installed or its complete workflow validated.

---

## Available Methods

| Method | Key | Requires | Annotations | Model assumption |
|---|---|---|---|---|
| SuSiE | `"susie"` | `susieR` (CRAN) | No | Multiple causal variants (sparse) |
| SuSiE-inf | `"susie_inf"` | `susieR >= 0.15.0` | No | Sparse + infinitesimal background |
| FINEMAP | `"finemap"` | External binary (auto-downloadable) | No | Multiple causal variants (stochastic search) |
| ABF | `"abf"` | Nothing | No | Single causal variant |
| Funmap | `"funmap"` | Python + `funmap` package (GitHub) | **Required** | Multiple causal variants + functional priors |
| PAINTOR | `"paintor"` | External binary (conda or source) | Optional | Multiple causal variants + annotation enrichment |
| BEATRICE | `"beatrice"` | Python + PyTorch (GitHub) | No | Variational Bayes with neural network PIP inference |
| CARMA | `"carma"` | `CARMA` R package (GitHub) | No | Spike-slab EM with LD-discrepancy outlier detection |
| marginal_z | `"marginal_z"` | Nothing | No | Model-free baseline (pip = \|z\| / sum\|z\|) |
| Enrichment oracle | `"polyfun_oracle"` | `susieR` (CRAN) | **Reads truth** | SuSiE with simulator enrichment weights; diagnostic comparator |
| Annotation regression surrogate | `"polyfun_est"` | `susieR` (CRAN) | Optional | SuSiE with empirical chi-squared/annotation regression priors |
| SparsePro | `"sparsepro"` | Python + cloned SparsePro repo | No (not used here) | Variational fine-mapping with a sparse-projection prior on effect groups |
| Functional BEATRICE | `"functional_beatrice"` | Python + PyTorch + source checkout | Required | BEATRICE with learned annotation priors |
| LD-score prior surrogate | `"polyfun_ldsc"` | `susieR` | Optional | NNLS LD-score prior with chromosome/region holdout |
| SBayesRC-style surrogate | `"sbayesrc"` | `nnet` | Optional | Experimental spike/slab Gibbs and empirical annotation refits |
| CAVIAR | `"caviar"` | External C++ binary | No | Configuration enumeration under an MVN summary likelihood |
| FiniMOM | `"finimom"` | `finimom` R package | No | Nonlocal effect prior and model-size prior |
| FINEMAP-inf | `"finemap_inf"` | Python + fine-mapping-inf checkout | No | FINEMAP sparse and infinitesimal extension |
| Pooled FB experiment | `"fb_pooled"` | Python + joint trainer | Required | Scenario-wide learned logistic annotation head |
| Cross-region FB experiment | `"fb_xregion"` | Python + joint trainer | Required | Scenario-wide learned LassoNet annotation head |


---


## Interpretation and comparison rules

`polyfun_oracle`, `polyfun_est`, `polyfun_ldsc`, and `sbayesrc` are local
comparators, not validated reproductions of upstream PolyFun or SBayesRC.
The legacy `sbayesrc` key is retained for compatibility; outputs explicitly
record `implementation = "sbayesrc_style_surrogate"` and
`upstream_equivalent = FALSE`. Running it across many regions does not make
it the official genome-wide SBayesRC analysis. It uses a fixed slab grid,
short single-chain Gibbs runs, regularized empirical annotation regression,
and a pilot/final-fit pipeline. No convergence or effective-sample-size
claim follows from the default 300 iterations.

The SBayesRC surrogate now uses actual sampled mixture assignments in its
pilot rather than assigning slab classes by ranks of posterior mean effects.
Pilots hold a sparse annotation-neutral prior fixed, and fit the shared
annotation mapping once after pooling draws. Repeated local annotation refits
in the pilots created a feedback loop that inflated inclusion at small loci.
For the simulator's OLS t statistics it uses
`beta_hat = t / sqrt(n - 2 + t^2)` on a standardized phenotype and
likelihood precision `(n - 1) / residual_variance`; `statistic_type =
"normal_score"` selects the approximate Gaussian score transform instead.
The default residual variance of 1 is still an approximation. Without
annotations the mixture weights remain at their fixed initialization; with
annotations, empirical refits are not full posterior sampling of the prior.
Its pooled fit learns from the scenario being analyzed and is a transductive
experiment, not an independently trained annotation predictor.

`polyfun_ldsc` regresses `(chi_squared - 1) / N` on the all-SNP LD score
plus annotation LD scores. The all-SNP genetic baseline is distinct from
the fixed null intercept of 1. It does not put the chi-squared regression
intercept into a genetic prior. NNLS fits must satisfy a projected-gradient
convergence check. This remains a restricted comparator: coefficients are
nonnegative, LD-score weights are simplified, the confounding intercept is
fixed, and it cannot represent arbitrary annotation depletion. Redundant
annotations can make coefficients nonidentifiable even when fitted values
are stable. When chromosome metadata exist, prior training excludes all
regions on the target chromosome; otherwise it excludes the target region.
`additional$prior_training_unit` records that distinction. A single training
group has no holdout and falls back to local estimation.

`run_methods()` removes phenotype truth before observed-data adapters and
scenario setup hooks; the explicitly labeled enrichment oracle receives it.
This is a guard against accidental use of labels, not a secure execution
sandbox. Region wrappers still receive individual genotypes and observed
phenotypes; a prospective challenge must separately enforce its declared
input track. `L`/maximum-causal settings should be prespecified or tuned on
training data, rather than set from a held-out region's true causal count.
Scenario setup errors now fail the affected fits instead of silently changing
an estimator, and malformed PIPs/credible-set indices are recorded as failures.

Normalized absolute-z scores support ranking comparisons, not PIP calibration.
PAINTOR and the SBayesRC surrogate do not return native posterior credible
sets through these wrappers. Their descriptive PIP/score mass selections live
in `additional$pip_mass_set`, with `credible_sets = list()` and
`credible_set_semantics = "not_provided"`. Exclude unsupported credible-set
metrics rather than treating their absence as failed statistical coverage.

## Workflow: from simulation to results

### Step 1 — Simulate data

```r
library(fmbenchmark)        # or, in a source checkout: pkgload::load_all()

sim <- run_simulation(
  n_regions = 3,
  n         = 500,
  p         = 200,
  S         = c(1, 2, 3),
  phi       = c(0.1, 0.2, 0.4),
  n_iter    = 5,
  seed      = 42
)
```

`sim` contains `sim$genotypes`, `sim$scenarios`, and `sim$params`.

---

### Step 2 — Run methods

All method wrappers and the orchestrator are available once the package is
loaded (Step 1). Run one or more methods with `run_methods()`:

```r
out <- run_methods(
  simulation  = sim,
  methods     = c("susie", "abf", "paintor"),
  method_args = list(
    susie   = list(L = 5, coverage = 0.9),
    abf     = list(prior_variance = 0.04),
    paintor = list(max_causal = 2)
  ),
  save        = FALSE,
  output_dir  = "results",
  verbose     = TRUE
)
```

**`methods`** — character vector of one or more method keys (case-insensitive).

**`method_args`** — named list of argument overrides, one entry per method.
Any argument not listed here uses the method's default value. Methods not
mentioned in `method_args` use all defaults.

**`save`** — if `TRUE`, each method's results are written to
`output_dir/<sim_label>/<method>.rds`, plus a `run_metadata.rds` file.

---

### Step 3 — Access results

Results are returned as a named list, one element per method:

```r
# PIPs for the first scenario × region combination
out$susie$results[[1]]$pip
out$abf$results[[1]]$pip

# Credible sets (list of integer index vectors)
out$susie$results[[1]]$credible_sets

# Method-specific additional outputs
out$susie$results[[1]]$additional$alpha      # L x p alpha matrix
out$susie$results[[1]]$additional$converged  # did IBSS converge?
out$abf$results[[1]]$additional$log10_abf   # per-variant log10 ABF

# Scenario metadata attached to each result
out$susie$results[[1]]$scenario_id
out$susie$results[[1]]$S
out$susie$results[[1]]$phi
out$susie$results[[1]]$iter
out$susie$results[[1]]$region_id
```

Results are indexed as a flat list over all `scenarios × regions`. For a
simulation with `n_scenarios` scenarios and `n_regions` regions, index `i`
corresponds to scenario `ceiling(i / n_regions)` and region
`((i - 1) %% n_regions) + 1`. Each result also carries `$scenario_id` and
`$region_id` fields for unambiguous lookup.

---

## Standard output format

Every method returns the same top-level fields:

| Field | Type | Description |
|---|---|---|
| `pip` | numeric vector (length p) | Posterior inclusion probabilities for probabilistic methods; legacy ranking-score field for `marginal_z` |
| `credible_sets` | list of integer vectors | Credible sets as variant indices. Empty list if none reported. |
| `method` | character | Method name |
| `input_type` | character | `"summary"` or `"individual"` |
| `params` | list | All hyperparameters actually used |
| `runtime_seconds` | numeric | Wall-clock time for this fit |
| `additional` | list | Method-specific outputs (see below) |
| `error` | character or NULL | Error message if the fit failed; otherwise absent |

Scenario metadata (`scenario_id`, `region_id`, `S`, `phi`, `p_causal`,
`iter`) is attached to each result by `run_methods()`.

---

## Method-specific details and arguments

---

### SuSiE

Sum of Single Effects regression. Fits L independent single-effect components
and returns one credible set per component (after a purity filter).

**Setup:** `setup_susie()` — installs `susieR` from CRAN if needed.

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `L` | `10` | Maximum number of causal components |
| `estimate_residual_variance` | `TRUE` | Estimate residual variance |
| `estimate_prior_variance` | `TRUE` | Estimate per-component prior variance |
| `prior_variance` | `NULL` | Prior effect variance (auto-set if NULL: `0.1 * var(y)` for individual-level, `0.1` for summary stats) |
| `coverage` | `0.95` | Credible set coverage |
| `min_abs_corr` | `0.5` | Minimum purity (min absolute correlation within a CS) |
| `max_iter` | `100` | Maximum IBSS iterations |
| `use_individual` | `FALSE` | Use individual-level data (X, y) instead of summary stats |

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `alpha` | L × p matrix of per-component posterior assignment probabilities |
| `posterior_mean` | Posterior mean effect size per variant (summed across components) |
| `lbf` | Log Bayes factor per component (length L) |
| `cs_purity` | Data frame of min/mean/median absolute correlation per CS |
| `converged` | Logical: did IBSS converge? |
| `elbo` | Final ELBO value |
| `n_iter_run` | Number of IBSS iterations run |

---

### SuSiE-inf

Extends SuSiE by adding a genome-wide infinitesimal (polygenic) component.
More robust when there is polygenic background signal within the region.

**Setup:** `setup_susie_inf()` — checks `susieR >= 0.15.0`.
If your version is older, it will print upgrade instructions.

> **Note:** `susieR >= 0.15.0` is required. The lockfile currently tracks
> v0.15.57 (installed from GitHub at `stephenslab/susieR`), which satisfies
> this requirement.

**Arguments:** identical to SuSiE (see above).

**Additional outputs (`$additional`):**

Same as SuSiE, plus:

| Field | Description |
|---|---|
| `tau2` | Estimated infinitesimal variance component |
| `theta` | Per-variant posterior means of the infinitesimal effects (BLUP), or NULL |

---

### FINEMAP

Stochastic shotgun search over causal configurations. Allows multiple causal
variants and reports credible sets for the most probable number of causal
variants k*.

The adapter reads model-size probabilities from the complete native log,
because the configuration table contains only the top retained configurations.
It reads only the `cred1`, `cred2`, … identifier columns as sets, retaining their
paired conditional probabilities separately. Invalid IDs, probabilities or
nonzero native exit status are explicit failures.
FINEMAP's sets condition on a selected causal configuration; see its
[official output specification](https://www.christianbenner.com/).

**Setup:** `setup_finemap()` — checks for the binary. If not found,
downloads it automatically (macOS and Linux only):

```r
fp <- setup_finemap()   # downloads if needed, returns path
```

Then pass the path when running:

```r
run_methods(sim, methods = "finemap",
            method_args = list(finemap = list(finemap_path = fp)))
```

If you already have FINEMAP on your PATH, `finemap_path` can be omitted.

To disable auto-download and see manual installation instructions:

```r
setup_finemap(download = FALSE)
```

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `finemap_path` | `"finemap"` | Path to binary, or name if on PATH |
| `n_causal` | `5` | Maximum number of causal variants |
| `n_iter` | `100000` | Number of stochastic shotgun search iterations |
| `prior_std` | `0.05` | Prior standard deviation on effect sizes |
| `coverage` | `0.95` | Credible set coverage |

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `log10bf` | Log10 Bayes factor per variant |
| `posterior_mean` | Posterior mean effect size per variant |
| `posterior_sd` | Posterior SD of effect size per variant |
| `k_posterior` | Named numeric vector: posterior probability that exactly k variants are causal |
| `best_k` | Most probable number of causal variants |
| `configs` | Data frame of top causal configurations and their posterior probabilities |

---

### ABF

Wakefield (2009) Approximate Bayes Factor. Computes a Bayes factor for each
variant under a single-causal-variant assumption and normalises to PIPs.
Returns one credible set.

No setup or installation needed.

> **Limitation:** ABF assumes exactly one causal variant per region. PIPs
> will be miscalibrated in scenarios where S > 1. It is best used as a
> simple baseline.

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `prior_variance` | `0.04` | Prior variance W on effect sizes under H1 (prior SD = 0.2) |
| `coverage` | `0.95` | Credible set coverage |

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `log10_abf` | Log10 approximate Bayes factor per variant |

### marginal_z

Model-free ranking baseline. The legacy `pip` field contains normalized
absolute marginal z scores, `|z_j| / sum(|z|)`, with a uniform fallback when
all scores are zero. These are not posterior inclusion probabilities.
`credible_sets` is empty; a descriptive score-mass selection using `coverage`
is returned in `additional$pip_mass_set`. `additional$is_posterior_probability`
is `FALSE`. The ranking is a useful comparator, but no method has a universal
finite-sample obligation to beat it under every architecture.

No installation or external dependencies are required. `coverage` defaults
to 0.95 and controls only the descriptive mass selection.

### polyfun_oracle

Information-advantaged SuSiE comparator using the simulator's annotation
enrichment coefficients. It forms normalized `exp(A %*% log(enrichment))`
weights and passes them to `susie_rss(prior_weights = ...)`. It reports
`oracle_access = "simulator_enrichment"` and
`strict_performance_ceiling = FALSE`.

These weights are not exact marginal inclusion probabilities under weighted
sampling without replacement. They also omit the genome-wide generator's
probability saturation and global calibration. Consequently this is an
annotation-enrichment oracle, not an exact posterior under the full simulator.
Approximate inference, model mismatch, LD mismatch, and the choice of loss
or metric can let another method outperform it. Bayes optimality for an
expected loss does not establish a pointwise ceiling on observed AUPRC.
Separate it from observed-data leaderboards; gaps mix prior information and
inference effects and do not by themselves measure room for improvement.

**Arguments:** all standard `susie_rss` arguments are forwarded — see
SuSiE section above. The wrapper itself takes no extra parameters.

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `prior_weights` | Normalized simulator enrichment weights, not exact marginal inclusion probabilities |
| `alpha`, `posterior_mean`, `lbf`, `converged`, `elbo`, `n_iter_run` | Forwarded from the underlying SuSiE fit (same as `susie`) |

### polyfun_est

Experimental annotation-regression comparator. It fits OLS of observed
chi-squared statistics on each variant's own annotations, clamps negative
coefficients, and converts fitted values to normalized SuSiE prior weights.
The scenario hook pools the same regions that are then analyzed; standalone
calls fit on the target region. The output identifies this as an
`annotation_regression_susie_surrogate` and records its training scope.

This is not S-LDSC: tagging effects are omitted, a null/confounding intercept
enters the weights, and clipping OLS coefficients does not solve constrained
least squares. It remains as a simplified comparator whose shortcomings can
be tested against `polyfun_ldsc`. Informative annotations alone do not
guarantee improved ranking or calibration after this empirical fit.
The upstream PolyFun pipeline requires its actual inputs and workflow; those
requirements do not make LD correction unnecessary on simulated data.

**Arguments:** all standard `susie_rss` arguments are forwarded — see
SuSiE section above. The wrapper itself takes no extra parameters.

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `prior_weights` | Per-SNP prior weights derived from the estimated `tau` |
| `tau` | Estimated per-annotation contributions (intercept + one entry per annotation column) |
| `alpha`, `posterior_mean`, `lbf`, `converged`, `elbo`, `n_iter_run` | Forwarded from the underlying SuSiE fit |

### Funmap

Annotation-informed fine-mapping that integrates functional annotations via a
random effects model. Uses a three-stage pipeline: baseline SuSiE, annotation
initialisation, then refinement with learned annotation weights. Addresses
overfitting when integrating high-dimensional annotations.

**Funmap requires functional annotations.** If no annotations were simulated
for a region, the method is skipped with an informative message. To produce
annotations, run `run_simulation()` with `annotations = "binary"` or
`annotations = "continuous"`.

**Funmap is a Python package** called from R via `reticulate`. See setup below.

**Setup:** `setup_funmap()` — checks reticulate, Python, and the funmap module
in one step. Install the funmap Python package with:

```bash
git clone https://github.com/LeeHITsz/Funmap.git
cd Funmap
pip install -r requirements.txt
pip install .
```

Then verify from R:

```r
setup_funmap()
```

If Python is not yet configured, `setup_funmap()` will print options for
pointing reticulate to the right environment.

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `L` | `10` | Maximum number of causal components |
| `max_iter` | `100` | Maximum iterations for convergence |
| `tol` | `5e-5` | Convergence tolerance |

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `alpha` | L × p matrix of per-component posterior inclusion indicators |
| `posterior_mean` | Posterior mean effect size per variant (summed across components) |
| `sigma2` | Estimated residual variance |
| `converged` | Logical |
| `cs_purity` | Data frame of min/mean/median absolute correlations per CS, or NULL |

**Note:** Funmap credible set indices are converted from Python's 0-based
indexing to R's 1-based indexing automatically.

---

### PAINTOR

Probabilistic Annotation INTegraTOR. Uses an EM algorithm to jointly estimate
annotation enrichment weights and posterior inclusion probabilities. Supports
multiple causal variants and optionally leverages functional annotations.

Unlike Funmap, annotations are optional — without them PAINTOR runs with a
uniform prior (intercept only), which is still a valid use case.

**PAINTOR is an external C++ binary.** It must be compiled from source or
installed via conda. No auto-download is provided.

**Setup:** `setup_paintor()` — checks for the binary and prints install
instructions if not found.

**macOS (Apple Silicon / arm64):** No pre-built conda package exists for
arm64. Compile from source with the patches below.

```bash
# --- Option 1: conda (Linux / macOS x86_64 only) ---
conda install -c bioconda paintor

# --- Option 2: compile from source (required on macOS arm64) ---

# Prerequisites
brew install eigen   # Eigen3 headers

# Clone and patch
git clone https://github.com/gkichaev/PAINTOR_V3.0
cd PAINTOR_V3.0

# Patch 1: replace the deprecated Eigen/Array stub (Eigen2 compat shim
#          that errors out; just forward to Core)
printf '#ifndef EIGEN_ARRAY_MODULE_H\n#define EIGEN_ARRAY_MODULE_H\n#include "Core"\n#endif\n' \
  > eigen/Eigen/Array

# Patch 2: update Makefile to use system Eigen3 and C++14
sed -i '' 's/-std=c++11/-std=c++14/g' Makefile
sed -i '' \
  's|-I/$(curr)/eigen/Eigen|-I/$(curr)/eigen_shim -I/opt/homebrew/include/eigen3|' \
  Makefile

# Patch 3: fix float-as-array-index (Eigen3 no longer allows this)
sed -i '' \
  's/causal_config_bit_vector(causal_index\[i\]-1) = 1/causal_config_bit_vector((int)(causal_index[i]-1)) = 1/' \
  Functions_model.cpp

# Create the shim directory so #include <Eigen> resolves to Eigen/Dense
mkdir -p eigen_shim
printf '#include <Eigen/Dense>\n' > eigen_shim/Eigen

# Build nlopt (bundled), then PAINTOR
bash install.sh

# Copy binary somewhere permanent
mkdir -p ~/tools && cp PAINTOR ~/tools/PAINTOR
```

Verify from R:

```r
setup_paintor("~/tools/PAINTOR")
```

Then run:

```r
run_methods(sim, methods = "paintor",
            method_args = list(paintor = list(paintor_path = "~/tools/PAINTOR",
                                              max_causal = 2)))
```

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `paintor_path` | `"PAINTOR"` | Path to binary, or name if on PATH |
| `max_causal` | `2` | Maximum number of causal variants (enumeration mode) |
| `mcmc` | `FALSE` | Use MCMC search instead of exact enumeration. Recommended for large regions (p > 500) or `max_causal > 3` |
| `coverage` | `0.95` | Descriptive PIP-mass selection threshold; does not establish credible-set coverage |

**Credible-set interpretation:** PAINTOR does not report native credible sets
through this adapter. `credible_sets` is empty. The previous greedy selection
is retained only as `additional$pip_mass_set`; a sum of marginal PIPs in a
multiple-causal model is not a posterior probability for a causal event.

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `annotations_used` | Integer. Number of annotation columns passed to PAINTOR (0 if none) |
| `log_bayes_factor` | Numeric vector (length p) or NULL. Log Bayes factor per variant, if present in PAINTOR output |

---

### BEATRICE

Bayesian fine-mapping using variational inference with a neural network.
BEATRICE trains a 3-layer network to infer posterior inclusion probabilities
via binary concrete (Gumbel-Softmax) sampling, with ABF-based likelihoods
and KL-divergence regularisation. No annotations required.

**BEATRICE is a Python script** called from R via `system2`. Requires PyTorch.

**Setup:** Clone the repository and install PyTorch (and other deps) into
your Python environment. On macOS the existing Anaconda base works if you add
`torch`:

```bash
# Clone
git clone https://github.com/sayangsep/Beatrice-Finemapping ~/Beatrice-Finemapping

# Install PyTorch (if not already present)
pip install torch

# Verify from R
```

```r
setup_beatrice(
  beatrice_dir = "~/Beatrice-Finemapping",
  python       = "/opt/anaconda3/bin/python3"
)
```

`setup_beatrice()` checks that `beatrice.py` exists and that `torch`, `numpy`,
`scipy`, and `pandas` are importable from the specified Python.

Then run:

```r
run_methods(sim, methods = "beatrice",
            method_args = list(beatrice = list(
              beatrice_dir = "~/Beatrice-Finemapping",
              python       = "/opt/anaconda3/bin/python3",
              max_iter     = 2000
            )))
```

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `beatrice_dir` | *(required)* | Path to cloned Beatrice-Finemapping repository |
| `python` | `"python"` | Path to Python executable |
| `max_iter` | `2000` | Training iterations. Minimum 500. Reduce for speed at the cost of noise. |
| `n_caus` | `5` | Expected number of causal variants |
| `sigma_sq` | `0.05` | Prior variance on effect sizes |
| `gamma_coverage` | `0.95` | Coverage threshold for credible sets |
| `sparse_concrete` | `50` | Sparsity parameter (top-K variants sampled per iteration) |

> **Note on runtime:** With `max_iter = 500`, BEATRICE takes ~10 s per region on a laptop CPU. At the default `max_iter = 2000` expect ~30–40 s per region. For large benchmarks consider reducing `max_iter` or running on a machine with a GPU.

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `cs_pip` | List of numeric vectors. Conditional inclusion probability of each variant within its credible set, as output by BEATRICE. NULL if the file was not produced. |

---

### SparsePro

Variational fine-mapping with a sparse-projection prior on effect groups
(Zhang et al. 2023, PLoS Genetics). A modern alternative to SuSiE that uses
a different variational family — useful in the benchmark as a fast modern
comparator that doesn't depend on annotations.

SparsePro is distributed as a Python script (`sparsepro_zld.py`) in its
upstream GitHub repository rather than as an importable Python package, so
this wrapper follows the same pattern as the BEATRICE / FINEMAP / PAINTOR
wrappers: write input files to a temp directory, invoke the script via
`system2()`, parse the output files, and return the standardised result list.

**Setup:** `setup_sparsepro(sparsepro_dir, python)` — verifies that
`sparsepro_zld.py` is present in `sparsepro_dir` and that the Python env has
`numpy`, `scipy`, and `pandas`. The SparsePro repo must be cloned manually
(no auto-clone):

```bash
git clone https://github.com/zhwm/SparsePro
cd SparsePro
pip install -r requirements.txt
```

Then pass the path via `method_args`:

```r
run_methods(sim, methods = "sparsepro",
            method_args = list(sparsepro = list(
              sparsepro_dir = "~/SparsePro",
              python        = "~/anaconda3/envs/beatrice/bin/python"
            )))
```

> **Limitation:** the wrapper currently does not pass functional annotations
> to SparsePro (it supports them upstream via the `--anno` flag). The
> annotation-aware comparators in this benchmark are Funmap, PAINTOR,
> Functional BEATRICE, polyfun_oracle, and polyfun_est; SparsePro's role
> here is the "modern variational baseline without annotations" slot.

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `sparsepro_dir` | (required) | Path to the cloned SparsePro repo containing `sparsepro_zld.py` |
| `python` | `"python"` | Python executable (full path if not on PATH) |
| `K` | `5` | Maximum number of effect groups (causal signals) SparsePro will consider |
| `cthres` | `0.95` | Coverage level for credible sets |

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `cs_pip` | List of numeric vectors. Per-variant inclusion probabilities within each credible set (the `pip` column of `.cs`). NULL if not parseable. |
| `cs_effect_size` | List of numeric vectors. Posterior effect sizes per variant per credible set (the `effect_size` column of `.cs`). NULL if not parseable. |
| `run_log` | Character vector. Captured stdout/stderr from the SparsePro invocation; useful for debugging unexpected failures. |

---

### CARMA

Contextual Adaptive Robust Marginal Analysis (Yang et al. 2023, Nature
Genetics). A spike-slab Bayesian fine-mapping method that jointly fits
causal configurations via EM, with an outlier-detection component that
identifies variants where the z-score is inconsistent with the LD reference
panel. This makes CARMA more robust than SuSiE or FINEMAP when the GWAS
sample and the LD reference are mismatched.

Unlike SuSiE, CARMA does not decompose signals into independent components.
It returns multiple LD-clustered sets, built using summed marginal PIPs and
`rho.index`. The adapter preserves each native group. Summed PIPs give an
expected causal count, so reaching this threshold does not establish a joint
coverage guarantee. Per-set hits and distinct-signal recovery remain empirical
diagnostics. See the [official construction](https://github.com/ZikunY/CARMA/blob/51f0fcf6a74ee74c75fd2ce70912fdd6cdd214fe/R/CARMA.R#L146-L200).

**Setup:** Pure R package, no external binary or Python needed.

```r
setup_carma()   # installs from GitHub (ZikunY/CARMA) if needed

# Or install manually:
# remotes::install_github("ZikunY/CARMA")
```

**Arguments:**

| Argument | Default | Description |
|---|---|---|
| `rho.index` | `0.95` | Credible set coverage threshold |
| `num.causal` | `10` | Maximum number of causal variants to consider |
| `tau` | `0.04` | Prior variance on effect sizes (prior SD = 0.2, matches ABF default) |
| `effect.size.prior` | `"Spike-slab"` | Prior type: `"Spike-slab"`, `"Cauchy"`, or `"Hyper-g"` |
| `outlier.switch` | `TRUE` | Enable LD-discrepancy outlier detection |
| `all.iter` | `3` | Number of outer EM iterations |

**Additional outputs (`$additional`):**

| Field | Description |
|---|---|
| `outliers` | Data frame of detected outlier variants (z-score inconsistent with LD). Zero rows if none detected or `outlier.switch = FALSE`. |

---

## Passing arguments — quick reference

All method-specific arguments are passed via the `method_args` list in
`run_methods()`. Each entry is named after the method key:

```r
run_methods(
  sim,
  methods     = c("susie", "susie_inf", "finemap", "abf", "funmap", "paintor", "beatrice", "carma"),
  method_args = list(
    susie     = list(L = 10, coverage = 0.95, min_abs_corr = 0.5),
    susie_inf = list(L = 10, coverage = 0.95),
    finemap   = list(n_causal = 3, prior_std = 0.05, finemap_path = "~/tools/finemap"),
    abf       = list(prior_variance = 0.04),
    funmap    = list(L = 10, max_iter = 100),
    paintor   = list(max_causal = 2, paintor_path = "~/tools/PAINTOR"),
    beatrice  = list(beatrice_dir = "~/Beatrice-Finemapping",
                     python = "/opt/anaconda3/bin/python3",
                     max_iter = 2000),
    carma     = list(rho.index = 0.95, num.causal = 10)
  )
)
# Notes:
# - Funmap requires sim to have been run with annotations = "binary" or "continuous"
# - PAINTOR uses annotations automatically if present; falls back to uniform prior if not
# - BEATRICE requires PyTorch; ~10 s/region at max_iter = 500, ~35 s at max_iter = 2000
# - CARMA preserves its separate LD-clustered PIP-mass sets
```

Arguments not listed use the method's default value. Methods not mentioned
in `method_args` at all use all defaults.

---

## Adding a new method

1. Create `R/wrapper_<method>.R` implementing:
   - `run_<method>(...)` — single-region runner returning the standard format
   - `run_<method>_region(region_geno, region_pheno, ...)` — adapter for `run_methods()`
2. Add one line to `.FM_REGISTRY` in `R/run_methods.R`:
   ```r
   .FM_REGISTRY <- list(
     ...,
     mymethod = "run_mymethod_region"
   )
   ```

The standard return format that `run_<method>()` must produce:

```r
list(
  pip             = <numeric vector, length p>,
  credible_sets   = <list of integer vectors, or list()>,
  method          = "<method name>",
  input_type      = "summary" or "individual",
  params          = list(<all hyperparameters used>),
  runtime_seconds = <numeric>,
  additional      = list(<method-specific outputs>)
)
```


## Simulation architecture changes

When annotation `enrichment` is omitted, the local simulator now samples one
fold vector per `run_simulation()` call and holds it fixed across regions
and scenario replicates. `params$enrichment_actual` records the sampled folds.
A direct `simulate_phenotypes()` call similarly shares one sampled vector
across its regions. Previously every region in every replicate silently
redrew its mapping, confounding shared-prior experiments with unlabelled
heterogeneity. This is an intentional change to the random sequence and
architecture; old seeds do not reproduce the old annotated scenarios.
The genome-wide generator already shared the mapping within a scenario.
