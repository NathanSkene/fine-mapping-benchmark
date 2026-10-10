# Validation record — 10 October 2026

This records executed checks, not a universal method ranking. Original upstream
source was `f053fdf318a3abe2c44eeb0ac3f3c99a03979b81`. The implementation and
development artifacts are frozen separately; run provenance records source
commit and individual R file hashes. The [review](REVIEW.md) and
[scientific contract](challenge-science.md) define interpretation limits.

## Executed correctness checks

- The full source R suite passed **947 expectations**, with zero failures or
  errors, 15 explicitly skipped dependency cases and one warning about the
  legacy single-region SBayesRC-style annotation refit. Four additional native
  scope assertions were added afterwards and passed separately. The release
  package check is recorded below after completion.
- The focused invariant suite checks hand-computed Gaussian Bayes factors,
  prior-predictive calibration, exact tied AP, raw scores, distinct set matching,
  intercept OLS, realised variance/covariance, annotation probability prevalence,
  singular LD, input recoding, known study overlap, unavailable causes, malformed
  predictions, release tampering, failure denominators and tuning metadata.
- **21 Python tests passed**, including real single-locus LassoNet and joint
  two-locus training, credible-set expansion and numerical failure paths.
  Independent Gaussian likelihood and final Bernoulli-prior calculations match
  the archived posterior, checkpoint priors and exported PIPs. Default mode
  runs 16 dependency-free tests and skips five opt-in Torch tests.
- The existing iteration-004 analysis fixture passed 37 checks. This establishes
  analysis compatibility, not correctness of historical simulated rankings.
- Tuning tests exercised CLI preflight, complete-grid enforcement, failed trials,
  actual Optuna journal initialisation and incompatible study reuse. No full
  hyperparameter campaign was rerun.

Native local environment: macOS ARM, R 4.6.1, susieR 0.15.57 at pinned commit
`732181c662d118930e4381a5d6dc44c89ce2cea0`, sim1000G 1.40 and hapsim 0.31.
The clean Python environment used Python 3.14.7, Torch 2.14.1, NumPy 2.5.3,
SciPy 1.18.1 and Optuna 5.0.0; exact session and dependency versions belong
with each run. No duplicate-OpenMP override was used.

```bash
Rscript --vanilla scripts/challenge/test_core.R
FMB_NEURAL_SMOKE=1 .python-env/bin/python -m unittest discover -s tests/python -v
```

The full R suite requires the documented genotype dependencies; a source-only
challenge does not. Archived `sim1000G`/`hapsim` and optional method binaries
cannot be assumed present in a generic CRAN environment. The CI invariant
workflow deliberately checks the generator/evaluator without those downloads.

## Real external executables on Imperial

| Software and PBS job | Executed check | Result and limit |
|---|---|---|
| FINEMAP 1.4.2, `4309937.pbs-7` | Native and wrapper fits on a 12-variant synthetic singleton locus; joint summary/LD permutation | Exit 0; native/wrapper PIPs and permuted aligned PIPs agree exactly. This tests integration, not multi-causal accuracy. |
| GCTB 2.5.5, `4309940.pbs-7` | Official `--sbayes RC`; 256 chromosome-1 SNPs in four identity-eigen-LD blocks; 300 iterations, 150 burn-in | Exit 0; native output imported and all frozen-file/ID/allele checks passed. This short chain establishes executable integration, not MCMC convergence or genome-wide performance. |

Both ran through PBS with R 4.5.2. FINEMAP's publisher HTTPS certificate was
rejected; the publisher's HTTP archive endpoint and archive/binary SHA-256
were recorded without bypassing certificate checks. GCTB came from its official
release. See [FINEMAP smoke](../scripts/challenge/native_finemap_smoke.R) and
[official SBayesRC contract](native-sbayesrc.md) for repeatable inputs.

The SBayesRC short-chain PIP sum was 61.86 with 16 planted effects. Successful
execution plainly does not establish calibration. The comparator uses standard
`sbayes RC`, not the newer upstream `gwfm RC` workflow. Genome-wide accuracy
requires a complete panel, adequate chains and mode-specific evaluation.

## Development challenge runs

The final local and Imperial release runs, their commits, summaries and paired
contrasts are recorded here after their frozen-source rerun. Earlier snapshots
were useful for finding defects but do not replace those final artifacts.
The initial Imperial run `4309900.pbs-7` passed 1,000 tasks across 20 mechanisms
with four predictors. Later tracks expanded the specification to 22 mechanisms.

Development runs are not blind tests. Their assessor seeds and truth are
available to the developer. A final release needs independent custody and
unseen genomic/experimental targets, not simply another random seed.

## Remaining execution gaps

The local suite skipped CARMA, PAINTOR, SparsePro, original BEATRICE and Funmap
runtime tests, and local FINEMAP tests; FINEMAP was subsequently exercised on
Imperial. Bundled Functional BEATRICE was exercised in its clean Python
environment. Mocked optional adapters and source loading are never counted as
real executions. This record does not certify every catalogued method.

Real ancestry-specific reference panels, dosage/imputation mechanisms, complete
genome-wide upstream comparisons and prospectively withheld experimental
causal evidence remain separate scientific validation requirements. The new
stress tracks establish explicit mathematical mechanisms; they do not supply
those empirical claims.
