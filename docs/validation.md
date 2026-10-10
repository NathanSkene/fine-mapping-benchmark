# Validation record — 10 October 2026

This records executed checks, not a universal method ranking. Original upstream
source was `f053fdf318a3abe2c44eeb0ac3f3c99a03979b81`. Validated implementation
is [`3ac28638a8e4249c2dec2b76eb131c6ff6bf27aa`](https://github.com/NathanSkene/fine-mapping-benchmark/tree/3ac28638a8e4249c2dec2b76eb131c6ff6bf27aa).
The implementation and
development artifacts are frozen separately; run provenance records source
commit and individual R file hashes. The [review](REVIEW.md) and
[scientific contract](challenge-science.md) define interpretation limits.

## Executed correctness checks

- The final installed-package suite passed **947 expectations**, with zero
  failures or errors, 15 explicitly skipped dependency cases and one warning
  about the legacy single-region SBayesRC-style annotation refit. The 17
  repository-only tuning assertions also passed: **964 R expectations** across
  these checks. The source-only invariant suite passed **363 expectations**.
- The focused invariant suite checks hand-computed Gaussian Bayes factors,
  prior-predictive calibration, exact tied AP, raw scores, distinct set matching,
  intercept OLS, realised variance/covariance, annotation probability prevalence,
  singular LD, input recoding, known study overlap, unavailable causes, malformed
  predictions, release tampering, failure denominators and tuning metadata.
- `R CMD build` and `R CMD check --no-manual` passed **Status: OK**, including
  vignette construction/rebuilding, examples, installed-package loading,
  documentation/API consistency and all package tests. The haplotype vignette
  leaves download-dependent chunks opt-in; its Gaussian example executes.
- GitHub's full package check and scientific-invariant jobs passed on both
  Linux and macOS for the validated implementation. See the
  [package-check run](https://github.com/NathanSkene/fine-mapping-benchmark/actions/runs/38074856951)
  and [invariant run](https://github.com/NathanSkene/fine-mapping-benchmark/actions/runs/38074857462).
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

Both final development runs used the validated implementation above, 22
mechanisms, seed `20261010`, and sample-size parameter `10000`:

| Run | Tasks and predictors | Evidence |
|---|---|---|
| Local ARM R 4.6.1 | 220 tasks, 24 variants, 10 independent synthetic draws per mechanism; six predictors including pinned SuSiE/SuSiE-inf | [Summary](validation-artifacts/local-development-final/summary.csv), [paired contrasts](validation-artifacts/local-development-final/paired-comparisons.csv), [PDF](validation-artifacts/local-development-final/report.pdf), [provenance](validation-artifacts/local-development-final/provenance.txt) |
| Imperial PBS `4309969.pbs-7`, R 4.5.2 | 1,100 tasks, 40 variants, 50 draws per mechanism; four dependency-light predictors | [Summary](validation-artifacts/imperial-development-final/summary.csv), [paired contrasts](validation-artifacts/imperial-development-final/paired-comparisons.csv), [PDF](validation-artifacts/imperial-development-final/report.pdf), [PBS status](validation-artifacts/imperial-development-final/pbs-status.txt) |

The final Imperial job exited 0 in 2m39s wall time, with 2m30s CPU time and
282,148 kB peak memory. Source SHA-256 files and R sessions accompany both
artifact directories. Their recorded R file hashes were independently checked
against the final implementation. Local Git state was clean; Imperial used an
exported snapshot with explicit source commit and unavailable Git dirty state.

Unsupported likelihoods remain visible: the exact/forced/annotation references
complete 1,050/1,100 tasks and reject missing-overlap metadata. ABF completes
900/1,100, rejecting variable-n and multi-cohort inputs. No other unexpected
failures remain in this run. Those counts describe scope, not an overall rank;
proper-score comparisons use eligible paired tasks and mechanism-specific
uncertainty. The local SuSiE wrappers have the same regional-input restriction.

The initial Imperial run `4309900.pbs-7` tested an earlier 20-mechanism snapshot.
Release attempt `4309955.pbs-7` failed before producing benchmark results:
minimal R lacked `nnet`, exposing an unnecessary dependency in fixed-prior
sampling. The corrected fixed-prior path runs in base R; annotation learning
fails explicitly without its dependency. Ubuntu CI also exposed missing libuv
headers while installing its test library; the workflow now uses a binary
repository/bundled libuv and verifies successful installation. These failures
were corrected before the final runs rather than counted as successes.

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
