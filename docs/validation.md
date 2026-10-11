# Validation record — 11 October 2026

This records executed checks, not a universal method ranking. Original upstream
source was `f053fdf318a3abe2c44eeb0ac3f3c99a03979b81`. The earlier development
implementation is [`3ac28638a8e4249c2dec2b76eb131c6ff6bf27aa`](https://github.com/NathanSkene/fine-mapping-benchmark/tree/3ac28638a8e4249c2dec2b76eb131c6ff6bf27aa).
The implementation and
development artifacts are frozen separately; run provenance records source
commit and individual R file hashes. The [review](REVIEW.md) and
[scientific contract](challenge-science.md) define interpretation limits.

## Integration failures and repairs — 11 October

The immutable `eb2d9bc34190c9fc89cfa25b43ecb2f070139fc8` integration finished
168/168 task receipts and 1,404 locus attempts. Its identity/completion census
passed: **1,158 successful attempts and 246 recorded failures**. Of these,
240 were declared input or algorithm absences: 152 non-scalar-N cases, 72
unsupported Gaussian INF cases, four exact/forced-one overlap cases and 12
annotation-free Funmap cases. Two FINEMAP probability-coherence rejections and
four SparsePro set-parser rejections required repairs. Those original errors
and outputs remain preserved. Completion is not equivalent to scientific
success.

The [repair ledger](validation-artifacts/integration-repairs-20261011/README.md)
and [source-bound records](validation-artifacts/integration-repairs-20261011/evidence.json)
distinguish the following evidence:

- A replay of one actual planned `p = 1000` known-overlap case diagnosed both
  zero-LD NaNs and subnormal division distortion. Direct normalized-Gram
  construction preserves every other public-input field and all truth fields,
  with unchanged validators and PSD thresholds.
- Two actual FINEMAP replays returned archived PIPs exactly. Inspected binary
  output precision explains the small union deficits without changing PIPs or
  scores. Public evaluation remains strict unless an assessor explicitly
  enables the exact-binary, token-bound contract. The count-log half-unit
  policy is separate and its active rounding mode is unverified.
- Five actual SparsePro native probes passed the v2 transport/parser replay
  with unchanged native PIPs. Its earlier harness failure is retained. The
  v3 parser-only harness failed before parsing on a Linux missing-path guard;
  the strict v4 replay then passed on all five sealed outputs, with original
  input hashes, exact copied bytes and independent PIP/set/effect alignment.
  It made five intercepted calls and zero new inference runs.

The final repaired-source v3 R package check passed **Status: OK**, with
**1,824 passing expectations, 24 explicit optional/repository-only skips and
zero errors, warnings or notes**. `NOT_CRAN=true` enabled the pipeline test.
Separate source-only native orchestration checks passed **218 expectations**
with no failures, errors, warnings or skips. The repair ledger binds the exact
archive, all 159 archive file hashes, validation driver and logs; 264 current
scientific/source files were independently verified against the frozen manifest.
These checks certify the tested software snapshot. Fresh immutable campaign
source/plans and execution gates remain pending here. Prior native integration quality gates
passed specified execution contracts with 20 retained draws per chain; this
cannot establish convergence or performance. Historical direct rescoring also
still needs Lucas's original raw simulations and per-method predictions, as
specified in the [original-results audit](original-results-audit.md).

A separate resource-probe failure exposed nominal-versus-realized genotype
size: monomorphic filtering can remove variants after a frozen upstream draw.
The exact-source audit records per-locus dimensions/fractions and all warnings,
preserves genotype draws/RNG and permits no resampling. Canonical VCF shortage
and bundled-capping warnings remain fatal, along with unexplained dimension
loss. Its 74 focused assertions passed independently; the source-bound
302-assertion proof also includes existing campaign and precision regressions.
This is not a regenerated full panel. An intermediate staged-v2 package
check passed for its archived source, before this audit and a FINEMAP
timer-accounting correction; it does not certify the final repaired source.

## Prior campaign hardening checks

The preceding source snapshot's complete R package check passed **Status: OK**, with **1,559
passing expectations, zero failures and zero warnings**. Its 24 explicit skips
comprise seven CARMA, one configured legacy BEATRICE, two FINEMAP, one Funmap,
two PAINTOR, two SparsePro and nine repository-only native orchestration tests.
Those nine require helpers deliberately excluded from the R package; they have
separate repository checks. Optional native execution is a separate server gate.
The [package proof](validation-artifacts/campaign-hardening-20261010/R-package-check.json)
records the exact archive, all 157 archived file hashes, check-log hashes and
skip reasons. All 155 frozen package-source files matched the working source
when that check completed. Controller and collector scripts are checked
separately and are outside this archive.

Actual GCTB execution established the [raw trace layout and all/core/no-text
output equivalence](../scripts/native/campaign/retention-proof.json). A separate
[R wrapper check](../scripts/native/campaign/wrapper-proof.json) imported
identical PIPs, effect means and effect SDs with all versus core output and
verified all seven exact text sinks. The [FINEMAP parser proof](../scripts/campaign/finemap-parser-proof.json)
preserves native rounded model-size probabilities without renormalising them.
These small execution fixtures establish interfaces and retention semantics;
they cannot establish full-size costs, accuracy or convergence.

The [BEATRICE event-import proof](validation-artifacts/campaign-hardening-20261010/beatrice-events-proof.json)
records actual 500-iteration vanilla and functional training on a tiny synthetic
locus, expansion, native output and strict R import of conditional joint events.
It binds all six production R/Python/contract files by SHA-256. These complete
interfaces supplement the analytical probability examples and real neural
tests; the tiny fixture cannot estimate comparative accuracy or campaign cost.

The [fresh command-line campaign check](validation-artifacts/campaign-hardening-20261010/campaign-cli-check.json)
ran ten actual ABF/marginal-score tasks with 34 dataset records and 148 locus
rows. It verified safe resume, exact-fit reuse, sealed bundle corruption refusal,
incomplete-work denominators and empty-analysis rebuilding. The bounded collector
matched the in-memory reference for typed raw metrics, all summaries, scoring
attribution and paired intervals. Separate stream tests exercised extreme double
values, infinite losses, missing values and quoted text; resource tests verified
packed-task overhead allocation. These are small correctness experiments.

Fifty-five scheduler regressions passed, including submission intent committed
before `qsub`, ambiguous outcomes, failed observations, retry cleanup and altered
immutable projections. Regional attempts count durable task-start intents;
unstarted batch tails remain eligible. A separate limit bounds failures before
the first start. Terminal observations trigger a final start-ledger read before
assignment release. Native jobs count confirmed unit submissions. Both
controllers report exhausted work explicitly and claim ownership by atomic
directory creation before accessing state. Unknown submissions, malformed
start ledgers and stale ownership claims require reconciliation.

The [local full-design ledger probe](validation-artifacts/campaign-hardening-20261010/sqlite-controller-scale.json)
records source-bound measurements for 649,220 method tasks and 240 generation
units. The schema-3 probe used approximately 32.05 MiB peak process memory and
produced a 1,280-byte status projection after 100,000 synthetic terminal
jobs/events. Its indexed check and synced empty retry-blocker report took
2.44 ms; this timing excludes a populated dependency report. The earlier
schema-2 measurements are preserved separately. These are local storage and
memory measurements. Real PBS execution and shared-filesystem ownership have
separate integration gates.

The [expanded Python proof](validation-artifacts/campaign-hardening-20261010/python-validation.json)
records **97 passing tests**, zero failures, errors or skips, with resource
warnings treated as errors. This fresh run includes all 55 scheduler regressions
and actual R/Torch campaign integrations. The earlier **30-test RL/Docker**
proof is retained with its unchanged relevant source hashes; it is not used to
claim coverage of the repaired controllers. The R, BEATRICE and RL scientific
source bytes remained unchanged through these scheduler repairs.

The [actual server filesystem record](validation-artifacts/campaign-hardening-20261010/server-filesystem-ownership.json)
preserves a failed cross-node advisory-lock check: PBS job `4311453.pbs-7`
acquired `flock` while the login host held the same shared file. This failure
motivated the atomic directory ownership guard. The actual frozen production
helper subsequently passed on login and compute hosts in `4311477.pbs-7`:
the second owner was refused, synced metadata was readable, claims released
cleanly and serialized SQLite commits were visible in both directions. The
[reproduction fixture](validation-artifacts/campaign-hardening-20261010/probe-controller-ownership-cross-node.py)
is bound to the tested helper hash. These checks establish observed coordination
behavior; they do not simulate storage failure or automate stale-owner recovery.

## Earlier development correctness checks

The checks below describe earlier source snapshots. The campaign-hardening
record above and its exact source hashes supersede their test totals; archived
development predictions retain their original interpretation limits.

- The earlier installed-package suite passed **947 expectations**, with zero
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

The local SuSiE/SuSiE-inf rows below are archived prototype results and have been
superseded for Gaussian-score comparisons. Their adapters used an OLS-t
reconstruction for normal scores. The current SuSiE adapter fits the typed
fixed-noise Gaussian likelihood; SuSiE-inf and FINEMAP-inf explicitly reject
Gaussian tasks without supported phenotype sufficient statistics. The exact,
forced-one, annotation reference and ABF rows are unaffected by this correction.
Use newly frozen campaign artifacts for current method comparisons.

Both final development runs used the validated implementation above, 22
scenario tracks including controls, seed `20261010`, and sample-size parameter `10000`:

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
