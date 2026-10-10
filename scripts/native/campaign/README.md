# Native SBayesRC comparison campaign

This is an executable experiment for official GCTB 2.5.5 SBayesRC. It is a
representative **subset** comparison, not a replication of dense genome-wide
claims. Nothing in these scripts submits jobs. The server controller owns
resource allocation and retries.

The frozen full design has five independent genotype draws for each of two
architectures, `sparse` and `sparse_inf` (sparse variance fraction .5). Each draw
has 128 empirical EUR source blocks, 200 retained variants per block and 5,000
individuals, with one phenotype shared across all blocks. The phenotype grid
uses polygenicity .001/.01, total heritability .1/.3/.5 and three independent
conditional phenotype replicates. Three binary annotations, with prevalences
.1/.2/.3 and enrichment 5, generate causal selection.

For each phenotype, expose no categories, the informative categories, or a
joint row permutation within each block. All three controls use identical
phenotypes, genotypes, summary statistics and LD; permutations preserve each
block's column burdens and between-category co-occurrence. There are 540
inputs and 2,160 native chains. The intercept is retained in the no-category
control. Regional methods receive the same phenotype summaries, block LD and
annotation exposures, with their input scope and learned-prior scope retained.
Twenty registered regional methods imply another 10,800 method-panel tasks;
they are separate from the native chain count.

The native files use verified retained source VCF coordinates and allele pairs.
`X_raw` counts ALT copies, so ALT is A1 and REF is A2. Both standardized marginal
OLS effects and SEs are divided by empirical raw-dosage SD, and A1 frequencies
are `mean(X_raw)/2`. IDs are `CHROM:POS:REF:ALT`; the variant mapping preserves
source IDs. Empirical LD eigenvalues/eigenvectors are written in GCTB's verified
little-endian float32 layout, with reconstruction checks. Numerical zero
eigenvalues are omitted and the retained variance fraction is recorded.

The joint LD representation treats blocks as independent. It omits finite
sample cross-block correlations even though summaries come from a shared
phenotype. This is an explicit working approximation, not exact joint
individual-data inference. The 128 blocks and five draws describe template and
simulation uncertainty; they do not establish transport to every human ancestry
or genomic region.

## Freeze, prepare, run and collect

Prepare all source VCFs/maps and install GCTB 2.5.5 and the R `posterior` package
before creating a plan. The plan freezes source code, VCF/map bytes, seeds,
method arguments and executable bytes. Use the same absolute plan path on
every worker. `dispatch-manifest.json` commits to the R plan and the four CSV
projections used by the scheduler; altered projections cannot supply new jobs
under an existing plan commitment.

`--regions-file` accepts a frozen coordinate override produced by deterministic
source-window repair. Full plans require 128 unique source regions spanning all
22 autosomes. Actual coordinates and the override file hash are recorded;
shifted windows are new corrected inputs, never historical data recovery.

```sh
Rscript --vanilla scripts/native/campaign/create_plan.R \
  --out /absolute/native-full --profile full \
  --vcf-dir /absolute/data/gwfm_vcf \
  --regions-file /absolute/data/native-regions.csv \
  --genetic-map-dir /absolute/data/genetic_maps \
  --gctb-path /absolute/gctb
Rscript --vanilla scripts/native/campaign/prepare.R --plan /absolute/native-full/plan.rds --index 1
Rscript --vanilla scripts/native/campaign/run_chain.R --plan /absolute/native-full/plan.rds --index 1
Rscript --vanilla scripts/native/campaign/run_regional.R --plan /absolute/native-full/plan.rds --index 1
Rscript --vanilla scripts/native/campaign/run_regional.R --plan /absolute/native-full/plan.rds --input-id input_00001
Rscript --vanilla scripts/native/campaign/collect.R --plan /absolute/native-full/plan.rds
```

Indices refer to `draws.csv`, `chains.csv` and `regional.csv`. The optional
`--input-id` mode runs every regional method for one input while loading the
individual genotype archive once, retaining separate method checkpoints.
Preparation is
memory intensive: the complete individual genotype panel and raw dosages are
kept for common-input regional comparisons. Use at least 16 GiB for full
preparation/regional workers and collection; measure actual requirements on
the canary before scaling concurrency. Native chain workers read the compact
summary/eigen files and do not repeatedly hash or load the individual genotype
archive. Input checks reject incomplete preparation and changed files. A worker
lock is removed only after the worker has stopped. Preserve failed output
directories; a failed native execution is a terminal attempted failure, while
a scheduler-killed worker without its receipt remains explicitly incomplete.
Resuming a completed worker repairs an absent terminal receipt without changing
the checkpoint. A present receipt must match the checkpoint bytes; changed
checkpoints and sealed native output files, including parameter traces, are
explicit errors. Collection accepts only receipted preparation and method
outputs.

Collection holds `analysis.lock`, moves the previous derived analysis into
`analysis-history/`, and rebuilds every table from current verified outputs.
`analysis-manifest.rds` seals the derived files and the readable JSON manifest,
and its terminal receipt commits to its bytes. The manifest records the frozen
plan, source/runtime fingerprints, explicit coverage, and collection status
(`complete`, `incomplete` or `failed`). A completed collection with pending
campaign work is `incomplete`; a validation error is retained as `failed`.

`integration` uses two blocks, n=200/p=24 and short chains. `canary` uses all 128
blocks at n=5,000/p=200 and 500-iteration chains. Both are installation/resource
experiments; neither supplies full scientific results or a convergence claim.

A worker interrupted before its terminal checkpoint can leave native output
and an inner lock. After the controller establishes that the worker has stopped
and clears its outer lock, a fresh worker acquires that lock and moves the old
output intact into `native-attempts/<attempt-id>/native-output`, recording an
archive manifest before starting again. A valid completed checkpoint is resumed
without rerunning or overwriting its output.

## Diagnostics and inference

Each input starts with four independent seeds, 25,000 iterations, 5,000 burn-in,
and thinning 10. The observed GCTB 2.5.5 release executable writes **one raw
CoreParameters row per iteration**, without an `Iter` column. Before its first
retained draw it can write uninitialised values, and between retained draws it
writes the cached previous draw. These rows are not independent samples.
The [frozen trace contract](trace-contract.json) binds the observed 14-column
header and layout to executable SHA-256
`649123c2fefacc23aa71c6537c9226b5a69f06f6045863672118009cbe680346`.
The parser requires exactly `n_iter` raw rows, derives one-based `Iter` from raw
row position, and selects only `Iter > burn_in` and `Iter %% thin == 0`. It
requires `floor(n_iter/thin) - floor(burn_in/thin)` finite retained draws; the
default full design retains 2,000 draws per chain. Both the raw core file and
`retained-core.tsv`, with explicit iteration numbers, are sealed with an
extraction record. Iteration numbers are excluded from diagnostics. Unexpected
headers, row counts or retained values fail validation, rather than triggering
a guessed alternate layout. Diagnostics use the official `posterior::rhat`, `ess_bulk` and
`ess_tail`: rank normalized/folded split R-hat must be at most 1.01, and both ESS
values at least 400. Missing traces, constant sampled parameters and nonfinite diagnostics
cannot establish convergence. These checks cover recorded core parameters.
The prespecified RC/default-gamma contract identifies only `Vg1 = 0` as the
structural zero-component genetic variance. Every retained value in all chains
must equal zero. It is reported as `declared_fixed_constant` with R-hat/ESS
missing and excluded explicitly from sampled-parameter pass criteria. This
contract is supported by the sampling code and the observed release log/draws.
Other constant columns are marked `constant_unassessable` and have their
within-chain variance recorded; constancy alone never creates an exclusion.
For example, the source's conditional residual-variance update can clamp
`ResVar` to trait variance; a constant trace there does not by itself prove
nonconvergence. Assessment status distinguishes unassessable parameters from
finite R-hat/ESS threshold failures. Neither status proves that the imported
posterior predictions are wrong. The extension table records the reason;
longer chains may not resolve a branch-clamped ancillary trace.
They cannot establish mixing of every SNP indicator. Across-chain PIP ranges
are a separate stability diagnostic: median, 95th percentile and maximum SNP
PIP range, the fraction of variants whose chains disagree on a 0.95 call, and
the range of each block's summed PIPs. These use predictions alone. Block sums
can agree while chains favour different LD-equivalent variants, so variant and
block summaries are both retained. No PIP-range or call threshold is treated as
proof of MCMC convergence or posterior precision. Per-block values are exported
in `pip-stability-by-block.csv`.

The installed-binary check used 150 iterations, 50 burn-in and thinning 5:
150 raw rows yielded 20 retained draws. Matched same-seed runs with and without
text output produced byte-identical final SNP summaries. A second same-seed
all/core/no-text retention check also gave identical SNP, parameter,
parameter-set and enrichment summaries, plus identical retained core draws
between all and core output. The
[actual Linux R-wrapper proof](wrapper-proof.json) additionally records
identical imported PIPs, effect means and effect standard deviations for all
versus core output on 256 variants, with all seven sink targets verified.
This establishes that
specific trace/output contract, not accuracy or convergence. The pinned source
commit is a source-audit reference; the compiled commit of the later release
executable has not been established. The source's robust-restart branch shadows
its original sample vector, which raises a possible final-summary/replacement-
trace coherence issue. This is not a demonstrated bug in the release binary.
Actual restart/disruption log markers therefore prevent ordinary trace
diagnostics and receive a coherence-unverified flag; execution results and logs
remain in the archive. See the [source restart branch](https://github.com/jianzeng/GCTB/blob/cc7fa7d765c83a89c6375946cf77fe50ba1a317e/scr/main.cpp#L255-L276).

The campaign explicitly requests `mcmc_scope = "core"` from the unmodified
native executable. Seven full-variant text arrays (`PIP`, `PEP`,
`DeltaPi1`–`DeltaPi5`) are redirected through exact `/dev/null` links; raw core,
small annotation/block traces and final SNP summaries are preserved. Each
checkpoint records discarded-array paths and the retention policy. The sink
links are excluded from byte manifests, with their exact targets validated.
No per-variant MCMC draws are saved under this policy, and it does not avoid
the executable's formatting CPU cost. All-output text retention would print
4.48 billion variant-array numbers per default full-sized chain, making it an
inappropriate default for 2,160 fits.

All initial outcomes stay in the archive. A prespecified optional extension
uses four fresh chains of 100,000 iterations/20,000 burn-in/thinning 10 for every
input failing core diagnostics, independent of accuracy. Extensions require
their own frozen plan and are reported separately; the current scripts do not
silently extend or replace initial results.

Collection records planned, attempted, completed and successful counts
separately, including missing scheduler outputs. A native consensus uses all
four successful equal-length chains; partial-chain averages do not become
primary results. Official PIPs target the event of **any nonzero mixture
component**. Both planted sparse membership and true any-nonzero-effect events
are scored, with event alignment flags. Sparse-member rankings for native
SBayesRC are operational diagnostics, not calibrated estimates of that same
event. In `sparse_inf`, the full nonzero target may have no negative variants;
its probability score has limited discrimination. Normalized absolute-z scores
are excluded from probability scores. Native credible sets are not fabricated.

`task-metrics.csv` contains event-typed locus outcomes and explicit completion
state. `draw-metrics.csv` aggregates each independent genotype draw;
`paired-draw-comparisons.csv` compares native/regional methods and annotation
controls on paired complete draw panels. Intervals resample genotype draws,
preserving shared phenotype and learned-prior dependence. Five draws yield
limited generalization precision, so inspect draw-level effects alongside the
intervals. Successful-only AP/Brier contrasts remain conditional; recall counts
all completed failures. Pending work is reported as pending and cannot support
a final ranking. Full planned chain coverage is exported even for draws not yet
prepared. Native chain logs, parameter draws and consensus PIPs remain archived.

Task rows distinguish native execution/import success from
`core_diagnostics_pass`, `native_core_assessment_status`,
`native_trace_coherence` and `native_accuracy_eligible`, and carry the PIP
stability fields directly. Native-specific chain diagnostics are missing for
regional comparators; they are not subject to a different method's MCMC test.
Raw per-locus probability scores remain available as measured operational
outcomes. Draw tables preserve `*_imported` summaries. Ordinary native accuracy
columns and their paired contrasts require every input in the draw/setting to
have four successful, coherent chains passing the core criterion; otherwise
those columns are missing. This prevents comparing a favourable native subset
against the regional comparator's complete setting. Passing this strict core
criterion still does not establish mixing of unrecorded variant trajectories.
Posterior conditioning and native null-support metadata accompany the scores;
a model forced to contain a cause can receive useful operational probability
scores on null tasks, but conditional PIPs cannot be described as calibrated
unconditional probabilities merely because their membership labels match.

Primary implementations and diagnostic definitions: [GCTB pinned source](https://github.com/jianzeng/GCTB/tree/cc7fa7d765c83a89c6375946cf77fe50ba1a317e),
[GCTB trace writer](https://github.com/jianzeng/GCTB/blob/cc7fa7d765c83a89c6375946cf77fe50ba1a317e/scr/mcmc.cpp),
[GCTB conditional residual-variance update](https://github.com/jianzeng/GCTB/blob/cc7fa7d765c83a89c6375946cf77fe50ba1a317e/scr/model.cpp#L2141-L2155),
[posterior R-hat](https://mc-stan.org/posterior/reference/rhat.html),
[posterior bulk ESS](https://mc-stan.org/posterior/reference/ess_bulk.html).
