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

## Diagnostics and inference

Each input starts with four independent seeds, 25,000 iterations, 5,000 burn-in,
and thinning 10. GCTB's `CoreParameters.mcmcsamples.txt` contains post-burn-in
draws already; the parser verifies the expected count and does not discard
burn-in again. Diagnostics use the official `posterior::rhat`, `ess_bulk` and
`ess_tail`: rank normalized/folded split R-hat must be at most 1.01, and both ESS
values at least 400. Missing traces, constant chains and nonfinite diagnostics
cannot establish convergence. These checks cover recorded core parameters.
Constant columns are marked `constant_unassessable` and have their within-chain
variance recorded; no column is inferred to be a known fixed model constant
from its observed constancy.
They cannot establish mixing of every SNP indicator. Across-chain PIP ranges
are a separate stability diagnostic. No PIP threshold is treated as proof of
MCMC convergence.

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

Primary implementations and diagnostic definitions: [GCTB pinned source](https://github.com/jianzeng/GCTB/tree/cc7fa7d765c83a89c6375946cf77fe50ba1a317e),
[GCTB trace writer](https://github.com/jianzeng/GCTB/blob/cc7fa7d765c83a89c6375946cf77fe50ba1a317e/scr/mcmc.cpp),
[posterior R-hat](https://mc-stan.org/posterior/reference/rhat.html),
[posterior bulk ESS](https://mc-stan.org/posterior/reference/ess_bulk.html).
