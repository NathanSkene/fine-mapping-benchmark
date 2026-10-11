# Full paired campaign

This runner preserves the complete original Iteration-004 grid and adds
Iteration-003-compatible controls and a separate stress family. A `full` plan
contains 45 original-design rows plus 135 control rows, 250 scenario draws per
row and ten loci per scenario. The control rows restore sample size 1,000,
variant counts 100–1,000, heritabilities .0075–.4, and independent LD references
of 200 and 500 individuals. The latest design retains sample size 5,000 and
500–2,000 nominal variants. Truncating an undersized VCF window fails generation.
The exact frozen upstream genotype source (SHA256 `651dbb1640b63187028c12406a82801a676f10be1f9a52b635aedad5a49be139`)
selects the requested variants before simulation and then removes columns
monomorphic in the GWAS or reference sample. Valid losses from that known
filter are retained with a per-locus audit of requested, realized and dropped
counts and fractions. Other generator sources must realize exact dimensions;
extra variants, invalid matrices and contradictory metadata fail generation.
Every original warning remains unchanged. The upstream warning that fewer
variants were available after VCF MAF filtering remains fatal even for the
known source. Bundled-VCF capping at 500 is also a fatal source shortage. There
is no padding, resampling, retry or change to genotype, phenotype or RNG draws.
The compact audit is stored in the simulation parameters
and generation checkpoint; metrics distinguish nominal design size from the
retained variant universe. Old failed generation attempts remain failed records
under their original source and plan.

With 20 registered methods and three paired stages this is **27 million locus
fits**, before the additional 13,200 Gaussian stress loci. The plan is the full
requested campaign, not evidence that these expensive jobs finished. Measure
actual preflight and canary runtimes, allocate resources, and inspect terminal
checkpoints before describing the campaign as complete. Native whole-panel
SBayesRC is a separate experiment with its own supplied-panel scope and chain
diagnostics; the registered `sbayesrc` remains the local surrogate.

## Attribution and historical evidence

Supply `--legacy-root` only for genuine preserved `job_*/sim.rds` plus
`scenario_*/results.rds` archives. Input and output hashes are frozen before
analysis. `archived_rescore` applies both old and corrected scoring to the same
saved predictions. Missing original predictions remain missing evidence.

When historical artifacts are unavailable, supply a frozen original source
checkout instead. The old generator produces new recorded-seed simulations;
these are **original-code reproductions with newly supplied haplotypes**, not
recovered historical results. The upstream and repaired inference methods then
receive identical observed inputs. Their saved PIPs are scored under both
evaluators. A third stage regenerates corrected phenotypes on the same frozen
genotypes. This final contrast changes the data-generating process and cannot
alone establish that an original ranking was wrong.

Unchanged implementations can reuse a realized fit only when their reachable
authored R functions/globals, exact observed inputs, arguments, paired RNG seed,
and frozen runtime dependencies have identical fingerprints. Dynamic lookup,
ordinary methods reading private truth, and mutable authored globals opt out.
External-source and binary hashes are frozen as well. Reuse is recorded in each
checkpoint and locus table; it establishes that the two arms executed the same
estimator, not an additional independent replication or compute advantage.

## Prepare and execute

Create an original source snapshot at the immutable original commit. Preserve
the old Python trainer as well as the old R wrappers. Supply widened,
checksum-verified VCFs and matching genetic maps explicitly; dataset provenance
must state population composition and reference-window changes.

```sh
FMB_TOOL_ROOT=/absolute/tools FMB_PYTHON=/absolute/env/bin/python \
  Rscript --vanilla scripts/campaign/create_plan.R \
  --out /absolute/frozen-campaign \
  --original-source /absolute/original-snapshot \
  --vcf-dir /absolute/verified-vcfs \
  --genetic-map-dir /absolute/maps

Rscript --vanilla scripts/campaign/preflight.R --plan /absolute/frozen-campaign/plan.rds --index 1
Rscript --vanilla scripts/campaign/generate.R --plan /absolute/frozen-campaign/plan.rds --index 1
Rscript --vanilla scripts/campaign/run_task.R --plan /absolute/frozen-campaign/plan.rds --index 1
Rscript --vanilla scripts/campaign/collect.R --plan /absolute/frozen-campaign/plan.rds
```

Submit these compute commands through the site's scheduler, using
`generation.csv` and `tasks.csv` for array dimensions. Generation and method
preflight must finish before downstream tasks. Default method task chunks hold
five whole scenarios, preserving every cross-region prior; CARMA holds one
because its measured large-locus cost is substantial. There is no method-wise
fallback or missing-dependency substitution. Preflight requires actual
successful execution on two annotated loci, including native external binaries.
An `integration` profile is solely a small installation test and is always
labelled accordingly.
SparsePro additionally runs a fixed coherent perfect-LD transport fixture and
must import a multi-member native set; singleton fits alone cannot exercise its
slash-separated set fields. This gate tests the interface, not localization of
indistinguishable variants. A later source change requires new source-bound
preflight receipts rather than reusing an earlier successful installation.
The `canary` profile preserves the full nominal 500–2,000-variant,
5,000-individual locus designs while selecting four sparse/background and annotation arms with one
scenario each. It measures real campaign resource costs and never replaces the
full scientific profile.
The `resource_probe` profile selects the same four full-size legacy arms and
one exact full stress batch: p=1,000, N=10,000, ten draws per mechanism, 220 mixed
public tasks. This probes pooled prior fitting and matrix storage at the actual
largest batch dimensions. With 20 regional methods it has five generation units
and 259 method task units, rather than the full scientific campaign.

An optional `--method-args` RDS overrides the exact frozen method arguments.
The default matches the original worker, including FINEMAP maximum five effects
and the original Functional BEATRICE tuning values. These capacities and tuning
choices are disclosed, not claimed to be a uniformly fair compute allocation.
Use a separate frozen plan for changed support or hyperparameters. A missing
annotation arm is out of scope for Funmap, not a failed dependency installation.

## Resume and scoring

Each dataset-method-stage artifact is atomically checkpointed, preserving both
predictions and evaluations. Resume requires identical plan, task and source
hashes. Data hashes are checked before loading. A lock protects live tasks;
after a scheduler confirms a terminal killed process, an operator can remove
that task's stale lock before resuming. A timeout is not successful completion.
New plans compact each completed task into one `completed.rds` bundle and its
byte-hash receipt. The bundle retains every ordered dataset prediction,
evaluation, session record, and started-attempt record. Per-dataset fingerprints
and the exact cleanup ledger are retained inside it. A task contains at most
five datasets; partial tasks keep their loose checkpoints for resume. Cleanup
starts only after the bundle and receipt are sealed and all remaining loose
files match their recorded hashes. A changed loose file remains available for
investigation. Interruptions between bundle, receipt and cleanup writes are
explicitly resumable. Collectors accept both the earlier loose layout and this
sealed bundle layout and never count both copies.

The default full plan has 649,220 task units and 2,701,220 dataset-method-stage
records. Keeping one attempt RDS, prediction RDS and inner receipt for every
record plus two outer files would create 9,402,100 regular files before caches
or scheduler logs. The former broad fit cache could add millions more. The new
terminal layout uses 1,298,440 bundle/receipt files and 649,220 task directories;
inputs, temporary partial-task files, scheduler scripts/logs and analysis
artifacts must still fit the quota. Exact fit reuse reads the immutable paired
original record after verifying source/input/settings/RNG/runtime equality;
there is no duplicate persistent fit cache for new plans. Old layouts are not
silently rewritten under their frozen plan IDs: changed scientific workers
require new plans.

The stress family covers all 22 tracks at p=40,200,1,000 and N=1,000,10,000,
100 independent draws per cell. Pooled methods see mixed-mechanism panels of
220 public tasks, with no mechanism labels or truth. The explicit small-model
reference is restricted to p≤200. Regional scalar-N methods record unsupported
multi-cohort and variant-specific-N tasks as attempted failures; the covariance-
aware reference can use the supported known-overlap inputs. Gaussian-summary
and real-haplotype results are reported separately.
SuSiE, PolyFun and the local SBayesRC surrogate receive explicitly typed normal
scores with specified residual variance. The installed SuSiE-inf path rejects
this fixed-noise summary experiment because its infinitesimal-variance estimator
needs additional phenotype sufficient statistics. FINEMAP-inf similarly rejects
normal scores because its author pipeline needs the phenotype second moment.
These unsupported events are retained; the legacy OLS experiments execute both.

`analysis/status.csv` is the authoritative completion census.
`task-metrics.csv` contains every checkpointed attempted locus;
`planned-coverage.csv` separately retains planned, checkpointed and missing
dataset denominators. Interrupted workers preserve started-attempt records and
scheduler outcomes, rather than turning absent outputs into successful fits.
Interim scientific summaries are incomplete until this census is complete.
The default collector streams one task bundle at a time, then aggregates one
stage/design stratum at a time. `--max-stratum-rows` sets an explicit bound
(100,000 by default); exceeding it fails collection and leaves a failed analysis
manifest rather than changing the denominator. The full profile rejects the
all-memory backend. The `--backend memory` path remains available for small
integration regression tests against the streamed summaries and paired intervals.

`task-metrics.csv` retains every metric column with 17 significant decimal
digits, logical values, infinities and an unambiguous plan-specific missing
marker. Read its column types and marker from `task-metrics-schema.rds` or JSON.
The exact typed values also remain in `task-metrics.blocks`: a bounded stream
of length-prefixed R-serialized data frames. In the streamed layout,
`task-metrics.rds` is a sealed index for that stream, not a materialized data
frame. `.campaign_each_metric_block()` iterates it without loading the census.
The analysis manifest commits to both files. Temporary stratum partitions use
hexadecimal doubles so R's decimal parser cannot perturb very small numbers;
they are deleted after aggregation. Correctly rounded decimal readers recover
the CSV's doubles; use the binary stream when exact R values matter.

`scoring-attribution.csv`
compares old and corrected scoring on identical predictions;
`regime-summary.csv` and `paired-comparisons.csv` report settings separately.
Reused genotypes are clustered by frozen panel and region. Comparisons involving
joint or pooled estimators use the whole panel as the uncertainty cluster,
because their learned prior couples loci. Empirical template generalization is
distinct from conditional phenotype-noise uncertainty. Normalized absolute-z
baseline scores are excluded from PIP probability scores. Availability, failures,
calibration, retrieval and runtime are retained rather than hidden in an overall
leaderboard.

Exact log loss retains infinity when a method assigns zero probability to the
realized class. Summaries report these counts and include infinity in the mean;
the challenge's clipped loss is retained as an explicitly named diagnostic.
AP eligibility requires an observed causal variant. The auxiliary failure-zero
AP composite uses that fixed eligibility for both successes and failures;
null and entirely missing-cause tasks remain in availability/false-call counts.
SBayesRC surrogate PIPs retain their any-nonzero-mixture-component event and an
explicit probability-target alignment flag. Against sparse membership their
Brier/log-loss contrasts are operational diagnostics. When total-effect truth
is available, separate any-nonzero-total-effect scores are also retained; older
inputs without that decomposition leave those scores missing.
CAVIAR's native causal superset is evaluated by all-causal containment, support
size and causal-variant recall. Its per-signal distinct-set recovery is undefined,
since one superset does not assert one signal. A native containment event is
reported separately from generic any-causal-hit diagnostics.
Posterior conditioning, native null-model support, reported null probabilities
and configuration counts are also retained. PAINTOR enumeration excludes the
null configuration; FINEMAP support is taken from the native logged prior.
Conditional nonnull model scores remain operational diagnostics against planted
sparse membership, not unconditional calibration claims. An unreported null
probability remains missing; a similarly named flag cannot fabricate one.

After real full-size canaries finish, `cost_report.R --canary-plan PATH
--full-plan PATH --workers 48 --out PATH` measures runtimes and compressed
retention bytes, including bundle headers, attempt ledgers and terminal receipts.
It reports adapter and end-to-end worker-hour extrapolations separately, and a
file-count census alongside byte estimates. A `resource_probe` plan also supplies stress costs from
the full 220-task batch. Largest-panel estimates extrapolated to smaller stress
dimensions are planning anchors, not measured runtimes at every dimension.
Its estimates disclose unresolved methods and exclude input archives, queue
delay and exact-fit reuse savings. The runner
currently retains complete adapter results, including available posterior means,
convergence and component diagnostics; plan storage from measured bytes before
dispatching the full campaign.
The worker timer spans the first R script expression through completed bundle
compaction; it excludes R process startup, the final timing-receipt write and
lock removal. Dataset processing, input IO, prelude and compaction are retained
separately in the existing terminal receipt. The estimate scales dataset work
by loci and startup/input/compaction overhead by planned worker tasks, so a
one-scenario canary does not multiply startup five times for a packed worker.
Interrupted timing writes and partial-resume timings remain unknown for this
extrapolation. PBS `resources_used` in the scheduler ledger measures the full
process/batch separately, including failed or timed-out attempts. Input archive,
temporary partition, retained metric-stream/CSV and prior-analysis archive bytes
are outside the primary-checkpoint storage estimate and need their own budget.
The cost input manifest records the exact sampled checkpoint/receipt bytes.

## Bounded PBS orchestration

`dispatch.py` fills a bounded queue rather than submitting the entire plan as a
huge array. Its defaults allow 64 single-CPU jobs, at most 16 simultaneous data
generation jobs, and four method tasks per worker (two for CARMA). It rotates
through method/stage groups so optional methods make progress early. Supply a
frozen bootstrap shell script with the actual R/Python module and library paths:

```sh
python3 scripts/campaign/dispatch.py \
  --plan /absolute/frozen-campaign/plan.rds \
  --repo /absolute/source-snapshot \
  --bootstrap /absolute/frozen-environment.sh \
  --max-workers 64 --max-generation 16 --once --dry-run
```

Remove `--dry-run` for a real dispatch cycle; remove `--once` for a persistent
controller. The schema-3 `dispatch/ledger.sqlite3` stores the frozen plan
projection, per-unit assignments and start-intent budgets, every job and its
exact script hash, and the append-only observation/retry history. The controller streams the
CSV projection into SQLite once and queries indexed ready work; it does not
load all tasks or historical jobs into Python on each poll. The small
`dispatch/state.json` contains the configuration hashes, counts and currently
active jobs. It is a status view rather than submission or completion authority.
Terminal job history remains in SQLite. The dispatch manifest commits to the
plan RDS bytes and the task,
dataset and generation CSV projections; editing a CSV cannot silently change
the controller's interpretation of the plan. An immutable-projection digest
also detects changed fixed rows in the ledger on resume. Terminal receipts
commit to the R checkpoint bytes and plan ID.

Before opening SQLite or reading scheduler state, the controller claims
`dispatch/controller.owner/` with an atomic directory creation. Its synced
`owner.json` records hostname, PID, a random ownership nonce and creation time.
Any existing claim blocks another controller, including a directory with
missing, malformed or stale metadata. A matching owner removes its own claim
on graceful normal or exceptional exit; a killed process leaves a claim for
explicit operator investigation. Never remove it until the recorded process
and scheduler state have been reconciled. A local `flock` remains a secondary
guard: an actual deployment probe showed it alone did not exclude a second
controller on another node of the shared filesystem.

The controller uses SQLite rollback journaling
(`journal_mode=DELETE`, `synchronous=FULL`). Submission intent and all reserved
unit assignments commit in one transaction **before** calling `qsub`. A returned
handle commits without charging task-start budgets; a timeout, malformed handle or
interrupted confirmation leaves a reserved unknown/pending submission and stops
the controller. Inspect its recorded unique PBS job name before reconciling it;
deleting the ledger or clearing reservations to obtain a retry can launch the
same computation twice. A temporary `qstat` failure retains the assignment and
counts toward concurrency. A newly issued missing handle has a 120-second
propagation grace period. An authoritative terminal/missing observation commits
the released assignments and a cleanup obligation before deleting abandoned
worker locks. Resume finishes that cleanup before a fresh attempt and preserves
the earlier scripts and observations.

Each job has a pre-created, synced `job-<token>.starts.jsonl`. Immediately
before each sequential R command, a small Python worker mode appends and syncs
an ordered task-start intent. `--start-python` selects an interpreter available
on compute nodes; its path and bytes are bound to the controller configuration.
The controller validates the ledger's exact job token, plan/source commitments,
task kind, assigned index and ordinal, then counts each start once in an atomic
SQLite transaction. It re-reads the final prefix after terminal scheduler
evidence before releasing assignments. An early batch failure leaves untouched
suffix tasks with their entire retry budgets. A start intent is conservative:
a process killed after syncing the intent but before invoking R consumes that
one task's budget. It does not prove that R executed.

A missing, malformed, truncated or rewritten start ledger stops dispatch and
retains assignments. A valid empty ledger differs from a missing ledger. Only
authoritative terminal/missing scheduler evidence with that valid empty ledger
charges a separate pre-start infrastructure failure to the batch's first unit;
it never charges the untouched suffix. `--max-attempts` bounds task-start intents
and `--max-prestart-failures` bounds these infrastructure failures, both default
to three. Transient observations and unresolved submissions consume neither
budget and never authorize duplicate launch.

The controller requires a filesystem with reliable atomic directory creation
and SQLite locking/sync semantics; validate those on deployment storage across
nodes, not only on one host. Workers never
open this database. WAL is deliberately avoided because SQLite documents its
[network-filesystem restriction](https://www.sqlite.org/wal.html). For a live
database backup, use SQLite's backup API. For a plain file copy, first stop the
controller cleanly and close database readers, then preserve the entire
`dispatch/` directory, including the ledger, job scripts and start ledgers. Do not remove a hot
rollback journal left by an interrupted process. Legacy schema-1/2 scheduler state is
explicitly rejected rather than silently migrated; a schema-3 status file with
a missing ledger is also an error. Preserve and inspect any existing history
before deliberately starting a new campaign ledger.

The controller revalidates all terminal receipts once at startup, including
previously cached completed IDs, rebuilds its status counters from those
verified hints and checks live assignment consistency before submitting work.
Loose dataset predictions and earlier cache
layouts have their own byte-hash receipts; a completed method task additionally seals
its full dataset/receipt ledger (or embeds its records and ledger in a sealed
bundle). A missing receipt is pending work. An invalid
receipt or changed sealed artifact is an integrity error. Unsealed dataset RDS
files are preserved outside the collection namespace and recomputed on resume.
Collection holds its own lock, archives the previous derived analysis, and
builds a fresh analysis with a sealed manifest of the plan, source, runtime,
completion status and derived-file hashes. An empty or failed collection cannot
leave previous scientific tables in the current analysis directory.

Runtime guards include installed package versions, external Python repository
source hashes, the exact BEATRICE credible-set contract JSON, configured
executable/interpreter bytes and threading settings.
They do not hash installed R libraries or installed Python package bytes.
Same-version reinstalls can evade that version guard: keep the installed
environment immutable for the campaign and retain its installation manifests.
This is not a byte-complete commitment to every installed dependency.
An explicit terminal/missing PBS handle permits retry of unfinished assigned
units; transient query errors and timeouts do not. An ambiguous submission
stops dispatch with its unique job name recorded, preventing duplicate launch.
The controller preserves every prior job outcome. `retry-blockers.csv` reports
unfinished units with exhausted start/infrastructure budgets, and method units
blocked by an exhausted generation or original-refit prerequisite. For paired
attribution, every corrected refit waits for its matching original refit even
when a method adapter would not require cached-output reuse. Independent ready
work continues; when no live or eligible work remains, exhausted work causes an
explicit incomplete-campaign error instead of an indefinite poll. The CSV and
bounded status counts identify the unresolved work. Runtime/timeout failures
still require inspection; exhausting retries never constitutes a finished benchmark.

## Prespecified priority order and planning limits

`build_priority_schedule.R` creates an external scheduling prefix from the
frozen design. It reads no predictions. The full Gaussian family runs first,
followed by a paired legacy audit: all 45 Iteration-004 rows, all five S levels,
phi .05/.2/.6 and iterations 1/2. The control audit retains sparse and
sparse-inf fractions .5/.9, no annotations and binary/continuous enrichment
5.4, all three LD references, all five S levels, phi .0075/.1/.4 and iterations
1/2. Every selected scenario retains all ten loci. Any overlapping packed task
is included whole, and boundary extras appear in the accompanying audit CSV.

For the default full plan this prefix contains 60 Gaussian generation units
and 72 legacy units, 32,324 method task units, 268,400 Gaussian method attempts
and 1,296,000 paired legacy locus fits. The remaining 616,896 method tasks and
108 generation units remain in the plan and are scheduled afterward. This
order gives earlier evidence across mechanisms and methods without changing
the full campaign's intended denominator or selecting tasks by performance.

```sh
Rscript --vanilla /absolute/new-controller/scripts/campaign/build_priority_schedule.R \
  --plan /absolute/frozen-full/plan.rds \
  --worker-repo /absolute/unchanged-scientific-source \
  --out /absolute/external-schedules/full-priority.json
python3 /absolute/new-controller/scripts/campaign/dispatch.py \
  --plan /absolute/frozen-full/plan.rds \
  --repo /absolute/unchanged-scientific-source \
  --bootstrap /absolute/frozen-environment.sh \
  --priority-manifest /absolute/external-schedules/full-priority.json \
  --max-workers 48 --max-generation 16 --once --dry-run
```

The builder validates the frozen workers' source and runtime, and the plan RDS
and CSV commitments, even when the scheduling code is newer. Its manifest binds
the exact plan ID and the ordered unique method/generation indices. The
controller records the manifest's byte hash and checks it before submissions;
editing, removing or replacing a schedule on resume is rejected. Ready priority
tasks precede ordinary tasks; blocked preflights, generation and paired-refit
prerequisites still apply. Ordinary ready work can fill otherwise unused slots.
Omitting the manifest preserves the original scheduling order. A newer external
controller cannot authorize a changed scientific worker against an old plan.

Twenty phenotype outcomes per exact legacy p/S/phi setting consist of ten
phenotype draws on each of **two genotype/LD loci**, or **one shared ten-locus
panel** for pooled methods. The draws measure conditional phenotype/effect
variability; they are not twenty independent empirical loci. Two-cluster
bootstrap intervals provide weak evidence for locus generalization, and one
pooled panel cannot identify between-panel uncertainty. Gaussian results have
ten independently generated pooled panels per p/N cell. Template provenance
and held-out empirical loci remain separate generalization questions.

The 27-million-fit legacy workload includes 6.75 million original-grid fits and
20.25 million control fits. It is a reproducibility and sensitivity audit, not
27 million independent replicates. A subsequent confirmatory campaign should
prespecify more independent panels, withheld locus templates, paired LD
references and fewer within-panel phenotype repetitions. Precision stopping
would require independent panel batches, fixed minimum/maximum sample sizes
and sequentially valid intervals; ordinary repeated confidence intervals do
not justify stopping when a favored ranking appears.

Canary estimates are resource anchors. The four large-p legacy arms do not
measure every S/phi/reference/annotation cell; extrapolating their cost to the
smaller controls is not a runtime bound. Adapter runtime omits orchestration,
loading, validation, evaluation and serialization; end-to-end checkpoint times
and file counts should also be inspected. The 220-task p1000 probe measures the
largest Gaussian pooled batch, while small-p exact references have their own
costs. Use the actual regional cap (48 in the shared 60-worker allocation),
retain unresolved methods, and measure tail cells before predicting a finish
date. Completion remains determined by sealed outputs for the full plan.
