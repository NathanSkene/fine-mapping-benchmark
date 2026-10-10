# Full paired campaign

This runner preserves the complete original Iteration-004 grid and adds
Iteration-003-compatible controls and a separate stress family. A `full` plan
contains 45 original-design rows plus 135 control rows, 250 scenario draws per
row and ten loci per scenario. The control rows restore sample size 1,000,
variant counts 100–1,000, heritabilities .0075–.4, and independent LD references
of 200 and 500 individuals. The latest design retains sample size 5,000 and
500–2,000 variants. Actual variant counts must equal the requested counts;
truncating an undersized VCF window fails generation.

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
The `canary` profile preserves the full 500–2,000-variant, 5,000-individual locus
sizes while selecting four sparse/background and annotation arms with one
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

After real full-size canaries finish, `cost_report.R --canary-plan PATH
--full-plan PATH --workers 64 --out PATH` measures runtimes and compressed
checkpoint/cache bytes. A `resource_probe` plan also supplies stress costs from
the full 220-task batch. Largest-panel estimates extrapolated to smaller stress
dimensions are planning anchors, not measured runtimes at every dimension.
Its estimates disclose unresolved methods and exclude input archives, queue
delay and exact-fit reuse savings. The runner
currently retains complete adapter results, including available posterior means,
convergence and component diagnostics; plan storage from measured bytes before
dispatching the full campaign.

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
controller. Each submitted job ID and exact task assignment is recorded before
observation. The dispatch manifest commits to the plan RDS bytes and the task,
dataset and generation CSV projections; editing a CSV cannot silently change
the controller's interpretation of the plan. Terminal receipts commit to the R
checkpoint bytes and plan ID.
The controller revalidates all terminal receipts once at startup, including
previously cached completed IDs. Dataset predictions and reusable fit caches
have their own byte-hash receipts; a completed method task additionally seals
its full dataset/receipt ledger. A missing receipt is pending work. An invalid
receipt or changed sealed artifact is an integrity error. Unsealed dataset RDS
files are preserved outside the collection namespace and recomputed on resume.
Collection holds its own lock, archives the previous derived analysis, and
builds a fresh analysis with a sealed manifest of the plan, source, runtime,
completion status and derived-file hashes. An empty or failed collection cannot
leave previous scientific tables in the current analysis directory.

Runtime guards include installed package versions, external Python repository
source hashes, configured executable/interpreter bytes and threading settings.
They do not hash installed R libraries or installed Python package bytes.
Same-version reinstalls can evade that version guard: keep the installed
environment immutable for the campaign and retain its installation manifests.
This is not a byte-complete commitment to every installed dependency.
An explicit terminal/missing PBS handle permits retry of unfinished assigned
units; transient query errors and timeouts do not. An ambiguous submission
stops dispatch with its unique job name recorded, preventing duplicate launch.
The controller caps automatic attempts at three while preserving every prior
job outcome. Runtime/timeout failures still require inspection; exhausting
retries is incomplete work, not a finished benchmark.
