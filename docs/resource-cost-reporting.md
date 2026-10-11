# Resource evidence and conditional campaign cost reports

`scripts/reporting/resource_cost_report.R` is an analysis utility. It reads the
unchanged frozen campaign source and sealed resource-probe evidence, plus any
explicitly selected full-plan runtime tails. It does not launch jobs, call a
fine-mapping method, regenerate inputs, rescore predictions, rewrite receipts or
change the campaign plan. The existing `scripts/campaign/cost_report.R` and its
reports remain preserved.

This separate route matters because the original empirical generator requests
a nominal variant count and then drops columns that are monomorphic in the
simulated GWAS sample or independent reference sample. A nominal 2,000-variant
locus can legitimately retain 1,998 variants. Requiring an observed maximum of
exactly 2,000 would reject that intended experiment. Conversely, a very quick
failed fit, an unsupported-input rejection or a small installation run cannot
establish that full-size inference is feasible.

## Immutable inputs and execution

Run the utility in the frozen runtime used to create the plans. Its explicit
`--source` must be the exact deployed source directory. The campaign's
`common.R` derives its repository from the top-level script path and does not
honour `FMB_REPO`, so this external utility loads the unchanged `R/` helpers into
an isolated environment using the supplied source path. It checks every listed
source byte against the frozen snapshot before loading helpers, then uses the
unchanged `.campaign_assert_plan`, runtime, checkpoint, receipt, bundle and
artifact validators. It validates the source, runtime and all tracked input
bytes again before sealing its completion manifest.
Its pre-source guard also requires the exact listed `R/` file inventory: an
extra unlisted helper is rejected before it can execute. A sealed method record
requires its own compatible sealed generation and input artifacts, for both
empirical and Gaussian inputs.
Probe and full plans must also share the exact method settings and original
generator byte snapshot. Identical installed libraries/executables alone do
not bind iteration counts, search limits, tolerances or priors.

The default read bound is 512 MiB per serialized plan, checkpoint or tail-selection
metadata file. A tail selection's size is checked before JSON parsing; at most
1,000 distinct selector files are admitted. The utility
loads only one selected task at a time; it never accumulates all fit archives.
Its full-plan pass reads design metadata and counts every planned worker,
dataset record and locus attempt, including unmeasured exact/forced-one workers.
The default bounds are 10,000 measured tasks, 100,000 output design groups,
1,000,000 planned tasks, 100,000 planned datasets and 1,000 generation units.
Increase a command-line read bound only after reviewing an actual artifact's
size. Exceeding a bound stops reporting rather than silently omitting evidence.

Choose a new output directory with a versioned `resource-cost-...` basename and
an existing parent outside the source and campaign roots:

```sh
Rscript --vanilla /path/to/reviewed/resource_cost_report.R \
  --source /path/to/source-0e4a7ead774b \
  --probe-plan /path/to/regional-resource_probe/plan.rds \
  --full-plan /path/to/regional-full/plan.rds \
  --tail-selections /path/to/regional-runtime-tail-selection.json,/path/to/regional-runtime-tail-finite-ld-addendum.json \
  --pbs-evidence /path/to/resource-pbs-bindings-v1.json \
  --workers 48 \
  --out /path/to/reports/resource-cost-20261011-v1
```

`--tail-selections` and `--pbs-evidence` are optional. Missing evidence is
reported as unresolved. A tail selection must bind the exact full-plan RDS,
plan digest, dispatch projection and legacy grid; an addendum must be supplied
with its exact parent selection. All selected datasets must retain their
declared generation dependencies. An original-scenario source hash is checked
when the selection contains it. No source hash, plan hash or old report is
updated to make incompatible evidence fit.

If critical custody fails, a partial output directory may remain for diagnosis.
It has no completed report manifest and must not be treated as a successful
report. Choose a fresh output directory after fixing the identified problem.

## What the report measures

The required empirical resource probe has ten loci, nominal requested sizes
`500|500|750|750|1000|1000|1500|1500|2000|2000` and 5,000 simulated individuals.
Its sealed `generation_dimension_audit` supplies the requested, realized and
dropped counts, per-locus simulated GWAS/reference Ns, fractions, source-file
hash, filter basis and preserved warnings. Method tables must match those exact
locus dimensions and the frozen scenario design.

Dimension loss is accepted only under the known original generator hash and
its explicit post-simulation monomorphic-column filter. An unknown filter
cannot explain loss. VCF undersupply after the pre-simulation MAF filter, and
the bundled example's 500-variant cap, remain errors even for the known
generator. The utility exports the entire realized vector and the ratios
`sum(realized p^k) / sum(requested p^k)` for `k = 1, 2, 3`. These describe size
sensitivity. They are not inferred runtime complexity, runtime bounds or memory
bounds for a method.

The required Gaussian resource probe is a 220-attempt, nominal p=1,000,
N=10,000 batch: ten replicates of each of 22 mechanisms. The missing-causal arm
intentionally removes one candidate, so its realized p is 999. That reduction
is reported separately and must appear only in that mechanism. Gaussian N is
a direct-score design parameter, not a count of observed individual phenotype
rows; variant-wise N and multicohort mechanisms remain separately labelled.
This utility selects no disease GWAS or published cohort release. Restricted
cohort and 23andMe membership are not applicable to these synthetic inputs.

Each locus remains in `locus-attempts.csv`, with its unchanged recorded runtime,
success, error and one of three execution classes:

- successful execution;
- a specific declared unsupported-input rejection;
- an unexpected execution failure.

Known unsupported cases use the exact method/input/error combination:
regional scalar-N methods on multicohort or variable-N Gaussian inputs;
Gaussian SuSiE-inf and FINEMAP-inf likelihood requirements; Funmap without
annotations; and exact/forced-one on variant-wise N without overlap metadata.
An unfamiliar error remains unexpected. Finite failed/rejected runtimes are
retained and counted. They do not qualify an inference cost anchor.

A qualifying new inference anchor requires successful execution at every
applicable largest nominal locus and no unexpected failures in the input.
An all-rejected input, a failed largest locus or an exactly reused fit cannot
clear that condition. Reused output remains in the measured census. Passing
this execution condition does not establish convergence, calibration or a
scientific winner.

## Strata and conditional end-to-end envelopes

The report keeps stage, nominal locus-size vector, N, reference N, architecture,
annotation, enrichment, annotation correlation, sparse/background proportion,
S pack, phi pack, replicate count and worker packing distinct. N=1,000 tails
and reference sizes of 500 or 200 are therefore valid separate strata.
Non-CARMA tail workers pack S=1/2/3/5/10; the selected CARMA worker contains only
S=10. Those workers are not pooled. Separate genotype realizations in the
in-sample and reference-panel tails do not constitute an isolated paired
LD-only intervention.

The end-to-end envelope uses only complete fresh-worker timing that reconciles
R script elapsed time with data processing and prelude/input-I/O/compaction
overhead. Its input must also pass the execution condition for every included
dataset. For an exactly matching nominal worker design, observed minimum,
median and maximum script times are multiplied by that design's planned task
count. These are conditional sampled anchors. They are neither confidence
intervals nor worst-case bounds. Mixed S batches retain their actual whole
worker time; it is not allocated to individual S values.

Unmeasured exact worker designs retain `NA` costs and an explicit reason. An S=2
probe does not establish the cost of an S=1/2/3/5/10 batch. No Gaussian stress
samples cannot become a zero-cost estimate. Exact/forced-one remain limited to
p<=200 and need their own smaller-p measurements; they cannot be presented as
measured p=1,000 stress inference. The report's aggregate hours cover only its
qualified measured subset. `full_R_cost_established` and
`inference_feasibility_cleared` remain false.

The R timing scope excludes process startup before the first script expression,
the final timing receipt write, lock cleanup, scheduler queue time and earlier
failed PBS attempts. Input generation and its retained bytes are reported
separately. Primary task-retention sizes include the sealed started-attempt
ledger; storage extrapolations exclude scheduler scripts/logs, old archives and
cache deduplication. They cannot be described as total required storage.
Gaussian input retention includes `panel.rds`, both release manifests and all
four public/assessor child payloads. The utility checks each bounded release
manifest's exact filename/hash schema, rejects unlisted or escaping paths, and
hashes/stats child files without deserializing assessor truth or generator
payloads. Generation checkpoint/receipt bytes are a separate column. These are
logical regular-file sizes; filesystem blocks and directory metadata remain
outside the measurement.

## Actual PBS resource evidence

The optional PBS file is a small mapping to preserved raw terminal `qstat -xf
-F json` output and submitted scripts. It must use this schema:

```json
{
  "schema": 1,
  "kind": "fmbenchmark-resource-pbs-bindings-v1",
  "attempts": [
    {
      "plan_sha256": "the frozen 64-character plan digest",
      "kind": "methods",
      "indices": [1, 2],
      "job_id": "1234567.pbs-7",
      "script_path": "/absolute/submitted/job-script.sh",
      "script_sha256": "the actual 64-character submitted-script hash",
      "qstat_path": "/absolute/preserved/qstat-observation.json",
      "qstat_sha256": "the actual 64-character observation-file hash"
    }
  ]
}
```

The raw observation can be direct PBS JSON or an existing observation wrapper
with `returncode: 0` and JSON in `stdout`. The utility verifies hashes, job ID,
terminal state, exit status, job name, submitted script route, exact frozen
source/plan worker commands and every declared unit index. It accepts terminal
failed jobs and preserves their actual wall, CPU, exit and peak `mem` resources.
PBS `kb` uses 1,024 bytes; requested memory and `vmem` are not substituted for
measured peak `mem`. Batch resource values stay whole-job values and are not
charged once per included method. Unknown quoting/control syntax fails closed;
the reviewed frozen launchers use ordinary paths without shell metacharacters.
Only complete anchored worker command lines using `Rscript` or the current
runtime's Rscript path are admitted. A commented, echoed or compound command
cannot establish that the declared worker ran.

The utility does not authenticate a scheduler response or contact PBS. Its
collector must preserve canonical observations and submitted-script/assignment
provenance. The report binds those supplied bytes and declares this limit. A
missing observation, missing peak measurement, pre-worker failure or unresolved
handle remains unresolved. The default report never claims that all PBS
attempt costs or resources have been established.

## Outputs and validation

Outputs include generation dimensions and costs, all selected task statuses,
all started attempts, locus attempts, measured records/tasks, actual PBS
attempts when supplied, the complete full-plan workload and the conditional
cost envelope. Input and output artifact manifests bind the report to file
bytes. Completed JSON/RDS manifests record scope, bounds and unresolved counts.
Read `conditional-cost-envelope.csv` alongside `task-status.csv` and
`generation-status.csv`; an observed R timing alone is not a sufficient gate.

Run the meaningful pure/fixture checks with the existing dependencies:

```sh
RENV_CONFIG_AUTOLOADER_ENABLED=false R_LIBS_USER=.local-lib \
  Rscript --vanilla scripts/reporting/test_resource_cost_report.R
```

The tests use disposable synthetic metadata, source hashes and sealed receipts.
They cover valid filtered dimensions, unexplained loss/undersupply, contradictory
audit/table metadata, missing stress, cheap failed/rejected inference, largest
locus success, reference/S packing distinctions, whole PBS job resources and
unchanged source/plan/receipt bytes. They also reject otherwise valid resealed
plans with different scalar method settings or original-source lineage,
unlisted source helpers before execution, missing Gaussian generation custody,
changed opaque release payloads, and commented/echoed/compound PBS worker
commands. They run no inference or scheduler command.
