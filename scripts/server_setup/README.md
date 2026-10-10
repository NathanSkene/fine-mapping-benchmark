# Imperial campaign operation

These scripts prepare isolated dependencies and execute the frozen paired
campaign through PBS. Downloads run on the login node; compilation, R planning,
generation and inference run on scheduled compute nodes. The SSH helper is
documented in Nathan's Imperial HPC skill. No other user's jobs are modified.

`environment.sh` specifies the verified module versions, isolated R/Python
libraries, external executables and single-thread BLAS settings. Run
`download.sh` and `download_data.sh` before `install_pbs.sh`. Compute nodes have
no internet, so installation consumes the downloaded source archives/wheels.
The recorded finimom C++17 adjustment and CARMA C++17/GSL-path adjustment are
compiler compatibility patches and remain separate from inference code. A clean
source build removes CARMA's bundled foreign-platform objects. Downloaded repository commits,
dependency versions and input bytes are recorded before plan creation.
PAINTOR's upstream installation builds its bundled NLopt 2.4.2 before the
PAINTOR executable; the runtime environment exposes that isolated library.
The Python environment uses a physical interpreter copy. Historical adapters
resolve executable symlinks with `normalizePath()`, which can otherwise escape
a virtual environment and lose its dependencies. This installation adjustment
keeps the unchanged historical source in the same isolated Python environment.

The legacy comparison uses new widened 2 Mb ALL-population source windows to
preserve the requested variant counts. The separate native campaign uses EUR
templates from 128 source blocks. Empty native windows shift forward in 1 Mb
steps, staying before the next source block. `native-regions.csv` and
`download-regions.tsv` record actual coordinates, original starts, sample
population and source URLs. These files are new input realizations; they do not
recover Lucas's historical predictions or data.

Upload a reviewed, immutable source snapshot to a new directory. Submit
`plans_pbs.sh` with `FMB_REPO_SOURCE`, `FMB_SOURCE_COMMIT` and `FMB_RUN_ROOT`
set explicitly. It creates small integration plans, resource probes and full
plans; creation is not execution. Never overwrite a plan or change its source
snapshot. Regional preflight jobs use `preflight_pbs.sh` and one method index,
requiring an actual successful fit on both annotated loci. A failed preflight
is retained; repair the cause and use a new plan.

For the standard server plans with all 20 methods, `submit_preflights.sh SOURCE
PLAN PREFIX` submits one job per method and records intent before each `qsub`.
Use a unique PBS name prefix. It refuses an existing submission ledger; an
interrupted or ambiguous submission requires inspecting the ledger and scheduler
before any further action. This preserves actual failed runs and prevents an
automatic rerun from duplicating jobs. Custom method subsets require their own
explicit indices.

Plan creation also freezes `regional-full/priority-prefix.json` from design
identities alone. The full regional controller uses this prefix to run the
Gaussian stress family and the paired audit earlier; every remaining task
stays in the full plan. No predictions enter scheduling selection. Small
integration and resource-probe controllers omit this full-plan-only prefix.

Run the controllers through `nohup controller.sh ...` only after the relevant
gates pass. Regional and native controllers are capped at48 and12 jobs,
respectively: at most60 campaign jobs together, one CPU each. Native
preparation, chain and regional jobs also have limits of2,6 and4. Shared
preflight and planning jobs finish before the full controllers start. Canary
measurements determine resource changes before the final launch; keep the
combined bound moderate if changing any limit. These are limits on submitted
and running assignments, not a request to occupy a fixed number of nodes.

Controllers omit a PBS queue directive by default, letting the scheduler choose
the account's permitted default queue. Pass `--queue QUEUE` to a Python
controller only after verifying access; a named queue can exist while the
account is denied access to it. Preflight and planning submissions likewise use
the scheduler default unless an accessible queue is deliberately selected.

A reviewed controller may run from a newer source directory with `--repo`
pointing to the unchanged worker snapshot used by an existing plan. Generated
jobs source that snapshot's bootstrap and execute its R workers; the controller
records its own file hash separately. The native controller also records the
hash of scheduler helpers loaded from the worker snapshot. This allows a
scheduler-only repair without rewriting the scientific plan or worker source.
For this handoff, invoke the newer Python controller by its absolute path;
`controller.sh` changes directory into its worker repository and chooses that
repository's controller. For example:

```sh
"$FMB_PYTHON" /absolute/new-controller/scripts/campaign/dispatch.py \
  --repo /absolute/unchanged-worker-source --plan /absolute/regional-plan/plan.rds \
  --bootstrap /absolute/unchanged-worker-source/scripts/server_setup/environment.sh \
  --max-workers 48 --max-generation 8 --batch-tasks 4
```

An existing controller state from a different controller hash requires an
explicitly audited handoff; do not delete live scheduler assignments to force a
resume. Resource reports for the regional wrapper should pass `--workers 48`,
matching its actual cap; the combined regional/native cap is 60.

```sh
bash scripts/server_setup/controller.sh /absolute/source regional /absolute/regional-full/plan.rds
bash scripts/server_setup/controller.sh /absolute/source native /absolute/native-full/plan.rds /absolute/regional-full/plan.rds
```

Native regional panels run all methods for one observed input after loading its
shared genotype archive once. Each method keeps its own checkpoint. A killed
panel resumes only incomplete methods. Four native chains remain independent
jobs and retain individual seeds, traces and output files. The controllers
verify plan/projection/bootstrap/checkpoint hashes. An uncertain `qstat` or
`qsub` result does not permit a duplicate execution. Authoritatively terminal
jobs without a receipt can retry up to three attempts, preserving scheduler
outcomes and clearing abandoned locks only after the process has stopped.
For native jobs, `--max-attempts` counts confirmed PBS submissions for the
whole preparation, chain or regional-panel unit, including bootstrap failures
before R starts. If all assigned handles have stopped and its required receipts
are still missing at that limit, the controller records each exhausted kind,
unit ID and missing receipt in `dispatch/state.json` and exits with an error,
including in `--once` mode. Its printed summary lists at most twenty examples;
the state retains the complete report and all prior assignments/observations.
Other live assignments keep their handles and locks. Inspect those records
before explicitly extending the retry budget or supplying valid completed
receipts; exhaustion does not mark any scientific result complete.

Both controllers claim `dispatch/controller.owner` by atomic directory creation
before reading or changing scheduler state. Its `owner.json` records the host,
PID and a unique ownership nonce; any existing claim blocks another controller,
including a directory with missing metadata. The shared scheduler helper keeps
the local advisory lock as a secondary guard. Graceful exit releases only its
own unchanged claim. A killed process can leave a claim behind: inspect its
owner and scheduler assignments before a deliberate manual release; controllers
never automatically clear stale or unknown controller claims. A newer native
controller requires worker scheduler helpers with this shared directory-lock
API and refuses older helpers that provide only advisory locking.

Run both integration campaigns and the full-size resource probes before
scientific interpretation. Estimate retained output bytes and compute from the
measured canaries. The individual ephemeral allocation currently provides
roughly11 TB; check current quota before launch, and keep compact final metrics,
provenance and paper artifacts in home storage. Ephemeral storage is not backed
up. Complete the collectors and their planned/attempted/successful census before
writing winners or calibration claims into the manuscript.
