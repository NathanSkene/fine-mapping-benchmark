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

Run the controllers through `nohup controller.sh ...` only after the relevant
gates pass. Regional and native controllers are capped at48 and12 jobs,
respectively: at most60 campaign jobs together, one CPU each. Native
preparation, chain and regional jobs also have limits of2,6 and4. Shared
preflight and planning jobs finish before the full controllers start. Canary
measurements determine resource changes before the final launch; keep the
combined bound moderate if changing any limit. These are limits on submitted
and running assignments, not a request to occupy a fixed number of nodes.

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

Run both integration campaigns and the full-size resource probes before
scientific interpretation. Estimate retained output bytes and compute from the
measured canaries. The individual ephemeral allocation currently provides
roughly11 TB; check current quota before launch, and keep compact final metrics,
provenance and paper artifacts in home storage. Ephemeral storage is not backed
up. Complete the collectors and their planned/attempted/successful census before
writing winners or calibration claims into the manuscript.
