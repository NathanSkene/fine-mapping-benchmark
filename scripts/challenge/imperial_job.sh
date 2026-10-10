#!/usr/bin/env bash
#PBS -N fmchallenge
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=00:45:00
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the repository root}"
source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/rds/general/apps/modules/modulegroups:/rds/easybuild/noarch/apps/modules/all:/etc/modulefiles:/usr/share/modulefiles
module purge 2>/dev/null || true
module load tools/prod R/4.5.2-gfbf-2025b
export RENV_CONFIG_AUTOLOADER_ENABLED=false
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p results
Rscript --vanilla scripts/challenge/test_core.R > results/server-invariants.log 2>&1
Rscript --vanilla scripts/challenge/run_demo.R --replicates 50 --variants 40 --out results/server-development-v1 > results/server-demo.log 2>&1
printf 'Completed scientific invariants and development run\n'
