#!/usr/bin/env bash
#PBS -N fmnative
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=00:05:00
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from isolated native smoke directory}"
source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/rds/general/apps/modules/modulegroups:/rds/easybuild/noarch/apps/modules/all:/etc/modulefiles:/usr/share/modulefiles
module purge 2>/dev/null || true
module load tools/prod R/4.5.2-gfbf-2025b
export RENV_CONFIG_AUTOLOADER_ENABLED=false
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript --vanilla sbayesrc_smoke.R wrapper_sbayesrc_native.R gctb_2.5.5_Linux/gctb smoke-output > native-smoke.log 2>&1
