#!/usr/bin/env bash
#PBS -N fm-native
#PBS -l select=1:ncpus=2:mem=1gb
#PBS -l walltime=00:10:00
set -euo pipefail
cd "${PBS_O_WORKDIR:?Submit from the separate finemap-smoke directory}"
source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/rds/general/apps/modules/modulegroups:/rds/easybuild/noarch/apps/modules/all:/etc/modulefiles:/usr/share/modulefiles
module purge 2>/dev/null || true
module load tools/prod R/4.5.2-gfbf-2025b
export RENV_CONFIG_AUTOLOADER_ENABLED=false
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
binary="$PWD/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64"
Rscript --vanilla native_finemap_smoke.R --binary "$binary" \
  --repo "$PWD/.." --out "$PWD/results" > smoke.log 2>&1
