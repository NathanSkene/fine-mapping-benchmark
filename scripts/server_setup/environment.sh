#!/usr/bin/env bash
# Source in campaign jobs. Dependencies are isolated from other projects.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/rds/general/apps/modules/modulegroups:/rds/easybuild/noarch/apps/modules/all:/etc/modulefiles:/usr/share/modulefiles
module purge 2>/dev/null || true
module load tools/prod R/4.5.2-gfbf-2025b GSL/2.8-GCC-14.3.0 Eigen/3.4.0-GCCcore-14.3.0 HTSlib/1.22.1-GCC-14.3.0 CMake/3.31.8-GCCcore-14.3.0
export FMB_CAMPAIGN_ROOT=${FMB_CAMPAIGN_ROOT:-/rds/general/user/nskene/home/fmbenchmark-campaign-20261010}
export FMB_TOOL_ROOT="$FMB_CAMPAIGN_ROOT/tools"
export R_LIBS_USER="$FMB_CAMPAIGN_ROOT/R-library"
export FMB_PYTHON="$FMB_TOOL_ROOT/python-env/bin/python"
export RETICULATE_PYTHON="$FMB_PYTHON"
export FMB_FINEMAP_PATH=/rds/general/user/nskene/home/fine-mapping-benchmark-review-20261010/finemap-smoke/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64
export FMB_PAINTOR_PATH="$FMB_TOOL_ROOT/PAINTOR_V3.0/PAINTOR"
export LD_LIBRARY_PATH="$FMB_TOOL_ROOT/PAINTOR_V3.0/lib:${LD_LIBRARY_PATH:-}"
export FMB_CAVIAR_PATH="$FMB_TOOL_ROOT/caviar/CAVIAR-C++/CAVIAR"
export FMB_SPARSEPRO_DIR="$FMB_TOOL_ROOT/SparsePro"
export FMB_FINEMAP_INF_DIR="$FMB_TOOL_ROOT/fine-mapping-inf"
export FMB_GCTB_PATH=/rds/general/user/nskene/home/fine-mapping-benchmark-review-20261010/native-smoke-isolated/gctb_2.5.5_Linux/gctb
export RENV_CONFIG_AUTOLOADER_ENABLED=false
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH="/opt/pbs/bin:$PATH"
