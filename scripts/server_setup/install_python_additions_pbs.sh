#!/usr/bin/env bash
#PBS -N fmpydeps
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -l walltime=00:30:00
# Add omissions found during import audit without recompiling finished R tools.
set -euo pipefail
repo=${FMB_BOOTSTRAP_REPO:-/rds/general/user/nskene/home/fmbenchmark-campaign-20261010/bootstrap}
source "$repo/scripts/server_setup/environment.sh"
test -f "$FMB_CAMPAIGN_ROOT/INSTALL_COMPLETE"
# Historical adapters normalizePath() the final executable, which otherwise
# follows a venv symlink out to the base interpreter. Preserve their source
# while using the same isolated packages: a venv-local physical interpreter.
if [[ -L "$FMB_PYTHON" ]]; then
  target=$(readlink -f "$FMB_PYTHON")
  replacement="$FMB_PYTHON.physical.$$"
  cp -- "$target" "$replacement"
  chmod --reference="$target" "$replacement"
  mv -Tf "$replacement" "$FMB_PYTHON"
fi
"$FMB_PYTHON" -m pip install --no-index --find-links "$FMB_CAMPAIGN_ROOT/downloads/python" absl-py imageio
"$FMB_PYTHON" -c 'import absl, imageio, numpy, scipy, torch, pandas, sys; print("Scientific Python imports verified", sys.executable, sys.prefix)'
Rscript --vanilla -e 'p <- Sys.getenv("FMB_PYTHON"); stopifnot(identical(normalizePath(p), p)); cat("Historical normalized interpreter remains inside the isolated venv\n")'
"$FMB_PYTHON" -m pip freeze > "$FMB_CAMPAIGN_ROOT/python-versions.txt"
touch "$FMB_CAMPAIGN_ROOT/PYTHON_IMPORTS_COMPLETE"
