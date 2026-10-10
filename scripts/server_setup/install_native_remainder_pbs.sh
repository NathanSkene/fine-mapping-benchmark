#!/usr/bin/env bash
#PBS -N fmnativebuild
#PBS -l select=1:ncpus=4:mem=12gb
#PBS -l walltime=01:00:00
# Resume the native/Python remainder after the completed isolated R install.
set -euo pipefail
repo=${FMB_BOOTSTRAP_REPO:-/rds/general/user/nskene/home/fmbenchmark-campaign-20261010/bootstrap}
source "$repo/scripts/server_setup/environment.sh"
test -x "$FMB_CAVIAR_PATH"
export MAKEFLAGS=-j4
cd "$FMB_TOOL_ROOT/PAINTOR_V3.0/nlopt-2.4.2"
./configure --prefix="$FMB_TOOL_ROOT/PAINTOR_V3.0"
make -j4
make install
cd ..
make -j4
test -x "$FMB_PAINTOR_PATH"
"$FMB_PYTHON" -m pip install --no-index --find-links "$FMB_CAMPAIGN_ROOT/downloads/python" absl-py imageio
"$FMB_PYTHON" -c 'import absl, imageio, numpy, scipy, torch, pandas, finemapinf, susieinf, funmap; print("Scientific Python imports verified")'
"$FMB_PYTHON" -m pip freeze > "$FMB_CAMPAIGN_ROOT/python-versions.txt"
touch "$FMB_CAMPAIGN_ROOT/INSTALL_COMPLETE" "$FMB_CAMPAIGN_ROOT/PYTHON_IMPORTS_COMPLETE"
