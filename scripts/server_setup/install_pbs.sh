#!/usr/bin/env bash
#PBS -N fmdeps
#PBS -l select=1:ncpus=4:mem=24gb
#PBS -l walltime=02:00:00
set -euo pipefail
repo=${FMB_BOOTSTRAP_REPO:-/rds/general/user/nskene/home/fmbenchmark-campaign-20261010/bootstrap}
source "$repo/scripts/server_setup/environment.sh"
test -f "$FMB_CAMPAIGN_ROOT/downloads/DOWNLOAD_COMPLETE"
mkdir -p "$FMB_CAMPAIGN_ROOT/build"
export MAKEFLAGS=-j4
# Record compiler/library compatibility adjustments separately from inference.
sed -i 's/^CXX_STD = CXX11/CXX_STD = CXX17/' "$FMB_TOOL_ROOT/finimom/src/Makevars"
git -C "$FMB_TOOL_ROOT/finimom" diff > "$FMB_CAMPAIGN_ROOT/build/finimom-build.patch"
# CARMA's source repository includes foreign compiled objects and a macOS GSL
# path. Rebuild all objects against the loaded compiler and official GSL.
cat > "$FMB_TOOL_ROOT/CARMA/src/Makevars" <<'MAKEVARS'
CXX_STD = CXX17
PKG_CPPFLAGS = $(shell gsl-config --cflags)
PKG_LIBS = $(shell gsl-config --libs)
MAKEVARS
git -C "$FMB_TOOL_ROOT/CARMA" diff > "$FMB_CAMPAIGN_ROOT/build/CARMA-build.patch"
Rscript --vanilla "$repo/scripts/server_setup/install.R" > "$FMB_CAMPAIGN_ROOT/logs/install-R.log" 2>&1
"$FMB_PYTHON" -m pip install --no-index --find-links "$FMB_CAMPAIGN_ROOT/downloads/python" \
  'numpy==2.2.6' 'scipy==1.15.3' 'pandas==2.2.3' 'setuptools==80.9.0' wheel \
  matplotlib seaborn scikit-learn tqdm cython optuna gymnasium build absl-py imageio 'torch==2.7.1+cpu'
"$FMB_PYTHON" -m pip install --no-index --no-build-isolation --no-deps "$FMB_CAMPAIGN_ROOT"/downloads/python/bgzip-*.tar.gz
"$FMB_PYTHON" -m pip install --no-index --no-build-isolation --no-deps "$FMB_TOOL_ROOT/Funmap"
"$FMB_PYTHON" -m pip install --no-index --no-build-isolation --no-deps "$FMB_TOOL_ROOT/fine-mapping-inf/susieinf"
"$FMB_PYTHON" -m pip install --no-index --no-build-isolation --no-deps "$FMB_TOOL_ROOT/fine-mapping-inf/finemapinf"
cd "$FMB_TOOL_ROOT/caviar/CAVIAR-C++"
g++ -std=gnu++11 caviar.cpp PostCal.cpp Util.cpp TopKSNP.cpp \
  -I armadillo/include/ -DARMA_DONT_USE_WRAPPER \
  -L"$EBROOTFLEXIBLAS/lib" -lflexiblas \
  -L"$EBROOTGSL/lib" -lgsl -lgslcblas -o CAVIAR
cd "$FMB_TOOL_ROOT/PAINTOR_V3.0"
# Upstream install.sh builds this bundled library before invoking make.
cd nlopt-2.4.2
./configure --prefix="$FMB_TOOL_ROOT/PAINTOR_V3.0"
make -j4
make install
cd ..
make -j4
"$FMB_PYTHON" -c 'import absl, imageio, numpy, scipy, torch, pandas, finemapinf, susieinf, funmap; print("Scientific Python imports verified")'
"$FMB_PYTHON" -m pip freeze > "$FMB_CAMPAIGN_ROOT/python-versions.txt"
touch "$FMB_CAMPAIGN_ROOT/INSTALL_COMPLETE" "$FMB_CAMPAIGN_ROOT/PYTHON_IMPORTS_COMPLETE"
