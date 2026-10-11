#!/usr/bin/env bash
# Run through nohup on a login node (network); installation is a separate PBS job.
set -euo pipefail
repo=${1:?Repository snapshot required}
source "$repo/scripts/server_setup/environment.sh"
mkdir -p "$FMB_TOOL_ROOT" "$FMB_CAMPAIGN_ROOT/downloads/python" "$R_LIBS_USER"
clone() {
  local name=$1 url=$2
  if [[ ! -d "$FMB_TOOL_ROOT/$name/.git" ]]; then git clone --depth 1 "$url" "$FMB_TOOL_ROOT/$name"; fi
}
clone CARMA https://github.com/ZikunY/CARMA.git
clone finimom https://github.com/vkarhune/finimom.git
clone SparsePro https://github.com/zhwm/SparsePro.git
clone fine-mapping-inf https://github.com/FinucaneLab/fine-mapping-inf.git
clone Funmap https://github.com/LeeHITsz/Funmap.git
clone Beatrice-Finemapping https://github.com/sayangsep/Beatrice-Finemapping.git
clone PAINTOR_V3.0 https://github.com/gkichaev/PAINTOR_V3.0.git
clone caviar https://github.com/fhormoz/caviar.git
clone susieR https://github.com/stephenslab/susieR.git
git -C "$FMB_TOOL_ROOT/susieR" fetch --depth 1 origin 732181c662d118930e4381a5d6dc44c89ce2cea0
git -C "$FMB_TOOL_ROOT/susieR" checkout --detach 732181c662d118930e4381a5d6dc44c89ce2cea0
Rscript --vanilla "$repo/scripts/server_setup/download_r_dependencies.R"
python -m venv --copies "$FMB_TOOL_ROOT/python-env"
"$FMB_PYTHON" -m pip download --only-binary=:all: --dest "$FMB_CAMPAIGN_ROOT/downloads/python" \
  'numpy==2.2.6' 'scipy==1.15.3' 'pandas==2.2.3' 'setuptools==80.9.0' wheel \
  matplotlib seaborn scikit-learn tqdm cython optuna gymnasium build absl-py imageio
"$FMB_PYTHON" -m pip download --no-deps --dest "$FMB_CAMPAIGN_ROOT/downloads/python" bgzip
"$FMB_PYTHON" -m pip download --only-binary=:all: --dest "$FMB_CAMPAIGN_ROOT/downloads/python" \
  --index-url https://download.pytorch.org/whl/cpu 'torch==2.7.1+cpu'
for dir in "$FMB_TOOL_ROOT"/*; do
  if [[ -d "$dir/.git" ]]; then printf '%s\t%s\n' "$(basename "$dir")" "$(git -C "$dir" rev-parse HEAD)"; fi
done > "$FMB_CAMPAIGN_ROOT/downloads/upstream-commits.tsv"
find "$FMB_CAMPAIGN_ROOT/downloads" -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > "$FMB_CAMPAIGN_ROOT/downloads/SHA256SUMS"
touch "$FMB_CAMPAIGN_ROOT/downloads/DOWNLOAD_COMPLETE"
