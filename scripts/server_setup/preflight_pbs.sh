#!/usr/bin/env bash
#PBS -N fmpreflight
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -l walltime=04:00:00
set -euo pipefail
repo=${FMB_REPO_SOURCE:?An immutable uploaded repository is required}
source "$repo/scripts/server_setup/environment.sh"
cd "$repo"
Rscript --vanilla scripts/campaign/preflight.R \
  --plan "${FMB_PREFLIGHT_PLAN:?A frozen regional plan is required}" \
  --index "${FMB_METHOD_INDEX:?One method index is required}"
