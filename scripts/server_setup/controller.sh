#!/usr/bin/env bash
# Lightweight login-node orchestration only; all R computation is in PBS jobs.
set -euo pipefail
repo=${1:?Immutable repository required}
kind=${2:?regional or native required}
plan=${3:?Frozen plan required}
source "$repo/scripts/server_setup/environment.sh"
cd "$repo"
case "$kind" in
  regional)
    exec "$FMB_PYTHON" scripts/campaign/dispatch.py --repo "$repo" --plan "$plan" \
      --bootstrap "$repo/scripts/server_setup/environment.sh" \
      --max-workers 48 --max-generation 8 --batch-tasks 4 --poll-seconds 60
    ;;
  native)
    preflight=${4:?Frozen regional preflight plan required}
    exec "$FMB_PYTHON" scripts/server_setup/native_dispatch.py --repo "$repo" --plan "$plan" \
      --bootstrap "$repo/scripts/server_setup/environment.sh" --preflight-plan "$preflight" \
      --max-workers 12 --max-prepare 2 --max-chains 6 --max-regional 4 --poll-seconds 60
    ;;
  *) printf 'Unknown controller kind: %s\n' "$kind" >&2; exit 2 ;;
esac
