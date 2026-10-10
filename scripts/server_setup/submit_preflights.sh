#!/usr/bin/env bash
# Lightweight login orchestration; every actual fit executes through PBS.
# An existing intent ledger requires inspection, never an automatic resubmit.
set -euo pipefail
repo=${1:?Immutable source directory required}
plan=${2:?Completed frozen regional plan required}
prefix=${3:?Unique job-name prefix required}
source "$repo/scripts/server_setup/environment.sh"
test -f "$plan"
root=$(dirname "$plan")
test -f "$root/dispatch-manifest.json"
test -f "$root/plan.sha256"
method_count=$("$FMB_PYTHON" - "$root/tasks.csv" <<'PY'
import csv
import sys
with open(sys.argv[1], newline="") as handle:
    methods = {row["method"] for row in csv.DictReader(handle)
               if row["stage"] != "challenge"}
print(len(methods))
PY
)
if [[ "$method_count" != 20 ]]; then
  printf 'Expected the standard 20-method server plan; found %s\n' "$method_count" >&2
  exit 1
fi
record="$root/preflight-submissions.tsv"
test ! -e "$record"
printf 'index\tjob_name\tstate\tjob_id\n' > "$record"
# The server plans register the complete ordered set of 20 methods. Never
# apply this helper to a custom method subset with different indices.
for index in $(seq 1 20); do
  name=$(printf '%s%02d' "$prefix" "$index")
  printf '%s\t%s\tpending\t\n' "$index" "$name" >> "$record"
  job=$(qsub -N "$name" \
    -o "$FMB_CAMPAIGN_ROOT/logs/$name.pbs.out" \
    -e "$FMB_CAMPAIGN_ROOT/logs/$name.pbs.err" \
    -v "FMB_REPO_SOURCE=$repo,FMB_PREFLIGHT_PLAN=$plan,FMB_METHOD_INDEX=$index" \
    "$repo/scripts/server_setup/preflight_pbs.sh")
  if [[ ! "$job" =~ ^[0-9]+([.][A-Za-z0-9._-]+)?$ ]]; then
    printf 'Unconfirmed submission for %s; inspect qselect -N %s\n' "$name" "$name" >&2
    exit 1
  fi
  printf '%s\t%s\tconfirmed\t%s\n' "$index" "$name" "$job" >> "$record"
  printf 'Method %s submitted as %s\n' "$index" "$job"
done
