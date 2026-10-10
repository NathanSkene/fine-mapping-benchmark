#!/usr/bin/env bash
#PBS -N fmplans
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -l walltime=02:00:00
# Freeze plans on a compute node after the immutable repository is uploaded.
set -euo pipefail
repo=${FMB_REPO_SOURCE:?An immutable uploaded repository is required}
source "$repo/scripts/server_setup/environment.sh"
test -f "$FMB_CAMPAIGN_ROOT/INSTALL_COMPLETE"
test -f "$FMB_CAMPAIGN_ROOT/PYTHON_IMPORTS_COMPLETE"
test -f "$FMB_CAMPAIGN_ROOT/data/DATA_COMPLETE"
export FMB_SOURCE_COMMIT=${FMB_SOURCE_COMMIT:?Pin the reviewed source commit}
campaign=${FMB_RUN_ROOT:?An unused campaign output root is required}
mkdir -p "$campaign"
cd "$repo"
for profile in integration resource_probe full; do
  out="$campaign/regional-$profile"
  if [[ -f "$out/plan.rds" ]]; then
    Rscript --vanilla scripts/server_setup/validate_plan.R --plan "$out/plan.rds" --kind regional
    printf 'Immutable existing plan retained: %s\n' "$out"
    continue
  fi
  Rscript --vanilla scripts/campaign/create_plan.R --profile "$profile" --out "$out" \
    --original-source "$FMB_CAMPAIGN_ROOT/original-f053" \
    --vcf-dir "$FMB_CAMPAIGN_ROOT/data/legacy-wide-vcf" \
    --genetic-map-dir "$FMB_CAMPAIGN_ROOT/data/genetic_maps" \
    --inputs-manifest "$FMB_CAMPAIGN_ROOT/data/inputs-manifest.csv"
  Rscript --vanilla scripts/server_setup/validate_plan.R --plan "$out/plan.rds" --kind regional
done
for profile in integration canary full; do
  out="$campaign/native-$profile"
  if [[ -f "$out/plan.rds" ]]; then
    Rscript --vanilla scripts/server_setup/validate_plan.R --plan "$out/plan.rds" --kind native
    printf 'Immutable existing plan retained: %s\n' "$out"
    continue
  fi
  Rscript --vanilla scripts/native/campaign/create_plan.R --profile "$profile" --out "$out" \
    --vcf-dir "$FMB_CAMPAIGN_ROOT/data/native-eur-vcf" \
    --regions-file "$FMB_CAMPAIGN_ROOT/data/native-regions.csv" \
    --genetic-map-dir "$FMB_CAMPAIGN_ROOT/data/genetic_maps" \
    --gctb-path "$FMB_GCTB_PATH"
  Rscript --vanilla scripts/server_setup/validate_plan.R --plan "$out/plan.rds" --kind native
done
touch "$campaign/PLANS_COMPLETE"
