#!/usr/bin/env Rscript
# A plan written before an interrupted CSV export is not a complete plan.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
source(file.path(repo, "scripts/campaign/common.R"))
opt <- .campaign_options(list(plan = "", kind = "regional"))
plan <- readRDS(opt$plan)
root <- dirname(normalizePath(opt$plan))
if (opt$kind == "regional") {
  .campaign_assert_plan(plan, repo)
  files <- c("tasks.csv", "datasets.csv", "generation.csv")
} else if (opt$kind == "native") {
  source(file.path(repo, "scripts/native/campaign/helpers.R"))
  .native_assert_plan(plan, repo, verify_source_data = TRUE)
  if (!identical(plan$executable_sha256, .challenge_sha256(plan$executable)))
    stop("Pinned native executable changed")
  files <- c("draws.csv", "inputs.csv", "chains.csv", "regional.csv")
} else stop("Unknown plan kind")
if (!identical(plan$runtime_fingerprint, .campaign_runtime_fingerprint(plan$method_args)))
  stop("Frozen runtime dependencies changed")
manifest <- jsonlite::fromJSON(file.path(root, "dispatch-manifest.json"), simplifyVector = FALSE)
if (!identical(manifest$schema, 1L) ||
    !identical(manifest$plan_sha256, plan$plan_sha256) ||
    !identical(readLines(file.path(root, "plan.sha256"), warn = FALSE), plan$plan_sha256) ||
    !identical(manifest$plan_rds_sha256, .challenge_sha256(opt$plan)) ||
    !setequal(names(manifest$files), files)) stop("Incomplete or changed plan commitment")
for (name in files) if (!identical(manifest$files[[name]], .challenge_sha256(file.path(root, name))))
  stop("Changed plan projection: ", name)
cat("Validated complete frozen", opt$kind, "plan:", root, "\n")
