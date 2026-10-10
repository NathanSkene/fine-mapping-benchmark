#!/usr/bin/env Rscript
# Assessor-side bridge from a trusted generate_challenge() panel to JSON releases.
# Never mount this process/panel/private directory in an agent worker.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: export_challenge.R panel.rds public_dir private_dir")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required")
panel <- readRDS(args[1])
if (!is.list(panel$inputs) || !is.list(panel$truth) ||
    !identical(names(panel$inputs), names(panel$truth))) stop("Invalid trusted panel")
pub <- normalizePath(args[2], mustWork = FALSE)
priv <- normalizePath(args[3], mustWork = FALSE)
if (pub == priv || startsWith(pub, paste0(priv, "/")) || startsWith(priv, paste0(pub, "/")))
  stop("Public/private paths overlap")
if (file.exists(pub) || file.exists(priv)) stop("Refusing to overwrite a release")
eligible <- vapply(panel$inputs, function(x) is.null(x$cohorts) && is.null(x$n_per_variant), TRUE)
if (!any(eligible)) stop("No single-cohort constant-n tasks available")
inputs <- lapply(panel$inputs[eligible], function(x) list(
  task_id = x$task_id, variant_id = unname(x$variant_id), z = unname(x$z),
  LD = unname(x$LD), annotations = unname(x$annotations), n = x$n,
  statistic = "known-noise-score"))
truth <- lapply(panel$truth[eligible], function(x) list(
  task_id = x$task_id, labels = as.numeric(x$observed_variant_id %in% x$causal_variant_id),
  null = x$null, scenario = x$scenario, group_id = x$cluster_id))
dir.create(pub, recursive = TRUE); dir.create(priv, recursive = TRUE, mode = "0700")
pub <- normalizePath(pub); priv <- normalizePath(priv)
if (pub == priv || startsWith(pub, paste0(priv, "/")) || startsWith(priv, paste0(pub, "/")))
  stop("Resolved public/private paths overlap")
write_json <- function(x, path) jsonlite::write_json(x, path, auto_unbox = TRUE,
                                                 digits = NA, pretty = FALSE, null = "null")
sha <- function(path) {
  out <- system2("openssl", c("dgst", "-sha256", shQuote(path)), stdout = TRUE)
  hit <- regmatches(out, regexpr("[a-fA-F0-9]{64}$", out))
  if (length(hit) != 1 || nchar(hit) != 64) stop("SHA-256 unavailable")
  tolower(hit)
}
for (pair in list(list(directory = pub, data = unname(inputs), file = "inputs.json"),
                  list(directory = priv, data = unname(truth), file = "truth.json"))) {
  file <- file.path(pair$directory, pair$file); write_json(pair$data, file)
  hashes <- list(); hashes[[pair$file]] <- sha(file)
  write_json(list(specification = "finemap-rl-gaussian-v1", n_tasks = length(inputs), files = hashes),
             file.path(pair$directory, "manifest.json"))
}
# Exclusion depends only on publicly declared unsupported input types. These
# tasks remain attempted in the full R campaign; this narrower bridge is not
# advertised as support for all mechanisms.
write_json(list(original_tasks = length(eligible), exported_tasks = sum(eligible),
                excluded_tasks = unname(names(panel$inputs)[!eligible]),
                excluded_reason = "cohort or variant-wise-n contract requires a future adapter"),
           file.path(pub, "bridge-report.json"))
cat("Exported", sum(eligible), "of", length(eligible), "tasks with separate assessor truth.\n")
