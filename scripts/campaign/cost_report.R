#!/usr/bin/env Rscript
# Resource extrapolation from actual full-size canary checkpoints. These are
# planning ranges, not confidence intervals or completed full-campaign evidence.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
opt <- .campaign_options(list(`canary-plan` = "", `full-plan` = "", workers = "64",
                              out = "results/campaign-cost-estimates.csv"))
if (!nzchar(opt$`canary-plan`) || !nzchar(opt$`full-plan`)) stop("Both plans are required")
canary <- readRDS(opt$`canary-plan`); full <- readRDS(opt$`full-plan`)
.campaign_assert_plan(canary, .campaign_repo); .campaign_assert_plan(full, .campaign_repo)
if (!canary$profile %in% c("canary", "resource_probe") || full$profile != "full")
  stop("Cost extrapolation requires canary/resource_probe and full profiles")
workers <- as.integer(opt$workers)
if (is.na(workers) || workers < 1L) stop("Invalid worker count")
root <- dirname(normalizePath(opt$`canary-plan`)); samples <- list()
for (task in split(canary$tasks, seq_len(nrow(canary$tasks)))) {
  directory <- file.path(root, "tasks", task$task_id)
  checkpoint <- file.path(directory, "completed.rds")
  .campaign_check_receipt(directory, checkpoint, canary$plan_sha256, task$task_id)
  .campaign_validate_checkpoint(checkpoint, canary, task)
  ids <- strsplit(task$dataset_ids, "|", fixed = TRUE)[[1L]]
  files <- list.files(directory, "\\.rds$", full.names = TRUE)
  files <- files[basename(files) != "completed.rds"]
  if (any(!basename(files) %in% paste0(ids, ".rds"))) stop("Unplanned canary checkpoint")
  for (file in files) {
    id <- sub("\\.rds$", "", basename(file))
    if (!.campaign_check_receipt(directory, file, canary$plan_sha256,
        paste(task$task_id, id, sep = "."), paste0(id, ".terminal.json"))) next
    d <- canary$datasets[match(id, canary$datasets$dataset_id), ]
    .campaign_validate_checkpoint(file, canary, list(task = task, dataset = d))
    z <- readRDS(file)
    if (!isTRUE(z$terminal) || !identical(z$plan_sha256, canary$plan_sha256)) stop("Foreign checkpoint")
    m <- z$method_task$method
    cache_path <- if (isTRUE(z$reuse$eligible)) file.path(root, "fit-cache", paste0(z$reuse$key, ".rds")) else ""
    if (nzchar(cache_path) && !.campaign_check_receipt(dirname(cache_path), cache_path,
          canary$plan_sha256, z$reuse$key, paste0(z$reuse$key, ".terminal.json")))
      stop("Measured fit cache is unsealed")
    challenge <- identical(z$metadata$family, "challenge")
    if (challenge && !m %in% c("exact", "forced_one") &&
        (nrow(z$table) != 220L || max(z$table$actual_p) != 1000L)) next
    samples[[length(samples) + 1L]] <- data.frame(method = m, stage = z$method_task$stage,
      family = z$metadata$family,
      seconds = if (challenge) z$elapsed_seconds else z$result[[m]]$total_runtime_seconds %||% z$elapsed_seconds,
      end_to_end_seconds = z$elapsed_seconds, n_loci = nrow(z$table),
      n_successful = sum(z$table$success), input_p_max = max(z$table$actual_p),
      n = if (challenge) z$table$n[1L] else z$table$n_gwas[1L], reused = isTRUE(z$reuse$reused),
      checkpoint_bytes = unname(file.info(file)$size),
      cache_bytes = if (nzchar(cache_path) && file.exists(cache_path)) unname(file.info(cache_path)$size) else 0,
      stringsAsFactors = FALSE)
  }
}
if (!length(samples)) stop("No actual full-size canary output; estimates cannot be manufactured")
samples <- do.call(rbind, samples)
legacy <- samples[samples$stage != "challenge", ]
stress <- samples[samples$stage == "challenge", ]
if (any(legacy$n != 5000L) || any(legacy$input_p_max != 2000L))
  stop("Canaries did not realize the full-size design")
if (canary$profile == "resource_probe" && nrow(stress) &&
    (any(stress$n != 10000L) || any(stress$input_p_max != 1000L) || any(stress$n_loci != 220L)))
  stop("Resource probe did not realize a full 220-task p1000/N10000 stress panel")
rows <- list()
for (key in unique(paste(full$tasks$stage, full$tasks$method, sep = "|"))) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1L]]; stage <- parts[1L]; method <- parts[2L]
  planned <- full$tasks[full$tasks$stage == stage & full$tasks$method == method, ]
  n_scenarios <- sum(vapply(strsplit(planned$dataset_ids, "|", fixed = TRUE), length, 1L))
  if (stage == "challenge") {
    ids <- unlist(strsplit(planned$dataset_ids, "|", fixed = TRUE), use.names = FALSE)
    plan_rows <- full$datasets[match(ids, full$datasets$dataset_id), ]
    draws <- full$generation$replicates[match(plan_rows$row, full$generation$generation_id)]
    planned_loci <- sum(draws) * nrow(challenge_scenarios())
  } else planned_loci <- 10L * n_scenarios
  measured <- samples[samples$stage == stage & samples$method == method, ]
  good <- is.finite(measured$seconds) & !measured$reused
  compute <- measured[good, ]
  lower <- middle <- upper <- NA_real_
  if (nrow(compute)) {
    per_locus <- compute$seconds / compute$n_loci
    lower <- min(per_locus) * planned_loci / 3600
    middle <- median(per_locus) * planned_loci / 3600
    upper <- max(per_locus) * planned_loci / 3600
  }
  storage <- cache_storage <- NA_real_
  if (nrow(measured)) {
    storage <- median(measured$checkpoint_bytes / measured$n_loci) * planned_loci / 1e9
    cache_storage <- median(measured$cache_bytes / measured$n_loci) * planned_loci / 1e9
  }
  rows[[length(rows) + 1L]] <- data.frame(stage = stage, method = method,
      planned_scenarios_or_panels = n_scenarios, planned_locus_fits = planned_loci,
      measured_nonreused_canary_units = nrow(compute),
      measured_units_with_any_success = sum(compute$n_successful > 0L),
      probe_success_rate = if (sum(compute$n_loci)) sum(compute$n_successful) / sum(compute$n_loci) else NA_real_,
      measured_storage_canary_units = nrow(measured),
      worker_hours_low = lower, worker_hours_median = middle, worker_hours_high = upper,
      ideal_days_at_worker_cap = middle / workers / 24,
      checkpoint_storage_gb_median = storage,
      fit_cache_storage_gb_before_dedup_median = cache_storage,
      scope = if (stage == "challenge")
        "largest observed stress panel extrapolated across dimensions; pooled panel-level cost; full retention; queue/inputs excluded" else
        "n5000 and 500-2000 variants; full-retention conditional extrapolation; cache before exact-fit dedup; inputs/queue excluded")
}
rows <- do.call(rbind, rows)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(rows, opt$out, row.names = FALSE)
write.csv(samples, sub("\\.csv$", "-measured.csv", opt$out), row.names = FALSE)
cat(sprintf("Measured %d full-size canary units. Estimates include %d rows with unresolved evidence.\n",
    nrow(samples), sum(is.na(rows$worker_hours_median))))
cat("Do not add the per-method ideal days as a promised completion date: queue delay, controls,\n")
cat("stress dimensions, convergence failures and legitimately reused fits change total cost.\n")
