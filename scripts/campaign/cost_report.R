#!/usr/bin/env Rscript
# Resource extrapolation from actual full-size canary checkpoints. These are
# planning ranges, not confidence intervals or completed full-campaign evidence.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
opt <- .campaign_options(list(`canary-plan` = "", `full-plan` = "", workers = "48",
                              out = "results/campaign-cost-estimates.csv"))
if (!nzchar(opt$`canary-plan`) || !nzchar(opt$`full-plan`)) stop("Both plans are required")
if (!grepl("[.]csv$", opt$out)) stop("--out must end in .csv")
canary <- readRDS(opt$`canary-plan`); full <- readRDS(opt$`full-plan`)
.campaign_assert_plan(canary, .campaign_repo); .campaign_assert_plan(full, .campaign_repo)
if (!canary$profile %in% c("canary", "resource_probe") || full$profile != "full")
  stop("Cost extrapolation requires canary/resource_probe and full profiles")
workers <- as.integer(opt$workers)
if (is.na(workers) || workers < 1L) stop("Invalid worker count")
root <- dirname(normalizePath(opt$`canary-plan`)); samples <- list()
cost_inputs <- c(normalizePath(opt$`canary-plan`), normalizePath(opt$`full-plan`))
for (task in split(canary$tasks, seq_len(nrow(canary$tasks)))) {
  directory <- file.path(root, "tasks", task$task_id)
  evidence <- .campaign_task_evidence(directory, canary, task)
  if (evidence$complete) cost_inputs <- c(cost_inputs, evidence$checkpoint, file.path(directory, "terminal.json"))
  resource <- if (evidence$complete) jsonlite::fromJSON(file.path(directory, "terminal.json"))$resource_measurement else NULL
  timed <- !is.null(resource) && identical(resource$schema, 1L) &&
       isTRUE(resource$n_new == resource$n_datasets) && isTRUE(resource$n_datasets == length(evidence$records))
  if (timed && !all(vapply(resource[c("script_elapsed_seconds", "data_processing_seconds", "overhead_seconds")],
       function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0, TRUE)))
    stop("Malformed worker resource measurement")
  for (id in names(evidence$records)) {
    z <- evidence$records[[id]]
    file <- if (evidence$bundled) evidence$checkpoint else file.path(directory, paste0(id, ".rds"))
    if (!evidence$bundled) cost_inputs <- c(cost_inputs, file, file.path(directory, paste0(id, ".terminal.json")))
    if (!isTRUE(z$terminal) || !identical(z$plan_sha256, canary$plan_sha256)) stop("Foreign checkpoint")
    m <- z$method_task$method
    cache_path <- if (isTRUE(z$reuse$eligible) && !is.null(z$reuse$key))
        file.path(root, "fit-cache", paste0(z$reuse$key, ".rds")) else ""
    if (nzchar(cache_path) && !.campaign_check_receipt(dirname(cache_path), cache_path,
          canary$plan_sha256, z$reuse$key, paste0(z$reuse$key, ".terminal.json")))
      stop("Measured fit cache is unsealed")
    challenge <- identical(z$metadata$family, "challenge")
    if (challenge && !m %in% c("exact", "forced_one") &&
        (nrow(z$table) != 220L || max(z$table$actual_p) != 1000L)) next
    allocation <- if (evidence$bundled) length(evidence$records) else 1L
    bytes <- function(paths) sum(file.info(paths)$size, na.rm = TRUE)
    extra_bytes <- if (evidence$bundled) bytes(file.path(directory, "terminal.json")) / allocation else
        bytes(file.path(directory, paste0(id, ".terminal.json"))) +
        bytes(c(file.path(directory, c("completed.rds", "terminal.json")),
          list.files(file.path(directory, "attempts"), "[.]rds$", full.names = TRUE))) / length(evidence$records)
    samples[[length(samples) + 1L]] <- data.frame(task_id = task$task_id, method = m, stage = z$method_task$stage,
      family = z$metadata$family,
      seconds = if (challenge) z$elapsed_seconds else z$result[[m]]$total_runtime_seconds %||% z$elapsed_seconds,
      dataset_processing_seconds = z$elapsed_seconds,
      end_to_end_seconds = if (timed) resource$script_elapsed_seconds / length(evidence$records) else NA_real_,
      worker_data_processing_seconds = if (timed) resource$data_processing_seconds / length(evidence$records) else NA_real_,
      worker_overhead_seconds = if (timed) resource$overhead_seconds else NA_real_,
      n_loci = nrow(z$table),
      n_successful = sum(z$table$success), input_p_max = max(z$table$actual_p),
      n = if (challenge) z$table$n[1L] else z$table$n_gwas[1L], reused = isTRUE(z$reuse$reused),
      checkpoint_bytes = unname(file.info(file)$size) / allocation,
      primary_retention_bytes = unname(file.info(file)$size) / allocation + extra_bytes,
      retention = if (evidence$bundled) "sealed_task_bundle; attempts included; bytes allocated equally by record" else "loose_dataset_checkpoint",
      cache_bytes = if (nzchar(cache_path) && file.exists(cache_path)) bytes(c(cache_path,
         file.path(dirname(cache_path), paste0(z$reuse$key, ".terminal.json")))) else 0,
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
  end_lower <- end_middle <- end_upper <- NA_real_
  if (nrow(compute)) {
    per_locus <- compute$seconds / compute$n_loci
    lower <- min(per_locus) * planned_loci / 3600
    middle <- median(per_locus) * planned_loci / 3600
    upper <- max(per_locus) * planned_loci / 3600
  }
  timed_compute <- measured[is.finite(measured$worker_data_processing_seconds) &
                             is.finite(measured$worker_overhead_seconds), ]
  if (nrow(timed_compute)) {
    # Data work scales with attempted loci; startup/input IO/compaction scales
    # with worker task units. A one-scenario canary must not charge its startup
    # five times to a five-scenario full worker.
    per_locus <- timed_compute$worker_data_processing_seconds / timed_compute$n_loci
    overhead <- timed_compute$worker_overhead_seconds[!duplicated(timed_compute$task_id)]
    end_lower <- (min(per_locus) * planned_loci + min(overhead) * nrow(planned)) / 3600
    end_middle <- (median(per_locus) * planned_loci + median(overhead) * nrow(planned)) / 3600
    end_upper <- (max(per_locus) * planned_loci + max(overhead) * nrow(planned)) / 3600
  }
  storage <- cache_storage <- NA_real_
  if (nrow(measured)) {
    storage <- median(measured$primary_retention_bytes / measured$n_loci) * planned_loci / 1e9
    cache_storage <- median(measured$cache_bytes / measured$n_loci) * planned_loci / 1e9
  }
  rows[[length(rows) + 1L]] <- data.frame(stage = stage, method = method,
      planned_scenarios_or_panels = n_scenarios, planned_locus_fits = planned_loci,
      measured_nonreused_canary_units = nrow(compute),
      measured_units_with_any_success = sum(compute$n_successful > 0L),
      probe_success_rate = if (sum(compute$n_loci)) sum(compute$n_successful) / sum(compute$n_loci) else NA_real_,
      measured_storage_canary_units = nrow(measured),
      measured_complete_worker_timing_units = length(unique(timed_compute$task_id)),
      worker_hours_low = lower, worker_hours_median = middle, worker_hours_high = upper,
      ideal_days_at_worker_cap = middle / workers / 24,
      end_to_end_worker_hours_low = end_lower, end_to_end_worker_hours_median = end_middle,
      end_to_end_worker_hours_high = end_upper,
      end_to_end_ideal_days_at_worker_cap = end_middle / workers / 24,
      end_to_end_scope = "R script entry through compaction; excludes R launch/final timing receipt/lock cleanup, PBS queue and failed prior attempts; dataset work scales by loci and overhead by task count",
      checkpoint_storage_gb_median = storage,
      fit_cache_storage_gb_before_dedup_median = cache_storage,
      scope = if (stage == "challenge")
        "largest observed stress panel extrapolated across dimensions; pooled panel-level cost; full retention; queue/inputs excluded" else
        "n5000 and 500-2000 variants; full-retention conditional extrapolation; cache before exact-fit dedup; inputs/queue excluded")
}
rows <- do.call(rbind, rows)
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.csv(rows, opt$out, row.names = FALSE)
.campaign_atomic_rds(list(schema = 1L, canary_plan_sha256 = canary$plan_sha256,
    full_plan_sha256 = full$plan_sha256, source = canary$corrected_source,
    runtime_fingerprint = canary$runtime_fingerprint,
    artifacts = .campaign_artifact_manifest(unique(cost_inputs))),
    sub("[.]csv$", "-inputs.rds", opt$out))
write.csv(samples, sub("\\.csv$", "-measured.csv", opt$out), row.names = FALSE)
packed <- lengths(strsplit(full$tasks$dataset_ids, "|", fixed = TRUE))
regular <- data.frame(retention = c("loose_checkpoints_with_one_attempt", "sealed_task_bundles"),
  planned_task_units = nrow(full$tasks), planned_dataset_records = sum(packed),
  primary_terminal_regular_files = c(3 * sum(packed) + 2 * nrow(full$tasks), 2 * nrow(full$tasks)),
  primary_terminal_directories = c(2 * nrow(full$tasks), nrow(full$tasks)),
  scope = "scheduler scripts/logs, inputs, retries and old caches excluded; file counts, not measured byte costs")
write.csv(regular, sub("\\.csv$", "-file-counts.csv", opt$out), row.names = FALSE)
cat(sprintf("Measured %d full-size canary units. Estimates include %d rows with unresolved evidence.\n",
    nrow(samples), sum(is.na(rows$worker_hours_median))))
cat("Do not add the per-method ideal days as a promised completion date: queue delay, controls,\n")
cat("stress dimensions, convergence failures and legitimately reused fits change total cost.\n")
