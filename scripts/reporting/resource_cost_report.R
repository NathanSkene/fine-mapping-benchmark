#!/usr/bin/env Rscript
# Analysis only: load immutable campaign helpers from --source. Do not source
# this file from a worker or put it inside scripts/campaign/.

.rc_known_filter_sha <- "651dbb1640b63187028c12406a82801a676f10be1f9a52b635aedad5a49be139"
.rc_filter_basis <- "known_f053_post_simulation_monomorphic_gwas_or_reference_filter"
.rc_scalar <- function(x, label, minimum = 0, integer = FALSE, missing = FALSE) {
  if (missing && is.numeric(x) && length(x) == 1L && is.na(x)) return(x)
  if (!is.numeric(x) || is.complex(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || (integer && x != floor(x))) stop("Invalid ", label)
  x
}
.rc_string <- function(x, label) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) stop("Invalid ", label)
  x
}
.rc_hash <- function(x, label = "SHA256") {
  .rc_string(x, label)
  if (!grepl("^[a-f0-9]{64}$", x)) stop("Invalid ", label)
  x
}
.rc_values <- function(x, label, integer = FALSE) {
  values <- suppressWarnings(as.numeric(strsplit(.rc_string(x, label), "|", fixed = TRUE)[[1L]]))
  if (!length(values) || any(!is.finite(values)) || any(values <= 0) ||
      (integer && any(values != floor(values)))) stop("Invalid ", label)
  values
}
.rc_key <- function(x) {
  # JSON escaping prevents delimiter collisions; NA is an explicit design value.
  jsonlite::toJSON(as.list(x), auto_unbox = TRUE, na = "null", digits = NA)
}
.rc_size_summary <- function(requested, realized) {
  data.frame(requested_p_vector = paste(requested, collapse = "|"),
    realized_p_vector = paste(realized, collapse = "|"),
    dropped_p_vector = paste(requested - realized, collapse = "|"),
    nominal_p_max = max(requested), actual_p_max = max(realized),
    requested_sum_p = sum(requested), realized_sum_p = sum(realized),
    requested_sum_p2 = sum(requested^2), realized_sum_p2 = sum(realized^2),
    requested_sum_p3 = sum(requested^3), realized_sum_p3 = sum(realized^3),
    realized_to_requested_sum_p = sum(realized) / sum(requested),
    realized_to_requested_sum_p2 = sum(realized^2) / sum(requested^2),
    realized_to_requested_sum_p3 = sum(realized^3) / sum(requested^3),
    scaling_interpretation = "dimension sensitivity only; these powers are not runtime or memory bounds",
    stringsAsFactors = FALSE)
}

.rc_dimension_census <- function(audit, requested, n, n_ref = NA_real_, full_probe = FALSE) {
  .rc_scalar(n, "simulated individual N", 2, TRUE)
  .rc_scalar(n_ref, "reference individual N", 2, TRUE, TRUE)
  if (!is.numeric(requested) || any(!is.finite(requested)) || any(requested < 1) ||
      any(requested != floor(requested))) stop("Invalid requested dimensions")
  if (full_probe && (!identical(as.numeric(requested), rep(c(500, 750, 1000, 1500, 2000), each = 2L)) ||
                     n != 5000 || length(requested) != 10L)) stop("Wrong nominal ten-locus full probe")
  if (!is.list(audit) || !identical(audit$schema, 1L) ||
      !identical(audit$source_file, "R/simulate_genotypes.R") ||
      !identical(audit$design_semantics, "nominal_requested_variant_count") ||
      !identical(audit$variant_universe_semantics, "retained_polymorphic_gwas_and_reference_if_present"))
    stop("Missing or unknown generation dimension audit")
  .rc_hash(audit$source_file_sha256, "generation source SHA256")
  known <- identical(audit$source_file_sha256, .rc_known_filter_sha)
  if (!identical(audit$known_post_simulation_filter, known)) stop("Contradictory filter origin")
  expected_basis <- if (known) .rc_filter_basis else "unverified_filter_origin_no_dimension_loss_permitted"
  if (!identical(audit$post_simulation_filter_basis, expected_basis)) stop("Contradictory filter basis")
  d <- audit$loci
  required <- c("region_id", "requested_p", "realized_p", "dropped_p", "retained_fraction",
                "dropped_fraction", "n_gwas", "n_reference", "post_simulation_filter_basis")
  if (!is.data.frame(d) || !all(required %in% names(d)) || nrow(d) != length(requested))
    stop("Malformed generation locus census")
  for (field in c("region_id", "requested_p", "realized_p", "dropped_p", "n_gwas"))
    if (!is.numeric(d[[field]]) || any(!is.finite(d[[field]])) ||
        any(d[[field]] != floor(d[[field]]))) stop("Invalid generation ", field)
  if (!identical(as.numeric(d$region_id), as.numeric(seq_along(requested))) ||
      !identical(as.numeric(d$requested_p), as.numeric(requested)) ||
      any(d$realized_p < 1 | d$realized_p > requested) ||
      any(d$dropped_p != requested - d$realized_p) || any(d$n_gwas != n) ||
      !identical(as.numeric(d$n_reference), rep(as.numeric(n_ref), length(requested))) ||
      !identical(d$post_simulation_filter_basis, rep(expected_basis, length(requested))))
    stop("Contradictory generation dimension audit")
  if (any(d$realized_p < requested) && !known) stop("Unknown origin of dimension reduction")
  for (field in c("retained_fraction", "dropped_fraction")) {
    expected <- if (field == "retained_fraction") d$realized_p / requested else d$dropped_p / requested
    if (!is.numeric(d[[field]]) || any(!is.finite(d[[field]])) ||
        any(abs(d[[field]] - expected) > 1e-14)) stop("Contradictory generation fraction")
  }
  w <- audit$warnings
  if (!is.data.frame(w) || !all(c("index", "message", "classification") %in% names(w)) ||
      !is.character(w$message) || anyNA(w$message)) stop("Malformed preserved generation warnings")
  undersupply <- grepl("requested[[:space:]]+[0-9]+[[:space:]]+SNPs.*only[[:space:]]+[0-9]+[[:space:]]+variants[[:space:]]+available.*MAF", w$message) |
    grepl("Capping p at 500 for regions where p > 500", w$message, fixed = TRUE)
  if (any(undersupply)) stop("VCF undersupply or bundled cap is not a legitimate dimension reduction")
  d$design_label <- if (full_probe) "nominal full-size filtered ten-locus design" else "nominal filtered ten-locus design"
  d$N_basis <- "actual simulated individual rows; sealed generator audit, not source VCF participant count"
  d
}

.rc_scenario_designs <- function(plan) {
  fields <- c("family", "nominal_p_vector", "n", "n_ref", "model", "annotation", "enrichment",
              "annotation_correlation", "p_causal", "S", "phi", "replicates", "n_loci")
  out <- vector("list", nrow(plan$datasets))
  grids <- vector("list", nrow(plan$grid))
  for (i in seq_len(nrow(plan$grid))) {
    g <- plan$grid[i, , drop = FALSE]
    S <- .rc_values(g$S_values, "S values", TRUE); phi <- .rc_values(g$phi_values, "phi values")
    if (g$model == "sparse") {
      sc <- expand.grid(S = S, phi = phi, iter = seq_len(g$n_iter), stringsAsFactors = FALSE)
      sc$p_causal <- NA_real_
    } else if (g$model == "sparse_inf") {
      sc <- expand.grid(S = S, phi = phi, p_causal = g$p_causal,
                        iter = seq_len(g$n_iter), stringsAsFactors = FALSE)
    } else stop("Unknown legacy architecture")
    grids[[i]] <- sc
  }
  for (i in seq_len(nrow(plan$datasets))) {
    d <- plan$datasets[i, , drop = FALSE]
    if (d$family == "challenge") {
      gen <- plan$generation[match(d$row, plan$generation$generation_id), , drop = FALSE]
      if (nrow(gen) != 1L || anyNA(gen$generation_id) || gen$family != "challenge") stop("Foreign Gaussian dataset")
      x <- list(family = "challenge", nominal_p_vector = as.character(gen$p), n = gen$n,
        n_ref = NA_real_, model = "conditional_Gaussian", annotation = "mixed_mechanisms",
        enrichment = NA_character_, annotation_correlation = NA_real_, p_causal = NA_real_,
        S = NA_real_, phi = NA_real_, replicates = gen$replicates, n_loci = 22L * gen$replicates)
    } else {
      if (is.na(d$row) || d$row < 1L || d$row > nrow(plan$grid)) stop("Only generated legacy plans are supported")
      g <- plan$grid[d$row, , drop = FALSE]; sc <- grids[[d$row]]
      if (is.na(d$scenario) || d$scenario < 1L || d$scenario > nrow(sc)) stop("Foreign legacy scenario")
      z <- sc[d$scenario, , drop = FALSE]
      x <- list(family = "legacy_reconstruction", nominal_p_vector = g$p_values, n = g$n,
        n_ref = g$n_ref, model = g$model, annotation = g$annotation_type,
        enrichment = g$enrichment_values, annotation_correlation = g$annotation_correlation,
        p_causal = z$p_causal, S = z$S, phi = z$phi, replicates = 1L, n_loci = g$n_regions)
    }
    if (!identical(names(x), fields)) stop("Internal design schema mismatch")
    out[[i]] <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  data <- do.call(rbind, out); data$dataset_id <- plan$datasets$dataset_id
  rownames(data) <- NULL; data
}

.rc_batch_design <- function(designs, indices, stage, method) {
  d <- designs[indices, , drop = FALSE]
  constant <- c("family", "nominal_p_vector", "n", "n_ref", "model", "annotation", "enrichment",
                "annotation_correlation", "p_causal", "replicates")
  for (field in constant) if (length(unique(d[[field]])) != 1L) stop("Mixed-design worker task: ", field)
  x <- d[1L, constant, drop = FALSE]
  x$S_pack <- paste(d$S, collapse = "|"); x$phi_pack <- paste(d$phi, collapse = "|")
  x$stage <- stage; x$method <- method
  x$dataset_records_per_task <- nrow(d); x$locus_attempts_per_task <- sum(d$n_loci)
  x$workload_key <- .rc_key(x[setdiff(names(x), c("workload_key", "locus_attempts_per_task"))])
  x
}

.rc_expected_rejection <- function(table, method, family, annotation) {
  error <- table$error
  expected <- rep(FALSE, nrow(table))
  if (family == "challenge") {
    scalar <- table$scenario %in% c("meta_missingness", "meta_known_overlap", "multi_ancestry_shared",
                                   "multi_ancestry_heterogeneous")
    if (!method %in% c("exact", "forced_one")) expected <- scalar &
      error == "unsupported_input: regional single-cohort scalar-N method"
    if (method %in% c("exact", "forced_one")) expected <- table$scenario == "meta_missingness" &
      error == "Variant-wise n requires sample-overlap metadata for this adapter"
    if (method == "susie_inf") expected <- expected | error == paste0(
      "unsupported_input: normal_score SuSiE-inf needs a supported fixed-noise infinitesimal variance estimator. ",
      "susieR estimates tau2 jointly with phenotype residual variance and requires observed phenotype sufficient statistics; ",
      "disabling this update silently removes the infinitesimal component.")
    if (method == "finemap_inf") expected <- expected | error == paste0(
      "unsupported_input: normal_score FINEMAP-inf jointly estimates phenotype residual and infinitesimal variances ",
      "using a phenotype second moment; the Gaussian score challenge does not supply those phenotype sufficient statistics.")
  }
  if (method == "funmap" && identical(annotation, "none")) expected <- expected |
    error == paste("No annotations available for this region.", "Funmap requires functional annotations.",
      "Re-run simulation with annotations = 'binary' or annotations = 'continuous'.")
  expected & !table$success
}

.rc_locus_census <- function(table, design, audit = NULL, method, stage) {
  if (!is.data.frame(table) || !all(c("actual_p", "success", "error") %in% names(table)) ||
      !is.logical(table$success) || anyNA(table$success) || !is.character(table$error) || anyNA(table$error) ||
      any(table$success & nzchar(table$error)) || any(!table$success & !nzchar(table$error)))
    stop("Malformed recorded execution outcomes")
  if (!all(c("method", "stage", "family", "dataset_id") %in% names(table)) ||
      anyNA(table[c("method", "stage", "family", "dataset_id")]) ||
      any(table$method != method) || any(table$stage != stage) ||
      any(table$family != design$family) || any(table$dataset_id != design$dataset_id))
    stop("Table execution identity disagrees with frozen task")
  for (field in c("actual_p")) if (!is.numeric(table[[field]]) || any(!is.finite(table[[field]])) ||
      any(table[[field]] < 1 | table[[field]] != floor(table[[field]]))) stop("Invalid table dimensions")
  if (nrow(design) != 1L || nrow(table) != design$n_loci) stop("Missing or extra measured locus attempts")
  family <- design$family
  if (family == "challenge") {
    if (method %in% c("exact", "forced_one") && as.numeric(design$nominal_p_vector) > 200L)
      stop("Exact/forced-one references are restricted to p<=200; no large-p stress measurement is permitted")
    needed <- c("task_id", "scenario", "p", "n")
    if (!all(needed %in% names(table)) || anyNA(table$task_id) || anyDuplicated(table$task_id) ||
        !all(table$p == as.numeric(design$nominal_p_vector)) || !all(table$n == design$n))
      stop("Contradictory Gaussian table dimensions or identity")
    mechanisms <- c("matched_prior", "null", "single_sparse", "multi_sparse", "many_sparse", "opposing_effects",
      "sparse_inf", "ld_reference_small", "ld_reference_shift", "annotation_informative", "annotation_irrelevant",
      "annotation_reversal", "missing_causal", "perfect_ld", "allele_flip", "allele_recoding", "selected_null",
      "selected_signal", "meta_missingness", "meta_known_overlap", "multi_ancestry_shared", "multi_ancestry_heterogeneous")
    counts <- table(factor(table$scenario, levels = mechanisms))
    if (anyNA(match(table$scenario, mechanisms)) || any(counts != design$replicates)) stop("Incomplete Gaussian mechanism batch")
    nominal <- rep(as.numeric(design$nominal_p_vector), nrow(table))
    expected <- nominal - as.integer(table$scenario == "missing_causal")
    if (any(table$actual_p != expected)) stop("Unexplained Gaussian candidate reduction")
    reduction <- ifelse(table$scenario == "missing_causal", "intentional one missing causal candidate", "none")
    ids <- table$task_id
    largest <- nominal == max(nominal) & table$scenario != "missing_causal"
    N_basis <- "direct Gaussian score design parameter; no individual phenotype rows; variable-N and cohort mechanisms distinct"
    nominal_N <- table$n
    reference_N <- rep(NA_real_, nrow(table))
  } else {
    needed <- c("region_id", "nominal_p", "n_gwas", "n_ref", "S", "phi", "p_causal", "annotation_type")
    if (is.null(audit) || !all(needed %in% names(table)) ||
        !identical(as.numeric(table$region_id), as.numeric(audit$region_id)) ||
        !identical(as.numeric(table$nominal_p), as.numeric(audit$requested_p)) ||
        !identical(as.numeric(table$actual_p), as.numeric(audit$realized_p)) ||
        !identical(as.numeric(table$n_gwas), as.numeric(audit$n_gwas)) ||
        !identical(as.numeric(table$n_ref), as.numeric(audit$n_reference)) ||
        any(table$S != design$S) || any(abs(table$phi - design$phi) > 1e-14) ||
        !identical(as.numeric(table$p_causal), rep(as.numeric(design$p_causal), nrow(table))) ||
        any(table$annotation_type != design$annotation)) stop("Table disagrees with sealed generator or scenario design")
    nominal <- audit$requested_p; reduction <- audit$post_simulation_filter_basis
    ids <- as.character(table$region_id); largest <- nominal == max(nominal)
    N_basis <- audit$N_basis[1L]; nominal_N <- table$n_gwas; reference_N <- table$n_ref
  }
  runtime <- if ("runtime_seconds" %in% names(table)) table$runtime_seconds else rep(NA_real_, nrow(table))
  if (!is.numeric(runtime) || any(is.infinite(runtime)) || any(runtime < 0, na.rm = TRUE)) stop("Malformed recorded fit runtime")
  rejection <- .rc_expected_rejection(table, method, family, design$annotation)
  classification <- ifelse(table$success, "successful_execution", ifelse(rejection,
    "expected_unsupported_rejection", "unexpected_execution_failure"))
  data.frame(locus_id = ids, scenario = if (family == "challenge") table$scenario else as.character(design$S),
    requested_p = nominal, actual_p = table$actual_p, dropped_p = nominal - table$actual_p,
    n = nominal_N, n_ref = reference_N, N_basis = N_basis, dimension_reduction_basis = reduction,
    largest_nominal_locus = largest, success = table$success, classification = classification,
    fit_runtime_seconds = runtime, finite_failed_fit_runtime = !table$success & is.finite(runtime),
    error = table$error, stringsAsFactors = FALSE)
}

.rc_execution_support <- function(loci, reused = FALSE) {
  applicable <- loci$classification != "expected_unsupported_rejection"
  large <- applicable & loci$largest_nominal_locus
  largest_success <- any(large) && all(loci$success[large])
  unexpected <- sum(loci$classification == "unexpected_execution_failure")
  supported <- !reused && largest_success && unexpected == 0L
  reason <- if (reused) "reused_fit_is_not_new_inference_measurement" else
    if (!any(applicable)) "all_attempts_are_expected_unsupported_rejections" else
    if (!any(large)) "no_applicable_largest_locus_execution" else
    if (!largest_success) "largest_applicable_locus_failed" else
    if (unexpected) "unexpected_execution_failure" else "observed_applicable_execution_supported"
  list(n_loci = nrow(loci), n_successful = sum(loci$success),
    n_expected_unsupported = sum(!applicable), n_unexpected_failure = unexpected,
    n_finite_failed_fit_runtime = sum(loci$finite_failed_fit_runtime),
    n_applicable_largest_loci = sum(large), n_successful_largest_loci = sum(loci$success & large),
    largest_locus_execution_supported = largest_success,
    new_inference_anchor_eligible = supported, inference_support_reason = reason)
}

.rc_worker_timing <- function(resource, n_records, n_new, elapsed_sum) {
  fields <- c("script_elapsed_seconds", "prelude_seconds", "input_io_seconds", "compaction_seconds",
              "data_processing_seconds", "overhead_seconds")
  empty <- stats::setNames(as.list(rep(NA_real_, length(fields))), fields)
  if (is.null(resource)) return(c(empty, list(complete_new_worker_timing = FALSE)))
  if (!is.list(resource) || !identical(resource$schema, 1L) ||
      resource$n_datasets != n_records || resource$n_new != n_new) stop("Foreign worker timing measurement")
  for (field in fields) .rc_scalar(resource[[field]], paste("worker timing", field))
  tolerance <- max(1e-5, resource$script_elapsed_seconds * 1e-7)
  if (abs(resource$overhead_seconds - sum(unlist(resource[c("prelude_seconds", "input_io_seconds", "compaction_seconds")]))) > tolerance ||
      abs(resource$script_elapsed_seconds - resource$data_processing_seconds - resource$overhead_seconds) > tolerance ||
      elapsed_sum > resource$data_processing_seconds + tolerance) stop("Contradictory worker timing decomposition")
  c(resource[fields], list(complete_new_worker_timing = n_new == n_records))
}

.rc_envelope <- function(workload, measured, workers = 48L) {
  .rc_scalar(workers, "worker cap", 1, TRUE)
  out <- workload
  out$observed_tasks <- out$qualified_timing_tasks <- integer(nrow(out))
  for (field in c("conditional_R_worker_hours_min", "conditional_R_worker_hours_median",
                  "conditional_R_worker_hours_max", "conditional_primary_retention_GB_median",
                  "ideal_R_days_at_worker_cap")) out[[field]] <- NA_real_
  out$unresolved <- "unmeasured_exact_worker_design"
  for (i in seq_len(nrow(out))) {
    m <- measured[measured$workload_key == out$workload_key[i], , drop = FALSE]
    out$observed_tasks[i] <- nrow(m)
    good <- m$task_complete & m$new_inference_anchor_eligible & m$complete_new_worker_timing &
      is.finite(m$script_elapsed_seconds)
    out$qualified_timing_tasks[i] <- sum(good)
    if (!any(good)) {
      if (nrow(m)) out$unresolved[i] <- paste(sort(unique(m$inference_support_reason)), collapse = ";")
      next
    }
    values <- m$script_elapsed_seconds[good] * out$planned_task_units[i] / 3600
    out$conditional_R_worker_hours_min[i] <- min(values)
    out$conditional_R_worker_hours_median[i] <- stats::median(values)
    out$conditional_R_worker_hours_max[i] <- max(values)
    out$conditional_primary_retention_GB_median[i] <- stats::median(m$primary_retention_bytes[good]) *
      out$planned_task_units[i] / 1e9
    out$ideal_R_days_at_worker_cap[i] <- out$conditional_R_worker_hours_median[i] / workers / 24
    out$unresolved[i] <- "conditional_sample_anchor_only; unmeasured_genotype_variation_tail_queue_retries_and_PBS_scope"
  }
  out$feasibility_cleared <- FALSE
  out$envelope_interpretation <- "observed min/median/max anchors for the exact nominal worker design, not confidence intervals or guaranteed bounds"
  out
}

.rc_workload_census <- function(plan, designs, max_groups = 100000L) {
  # Decode each exact dataset pack once, rather than rebuilding a data frame
  # and JSON for every repeated method/stage/phenotype-iteration worker.
  packs <- unique(plan$tasks$dataset_ids)
  if (length(packs) > max_groups) stop("Exact dataset-pack cache exceeds report bound")
  lookup <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_len(nrow(plan$datasets))) assign(plan$datasets$dataset_id[i], i, lookup)
  class_lookup <- new.env(hash = TRUE, parent = emptyenv())
  classes <- list(); pack_class <- integer(length(packs))
  for (i in seq_along(packs)) {
    ids <- strsplit(packs[i], "|", fixed = TRUE)[[1L]]
    if (length(ids) > 5L || anyDuplicated(ids)) stop("Unbounded or duplicated full workload pack")
    di <- vapply(ids, function(id) {
      if (!exists(id, lookup, inherits = FALSE)) stop("Foreign full workload dataset")
      get(id, lookup)
    }, 1L)
    batch <- .rc_batch_design(designs, di, "stage_placeholder", "method_placeholder")
    key <- batch$workload_key
    if (!exists(key, class_lookup, inherits = FALSE)) {
      j <- length(classes) + 1L
      if (j > max_groups) stop("Design cache exceeds report bound")
      classes[[j]] <- batch; assign(key, j, class_lookup)
    }
    pack_class[i] <- get(key, class_lookup)
  }
  for (field in c("stage", "method"))
    if (!is.character(plan$tasks[[field]]) || anyNA(plan$tasks[[field]]) ||
        any(!grepl("^[A-Za-z0-9_]+$", plan$tasks[[field]]))) stop("Malformed workload stage or method")
  index <- pack_class[match(plan$tasks$dataset_ids, packs)]
  group_keys <- paste(index, plan$tasks$stage, plan$tasks$method, sep = ":")
  first <- which(!duplicated(group_keys))
  if (length(first) > max_groups) stop("Workload group bound exceeded")
  counts <- table(group_keys)
  rows <- vector("list", length(first))
  for (j in seq_along(first)) {
    i <- first[j]; batch <- classes[[index[i]]]
    batch$stage <- plan$tasks$stage[i]; batch$method <- plan$tasks$method[i]
    batch$workload_key <- .rc_key(batch[setdiff(names(batch), c("workload_key", "locus_attempts_per_task"))])
    batch$planned_task_units <- as.integer(counts[[group_keys[i]]])
    batch$planned_dataset_records <- batch$planned_task_units * batch$dataset_records_per_task
    batch$planned_locus_attempts <- batch$planned_task_units * batch$locus_attempts_per_task
    rows[[j]] <- batch
  }
  result <- do.call(rbind, rows)
  result <- result[order(result$workload_key, method = "radix"), , drop = FALSE]
  rownames(result) <- NULL
  if (sum(result$planned_task_units) != nrow(plan$tasks) ||
      sum(result$planned_dataset_records) != sum(lengths(strsplit(plan$tasks$dataset_ids, "|", fixed = TRUE))))
    stop("Workload census lost planned attempts")
  result
}

.rc_load_helpers <- function(source, expected_snapshot = NULL) {
  source <- normalizePath(source, mustWork = TRUE)
  if (!dir.exists(file.path(source, "R")) || !file.exists(file.path(source, "scripts/campaign/common.R")))
    stop("--source must be an exact campaign source directory")
  if (nzchar(Sys.getenv("FMB_R_LIB"))) .libPaths(c(Sys.getenv("FMB_R_LIB"), .libPaths()))
  if (dir.exists(file.path(source, ".local-lib"))) .libPaths(c(file.path(source, ".local-lib"), .libPaths()))
  if (!is.null(expected_snapshot)) {
    if (!is.data.frame(expected_snapshot) || !identical(names(expected_snapshot), c("path", "sha256")) ||
        !is.character(expected_snapshot$path) || !is.character(expected_snapshot$sha256) ||
        anyNA(expected_snapshot) || anyDuplicated(expected_snapshot$path) ||
        any(startsWith(expected_snapshot$path, "/")) ||
        any(vapply(strsplit(expected_snapshot$path, "/", fixed = TRUE), function(x) ".." %in% x, TRUE)) ||
        any(!grepl("^[a-f0-9]{64}$", expected_snapshot$sha256))) stop("Malformed frozen source snapshot")
    expected_R <- sort(expected_snapshot$path[grepl("^R/[^/]+[.]R$", expected_snapshot$path)], method = "radix")
    actual_R <- sort(paste0("R/", list.files(file.path(source, "R"), "\\.R$")), method = "radix")
    if (!identical(expected_R, actual_R)) stop("Unlisted or missing R helper file before helper loading")
    for (i in seq_len(nrow(expected_snapshot))) {
      path <- file.path(source, expected_snapshot$path[i])
      if (!file.exists(path) || !identical(.rc_file_sha(path), expected_snapshot$sha256[i]))
        stop("Frozen source bytes changed before helper loading: ", path)
    }
  }
  env <- new.env(parent = globalenv())
  # Equivalent to common.R's helper loading, with an explicit immutable source
  # path; common.R derives its path from --file and cannot route an external script.
  for (path in list.files(file.path(source, "R"), "\\.R$", full.names = TRUE)) sys.source(path, envir = env)
  env$.campaign_repo <- source
  env
}

.rc_check_plan_shape <- function(plan, max_tasks = 1000000L, max_datasets = 100000L) {
  required <- c("plan_sha256", "profile", "corrected_source", "runtime_fingerprint", "method_args",
                "tasks", "datasets", "generation", "grid")
  if (!is.list(plan) || !all(required %in% names(plan))) stop("Malformed plan")
  if (!is.data.frame(plan$tasks) || nrow(plan$tasks) > max_tasks ||
      !is.data.frame(plan$datasets) || nrow(plan$datasets) > max_datasets ||
      !is.data.frame(plan$generation) || nrow(plan$generation) > 1000L || !is.data.frame(plan$grid))
    stop("Plan exceeds report resource bounds")
  for (spec in list(c("tasks", "task_id", "index", "dataset_ids", "stage", "method"),
                    c("datasets", "dataset_id", "family", "row", "scenario", "sim_file", "sim_sha256"),
                    c("generation", "generation_id", "family", "row", "p", "n", "replicates"))) {
    d <- plan[[spec[1L]]]
    if (!all(spec[-1L] %in% names(d)) || anyNA(d[[spec[2L]]]) || anyDuplicated(d[[spec[2L]]])) stop("Malformed plan identities")
  }
  if (!identical(as.numeric(plan$tasks$index), as.numeric(seq_len(nrow(plan$tasks)))) ||
      !identical(as.numeric(plan$generation$generation_id), as.numeric(seq_len(nrow(plan$generation)))))
    stop("Unsupported plan index ordering")
  invisible(TRUE)
}

.rc_generation_retention <- function(z, directory, family, track, helper, max_manifest_bytes = 1024^2) {
  directory <- normalizePath(directory, mustWork = TRUE)
  expected <- if (family == "challenge") file.path(directory, c("panel.rds", "public/manifest.csv", "assessor/manifest.csv")) else
    file.path(directory, "simulation.rds")
  if (!is.data.frame(z$artifacts) || !identical(names(z$artifacts), c("path", "sha256")) ||
      nrow(z$artifacts) != length(expected)) stop("Missing declared generation input artifacts")
  helper$.campaign_validate_artifacts(z$artifacts)
  actual <- normalizePath(z$artifacts$path, mustWork = TRUE)
  if (!setequal(actual, expected)) stop("Foreign generation artifact route")
  paths <- actual
  for (path in paths) track(path)
  if (family == "challenge") for (scope in c("public", "assessor")) {
    release <- file.path(directory, scope); manifest <- file.path(release, "manifest.csv")
    if (file.info(manifest)$size > max_manifest_bytes) stop("Release manifest exceeds report read bound")
    files <- if (scope == "public") c("inputs.rds", "metadata.rds") else c("truth.rds", "generator.rds")
    m <- utils::read.csv(manifest, stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE)
    if (!identical(names(m), c("file", "sha256")) || nrow(m) != 2L || anyNA(m) ||
        anyDuplicated(m$file) || !setequal(m$file, files) || any(!grepl("^[a-f0-9]{64}$", m$sha256)))
      stop("Malformed or foreign release child manifest")
    inventory <- list.files(release, all.files = TRUE, no.. = TRUE)
    if (!setequal(inventory, c(files, "manifest.csv"))) stop("Unlisted retained release payload")
    children <- file.path(release, m$file)
    resolved <- normalizePath(children, mustWork = TRUE)
    if (any(dirname(resolved) != normalizePath(release, mustWork = TRUE)) || any(dir.exists(resolved)))
      stop("Release child escapes its frozen namespace")
    helper$.campaign_validate_artifacts(data.frame(path = resolved, sha256 = m$sha256, stringsAsFactors = FALSE))
    for (i in seq_along(resolved)) if (!identical(track(resolved[i]), m$sha256[i])) stop("Changed release child while reporting")
    paths <- c(paths, resolved)
  }
  paths <- unique(paths)
  sizes <- file.info(paths)$size
  if (anyNA(sizes) || any(sizes < 0)) stop("Cannot measure retained input file size")
  list(paths = paths, retained_input_bytes = sum(sizes), retained_input_files = length(paths),
    input_retention_scope = "all declared generation and release payload logical file bytes; assessor payloads hashed/stat'ed only; filesystem blocks and directory metadata excluded")
}

.rc_generation_index <- function(plan, dataset_index, root) {
  d <- plan$datasets[dataset_index, , drop = FALSE]
  if (!identical(d$sim_sha256, "generated")) stop("Resource reporting requires generated sealed inputs")
  paths <- file.path(root, "datasets", sprintf("generation_%05d", plan$generation$generation_id))
  index <- match(dirname(d$sim_file), paths)
  if (is.na(index)) stop("Dataset generation route is not the frozen plan route")
  index
}

# PBS evidence uses a small declared mapping, but trusts resource values only
# after hashing raw qstat output and binding its Job_Name/Submit_arguments to
# the actual submitted script and its exact worker commands. It does not poll.
.rc_pbs_observations <- function(path, plans, source, track, max_bytes = 64 * 1024^2) {
  empty <- data.frame(plan_sha256 = character(), kind = character(), indices = character(),
    job_id = character(), job_state = character(), exit_status = integer(),
    wall_seconds = numeric(), cpu_seconds = numeric(), peak_memory_bytes = numeric(),
    resource_scope = character(), stringsAsFactors = FALSE)
  if (!nzchar(path)) return(empty)
  if (file.info(path)$size > max_bytes) stop("PBS evidence exceeds resource bound")
  track(path); v <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (!identical(v$schema, 1L) || !identical(v$kind, "fmbenchmark-resource-pbs-bindings-v1") ||
      !is.list(v$attempts) || length(v$attempts) > 10000L) stop("Unknown canonical PBS evidence format")
  duration <- function(x) {
    if (is.null(x)) return(NA_real_)
    .rc_string(x, "PBS duration")
    if (!grepl("^[0-9]+:[0-5][0-9]:[0-5][0-9]$", x)) stop("Invalid PBS duration")
    z <- as.numeric(strsplit(x, ":", fixed = TRUE)[[1L]]); sum(z * c(3600, 60, 1))
  }
  memory <- function(x) {
    if (is.null(x)) return(NA_real_)
    .rc_string(x, "PBS memory")
    if (!grepl("^[0-9]+(kb|mb|gb|b)$", x)) stop("Invalid PBS memory units")
    unit <- sub("^[0-9]+", "", x); as.numeric(sub("(kb|mb|gb|b)$", "", x)) *
      c(b = 1, kb = 1024, mb = 1024^2, gb = 1024^3)[[unit]]
  }
  rows <- vector("list", length(v$attempts)); seen <- character()
  for (i in seq_along(v$attempts)) {
    a <- v$attempts[[i]]; .rc_hash(a$plan_sha256, "PBS plan SHA")
    if (!a$plan_sha256 %in% names(plans)) stop("Foreign PBS plan")
    .rc_string(a$job_id, "PBS identity")
    if (!grepl("^[0-9]+[.][A-Za-z0-9_.-]+$", a$job_id) || a$job_id %in% seen) stop("Duplicate or malformed PBS identity")
    seen <- c(seen, a$job_id)
    if (!a$kind %in% c("generation", "methods")) stop("Foreign PBS unit kind")
    indices <- unlist(a$indices, use.names = FALSE)
    bound <- if (a$kind == "generation") nrow(plans[[a$plan_sha256]]$plan$generation) else nrow(plans[[a$plan_sha256]]$plan$tasks)
    if (!is.numeric(indices) || !length(indices) || any(!is.finite(indices)) || any(indices != floor(indices)) ||
        any(indices < 1 | indices > bound) || anyDuplicated(indices)) stop("Malformed PBS unit binding")
    for (field in c("script", "qstat")) {
      p <- a[[paste0(field, "_path")]]; expected <- a[[paste0(field, "_sha256")]]
      .rc_hash(expected, paste("PBS", field, "SHA")); track(p)
      if (!identical(.rc_file_sha(p), expected)) stop("Changed PBS evidence artifact")
      if (file.info(p)$size > max_bytes) stop("PBS artifact exceeds resource bound")
    }
    raw <- jsonlite::fromJSON(a$qstat_path, simplifyVector = FALSE)
    if (!is.null(raw$stdout)) {
      if (!identical(raw$returncode, 0L)) stop("PBS query failed")
      raw <- jsonlite::fromJSON(raw$stdout, simplifyVector = FALSE)
    }
    job <- raw$Jobs[[a$job_id]]
    if (is.null(job) || !identical(job$job_state, "F") || !is.numeric(job$Exit_status) ||
        length(job$Exit_status) != 1L || job$Exit_status != floor(job$Exit_status))
      stop("PBS evidence does not establish an actual terminal job")
    script <- normalizePath(a$script_path, mustWork = TRUE)
    if (!identical(normalizePath(job$Submit_arguments, mustWork = TRUE), script)) stop("PBS submitted script route mismatch")
    lines <- readLines(script, warn = FALSE)
    if (sum(grepl("^#PBS -N ", lines)) != 1L ||
        !identical(sub("^#PBS -N ", "", lines[grepl("^#PBS -N ", lines)]), job$Job_Name)) stop("PBS job name mismatch")
    # A deliberately restrictive parser: frozen launchers use paths without
    # shell metacharacters. Unknown quoting/control syntax is not executed.
    plan_path <- plans[[a$plan_sha256]]$path
    worker <- file.path(source, "scripts/campaign", if (a$kind == "generation") "generate.R" else "run_task.R")
    if (any(grepl("[[:space:]'\";$`\\\\]", c(worker, plan_path)))) stop("Unsupported PBS command path quoting")
    commands <- lines[grepl("scripts/campaign/(run_task|generate)[.]R --plan ", lines)]
    executable <- unique(c("Rscript", normalizePath(file.path(R.home("bin"), "Rscript"), mustWork = TRUE),
                          unname(Sys.which("Rscript"))))
    executable <- executable[nzchar(executable)]
    allowed <- lapply(indices, function(index) paste(executable, "--vanilla", worker, "--plan", plan_path, "--index", index))
    if (length(commands) != length(indices) ||
        !all(vapply(seq_along(indices), function(i) commands[i] %in% allowed[[i]], TRUE)))
      stop("PBS script is not bound to the exact source, plan and declared units")
    r <- job$resources_used %||rc% list()
    rows[[i]] <- data.frame(plan_sha256 = a$plan_sha256, kind = a$kind,
      indices = paste(indices, collapse = "|"), job_id = a$job_id, job_state = job$job_state,
      exit_status = job$Exit_status, wall_seconds = duration(r$walltime), cpu_seconds = duration(r$cput),
      peak_memory_bytes = memory(r$mem), resource_scope = "actual whole PBS job; batch resources are not allocated to individual methods; qstat and script hashes bound, scheduler authenticity supplied by collector",
      stringsAsFactors = FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else empty
}
`%||rc%` <- function(x, y) if (is.null(x)) y else x
.rc_file_sha <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required for SHA256")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
.rc_bounded_file <- function(path, max_bytes, label) {
  .rc_string(path, label); .rc_scalar(max_bytes, paste(label, "read bound"), 1, TRUE)
  if (!file.exists(path) || dir.exists(path)) stop(label, " must be an existing file")
  size <- file.info(path)$size
  if (is.na(size) || size > max_bytes) stop(label, " exceeds explicit read bound")
  normalizePath(path, mustWork = TRUE)
}

.rc_main <- function(args = commandArgs(TRUE), script_path = NULL) {
  defaults <- list(source = "", `probe-plan` = "", `full-plan` = "", out = "",
    `tail-selections` = "", `pbs-evidence` = "", workers = "48", `max-measured-tasks` = "10000",
    `max-checkpoint-bytes` = as.character(512 * 1024^2), `max-report-groups` = "100000")
  if (length(args) %% 2L) stop("Use --option value pairs")
  if (length(args)) for (i in seq(1L, length(args), 2L)) {
    k <- sub("^--", "", args[i]); if (!k %in% names(defaults)) stop("Unknown option: ", k)
    defaults[[k]] <- args[i + 1L]
  }
  o <- defaults
  for (k in c("source", "probe-plan", "full-plan", "out")) .rc_string(o[[k]], paste("--", k, sep = ""))
  for (k in c("workers", "max-measured-tasks", "max-checkpoint-bytes", "max-report-groups")) {
    o[[k]] <- suppressWarnings(as.numeric(o[[k]])); .rc_scalar(o[[k]], paste("--", k, sep = ""), 1, TRUE)
  }
  o$source <- normalizePath(o$source, mustWork = TRUE)
  o$`probe-plan` <- normalizePath(o$`probe-plan`, mustWork = TRUE)
  o$`full-plan` <- normalizePath(o$`full-plan`, mustWork = TRUE)
  if (file.exists(o$out) || dir.exists(o$out) || !grepl("^resource-cost-[A-Za-z0-9][A-Za-z0-9_.-]*$", basename(o$out)))
    stop("Choose a new versioned --out directory named resource-cost-...")
  parent <- normalizePath(dirname(o$out), mustWork = TRUE); o$out <- file.path(parent, basename(o$out))
  protected <- c(o$source, dirname(o$`probe-plan`), dirname(o$`full-plan`))
  if (any(startsWith(paste0(o$out, "/"), paste0(protected, "/")))) stop("Report output may not be inside source or frozen campaign roots")
  if (nzchar(Sys.getenv("FMB_R_LIB"))) .libPaths(c(Sys.getenv("FMB_R_LIB"), .libPaths()))
  if (dir.exists(file.path(o$source, ".local-lib"))) .libPaths(c(file.path(o$source, ".local-lib"), .libPaths()))
  for (path in c(o$`probe-plan`, o$`full-plan`))
    if (file.info(path)$size > o$`max-checkpoint-bytes`) stop("Plan exceeds explicit read bound")
  plan_bytes <- vapply(c(o$`probe-plan`, o$`full-plan`), .rc_file_sha, "")
  probe <- readRDS(o$`probe-plan`); full <- readRDS(o$`full-plan`)
  .rc_check_plan_shape(probe); .rc_check_plan_shape(full)
  helper <- .rc_load_helpers(o$source, probe$corrected_source)
  # Helpers are private to the exact-source environment. The parent global
  # environment supplies base/package functions, not scientific replacements.
  for (p in list(probe, full)) { .rc_check_plan_shape(p); helper$.campaign_assert_plan(p, o$source) }
  if (probe$profile != "resource_probe" || full$profile != "full") stop("Require resource_probe and full profiles")
  if (!identical(probe$corrected_source, full$corrected_source) ||
      !identical(probe$runtime_fingerprint, full$runtime_fingerprint)) stop("Probe and full lineage differ")
  if (!identical(probe$method_args, full$method_args) ||
      !identical(probe$original_source_hashes, full$original_source_hashes))
    stop("Probe/full method settings or original-generator byte lineage differ")
  if (!identical(probe$challenge$variants, 1000L) || !identical(probe$challenge$n, 10000L) ||
      probe$challenge$mixed_mechanism_batch_replicates != 10L) stop("Wrong required Gaussian resource-probe design")
  probe_design <- .rc_scenario_designs(probe); full_design <- .rc_scenario_designs(full)
  for (i in which(probe_design$family != "challenge")) {
    .rc_scalar(probe_design$n[i], "probe N", 2, TRUE)
    if (probe_design$n[i] != 5000 || !identical(.rc_values(probe_design$nominal_p_vector[i], "probe p", TRUE),
                  rep(c(500, 750, 1000, 1500, 2000), each = 2L)) || probe_design$n_loci[i] != 10L)
      stop("Resource probe lacks the nominal N5000 ten-locus full-size design")
  }
  tails <- if (nzchar(o$`tail-selections`)) strsplit(o$`tail-selections`, ",", fixed = TRUE)[[1L]] else character()
  if (length(tails) > 1000L || anyDuplicated(tails)) stop("Tail selection count exceeds report bound or repeats a file")
  if (length(tails)) tails <- vapply(tails, .rc_bounded_file, "", max_bytes = o$`max-checkpoint-bytes`, label = "Tail selection")
  if (anyDuplicated(tails)) stop("Tail selection repeats the same canonical file")
  tail_indices <- integer(); tail_generations <- integer()
  for (path in tails) {
    v <- jsonlite::fromJSON(path)
    if (!identical(v$schema, 1L) || !identical(v$plan_sha256, full$plan_sha256) ||
        !identical(v$plan_rds_sha256, helper$.challenge_sha256(o$`full-plan`)) ||
        !is.character(v$source_commit) || !grepl("^[a-f0-9]{40}$", v$source_commit)) stop("Foreign tail selection")
    for (field in c("projection_manifest_sha256", "legacy_grid_csv_sha256",
                   if (!is.null(v$original_scenario_source_sha256)) "original_scenario_source_sha256")) {
      .rc_hash(v[[field]], paste("tail", field))
      route <- switch(field,
        projection_manifest_sha256 = file.path(dirname(o$`full-plan`), "dispatch-manifest.json"),
        legacy_grid_csv_sha256 = file.path(dirname(o$`full-plan`), "legacy-grid.csv"),
        original_scenario_source_sha256 = file.path(full$original_source, "R/run_simulation.R"))
      if (!file.exists(route) || !identical(helper$.challenge_sha256(route), v[[field]])) stop("Changed frozen tail binding: ", field)
    }
    if (!is.null(v$parent_selection_sha256) && !v$parent_selection_sha256 %in% vapply(tails, .rc_file_sha, ""))
      stop("Tail addendum requires its exact parent selection")
    idx <- v$method_task_indices
    if (!is.numeric(idx) || !length(idx) || any(!is.finite(idx)) || any(idx != floor(idx)) ||
        any(idx < 1L | idx > nrow(full$tasks)) || anyDuplicated(idx)) stop("Malformed selected tail indices")
    gen <- v$generation_indices %||rc% v$generation_index
    if (!is.numeric(gen) || !length(gen) || any(!is.finite(gen)) || any(gen != floor(gen)) ||
        any(gen < 1L | gen > nrow(full$generation)) || anyDuplicated(gen)) stop("Malformed selected generation indices")
    tail_indices <- union(tail_indices, idx); tail_generations <- union(tail_generations, gen)
  }
  if (nrow(probe$tasks) + length(tail_indices) > o$`max-measured-tasks`) stop("Measured task selection exceeds explicit resource bound")
  if (!dir.create(o$out)) stop("Cannot exclusively create report directory")
  artifacts <- new.env(hash = TRUE, parent = emptyenv()); n_artifacts <- 0L
  track <- function(path) {
    path <- normalizePath(path, mustWork = TRUE)
    sha <- helper$.challenge_sha256(path)
    if (exists(path, artifacts, inherits = FALSE)) {
      if (!identical(get(path, artifacts), sha)) stop("Input changed while reporting: ", path)
    } else {
      n_artifacts <<- n_artifacts + 1L
      if (n_artifacts > 100000L) stop("Artifact manifest exceeds report bound")
      assign(path, sha, artifacts)
    }
    invisible(sha)
  }
  for (path in c(o$`probe-plan`, o$`full-plan`, tails, script_path)) if (!is.null(path) && nzchar(path)) track(path)
  if (!identical(unname(plan_bytes), unname(vapply(c(o$`probe-plan`, o$`full-plan`), track, ""))))
    stop("Plan bytes changed while loading")
  write_part <- function(data, name) {
    if (!nrow(data)) return(invisible(NULL))
    path <- file.path(o$out, name); new <- !file.exists(path)
    utils::write.table(data, path, sep = ",", quote = TRUE, row.names = FALSE, col.names = new,
                       append = !new, na = "NA", qmethod = "double")
  }
  plans <- stats::setNames(list(list(plan = probe, path = o$`probe-plan`, design = probe_design,
                  indices = seq_len(nrow(probe$tasks)), generations = seq_len(nrow(probe$generation)), scope = "resource_probe")), probe$plan_sha256)
  if (length(tail_indices)) plans[[full$plan_sha256]] <- list(plan = full, path = o$`full-plan`, design = full_design,
        indices = sort(tail_indices), generations = sort(tail_generations), scope = "prespecified_runtime_tails")
  # Full is also eligible for PBS mapping even when no tail selection is supplied.
  pbs_plans <- plans; pbs_plans[[full$plan_sha256]] <- list(plan = full, path = o$`full-plan`)
  pbs <- .rc_pbs_observations(o$`pbs-evidence`, pbs_plans, o$source, track)
  if (nrow(pbs)) write_part(pbs, "pbs-attempts.csv")
  measured <- list(); task_status <- list(); generation_status <- list()
  for (binding in plans) {
    p <- binding$plan; root <- dirname(binding$path); dim_cache <- new.env(hash = TRUE, parent = emptyenv())
    verified_generations <- new.env(hash = TRUE, parent = emptyenv())
    for (gi in binding$generations) {
      g <- p$generation[gi, , drop = FALSE]
      directory <- file.path(root, "datasets", sprintf("generation_%05d", gi)); checkpoint <- file.path(directory, "generated.rds")
      receipt <- file.path(directory, "terminal.json")
      complete <- helper$.campaign_check_receipt(directory, checkpoint, p$plan_sha256, paste0("generation_", gi))
      if (complete) {
        if (file.info(checkpoint)$size > o$`max-checkpoint-bytes`) stop("Generation checkpoint exceeds read bound")
        track(checkpoint); track(receipt)
        helper$.campaign_validate_checkpoint(checkpoint, p, g)
        z <- readRDS(checkpoint)
        retention <- .rc_generation_retention(z, directory, g$family, track, helper)
        assign(as.character(gi), TRUE, verified_generations)
        if (g$family != "challenge") {
          grid <- p$grid[g$row, , drop = FALSE]
          audit <- .rc_dimension_census(z$generation_dimension_audit,
            .rc_values(grid$p_values, "requested p", TRUE), grid$n, grid$n_ref,
            full_probe = binding$scope == "resource_probe")
          source_sha <- helper$.challenge_sha256(file.path(p$original_source, "R/simulate_genotypes.R"))
          if (!identical(source_sha, z$generation_dimension_audit$source_file_sha256)) stop("Audit generator origin disagrees with frozen original source")
          audit$generation_index <- gi; audit$plan_sha256 <- p$plan_sha256; audit$scope <- binding$scope
          write_part(audit, "generation-dimensions.csv"); assign(as.character(gi), audit, dim_cache)
          sizes <- .rc_size_summary(audit$requested_p, audit$realized_p)
        } else sizes <- .rc_size_summary(g$p, g$p)
        seconds <- as.numeric(difftime(z$completed, z$started, units = "secs"))
        .rc_scalar(seconds, "generation elapsed seconds")
        row <- cbind(data.frame(plan_sha256 = p$plan_sha256, scope = binding$scope, generation_index = gi,
          family = g$family, complete = TRUE, generation_elapsed_seconds = seconds,
          retained_input_bytes = retention$retained_input_bytes, retained_input_files = retention$retained_input_files,
          generation_checkpoint_receipt_bytes = sum(file.info(c(checkpoint, receipt))$size),
          input_retention_scope = retention$input_retention_scope, stringsAsFactors = FALSE), sizes)
        write_part(row, "generation-costs.csv")
        rm(z)
      }
      generation_status[[length(generation_status) + 1L]] <- data.frame(plan_sha256 = p$plan_sha256,
        scope = binding$scope, generation_index = gi, family = g$family, complete = complete,
        interpretation = if (complete) "sealed generation metadata and artifacts verified" else "missing sealed generation evidence; unresolved",
        stringsAsFactors = FALSE)
    }
    for (ti in binding$indices) {
      task <- p$tasks[ti, , drop = FALSE]; ids <- strsplit(task$dataset_ids, "|", fixed = TRUE)[[1L]]
      di <- match(ids, p$datasets$dataset_id)
      if (anyNA(di) || anyDuplicated(ids) || length(ids) > 5L) stop("Foreign or unbounded task dataset binding")
      directory <- file.path(root, "tasks", task$task_id)
      checkpoint <- file.path(directory, "completed.rds")
      if (file.exists(checkpoint) && file.info(checkpoint)$size > o$`max-checkpoint-bytes`) stop("Task bundle exceeds explicit read bound")
      for (path in file.path(directory, paste0(ids, ".rds")))
        if (file.exists(path) && file.info(path)$size > o$`max-checkpoint-bytes`) stop("Loose checkpoint exceeds explicit read bound")
      evidence <- helper$.campaign_task_evidence(directory, p, task)
      task_status[[length(task_status) + 1L]] <- data.frame(plan_sha256 = p$plan_sha256, scope = binding$scope,
        task_id = task$task_id, task_index = ti, stage = task$stage, method = task$method,
        complete = evidence$complete, expected_dataset_records = length(ids), sealed_dataset_records = length(evidence$records),
        started_attempts = evidence$started_attempts, stringsAsFactors = FALSE)
      batch <- .rc_batch_design(binding$design, di, task$stage, task$method)
      outcomes <- list(); elapsed <- numeric(); reused <- logical(); retained <- 0; n_new <- 0L; resource <- NULL
      if (evidence$complete) {
        track(checkpoint); track(file.path(directory, "terminal.json"))
        terminal <- jsonlite::fromJSON(file.path(directory, "terminal.json")); resource <- terminal$resource_measurement
        retained <- sum(file.info(c(checkpoint, file.path(directory, "terminal.json")))$size)
      }
      if (evidence$bundled) {
        bundle <- readRDS(checkpoint)
        helper$.campaign_validate_bundle(bundle, p, task, verify_record_hashes = TRUE)
        n_new <- bundle$n_new %||rc% 0L
        attempts <- bundle$attempt_records
      } else {
        attempt_paths <- list.files(file.path(directory, "attempts"), "[.]rds$", full.names = TRUE)
        if (length(attempt_paths) > 1000L) stop("Attempt ledger exceeds report read bound")
        if (length(attempt_paths) && (any(file.info(attempt_paths)$size > o$`max-checkpoint-bytes`) ||
            sum(file.info(attempt_paths)$size) > o$`max-checkpoint-bytes`))
          stop("Loose attempt ledger exceeds explicit checkpoint read bound")
        attempts <- stats::setNames(lapply(attempt_paths, function(path) { track(path); readRDS(path) }), basename(attempt_paths))
      }
      for (nm in names(attempts)) {
        a <- attempts[[nm]]
        if (!identical(a$plan_sha256, p$plan_sha256) || !identical(a$task, task) || !a$dataset$dataset_id %in% ids ||
            !identical(a$dataset, p$datasets[match(a$dataset$dataset_id, p$datasets$dataset_id), , drop = FALSE]))
          stop("Foreign started attempt")
        write_part(data.frame(plan_sha256 = p$plan_sha256, scope = binding$scope, task_id = task$task_id,
          dataset_id = a$dataset$dataset_id, attempt_name = nm, started = as.character(a$started),
          pid = a$pid, host = a$host, duration_available = FALSE,
          interpretation = "all started attempts retained; this ledger records starts, not duration or completion",
          stringsAsFactors = FALSE), "started-attempts.csv")
      }
      for (id in names(evidence$records)) {
        z <- evidence$records[[id]]; dindex <- match(id, binding$design$dataset_id)
        if (is.na(dindex)) stop("Foreign sealed dataset")
        design <- binding$design[dindex, , drop = FALSE]
        gi <- .rc_generation_index(p, match(id, p$datasets$dataset_id), root)
        if (!gi %in% binding$generations) stop("Tail selection omits a required generation unit")
        if (!exists(as.character(gi), verified_generations, inherits = FALSE))
          stop("Sealed ", if (design$family == "challenge") "Gaussian" else "legacy", " fit lacks compatible sealed generation/input custody")
        if (design$family != "challenge" && !exists(as.character(gi), dim_cache, inherits = FALSE))
          stop("Sealed legacy fit lacks compatible sealed generation audit")
        audit <- if (design$family == "challenge") NULL else get(as.character(gi), dim_cache)
        loci <- .rc_locus_census(z$table, design, audit, task$method, task$stage)
        loci$plan_sha256 <- p$plan_sha256; loci$scope <- binding$scope
        loci$task_id <- task$task_id; loci$dataset_id <- id; loci$stage <- task$stage; loci$method <- task$method
        loci$reused <- isTRUE(z$reuse$reused)
        write_part(loci, "locus-attempts.csv")
        support <- .rc_execution_support(loci, isTRUE(z$reuse$reused)); outcomes[[length(outcomes) + 1L]] <- support
        .rc_scalar(z$elapsed_seconds, "dataset processing seconds"); elapsed <- c(elapsed, z$elapsed_seconds)
        reused <- c(reused, isTRUE(z$reuse$reused))
        if (!evidence$bundled) {
          file <- file.path(directory, paste0(id, ".rds")); receipt <- file.path(directory, paste0(id, ".terminal.json"))
          track(file); track(receipt); retained <- retained + sum(file.info(c(file, receipt))$size)
        }
        r <- cbind(data.frame(plan_sha256 = p$plan_sha256, scope = binding$scope, task_id = task$task_id,
           dataset_id = id, stage = task$stage, method = task$method, task_complete = evidence$complete,
           dataset_processing_seconds = z$elapsed_seconds, reused = isTRUE(z$reuse$reused), stringsAsFactors = FALSE),
           design[setdiff(names(design), c("dataset_id", "n_loci"))], as.data.frame(support, stringsAsFactors = FALSE),
           .rc_size_summary(loci$requested_p, loci$actual_p))
        write_part(r, "measured-records.csv")
      }
      if (length(outcomes)) {
        timing <- .rc_worker_timing(resource, length(ids), n_new, sum(elapsed))
        eligible <- evidence$complete && length(outcomes) == length(ids) && all(vapply(outcomes, `[[`, TRUE, "new_inference_anchor_eligible"))
        reasons <- paste(sort(unique(vapply(outcomes, `[[`, "", "inference_support_reason"))), collapse = ";")
        m <- cbind(data.frame(plan_sha256 = p$plan_sha256, scope = binding$scope, task_id = task$task_id,
          task_complete = evidence$complete, started_attempts = evidence$started_attempts,
          observed_dataset_records = length(outcomes), n_reused_records = sum(reused),
          dataset_processing_seconds = sum(elapsed), primary_retention_bytes = retained,
          new_inference_anchor_eligible = eligible, inference_support_reason = reasons,
          stringsAsFactors = FALSE), batch, as.data.frame(timing, stringsAsFactors = FALSE))
        measured[[length(measured) + 1L]] <- m; write_part(m, "measured-tasks.csv")
      }
      rm(evidence); if (exists("bundle", inherits = FALSE)) rm(bundle)
    }
  }
  write_part(do.call(rbind, task_status), "task-status.csv")
  write_part(do.call(rbind, generation_status), "generation-status.csv")
  measured <- if (length(measured)) do.call(rbind, measured) else data.frame(workload_key = character(),
    task_complete = logical(), new_inference_anchor_eligible = logical(), complete_new_worker_timing = logical(),
    script_elapsed_seconds = numeric(), primary_retention_bytes = numeric(), inference_support_reason = character())
  # Plan metadata alone determines the complete workload, including unmeasured
  # exact/forced-one workers. No full fit checkpoint or truth is opened here.
  workload <- .rc_workload_census(full, full_design, o$`max-report-groups`)
  write_part(workload, "full-workload.csv")
  envelope <- .rc_envelope(workload, measured, o$workers); write_part(envelope, "conditional-cost-envelope.csv")
  paths <- sort(ls(artifacts, all.names = TRUE), method = "radix")
  expected <- unname(vapply(paths, function(path) get(path, artifacts), ""))
  # Revalidate all input bytes and frozen helpers/runtime after the complete
  # read, before producing a completion manifest. No receipt is rewritten.
  helper$.campaign_validate_artifacts(data.frame(path = paths, sha256 = expected, stringsAsFactors = FALSE))
  helper$.campaign_assert_plan(probe, o$source); helper$.campaign_assert_plan(full, o$source)
  write_part(helper$.campaign_artifact_manifest(paths), "input-artifacts.csv")
  totals <- list(schema = 1L, kind = "fmbenchmark-analysis-only-resource-report-v1",
    probe_plan_sha256 = probe$plan_sha256, full_plan_sha256 = full$plan_sha256,
    source_directory = o$source, corrected_source = probe$corrected_source,
    runtime_fingerprint = probe$runtime_fingerprint, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    interpretation = "conditional planning anchors and unresolved evidence, not worst-case guarantees, convergence or scientific winners",
    complete = TRUE, source_plans_and_receipts_unchanged = TRUE,
    planned_method_task_units = nrow(full$tasks), planned_dataset_records = sum(workload$planned_dataset_records),
    planned_locus_attempts = sum(workload$planned_locus_attempts), measured_task_rows = nrow(measured),
    qualified_exact_design_groups = sum(envelope$qualified_timing_tasks > 0L),
    unresolved_exact_design_groups = sum(envelope$qualified_timing_tasks == 0L),
    conditional_qualified_subset_R_hours_median = if (any(envelope$qualified_timing_tasks > 0L))
      sum(envelope$conditional_R_worker_hours_median, na.rm = TRUE) else NA_real_,
    full_R_cost_established = FALSE, inference_feasibility_cleared = FALSE,
    PBS_actual_attempts_observed = nrow(pbs), PBS_resources_complete_for_all_attempts = FALSE,
    input_artifacts = length(paths), max_checkpoint_bytes = o$`max-checkpoint-bytes`,
    generation_and_PBS_attempt_costs = "reported separately; prior failed attempts and input retention not erased; unavailable duration remains unknown",
    storage_scope = "measured task checkpoints/receipts, sealed attempt ledger included; generation inputs separate; scheduler/logs/retries/archives/fit-cache dedup excluded from extrapolation",
    stage_scope = "original_refit, corrected_refit and corrected_resimulation kept separate",
    tail_scope = "N1000 and reference sizes including500/200 are separate; packed S1/2/3/5/10 and CARMA S10-only are distinct worker designs; no paired LD-only interpretation",
    Gaussian_stress_scope = "required p1000/N10000/220-attempt mixed batch; intentional missing-cause reduction is explicit; exact/forced_one small-p scope remains unresolved without its own measurements")
  helper$.campaign_atomic_rds(totals, file.path(o$out, "report-manifest.rds"))
  jsonlite::write_json(totals, file.path(o$out, "report-manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  outputs <- list.files(o$out, full.names = TRUE)
  output_manifest <- helper$.campaign_artifact_manifest(outputs)
  utils::write.csv(output_manifest, file.path(o$out, "output-artifacts.csv"), row.names = FALSE)
  cat(sprintf("Report complete: %d measured tasks; %d exact worker designs have timing anchors; %d remain unmeasured. Full feasibility remains unresolved.\n",
    nrow(measured), totals$qualified_exact_design_groups, totals$unresolved_exact_design_groups))
  invisible(totals)
}

if (sys.nframe() == 0L) {
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
  .rc_main(script_path = normalizePath(file, mustWork = TRUE))
}
