#!/usr/bin/env Rscript
# Bounded collection: one task's fits, then one design stratum, in memory.
.campaign_csv_block <- function(x, path, append = file.exists(path), missing = "NA",
                                double_format = c("decimal", "hexadecimal")) {
  double_format <- match.arg(double_format)
  encoded <- lapply(x, function(v) {
    if (!is.atomic(v) || !is.null(dim(v))) stop("Unsupported metric column type")
    value <- if (is.double(v)) sprintf(if (double_format == "decimal") "%.17g" else "%a", v) else as.character(v)
    absent <- is.na(v)
    if (is.double(v)) absent <- absent & !is.nan(v)
    if (is.character(v) && any(v == missing, na.rm = TRUE))
      stop("A literal metric string collides with the declared missing-value marker")
    value[absent] <- missing
    value
  })
  encoded <- as.data.frame(encoded, stringsAsFactors = FALSE, optional = TRUE)
  names(encoded) <- names(x)
  utils::write.table(encoded, path, append = append, sep = ",", quote = which(vapply(x, is.character, TRUE)),
      row.names = FALSE, col.names = !append, qmethod = "double", fileEncoding = "UTF-8")
}

.campaign_write_metric_block <- function(table, connection) {
  value <- serialize(table, NULL, version = 2L)
  writeBin(as.double(length(value)), connection, size = 8L, endian = "little")
  writeBin(value, connection)
}

.campaign_each_metric_block <- function(path, callback) {
  connection <- file(path, "rb"); on.exit(close(connection), add = TRUE)
  count <- rows <- 0L
  repeat {
    header <- readBin(connection, "raw", n = 8L)
    if (!length(header)) break
    if (length(header) != 8L) stop("Truncated metric block header")
    size <- readBin(header, "double", n = 1L, size = 8L, endian = "little")
    if (!is.finite(size) || size < 1 || size != floor(size) || size > .Machine$integer.max)
      stop("Malformed metric block size")
    bytes <- readBin(connection, "raw", n = as.integer(size))
    if (length(bytes) != size) stop("Truncated metric block stream")
    table <- unserialize(bytes)
    if (!is.data.frame(table)) stop("Metric block is not a data frame")
    callback(table)
    count <- count + 1L; rows <- rows + nrow(table)
  }
  list(blocks = count, rows = rows)
}

.campaign_metric_classes <- function(table) {
  classes <- vapply(table, function(v) {
    if (!is.atomic(v) || !is.null(dim(v)) || !typeof(v) %in% c("logical", "integer", "double", "character"))
      stop("Unsupported metric class; collection must not silently coerce structured columns")
    typeof(v)
  }, "")
  classes
}

.campaign_normalize_metric_block <- function(table, classes) {
  for (name in setdiff(names(classes), names(table))) table[[name]] <- NA
  table <- table[names(classes)]
  for (name in names(classes)) table[[name]] <- switch(classes[[name]],
      logical = as.logical(table[[name]]), integer = as.integer(table[[name]]),
      double = as.double(table[[name]]), character = as.character(table[[name]]))
  table
}

.campaign_partition_writer <- function(max_open = 64L) {
  if (length(max_open) != 1L || is.na(max_open) || max_open < 1L)
    stop("Invalid open-partition bound")
  handles <- list(); clock <- 0L; peak <- 0L
  close_all <- function() {
    for (entry in handles) close(entry$connection)
    handles <<- list()
    invisible(NULL)
  }
  write <- function(table, path, missing) {
    clock <<- clock + 1L
    if (is.null(handles[[path]])) {
      if (length(handles) >= max_open) {
        oldest <- which.min(vapply(handles, `[[`, 1L, "last"))
        close(handles[[oldest]]$connection); handles[[oldest]] <<- NULL
      }
      append <- file.exists(path)
      handles[[path]] <<- list(connection = file(path, "at", encoding = "UTF-8"),
                               last = clock, append = append)
      peak <<- max(peak, length(handles))
    }
    .campaign_csv_block(table, handles[[path]]$connection, append = handles[[path]]$append,
                        missing = missing, double_format = "hexadecimal")
    handles[[path]]$append <<- TRUE; handles[[path]]$last <<- clock
    invisible(NULL)
  }
  list(write = write, close = close_all, peak = function() peak)
}

.campaign_streamed_group <- function(x) {
  mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  data.frame(stage = x$stage[1L], method = x$method[1L], stratum = x$stratum[1L],
    n_attempted = nrow(x), n_successful = sum(x$success), success_rate = mean(x$success),
    ap_successful = mean_na(x$ap), n_ap_eligible_attempts = sum(x$n_observed_causal > 0),
    n_ap_eligible_successful = sum(x$n_observed_causal > 0 & x$success),
    ap_observed_positive_attempts_zero_for_failure = mean_na(
        .campaign_ap_failure_zero(x$ap, x$success, x$n_observed_causal)),
    brier_successful = mean_na(x$brier), recall_all_attempts = mean_na(x$recall),
    log_loss_successful = mean_na(x$log_loss),
    pip_event = paste(unique(x$pip_event), collapse = ";"),
    posterior_conditioning = paste(unique(x$posterior_conditioning), collapse = ";"),
    native_null_support = paste(unique(x$native_null_support), collapse = ";"),
    n_successful_null_model_in_support = sum(x$success & x$null_model_in_support, na.rm = TRUE),
    n_successful_null_model_excluded = sum(x$success & !x$null_model_in_support, na.rm = TRUE),
    n_native_null_probabilities_reported = sum(x$success & x$null_probability_reported),
    score_target = paste(unique(x$score_target), collapse = ";"),
    n_probability_target_aligned_successful = sum(x$success & x$probability_target_aligned, na.rm = TRUE),
    n_probability_target_misaligned_successful = sum(x$success & !x$probability_target_aligned, na.rm = TRUE),
    n_probability_target_unknown_successful = sum(x$success & is.na(x$probability_target_aligned)),
    probability_score_interpretation = paste(unique(x$probability_score_interpretation), collapse = ";"),
    brier_any_nonzero_total_effect_successful = mean_na(x$brier_any_nonzero_total_effect),
    log_loss_any_nonzero_total_effect_successful = mean_na(x$log_loss_any_nonzero_total_effect),
    n_log_loss_any_nonzero_infinite = sum(is.infinite(x$log_loss_any_nonzero_total_effect)),
    n_log_loss_observed = sum(!is.na(x$log_loss)), n_log_loss_infinite = sum(is.infinite(x$log_loss)),
    fdr_micro = if (sum(x$tp + x$fp)) sum(x$fp) / sum(x$tp + x$fp) else NA_real_,
    cs_distinct_recall_all_attempts = mean_na(x$cs_distinct_recall),
    native_all_causal_coverage_reported = mean_na(x$native_all_causal_coverage_reported),
    native_all_causal_recovery_all_attempts = mean_na(x$native_all_causal_recovery_all_attempts),
    cs_causal_variant_recall_all_attempts = mean_na(x$cs_causal_variant_recall_all_attempts),
    mean_cs_support_size = mean_na(x$cs_support_size),
    median_runtime_seconds = if (any(is.finite(x$runtime_seconds))) median(x$runtime_seconds, na.rm = TRUE) else NA_real_,
    n_genotype_clusters = length(unique(x$cluster_id)), n_panel_clusters = length(unique(x$panel_cluster_id)),
    n_source_locus_clusters = length(unique(stats::na.omit(x$source_locus_cluster))), stringsAsFactors = FALSE)
}

.campaign_streamed_pairs <- function(d, plan, n_boot, ci, path) {
  methods <- sort(unique(d$method)); if (length(methods) < 2L) return(ci)
  for (pair in combn(methods, 2L, simplify = FALSE)) for (metric in c("ap", "brier", "recall", "cs_distinct_recall")) {
    a <- d[d$method == pair[1L], ]; b <- d[d$method == pair[2L], ]
    keys <- c("dataset_id", "region_id", "task_id")
    keys <- keys[keys %in% names(a) & vapply(keys, function(k) any(!is.na(a[[k]])), TRUE)]
    pooled <- c("fb_pooled", "fb_xregion", "polyfun_est", "polyfun_ldsc", "sbayesrc", "paintor")
    cluster_field <- if (any(pair %in% pooled)) "panel_cluster_id" else "cluster_id"
    alignment <- c(a$probability_target_aligned[a$success], b$probability_target_aligned[b$success])
    aligned <- length(alignment) > 0L && all(alignment %in% TRUE)
    paired <- merge(a[c(keys, cluster_field, metric)], b[c(keys, metric)], by = keys, suffixes = c("_a", "_b"))
    differences <- paired[[paste0(metric, "_a")]] - paired[[paste0(metric, "_b")]]
    good <- is.finite(differences); ds <- differences[good]
    clusters <- split(ds, paired[[cluster_field]][good]); interval <- c(NA_real_, NA_real_)
    if (length(clusters) >= 2L) {
      set.seed(plan$seed + ci)
      sampled <- replicate(n_boot, mean(unlist(clusters[sample(seq_along(clusters), length(clusters), replace = TRUE)], use.names = FALSE)))
      interval <- quantile(sampled, c(.025, .975), names = FALSE)
    }
    ci <- ci + 1L
    .campaign_csv_block(data.frame(stage = d$stage[1L], stratum = d$stratum[1L], method_a = pair[1L],
      method_b = pair[2L], metric = metric, difference = if (length(ds)) mean(ds) else NA_real_,
      lower = interval[1L], upper = interval[2L], n_pairs = sum(good), n_attempted_pairs = nrow(paired),
      n_clusters = length(clusters), probability_targets_aligned = if (metric == "brier") aligned else NA,
      interpretation = if (metric == "brier" && !aligned)
        "operational_sparse_membership_diagnostic_difference; event or posterior-conditioning mismatch; not calibration ranking" else
        "prespecified operational metric difference"), path, append = ci > 1L)
  }
  ci
}

.campaign_collect_streamed <- function(plan, root, analysis, n_boot, max_stratum_rows = 100000L) {
  if (is.na(n_boot) || n_boot < 20L || is.na(max_stratum_rows) || max_stratum_rows < 1L)
    stop("Invalid bootstrap or per-stratum memory bound")
  blocks_path <- file.path(analysis, "task-metrics.blocks")
  connection <- file(blocks_path, "wb")
  connection_open <- TRUE
  global_open <- summary_open <- pairs_open <- score_open <- FALSE
  partition_writer <- NULL
  on.exit({
    if (connection_open) close(connection)
    if (global_open) close(global_connection)
    if (summary_open) close(summary_connection)
    if (pairs_open) close(pairs_connection)
    if (score_open) close(score_connection)
    if (!is.null(partition_writer)) partition_writer$close()
  }, add = TRUE)
  classes <- character(); count <- rows <- complete_count <- checkpointed <- planned <- 0L
  score_count <- 0L; status_buffer <- list()
  coverage <- list(); rank <- c(logical = 1L, integer = 2L, double = 3L, character = 4L)
  for (i in seq_len(nrow(plan$tasks))) {
    task <- plan$tasks[i, ]; directory <- file.path(root, "tasks", task$task_id)
    evidence <- .campaign_task_evidence(directory, plan, task)
    expected <- length(strsplit(task$dataset_ids, "|", fixed = TRUE)[[1L]])
    status <- data.frame(task_id = task$task_id, stage = task$stage, method = task$method,
      complete = evidence$complete, datasets_checkpointed = length(evidence$records), datasets_expected = expected,
      started_attempts = evidence$started_attempts,
      storage = if (evidence$bundled) "sealed_task_bundle" else "loose_dataset_checkpoints",
      evidence_scope = if (evidence$complete) "checkpointed_terminal" else "incomplete; missing units excluded from scientific summaries")
    status_buffer[[length(status_buffer) + 1L]] <- status
    if (length(status_buffer) == 256L) {
      .campaign_csv_block(do.call(rbind, status_buffer), file.path(analysis, "status.csv"))
      status_buffer <- list()
    }
    complete_count <- complete_count + as.integer(evidence$complete)
    checkpointed <- checkpointed + length(evidence$records); planned <- planned + expected
    key <- paste(task$stage, task$method, sep = "|")
    if (is.null(coverage[[key]])) coverage[[key]] <- data.frame(stage = task$stage, method = task$method,
      datasets_checkpointed = 0L, datasets_expected = 0L, started_attempts = 0L)
    coverage[[key]]$datasets_checkpointed <- coverage[[key]]$datasets_checkpointed + length(evidence$records)
    coverage[[key]]$datasets_expected <- coverage[[key]]$datasets_expected + expected
    coverage[[key]]$started_attempts <- coverage[[key]]$started_attempts + evidence$started_attempts
    for (z in evidence$records) {
      table <- z$table
      table$stratum <- ifelse(table$family == "challenge", paste("challenge", table$scenario, table$p, table$n, sep = ":"),
        paste(table$family, table$n_gwas, table$n_ref, table$annotation_type, table$enrichment_fold,
              table$S, table$phi, table$p_causal, table$nominal_p, sep = ":"))
      next_classes <- .campaign_metric_classes(table)
      for (name in names(next_classes)) if (!name %in% names(classes) || rank[next_classes[[name]]] > rank[classes[[name]]])
        classes[name] <- next_classes[[name]]
      .campaign_write_metric_block(table, connection)
      count <- count + 1L; rows <- rows + nrow(table)
      if (!is.null(z$original_evaluation) && is.null(z$original_evaluation$error)) {
        m <- z$method_task$method; old <- z$original_evaluation[[m]]$global
        if (!score_open) {
          score_connection <- file(file.path(analysis, "scoring-attribution.csv"), "wt", encoding = "UTF-8")
          score_open <- TRUE
        }
        .campaign_csv_block(data.frame(dataset_id = z$metadata$dataset_id, stage = z$method_task$stage, method = m,
          original_grid_auprc = old$auprc %||% NA_real_, original_n_fits = old$n_fits %||% NA_integer_,
          original_n_failed = old$n_failed %||% NA_integer_, corrected_exact_auprc = z$corrected_raw_evaluation[[m]]$global$auprc,
          corrected_accepted_auprc = z$evaluation[[m]]$global$auprc,
          original_accepted_auprc = z$original_accepted_evaluation[[m]]$global$auprc %||% NA_real_,
          corrected_n_fits = z$evaluation[[m]]$global$n_fits, corrected_n_failed = z$evaluation[[m]]$global$n_failed),
          score_connection, append = score_count > 0L)
        score_count <- score_count + 1L
      }
    }
  }
  close(connection); connection_open <- FALSE
  if (score_open) { close(score_connection); score_open <- FALSE }
  if (length(status_buffer)) .campaign_csv_block(do.call(rbind, status_buffer), file.path(analysis, "status.csv"))
  census <- do.call(rbind, coverage); rownames(census) <- NULL
  census$missing_datasets <- census$datasets_expected - census$datasets_checkpointed
  census$campaign_complete <- census$missing_datasets == 0L
  .campaign_csv_block(census, file.path(analysis, "planned-coverage.csv"))
  if (!count) {
    unlink(blocks_path)
    .campaign_analysis_manifest(analysis, plan, "no_scientific_output", list(tasks_complete = complete_count,
      tasks_planned = nrow(plan$tasks), datasets_checkpointed = checkpointed, datasets_planned = planned, locus_rows = 0L))
    return(invisible(list(complete = complete_count, rows = 0L)))
  }
  missing <- paste0("<FMB_MISSING:", plan$plan_sha256, ">")
  schema <- list(schema = 1L, columns = names(classes), storage_classes = classes, missing_marker = missing,
      numeric_precision = "CSV: 17 significant decimal digits; binary metric blocks: exact IEEE doubles", encoding = "UTF-8",
      partition_double_encoding = "hexadecimal (%a), exact R parsing; discarded after aggregation",
      logical_values = c("TRUE", "FALSE"), block_stream = "8-byte little-endian double length + R serialization v2")
  .campaign_atomic_rds(schema, file.path(analysis, "task-metrics-schema.rds"))
  jsonlite::write_json(schema, file.path(analysis, "task-metrics-schema.json"), auto_unbox = TRUE, pretty = TRUE)
  partitions <- file.path(analysis, "stratum-work"); dir.create(partitions)
  partition_writer <- .campaign_partition_writer()
  global_connection <- file(file.path(analysis, "task-metrics.csv"), "wt", encoding = "UTF-8")
  global_open <- TRUE; global_rows <- 0L
  index <- new.env(hash = TRUE, parent = emptyenv()); descriptions <- list(); stages <- character()
  observed <- .campaign_each_metric_block(blocks_path, function(table) {
    table <- .campaign_normalize_metric_block(table, classes)
    .campaign_csv_block(table, global_connection, append = global_rows > 0L, missing = missing)
    global_rows <<- global_rows + nrow(table)
    for (key in unique(paste(table$stage, table$stratum, sep = "|"))) {
      select <- which(paste(table$stage, table$stratum, sep = "|") == key)
      if (!exists(key, index, inherits = FALSE)) {
        id <- length(descriptions) + 1L; assign(key, id, index)
        descriptions[[id]] <<- list(stage = table$stage[select[1L]], stratum = table$stratum[select[1L]], rows = 0L,
                                   file = file.path(partitions, sprintf("stratum_%06d.csv", id)))
        if (!table$stage[select[1L]] %in% stages) stages <<- c(stages, table$stage[select[1L]])
      }
      id <- get(key, index, inherits = FALSE)
      descriptions[[id]]$rows <<- descriptions[[id]]$rows + length(select)
      if (descriptions[[id]]$rows > max_stratum_rows) stop("Design stratum exceeds the explicit collection memory bound: ", key)
      partition_writer$write(table[select, , drop = FALSE], descriptions[[id]]$file, missing)
    }
  })
  close(global_connection); global_open <- FALSE; partition_writer$close()
  if (observed$blocks != count || observed$rows != rows) stop("Metric stream census changed")
  .campaign_atomic_rds(list(schema = 2L, format = "metric_block_stream", plan_sha256 = plan$plan_sha256,
    file = blocks_path, sha256 = .challenge_sha256(blocks_path), blocks = count, rows = rows,
    column_schema = schema), file.path(analysis, "task-metrics.rds"))
  order <- order(match(vapply(descriptions, `[[`, "", "stage"), stages), seq_along(descriptions), method = "radix")
  col_classes <- unname(ifelse(classes == "double", "numeric", classes)); ci <- 0L
  summary_connection <- file(file.path(analysis, "regime-summary.csv"), "wt", encoding = "UTF-8")
  pairs_connection <- file(file.path(analysis, "paired-comparisons.csv"), "wt", encoding = "UTF-8")
  summary_open <- pairs_open <- TRUE; summary_count <- 0L
  for (id in order) {
    d <- utils::read.csv(descriptions[[id]]$file, stringsAsFactors = FALSE, check.names = FALSE,
        colClasses = col_classes, na.strings = missing, fileEncoding = "UTF-8")
    if (nrow(d) != descriptions[[id]]$rows) stop("Stratum row census changed")
    for (method in sort(unique(d$method))) {
      .campaign_csv_block(.campaign_streamed_group(d[d$method == method, ]), summary_connection,
                          append = summary_count > 0L)
      summary_count <- summary_count + 1L
    }
    ci <- .campaign_streamed_pairs(d, plan, n_boot, ci, pairs_connection)
    unlink(descriptions[[id]]$file)
  }
  close(summary_connection); summary_open <- FALSE
  close(pairs_connection); pairs_open <- FALSE
  if (!ci) unlink(file.path(analysis, "paired-comparisons.csv"))
  unlink(partitions, recursive = TRUE)
  .campaign_analysis_manifest(analysis, plan, if (complete_count == nrow(plan$tasks)) "complete" else "incomplete", list(
     tasks_complete = complete_count, tasks_planned = nrow(plan$tasks), datasets_checkpointed = checkpointed,
     datasets_planned = planned, locus_rows = rows, collector = "streamed_strata_v1",
     max_stratum_rows = max_stratum_rows, observed_max_stratum_rows = max(vapply(descriptions, `[[`, 1L, "rows")),
     max_open_partition_files = 64L, observed_open_partition_files = partition_writer$peak()))
  invisible(list(complete = complete_count, rows = rows))
}

if (sys.nframe() == 0L) source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "collect.R"))
