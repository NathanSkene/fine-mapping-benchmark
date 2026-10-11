#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
main <- function() {
opt <- .campaign_options(list(plan = "results/full-campaign/plan.rds", `bootstrap-draws` = "1000",
                              backend = "auto", `max-stratum-rows` = "100000"))
plan <- readRDS(opt$plan); .campaign_assert_plan(plan, .campaign_repo)
root <- dirname(normalizePath(opt$plan))
.campaign_claim(file.path(root, "analysis.lock"))
on.exit(unlink(file.path(root, "analysis.lock"), recursive = TRUE), add = TRUE)
analysis <- file.path(root, "analysis")
if (dir.exists(analysis)) {
  archive <- file.path(root, "analysis-archive", paste0(format(Sys.time(), "%Y%m%dT%H%M%OS6"), "-", Sys.getpid()))
  dir.create(dirname(archive), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(analysis, archive)) stop("Cannot archive previous derived analysis")
}
dir.create(analysis, showWarnings = FALSE)
collection_finished <- FALSE
.campaign_analysis_manifest(analysis, plan, "collecting")
on.exit(if (!collection_finished)
  .campaign_analysis_manifest(analysis, plan, "failed_or_interrupted"), add = TRUE, after = FALSE)
backend <- if (opt$backend == "auto") "streamed" else opt$backend
if (!backend %in% c("streamed", "memory")) stop("Unknown collector backend")
if (plan$profile == "full" && backend == "memory") stop("Full campaigns require bounded streamed collection")
if (backend == "streamed") {
  source(file.path(.campaign_repo, "scripts/campaign/collect_streamed.R"))
  census <- .campaign_collect_streamed(plan, root, analysis, as.integer(opt$`bootstrap-draws`),
                                       as.integer(opt$`max-stratum-rows`))
  collection_finished <- TRUE
  cat(sprintf("%d/%d tasks complete; %d attempted locus fits collected with bounded strata\n",
              census$complete, nrow(plan$tasks), census$rows))
  return(invisible(NULL))
}
tables <- old_scores <- status <- list(); j <- 0L
for (i in seq_len(nrow(plan$tasks))) {
  task <- plan$tasks[i, ]; directory <- file.path(root, "tasks", task$task_id)
  evidence <- .campaign_task_evidence(directory, plan, task)
  complete <- evidence$complete
  n_read <- 0L
  for (z in evidence$records) {
    j <- j + 1L; tables[[j]] <- z$table
    n_read <- n_read + 1L
    if (!is.null(z$original_evaluation) && is.null(z$original_evaluation$error)) {
      m <- z$method_task$method; old <- z$original_evaluation[[m]]$global
      old_scores[[length(old_scores) + 1L]] <- data.frame(
        dataset_id = z$metadata$dataset_id, stage = z$method_task$stage, method = m,
        original_grid_auprc = old$auprc %||% NA_real_,
        original_n_fits = old$n_fits %||% NA_integer_,
        original_n_failed = old$n_failed %||% NA_integer_,
        corrected_exact_auprc = z$corrected_raw_evaluation[[m]]$global$auprc,
        corrected_accepted_auprc = z$evaluation[[m]]$global$auprc,
        original_accepted_auprc = z$original_accepted_evaluation[[m]]$global$auprc %||% NA_real_,
        corrected_n_fits = z$evaluation[[m]]$global$n_fits,
        corrected_n_failed = z$evaluation[[m]]$global$n_failed)
    }
  }
  status[[i]] <- data.frame(task_id = task$task_id, stage = task$stage, method = task$method,
       complete = complete, datasets_checkpointed = n_read,
       datasets_expected = length(strsplit(task$dataset_ids, "|", fixed = TRUE)[[1L]]),
       started_attempts = evidence$started_attempts,
       storage = if (evidence$bundled) "sealed_task_bundle" else "loose_dataset_checkpoints",
       evidence_scope = if (complete) "checkpointed_terminal" else "incomplete; missing units excluded from scientific summaries")
}
status <- do.call(rbind, status)
write.csv(status, file.path(root, "analysis/status.csv"), row.names = FALSE)
coverage <- aggregate(status[c("datasets_checkpointed", "datasets_expected", "started_attempts")],
                       status[c("stage", "method")], sum)
coverage$missing_datasets <- coverage$datasets_expected - coverage$datasets_checkpointed
coverage$campaign_complete <- coverage$missing_datasets == 0L
write.csv(coverage, file.path(root, "analysis/planned-coverage.csv"), row.names = FALSE)
if (!length(tables)) {
  .campaign_analysis_manifest(analysis, plan, "no_scientific_output", list(
       tasks_complete = sum(status$complete), tasks_planned = nrow(status),
       datasets_checkpointed = sum(status$datasets_checkpointed),
       datasets_planned = sum(status$datasets_expected), locus_rows = 0L))
  collection_finished <- TRUE
  cat(sprintf("%d/%d tasks complete; no scientific output yet\n", sum(status$complete), nrow(status)))
  return(invisible(NULL))
}
bind_different <- function(xs) {
  columns <- unique(unlist(lapply(xs, names)))
  do.call(rbind, lapply(xs, function(x) {
    for (column in setdiff(columns, names(x))) x[[column]] <- NA
    x[columns]
  }))
}
data <- bind_different(tables)
.campaign_atomic_rds(data, file.path(root, "analysis/task-metrics.rds"))
write.csv(data, file.path(root, "analysis/task-metrics.csv"), row.names = FALSE)
if (length(old_scores)) write.csv(do.call(rbind, old_scores),
    file.path(root, "analysis/scoring-attribution.csv"), row.names = FALSE)
data$stratum <- ifelse(data$family == "challenge",
    paste("challenge", data$scenario, data$p, data$n, sep = ":"),
    paste(data$family, data$n_gwas, data$n_ref, data$annotation_type, data$enrichment_fold,
          data$S, data$phi, data$p_causal, data$nominal_p, sep = ":"))
mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
groups <- split(seq_len(nrow(data)), paste(data$stage, data$method, data$stratum, sep = "|"))
summary <- do.call(rbind, lapply(groups, function(indices) {
  x <- data[indices, ]
  data.frame(stage = x$stage[1L], method = x$method[1L], stratum = x$stratum[1L],
      n_attempted = nrow(x), n_successful = sum(x$success), success_rate = mean(x$success),
      ap_successful = mean_na(x$ap),
      n_ap_eligible_attempts = sum(x$n_observed_causal > 0),
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
      n_log_loss_observed = sum(!is.na(x$log_loss)),
      n_log_loss_infinite = sum(is.infinite(x$log_loss)),
      fdr_micro = if (sum(x$tp + x$fp)) sum(x$fp) / sum(x$tp + x$fp) else NA_real_,
      cs_distinct_recall_all_attempts = mean_na(x$cs_distinct_recall),
      native_all_causal_coverage_reported = mean_na(x$native_all_causal_coverage_reported),
      native_all_causal_recovery_all_attempts = mean_na(x$native_all_causal_recovery_all_attempts),
      cs_causal_variant_recall_all_attempts = mean_na(x$cs_causal_variant_recall_all_attempts),
      mean_cs_support_size = mean_na(x$cs_support_size),
      median_runtime_seconds = if (any(is.finite(x$runtime_seconds)))
        median(x$runtime_seconds, na.rm = TRUE) else NA_real_,
      n_genotype_clusters = length(unique(x$cluster_id)),
      n_panel_clusters = length(unique(x$panel_cluster_id)),
      n_source_locus_clusters = length(unique(stats::na.omit(x$source_locus_cluster))),
      stringsAsFactors = FALSE)
}))
write.csv(summary, file.path(root, "analysis/regime-summary.csv"), row.names = FALSE)

# Comparisons use identical attempted units and retained cluster membership.
# No single all-mechanism leaderboard hides availability, calibration or scope.
contrasts <- list(); ci <- 0L
n_boot <- as.integer(opt$`bootstrap-draws`)
if (is.na(n_boot) || n_boot < 20L) stop("Invalid bootstrap draw count")
for (stage in unique(data$stage)) for (stratum in unique(data$stratum[data$stage == stage])) {
  d <- data[data$stage == stage & data$stratum == stratum, ]
  methods <- sort(unique(d$method)); if (length(methods) < 2L) next
  for (pair in combn(methods, 2L, simplify = FALSE)) for (metric in c("ap", "brier", "recall", "cs_distinct_recall")) {
    a <- d[d$method == pair[1L], ]; b <- d[d$method == pair[2L], ]
    keys <- c("dataset_id", "region_id", "task_id")
    keys <- keys[keys %in% names(a) & vapply(keys, function(k) any(!is.na(a[[k]])), TRUE)]
    pooled <- c("fb_pooled", "fb_xregion", "polyfun_est", "polyfun_ldsc", "sbayesrc", "paintor")
    cluster_field <- if (any(pair %in% pooled)) "panel_cluster_id" else "cluster_id"
    alignment <- c(a$probability_target_aligned[a$success], b$probability_target_aligned[b$success])
    aligned <- length(alignment) > 0L && all(alignment %in% TRUE)
    paired <- merge(a[c(keys, cluster_field, metric)], b[c(keys, metric)], by = keys,
                     suffixes = c("_a", "_b"))
    differences <- paired[[paste0(metric, "_a")]] - paired[[paste0(metric, "_b")]]
    good <- is.finite(differences)
    ds <- differences[good]; clusters <- split(ds, paired[[cluster_field]][good])
    interval <- c(NA_real_, NA_real_)
    if (length(clusters) >= 2L) {
      set.seed(plan$seed + ci)
      sampled <- replicate(n_boot, mean(unlist(clusters[sample(seq_along(clusters),
                                    length(clusters), replace = TRUE)], use.names = FALSE)))
      interval <- quantile(sampled, c(.025, .975), names = FALSE)
    }
    ci <- ci + 1L
    contrasts[[ci]] <- data.frame(stage = stage, stratum = stratum, method_a = pair[1L],
      method_b = pair[2L], metric = metric, difference = if (length(ds)) mean(ds) else NA_real_,
      lower = interval[1L], upper = interval[2L], n_pairs = sum(good),
      n_attempted_pairs = nrow(paired), n_clusters = length(clusters),
      probability_targets_aligned = if (metric == "brier") aligned else NA,
      interpretation = if (metric == "brier" && !aligned)
        "operational_sparse_membership_diagnostic_difference; event or posterior-conditioning mismatch; not calibration ranking" else
        "prespecified operational metric difference")
  }
}
if (length(contrasts)) write.csv(do.call(rbind, contrasts),
     file.path(root, "analysis/paired-comparisons.csv"), row.names = FALSE)
.campaign_analysis_manifest(analysis, plan, if (all(status$complete)) "complete" else "incomplete", list(
       tasks_complete = sum(status$complete), tasks_planned = nrow(status),
       datasets_checkpointed = sum(status$datasets_checkpointed),
       datasets_planned = sum(status$datasets_expected), locus_rows = nrow(data)))
collection_finished <- TRUE
cat(sprintf("%d/%d tasks complete; %d attempted locus fits collected\n",
    sum(status$complete), nrow(status), nrow(data)))
}
main()
