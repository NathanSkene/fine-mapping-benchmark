#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
main <- function() {
opt <- .native_options(list(plan = "results/native-campaign/plan.rds", `bootstrap-draws` = "1000"))
plan <- readRDS(opt$plan); .native_assert_plan(plan, .native_repo)
if (!requireNamespace("posterior", quietly = TRUE) ||
    !identical(plan$posterior_version, as.character(utils::packageVersion("posterior"))))
  stop("Frozen official posterior diagnostic dependency changed or is missing")
runtime <- .campaign_runtime_fingerprint(plan$method_args)
if (!identical(runtime, plan$runtime_fingerprint)) stop("Frozen collection tool/runtime dependencies changed")
root <- dirname(normalizePath(opt$plan)); analysis <- file.path(root, "analysis")
lock <- file.path(root, "analysis.lock")
.campaign_claim(lock)
on.exit(unlink(lock, recursive = TRUE), add = TRUE)
.native_archive_analysis(root)
dir.create(analysis, showWarnings = FALSE)
started <- Sys.time(); collection_status <- "failed"; collection_error <- NULL
coverage <- list(n_prepared = 0L, n_planned_draws = nrow(plan$draws),
  n_checkpointed_chains = 0L, n_planned_chains = nrow(plan$chains),
  n_checkpointed_regional = 0L, n_planned_regional = nrow(plan$regional))
on.exit(.native_seal_analysis(analysis, plan, runtime, collection_status, coverage,
                              started, collection_error), add = TRUE, after = FALSE)
tryCatch({
metrics_file <- file.path(analysis, "task-metrics.csv")
chain_status <- regional_status <- core_diagnostics <- input_diagnostics <- block_stability <- aggregates <- list()
first <- TRUE; ai <- ci <- ri <- di <- ii <- bi <- 0L; verified_prepared <- character()
for (draw_id in plan$draws$draw_id) {
  prepared_path <- file.path(root, "prepared", draw_id, "completed.rds")
  if (!.campaign_check_receipt(dirname(prepared_path), prepared_path, plan$plan_sha256, draw_id)) next
  prepared <- .native_read_prepared(plan, root, draw_id)
  verified_prepared <- c(verified_prepared, draw_id)
  sim_all <- readRDS(file.path(root, "prepared", draw_id, "simulation.rds"))
  controls <- readRDS(file.path(root, "prepared", draw_id, "annotation-controls.rds"))
  rows <- list(); j <- 0L
  for (input_id in plan$inputs$input_id[plan$inputs$draw_id == draw_id]) {
    input <- plan$inputs[match(input_id, plan$inputs$input_id), ]
    info <- prepared$inputs[match(input_id, prepared$inputs$input_id), ]
    sim <- .native_expose_annotations(sim_all, info$scenario, input$exposure, controls)
    tasks <- plan$chains[plan$chains$input_id == input_id, ]; fits <- vector("list", 4L)
    exists <- successful <- rep(FALSE, 4L); traces <- character(4L); trace_errors <- character(4L)
    for (k in 1:4) {
      task <- tasks[k, ]; path <- file.path(root, "chains", task$chain_id, "completed.rds")
      exists[k] <- .native_completed_checkpoint(path, plan, task, task$chain_id)
      z <- if (exists[k]) readRDS(path) else NULL
      successful[k] <- exists[k] && isTRUE(z$success)
      trace_errors[k] <- if (exists[k]) z$trace_error %||% "" else "missing_checkpoint"
      if (successful[k]) { fits[[k]] <- z$fit; traces[k] <- z$core_trace }
      ci <- ci + 1L
      chain_status[[ci]] <- cbind(task, data.frame(checkpointed = exists[k], success = successful[k],
        error = if (!exists[k]) "missing_checkpoint" else z$error %||% "",
        trace_error = if (exists[k]) z$trace_error %||% "" else "missing_checkpoint"))
    }
    diagnostic_error <- NULL; diagnostics <- NULL
    if (all(successful)) diagnostics <- tryCatch({
      if (any(nzchar(trace_errors))) stop(paste(unique(trace_errors[nzchar(trace_errors)]), collapse = "; "))
      .native_diagnostics(traces, plan$config, plan$trace_contract, plan$thresholds)
    }, error = function(e) { diagnostic_error <<- conditionMessage(e); NULL })
    diagnostic_pass <- .native_diagnostics_pass(diagnostics)
    core_assessment <- .native_core_assessment(diagnostics)
    coherence <- if (!all(exists)) "incomplete" else if (all(successful) && !any(nzchar(trace_errors)))
      "verified_ordinary_iteration_layout" else "unverified_restart_or_trace_error"
    accuracy_eligible <- diagnostic_pass && identical(coherence, "verified_ordinary_iteration_layout")
    if (!is.null(diagnostics)) {
      di <- di + 1L; core_diagnostics[[di]] <- cbind(input, diagnostics)
    }
    pip <- if (all(successful)) rowMeans(do.call(cbind, lapply(fits, `[[`, "pip"))) else NULL
    if (all(successful)) {
      ids <- lapply(fits, function(f) f$additional$variant_ids)
      map <- read.csv(file.path(root, "prepared", draw_id, "variant-map.csv"), stringsAsFactors = FALSE)
      if (!all(vapply(ids, identical, TRUE, map$SNP))) stop("Native chain PIP variant IDs differ")
      pip_matrix <- do.call(cbind, lapply(fits, `[[`, "pip"))
      stability <- .native_pip_stability(pip_matrix, map$region_id)
      bi <- bi + 1L; block_stability[[bi]] <- cbind(input, stability$blocks)
      .campaign_atomic_rds(list(input = input, variant_ids = map$SNP, pip = pip,
        pip_chains = pip_matrix, core_diagnostics_pass = diagnostic_pass,
        native_core_assessment_status = core_assessment,
        native_trace_coherence = coherence, native_accuracy_eligible = accuracy_eligible,
        pip_stability = stability,
        convergence_scope = "CoreParameters only; marginal PIP variation is a stability diagnostic, not SNP-level ESS"),
        file.path(analysis, "consensus", paste0(input_id, ".rds")))
    } else stability <- list(summary = data.frame(median_pip_chain_range = NA_real_,
      p95_pip_chain_range = NA_real_, maximum_pip_chain_range = NA_real_,
      mean_pip_chain_range = NA_real_, fraction_call_095_disagreement = NA_real_))
    ii <- ii + 1L
    input_diagnostics[[ii]] <- cbind(input, data.frame(n_chains_checkpointed = sum(exists),
      n_chains_successful = sum(successful), core_diagnostics_pass = diagnostic_pass,
      native_core_assessment_status = core_assessment,
      native_trace_coherence = coherence, native_accuracy_eligible = accuracy_eligible,
      diagnostic_error = diagnostic_error %||% if (all(successful)) "" else "fewer_than_four_successful_chains",
      extension_eligible = !diagnostic_pass,
      extension_reason = if (diagnostic_pass) "not_eligible" else if (core_assessment != "unavailable")
        core_assessment else coherence), stability$summary)
    j <- j + 1L
    native_flags <- data.frame(core_diagnostics_pass = diagnostic_pass,
      native_core_assessment_status = core_assessment,
      native_trace_coherence = coherence, native_accuracy_eligible = accuracy_eligible)
    native_stability <- stability$summary; names(native_stability) <- paste0("native_", names(native_stability))
    rows[[j]] <- .native_metric_rows(sim, input, "sbayesrc_native", pip = pip,
      success = all(successful), attempted = any(exists), complete = all(exists),
      error = if (all(successful)) "" else "not_four_successful_chains",
      runtime = if (all(successful)) sum(vapply(fits, `[[`, 1.0, "runtime_seconds")) else NA_real_,
      pip_event = "any_nonzero_mixture_component", native_diagnostics = cbind(native_flags, native_stability),
      posterior_support = data.frame(posterior_conditioning = "joint_whole_panel_mixture_model",
        native_null_model_in_support = NA, null_probability_not_provided = TRUE))
    regional_tasks <- plan$regional[plan$regional$input_id == input_id, ]
    for (k in seq_len(nrow(regional_tasks))) {
      task <- regional_tasks[k, ]; path <- file.path(root, "regional", task$regional_id, "completed.rds")
      complete <- .native_completed_checkpoint(path, plan, task, task$regional_id)
      z <- if (complete) readRDS(path) else NULL
      ri <- ri + 1L; regional_status[[ri]] <- cbind(task, input[setdiff(names(input), "input_id")],
        data.frame(checkpointed = complete, error = if (complete) "" else "missing_checkpoint"))
      j <- j + 1L
      rows[[j]] <- .native_regional_rows(sim, input, task$method,
                      if (complete) z$result else NULL, attempted = complete)
    }
  }
  rows <- do.call(rbind, rows)
  utils::write.table(rows, metrics_file, sep = ",", row.names = FALSE, col.names = first,
                     append = !first, qmethod = "double")
  first <- FALSE
  ai <- ai + 1L; aggregates[[ai]] <- .native_aggregate_draw(rows)
}
# Explicit plan coverage survives unprepared draws and scheduler-killed workers.
status <- plan$chains
status$checkpointed <- vapply(seq_len(nrow(status)), function(i)
  .native_completed_checkpoint(file.path(root, "chains", status$chain_id[i], "completed.rds"),
                               plan, plan$chains[i, ], status$chain_id[i]), TRUE)
status$prepared <- status$draw_id %in% verified_prepared
write.csv(status, file.path(analysis, "planned-chain-coverage.csv"), row.names = FALSE)
regional_coverage <- plan$regional
regional_coverage$checkpointed <- vapply(seq_len(nrow(regional_coverage)), function(i)
  .native_completed_checkpoint(file.path(root, "regional", regional_coverage$regional_id[i], "completed.rds"),
                               plan, plan$regional[i, ], regional_coverage$regional_id[i]), TRUE)
write.csv(regional_coverage, file.path(analysis, "planned-regional-coverage.csv"), row.names = FALSE)
coverage$n_prepared <- length(verified_prepared)
coverage$n_checkpointed_chains <- sum(status$checkpointed)
coverage$n_checkpointed_regional <- sum(regional_coverage$checkpointed)
collection_status <- if (coverage$n_prepared == coverage$n_planned_draws &&
  coverage$n_checkpointed_chains == coverage$n_planned_chains &&
  coverage$n_checkpointed_regional == coverage$n_planned_regional) "complete" else "incomplete"
if (length(chain_status)) write.csv(do.call(rbind, chain_status), file.path(analysis, "chain-status.csv"), row.names = FALSE)
if (length(regional_status)) write.csv(do.call(rbind, regional_status), file.path(analysis, "regional-status.csv"), row.names = FALSE)
if (length(core_diagnostics)) write.csv(do.call(rbind, core_diagnostics), file.path(analysis, "core-diagnostics.csv"), row.names = FALSE)
if (length(block_stability)) write.csv(do.call(rbind, block_stability), file.path(analysis, "pip-stability-by-block.csv"), row.names = FALSE)
if (length(input_diagnostics)) {
  ds <- do.call(rbind, input_diagnostics)
  write.csv(ds, file.path(analysis, "input-diagnostics.csv"), row.names = FALSE)
  write.csv(ds[ds$extension_eligible, ], file.path(analysis, "extension-candidates.csv"), row.names = FALSE)
}
if (!length(aggregates)) { message("No receipted prepared scientific inputs; plan coverage retained"); return(invisible(NULL)) }
draw_metrics <- do.call(rbind, aggregates)
write.csv(draw_metrics, file.path(analysis, "draw-metrics.csv"), row.names = FALSE)
comparisons <- list(); pi <- 0L
stratum <- interaction(draw_metrics$model, draw_metrics$pi, draw_metrics$h2, draw_metrics$event, drop = TRUE)
for (indices in split(seq_len(nrow(draw_metrics)), stratum)) {
  d <- draw_metrics[indices, ]
  for (metric in c("ap_successful", "brier_successful", "recall_all_completed")) {
    for (exposure in c("none", "informative", "permuted")) {
      a <- d[d$method == "sbayesrc_native" & d$exposure == exposure, ]
      for (method in setdiff(unique(d$method), "sbayesrc_native")) {
        b <- d[d$method == method & d$exposure == exposure, ]
        pi <- pi + 1L
        comparisons[[pi]] <- cbind(d[1L, c("model", "pi", "h2", "event")],
          data.frame(method_a = "sbayesrc_native", method_b = method, exposure_a = exposure, exposure_b = exposure),
          .native_paired_draws(a, b, metric, as.integer(opt$`bootstrap-draws`), plan$seed + pi))
      }
    }
    for (method in unique(d$method)) for (control in c("none", "permuted")) {
      a <- d[d$method == method & d$exposure == "informative", ]; b <- d[d$method == method & d$exposure == control, ]
      pi <- pi + 1L
      comparisons[[pi]] <- cbind(d[1L, c("model", "pi", "h2", "event")],
        data.frame(method_a = method, method_b = method, exposure_a = "informative", exposure_b = control),
        .native_paired_draws(a, b, metric, as.integer(opt$`bootstrap-draws`), plan$seed + pi))
    }
  }
}
if (length(comparisons)) write.csv(do.call(rbind, comparisons), file.path(analysis, "paired-draw-comparisons.csv"), row.names = FALSE)
message("Collected explicit coverage, native chain diagnostics, event-typed scores and genotype-draw contrasts. ",
        "Core trace diagnostics do not establish convergence of unrecorded SNP indicator trajectories.")

}, error = function(e) {
  collection_status <<- "failed"; collection_error <<- conditionMessage(e)
  stop(e)
})
}
main()
