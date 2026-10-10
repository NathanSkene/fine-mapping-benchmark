#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
main <- function() {
opt <- .campaign_options(list(plan = "results/full-campaign/plan.rds", index = "1",
                              `require-preflight` = "true"))
plan <- readRDS(opt$plan)
.campaign_assert_plan(plan, .campaign_repo)
index <- as.integer(opt$index)
if (is.na(index) || index < 1L || index > nrow(plan$tasks)) stop("Invalid task index")
task <- plan$tasks[index, ]
root <- dirname(normalizePath(opt$plan))
directory <- file.path(root, "tasks", task$task_id)
checkpoint <- file.path(directory, "completed.rds")
if (.campaign_validate_checkpoint(checkpoint, plan, task)) {
  .campaign_check_receipt(directory, checkpoint, plan$plan_sha256, task$task_id)
  if (is.null(readRDS(checkpoint)$artifacts)) stop("Completed task lacks sealed dataset artifacts")
  .campaign_receipt(directory, checkpoint, plan$plan_sha256, task$task_id)
  cat("Method task already complete and compatible\n"); quit(save = "no", status = 0L)
}
if (identical(opt$`require-preflight`, "true") && !task$method %in% c("exact", "forced_one")) {
  preflight_path <- file.path(root, "preflight", task$method, "completed.rds")
  if (!file.exists(preflight_path)) stop("Required preflight is missing for ", task$method)
  if (!.campaign_check_receipt(dirname(preflight_path), preflight_path, plan$plan_sha256, task$method))
    stop("Required preflight receipt is missing for ", task$method)
  preflight <- readRDS(preflight_path)
  if (!identical(preflight$plan_sha256, plan$plan_sha256) || !isTRUE(preflight$success))
    stop("Required native/dependency preflight did not succeed for ", task$method)
}
.campaign_claim(paste0(directory, ".lock"))
on.exit(unlink(paste0(directory, ".lock"), recursive = TRUE), add = TRUE)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
ids <- strsplit(task$dataset_ids, "|", fixed = TRUE)[[1L]]
datasets <- plan$datasets[match(ids, plan$datasets$dataset_id), ]
if (anyNA(datasets$dataset_id)) stop("Task references an unknown dataset")
original <- if (nzchar(plan$original_source))
    .campaign_load_source(plan$original_source) else NULL
args <- plan$method_args[[task$method]] %||% list()
if (task$stage == "original_refit" && task$method %in%
       c("beatrice", "functional_beatrice", "fb_pooled", "fb_xregion"))
  args$beatrice_dir <- file.path(plan$original_source, "BEATRICE_annot_sparse")
started <- Sys.time(); n_new <- 0L
cached_file <- ""; cached_sim <- NULL; cached_hash <- NULL; cached_genotype_hash <- NULL
for (i in seq_len(nrow(datasets))) {
  d <- datasets[i, ]; artifact <- file.path(directory, paste0(d$dataset_id, ".rds"))
  receipt_name <- paste0(d$dataset_id, ".terminal.json")
  receipt_id <- paste(task$task_id, d$dataset_id, sep = ".")
  identity <- list(task = task, dataset = d)
  if (.campaign_validate_checkpoint(artifact, plan, identity)) {
    if (.campaign_check_receipt(directory, artifact, plan$plan_sha256, receipt_id, receipt_name)) next
    # An RDS without its commit receipt is an interrupted write, not reusable
    # evidence. Preserve it outside the collection namespace and recompute.
    unsealed <- file.path(directory, "unsealed")
    dir.create(unsealed, showWarnings = FALSE)
    if (!file.rename(artifact, file.path(unsealed, paste0(basename(artifact), ".", Sys.getpid(),
          ".", format(Sys.time(), "%Y%m%d%H%M%OS6"))))) stop("Cannot preserve unsealed dataset checkpoint")
  }
  .campaign_start_attempt(directory, plan, task, d)
  if (d$sim_file != cached_file) {
    if (d$sim_sha256 == "generated") {
      generation_path <- file.path(dirname(d$sim_file), "generated.rds")
      generation_id <- plan$generation$generation_id[match(dirname(d$sim_file),
          file.path(root, "datasets", sprintf("generation_%05d", plan$generation$generation_id)))]
      if (is.na(generation_id) || !.campaign_check_receipt(dirname(d$sim_file), generation_path,
             plan$plan_sha256, paste0("generation_", generation_id)))
        stop("Dataset generation receipt is missing")
      generated <- readRDS(generation_path)
      if (!isTRUE(generated$terminal) || !identical(generated$plan_sha256, plan$plan_sha256))
        stop("Dataset generation is incomplete or belongs to another plan")
      expected <- generated$artifacts$sha256[match(d$sim_file, generated$artifacts$path)]
      cached_sim <- .campaign_read_frozen(d$sim_file, expected)
      cached_genotype_hash <- generated$genotypes_sha256
    } else cached_sim <- .campaign_read_frozen(d$sim_file, d$sim_sha256)
    if (d$sim_sha256 != "generated") cached_genotype_hash <- d$genotypes_sha256
    cached_file <- d$sim_file; cached_hash <- .challenge_sha256(d$sim_file)
  }
  began <- Sys.time(); metadata <- c(as.list(d), list(stage = task$stage,
                          input_panel_id = cached_hash))
  if (task$stage == "challenge") {
    panel <- cached_sim
    pred <- .campaign_challenge_predict(panel, seq_along(panel$inputs), task$method, args)
    ev <- evaluate_challenge(panel, stats::setNames(list(pred), task$method))
    table <- .campaign_challenge_task_metrics(ev$tasks, pred, panel, task$method)
    ev$tasks <- table
    for (si in seq_len(nrow(ev$summary))) {
      selected <- table[table$scenario == ev$summary$scenario[si], ]
      for (field in c("native_all_causal_coverage_reported", "native_all_causal_recovery_all_attempts",
                       "cs_causal_variant_recall_all_attempts", "cs_support_size"))
        ev$summary[[field]][si] <- .challenge_safe_mean(selected[[field]])
      ev$summary$cs_distinct_recall_all_attempts[si] <- .challenge_safe_mean(selected$cs_distinct_recall)
      ev$summary$log_loss_successful[si] <- .challenge_safe_mean(selected$log_loss)
      ev$summary$n_log_loss_infinite[si] <- sum(is.infinite(selected$log_loss))
    }
    if (task$method == "marginal_z") {
      for (field in c("brier", "log_loss_clipped", "null_brier")) table[[field]] <- NA_real_
      for (field in c("brier_successful", "log_loss_clipped_successful", "null_brier_successful"))
        ev$summary[[field]] <- NA_real_
      ev$tasks <- table
    }
    table$dataset_id <- d$dataset_id; table$family <- "challenge"; table$stage <- "challenge"
    table$panel_cluster_id <- cached_hash
    table$cluster_id <- paste(d$dataset_id, table$cluster_id, sep = ":")
    table$p <- panel$config$p; table$n <- panel$config$n
    table$actual_p <- vapply(panel$inputs[match(table$task_id, names(panel$inputs))],
                              function(x) length(x$z), 1L)
    fit <- pred; scored_original <- NULL
  } else {
    sim <- .campaign_one_scenario(cached_sim, d$scenario)
    if (task$stage == "corrected_resimulation") {
      seed <- plan$seed + 10000000L + match(d$dataset_id, plan$datasets$dataset_id)
      sim <- .campaign_corrected_phenotype(cached_sim, d$scenario, seed)
    }
    cache_manifest <- list(reused = FALSE, eligible = FALSE)
    if (task$stage == "archived_rescore") {
      archived <- .campaign_merge_archived(d$result_files, d$result_sha256)
      if (is.null(archived) || !task$method %in% archived$methods_run) {
        # Absent archive is missing evidence, never reconstructed old output.
        fit <- list(); fit[[task$method]] <- list(results = list())
        fit$methods_run <- task$method
      } else {
        fit <- archived[c(task$method, "methods_run")]; fit$methods_run <- task$method
        # Archives may contain a full row or a one-scenario unit. Normalize IDs.
        selected <- fit[[task$method]]$results
        if (length(unique(vapply(selected, function(f) f$scenario_id, 1L))) > 1L)
          selected <- Filter(function(f) identical(as.integer(f$scenario_id), as.integer(d$scenario)), selected)
        fit[[task$method]]$results <- lapply(selected, function(f) { f$scenario_id <- 1L; f })
      }
      scored_original <- tryCatch(original$evaluate_methods(sim, fit, verbose = FALSE),
                           error = function(e) list(error = conditionMessage(e)))
    } else {
      runner <- if (task$stage == "original_refit") original$run_methods else run_methods
      inference_seed <- plan$seed + 2000000L + 100L * match(d$dataset_id, plan$datasets$dataset_id) +
                           match(task$method, plan$methods)
      source_env <- if (task$stage == "original_refit") original else globalenv()
      cache_manifest <- .campaign_fit_cache_identity(sim, task$method, args, inference_seed,
                           source_env, plan$runtime_fingerprint, genotypes_sha256 = cached_genotype_hash)
      cached_fit_path <- if (isTRUE(cache_manifest$eligible)) file.path(root, "fit-cache",
                              paste0(cache_manifest$key, ".rds")) else ""
      cached_fit_receipt <- paste0(cache_manifest$key %||% "ineligible", ".terminal.json")
      cache_sealed <- nzchar(cached_fit_path) &&
        .campaign_check_receipt(dirname(cached_fit_path), cached_fit_path, plan$plan_sha256,
              cache_manifest$key, cached_fit_receipt)
      if (cache_sealed) {
        cached_fit <- readRDS(cached_fit_path)
        if (!identical(cached_fit$identity, cache_manifest) ||
            !identical(cached_fit$plan_sha256, plan$plan_sha256)) stop("Incompatible fit cache")
        fit <- cached_fit$fit; cache_manifest$reused <- TRUE
        cache_manifest$original_compute_source <- cached_fit$metadata
      } else {
        set.seed(inference_seed)
        fit <- runner(sim, methods = task$method,
               method_args = stats::setNames(list(args), task$method), verbose = FALSE)
        if (nzchar(cached_fit_path)) {
          .campaign_atomic_rds(list(identity = cache_manifest, fit = fit,
             plan_sha256 = plan$plan_sha256, metadata = metadata), cached_fit_path)
          .campaign_receipt(dirname(cached_fit_path), cached_fit_path, plan$plan_sha256,
                cache_manifest$key, cached_fit_receipt)
        }
        cache_manifest$reused <- FALSE
      }
    }
    # Isolate metric changes on identical raw predictions. The acceptance
    # policy contrast has separate old/new evaluations after rejection.
    scored_original <- if (!is.null(original)) tryCatch(original$evaluate_methods(sim, fit, verbose = FALSE),
                         error = function(e) list(error = conditionMessage(e))) else NULL
    scored_corrected_raw <- evaluate_methods(sim, fit, verbose = FALSE)
    fit <- .campaign_reject_substitution(fit, sim, task$method)
    ev <- evaluate_methods(sim, fit, verbose = FALSE)
    scored_original_accepted <- if (!is.null(original))
       tryCatch(original$evaluate_methods(sim, fit, verbose = FALSE),
                error = function(e) list(error = conditionMessage(e))) else NULL
    table <- .campaign_task_table(fit, sim, task$method, metadata)
    table$fit_reused <- isTRUE(cache_manifest$reused)
  }
  if (task$stage == "challenge") {
    cache_manifest <- NULL; scored_corrected_raw <- ev; scored_original_accepted <- NULL
    table$likelihood_interpretation <- if (task$method == "finemap_inf")
       "normal_score_requested; unsupported phenotype second moment required" else
       if (task$method == "susie_inf")
       "normal_score_requested; unsupported fixed-noise infinitesimal variance estimator" else
       if (task$method %in% c("susie", "polyfun_est", "polyfun_ldsc", "sbayesrc"))
       "normal_score; specified residual_variance=1; method-specific prior" else "declared_method_working_model"
  }
  output <- list(terminal = TRUE, plan_sha256 = plan$plan_sha256, task = identity,
      method_task = task, result = fit, evaluation = ev, original_evaluation = scored_original,
      corrected_raw_evaluation = scored_corrected_raw, original_accepted_evaluation = scored_original_accepted,
      reuse = cache_manifest,
      table = table, metadata = metadata, started = began, completed = Sys.time(),
      elapsed_seconds = as.numeric(difftime(Sys.time(), began, units = "secs")),
      session = utils::sessionInfo())
  .campaign_atomic_rds(output, artifact)
  .campaign_receipt(directory, artifact, plan$plan_sha256, receipt_id, receipt_name)
  n_new <- n_new + 1L
  cat(sprintf("%s %s %s: %d/%d successful in %.1fs\n", task$stage, task$method,
      d$dataset_id, sum(table$success), nrow(table), output$elapsed_seconds))
}
.campaign_atomic_rds(list(terminal = TRUE, plan_sha256 = plan$plan_sha256, task = task,
    n_datasets = nrow(datasets), n_new = n_new, completed = Sys.time(),
    artifacts = .campaign_artifact_manifest(c(file.path(directory, paste0(ids, ".rds")),
                      file.path(directory, paste0(ids, ".terminal.json")))),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))), checkpoint)
.campaign_receipt(directory, checkpoint, plan$plan_sha256, task$task_id)
cat(sprintf("Task %s terminal; all %d datasets checkpointed\n", task$task_id, nrow(datasets)))
}
main()
