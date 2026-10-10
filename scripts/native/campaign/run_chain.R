#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
main <- function() {
  opt <- .native_options(list(plan = "results/native-campaign/plan.rds", index = "1"))
  plan <- readRDS(opt$plan); .native_assert_plan(plan, .native_repo)
  if (!identical(plan$executable_sha256, .challenge_sha256(plan$executable))) stop("Pinned executable changed")
  i <- as.integer(opt$index); if (is.na(i) || !i %in% seq_len(nrow(plan$chains))) stop("Invalid chain index")
  task <- plan$chains[i, ]; root <- dirname(normalizePath(opt$plan))
  directory <- file.path(root, "chains", task$chain_id); checkpoint <- file.path(directory, "completed.rds")
  if (.campaign_validate_checkpoint(checkpoint, plan, task)) {
    .campaign_check_receipt(directory, checkpoint, plan$plan_sha256, task$chain_id)
    .campaign_receipt(directory, checkpoint, plan$plan_sha256, task$chain_id)
    return(invisible(NULL))
  }
  prepared <- .native_read_prepared(plan, root, task$draw_id, include_simulation = FALSE)
  input <- prepared$inputs[match(task$input_id, prepared$inputs$input_id), ]
  .campaign_claim(paste0(directory, ".lock")); on.exit(unlink(paste0(directory, ".lock"), recursive = TRUE))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  # A scheduler interruption can leave upstream logs and its inner lock.
  # The controller removes the outer lock only after that worker has stopped.
  # Once this fresh worker owns it, preserve the old attempt before rerunning.
  .native_archive_interrupted_output(directory, plan, task)
  provenance <- list(gwas_source = "Corrected simulate_gwfm_data; one global phenotype; exact marginal OLS effects",
    ld_source = "Empirical within-block cor(X); full positive eigenspectrum; source-verified float32 GCTB layout; finite cross-block sample LD omitted",
    annotation_source = paste("Truth generated with three binary annotations enrichment5; exposed control", task$exposure),
    genome_build = "GRCh37 retained source VCF coordinates", ancestry = "EUR empirical source haplotype templates simulated by sim1000G",
    panel_id = paste(task$draw_id, task$input_id, plan$plan_sha256, sep = ":"),
    preprocessing = "Stable CHROM:POS:REF:ALT IDs; ALT dosage is A1; standardized effects and SE divided by empirical raw dosage SD; no SNP imputation")
  started <- Sys.time(); error <- NULL
  fit <- tryCatch(run_sbayesrc_native(input$gwas_file, input$ld_dir, input$annotation_file,
    file.path(directory, "native-output"), scale = "subset", provenance = provenance,
    gctb_path = plan$executable, expected_version = plan$expected_version,
    n_iter = plan$config$n_iter, burn_in = plan$config$burn_in, thin = plan$config$thin,
    seed = task$seed, threads = 1L, save_mcmc = TRUE, mcmc_scope = "core"),
    error = function(e) { error <<- conditionMessage(e); NULL })
  core <- file.path(directory, "native-output", "sbayesrc_native.mcmcsamples", "CoreParameters.mcmcsamples.txt")
  trace_error <- NULL; extraction <- NULL
  retained_path <- file.path(directory, "native-output", "retained-core.tsv")
  restart_markers <- .native_trace_restart_markers(file.path(directory, "native-output"))
  if (!is.null(fit)) tryCatch({
    if (nrow(restart_markers) %||% 0L)
      stop("Native restart/disruption: final summary and core trace coherence is unverified")
    retained <- .native_read_core(core, plan$config, plan$trace_contract)
    utils::write.table(retained, retained_path, sep = "\t", row.names = FALSE, quote = FALSE)
    extraction <- attr(retained, "extraction")
    extraction$retained_core_sha256 <- .challenge_sha256(retained_path)
    jsonlite::write_json(extraction, file.path(directory, "native-output", "trace-extraction.json"),
                           auto_unbox = TRUE, pretty = TRUE)
  }, error = function(e) { trace_error <<- conditionMessage(e) })
  artifacts <- sort(list.files(file.path(directory, "native-output"), recursive = TRUE,
                               full.names = TRUE), method = "radix")
  artifacts <- artifacts[!dir.exists(artifacts)]
  sinks <- artifacts[nzchar(Sys.readlink(artifacts))]
  if (length(sinks) && (any(Sys.readlink(sinks) != "/dev/null") ||
      any(!basename(sinks) %in% .sbayesrc_native_text_sinks())))
    stop("Unexpected native-output symbolic link")
  artifacts <- setdiff(artifacts, sinks)
  result <- list(terminal = TRUE, task = task, plan_sha256 = plan$plan_sha256,
    success = !is.null(fit), fit = fit, error = error, core_trace = core,
    trace_error = trace_error, started = started, completed = Sys.time(),
    retained_core_trace = if (is.null(extraction)) NULL else retained_path,
    trace_extraction = extraction, restart_markers = restart_markers,
    text_retention = if (is.null(fit)) NULL else fit$additional$manifest$text_retention,
    discarded_array_links = data.frame(path = sinks, target = rep("/dev/null", length(sinks))),
    artifacts = data.frame(path = unname(artifacts),
      sha256 = unname(vapply(artifacts, .challenge_sha256, "")), row.names = NULL),
    session = utils::sessionInfo(), convergence_claim = "not established by execution")
  .campaign_atomic_rds(result, checkpoint)
  .campaign_receipt(directory, checkpoint, plan$plan_sha256, task$chain_id)
  if (!is.null(error)) message("Native attempted failure retained: ", error)
}
main()
