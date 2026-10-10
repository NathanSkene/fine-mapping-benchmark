#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
main <- function() {
  opt <- .native_options(list(plan = "results/native-campaign/plan.rds", index = "1", `input-id` = ""))
  plan <- readRDS(opt$plan); .native_assert_plan(plan, .native_repo)
  if (!identical(plan$runtime_fingerprint, .campaign_runtime_fingerprint(plan$method_args)))
    stop("Frozen regional tool/runtime dependencies changed")
  if (nzchar(opt$`input-id`)) {
    indices <- which(plan$regional$input_id == opt$`input-id`)
    if (!length(indices)) stop("Unknown input ID")
  } else {
    indices <- as.integer(opt$index)
    if (is.na(indices) || !indices %in% seq_len(nrow(plan$regional))) stop("Invalid regional index")
  }
  first_task <- plan$regional[indices[1L], ]; input <- plan$inputs[match(first_task$input_id, plan$inputs$input_id), ]
  root <- dirname(normalizePath(opt$plan))
  prepared <- .native_read_prepared(plan, root, input$draw_id)
  source_input <- prepared$inputs[match(input$input_id, prepared$inputs$input_id), ]
  sim <- readRDS(file.path(root, "prepared", input$draw_id, "simulation.rds"))
  controls <- readRDS(file.path(root, "prepared", input$draw_id, "annotation-controls.rds"))
  sim <- .native_expose_annotations(sim, source_input$scenario, input$exposure, controls)
  for (i in indices) {
    task <- plan$regional[i, ]; directory <- file.path(root, "regional", task$regional_id)
    checkpoint <- file.path(directory, "completed.rds")
    if (.campaign_validate_checkpoint(checkpoint, plan, task)) {
      .campaign_check_receipt(directory, checkpoint, plan$plan_sha256, task$regional_id)
      .campaign_receipt(directory, checkpoint, plan$plan_sha256, task$regional_id)
      next
    }
    .campaign_claim(paste0(directory, ".lock"))
    tryCatch({
      set.seed(plan$seed + 3000000L + i)
      started <- Sys.time()
      result <- run_methods(sim, task$method, stats::setNames(list(plan$method_args[[task$method]]), task$method), verbose = FALSE)
      result <- .campaign_reject_substitution(result, sim, task$method)
      .campaign_atomic_rds(list(terminal = TRUE, task = task, plan_sha256 = plan$plan_sha256,
        result = result, input = input, started = started, completed = Sys.time(), session = utils::sessionInfo()), checkpoint)
      .campaign_receipt(directory, checkpoint, plan$plan_sha256, task$regional_id)
    }, finally = unlink(paste0(directory, ".lock"), recursive = TRUE))
  }
}
main()
