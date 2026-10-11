#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
main <- function() {
  opt <- .native_options(list(plan = "results/native-campaign/plan.rds", index = "1"))
  plan <- readRDS(opt$plan); .native_assert_plan(plan, .native_repo, verify_source_data = TRUE)
  i <- as.integer(opt$index); if (is.na(i) || !i %in% seq_len(nrow(plan$draws))) stop("Invalid draw index")
  draw <- plan$draws[i, ]; root <- dirname(normalizePath(opt$plan))
  directory <- file.path(root, "prepared", draw$draw_id); checkpoint <- file.path(directory, "completed.rds")
  if (file.exists(checkpoint)) {
    .native_read_prepared(plan, root, draw$draw_id, require_receipt = FALSE)
    .campaign_receipt(directory, checkpoint, plan$plan_sha256, draw$draw_id)
    return(invisible(NULL))
  }
  .campaign_claim(paste0(directory, ".lock")); on.exit(unlink(paste0(directory, ".lock"), recursive = TRUE))
  sim <- simulate_gwfm_data(n = plan$config$n, p = plan$config$p, n_iter = plan$config$replicates,
    pi = plan$config$pi, h2 = plan$config$h2, model = draw$model, p_causal = .5,
    inf_model = "beatrice", annotations = "binary", n_annotations = 3L,
    annotation_proportions = c(.1, .2, .3), enrichment = rep(5, 3L),
    regions = plan$regions, vcf_dir = plan$vcf_dir, genetic_map_dir = plan$map_dir,
    seed = draw$seed, verbose = TRUE)
  z <- .native_prepare_draw(sim, draw, plan, directory)
  .campaign_atomic_rds(z, checkpoint)
  .campaign_receipt(directory, checkpoint, plan$plan_sha256, draw$draw_id)
}
main()
