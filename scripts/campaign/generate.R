#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
main <- function() {
opt <- .campaign_options(list(plan = "results/full-campaign/plan.rds", index = "1"))
plan <- readRDS(opt$plan)
.campaign_assert_plan(plan, .campaign_repo)
index <- as.integer(opt$index)
if (is.na(index) || index < 1L || index > nrow(plan$generation)) stop("Invalid generation index")
row <- plan$generation[index, ]
directory <- file.path(dirname(normalizePath(opt$plan)), "datasets",
                       sprintf("generation_%05d", row$generation_id))
checkpoint <- file.path(directory, "generated.rds")
if (.campaign_validate_checkpoint(checkpoint, plan, row)) {
  .campaign_check_receipt(directory, checkpoint, plan$plan_sha256, paste0("generation_", index))
  .campaign_receipt(directory, checkpoint, plan$plan_sha256, paste0("generation_", index))
  cat("Generation already complete and verified\n"); quit(save = "no", status = 0L)
}
.campaign_claim(paste0(directory, ".lock"))
on.exit(unlink(paste0(directory, ".lock"), recursive = TRUE), add = TRUE)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
started <- Sys.time()
generation_dimension_audit <- NULL
if (row$family == "challenge") {
  panel <- generate_challenge(n_replicates = row$replicates, p = row$p, n = row$n,
                               seed = row$seed, split = "final")
  .campaign_atomic_rds(panel, file.path(directory, "panel.rds"))
  write_challenge_release(panel, file.path(directory, "public"), file.path(directory, "assessor"),
                           release_id = paste0("campaign-", row$generation_id))
  artifacts <- c(file.path(directory, "panel.rds"),
                 file.path(directory, "public/manifest.csv"),
                 file.path(directory, "assessor/manifest.csv"))
} else {
  original <- .campaign_load_source(plan$original_source)
  g <- plan$grid[row$row, ]
  split_integer <- function(s) as.integer(strsplit(s, "|", fixed = TRUE)[[1L]])
  split_numeric <- function(s) as.numeric(strsplit(s, "|", fixed = TRUE)[[1L]])
  args <- list(n_regions = g$n_regions, n = g$n, p = split_integer(g$p_values),
      n_iter = g$n_iter, S = split_integer(g$S_values), phi = split_numeric(g$phi_values),
      model = g$model, annotations = g$annotation_type,
      n_annotations = g$n_annotations, enrichment = if (g$annotation_type == "none") NULL else
        split_numeric(g$enrichment_values),
      annotation_correlation = if (is.na(g$annotation_correlation)) 0 else g$annotation_correlation,
      seed = row$seed, vcf_dir = if (nzchar(plan$data_paths$vcf_dir)) plan$data_paths$vcf_dir else NULL,
      genetic_map_dir = plan$data_paths$genetic_map_dir, verbose = FALSE)
  if (!is.na(g$p_causal)) args$p_causal <- g$p_causal
  if (!is.na(g$n_ref)) args$n_ref <- g$n_ref
  # Extdata lookup is relative to the snapshot root in source-only mode.
  old <- getwd(); setwd(plan$original_source)
  captured <- tryCatch(.campaign_capture_generation(function()
    do.call(original$run_simulation, args)), finally = setwd(old))
  sim <- captured$value
  used <- sim$params$vcf_files_used
  if (length(used) && nzchar(plan$data_paths$vcf_dir))
    .campaign_check_inputs(used, plan$input_file_manifest)
  maps <- list.files(plan$data_paths$genetic_map_dir, full.names = TRUE, recursive = TRUE)
  if (length(maps)) .campaign_check_inputs(maps, plan$input_file_manifest)
  requested <- split_integer(g$p_values)
  generation_dimension_audit <- .campaign_generation_dimension_audit(sim, requested, g$n,
    n_ref = if (is.na(g$n_ref)) NULL else g$n_ref,
    source_file_sha256 = .challenge_sha256(file.path(plan$original_source, "R/simulate_genotypes.R")),
    warnings = captured$warnings)
  sim$params$campaign_generation_source <- "original code reconstruction"
  sim$params$generation_dimension_audit <- generation_dimension_audit
  .campaign_atomic_rds(sim, file.path(directory, "simulation.rds"))
  artifacts <- file.path(directory, "simulation.rds")
}
provenance <- list(terminal = TRUE, plan_sha256 = plan$plan_sha256, task = row,
  started = started, completed = Sys.time(), session = utils::sessionInfo(),
  genotypes_sha256 = if (row$family == "challenge") NA_character_ else .campaign_object_sha(sim$genotypes),
  generation_dimension_audit = generation_dimension_audit,
  artifacts = data.frame(path = artifacts, sha256 = vapply(artifacts, .challenge_sha256, "")))
.campaign_atomic_rds(provenance, checkpoint)
.campaign_receipt(directory, checkpoint, plan$plan_sha256, paste0("generation_", index))
cat(sprintf("Generated %s in %.1f seconds\n", row$family,
            as.numeric(difftime(Sys.time(), started, units = "secs"))))
}
main()
