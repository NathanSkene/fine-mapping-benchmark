#!/usr/bin/env Rscript
# Synthetic software fixtures only; stage bootstrap checks are not ranking evidence.
source("paper/regional_publication.R")
checks <- 0L
check <- function(ok) { stopifnot(ok); checks <<- checks + 1L }
reject <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  check(inherits(error, "error") && grepl(pattern, conditionMessage(error), fixed = TRUE))
}
workspace <- Sys.getenv("FMB_PAPER_FIXTURE_DIR", unset = tempfile("paper-regional-publication-"))
if (dir.exists(workspace)) stop("Use a new fixture directory")
dir.create(workspace)
if (dir.exists(".local-lib")) .libPaths(c(".local-lib", .libPaths()))
for (path in list.files("R", full.names = TRUE)) if (endsWith(path, ".R")) source(path)
source("scripts/campaign/collect_streamed.R")

base <- data.frame(dataset_id = paste0("gaussian_fixture_", 1:5), family = "challenge", stage = "challenge",
  region_id = NA_integer_, task_id = paste0("task_", 1:5), scenario = "mixed_fixture", p = 5L, n = 50L,
  n_gwas = NA_integer_, n_ref = NA_integer_, nominal_p = NA_integer_, actual_p = c(4L, 5L, 4L, 5L, 5L),
  S = NA_integer_, phi = NA_real_, p_causal = NA_real_, annotation_type = NA_character_, enrichment_fold = NA_real_,
  n_causal = c(2L, 0L, 1L, 5L, 2L), n_observed_causal = c(2L, 0L, 0L, 5L, 2L),
  n_missing_causal = c(0L, 0L, 1L, 0L, 0L), null = c(FALSE, TRUE, FALSE, FALSE, FALSE),
  n_any_nonzero_total_effect = c(2L, 0L, 4L, 5L, 5L), success = TRUE,
  ap = c(.7, NA, NA, 1, .6), brier = c(.1, 1, .2, .05, .2),
  log_loss = c(.1, Inf, .2, .05, .2), recall = c(.5, NA, 0, 1, .5),
  cs_distinct_recall = c(.5, NA, 0, 1, .5), score_target = "sparse_component_membership",
  posterior_conditioning = "unspecified", cluster_id = paste0("c", 1:5), panel_cluster_id = "panel")
base$stratum <- paste("challenge", base$scenario, base$p, base$n, sep = ":")
legacy <- base[1:3, ]
legacy$dataset_id <- "legacy_fixture"; legacy$family <- "legacy_reconstruction"; legacy$stage <- "corrected_refit"
legacy$region_id <- 1:3; legacy$task_id <- legacy$scenario <- NA_character_
legacy$p <- legacy$n <- NA_integer_; legacy$n_gwas <- 100L; legacy$n_ref <- 20L
legacy$nominal_p <- 10L; legacy$actual_p <- c(8L, 9L, 10L)
legacy$S <- 2L; legacy$phi <- .2; legacy$p_causal <- 1
legacy$annotation_type <- "binary"; legacy$enrichment_fold <- 5
legacy$n_causal <- legacy$n_observed_causal <- legacy$n_any_nonzero_total_effect <- 2L
legacy$n_missing_causal <- 0L; legacy$null <- FALSE
legacy$ap <- c(.12345678901234566, .5, .4); legacy$brier <- c(1e-300, .2, .3)
legacy$log_loss <- .3; legacy$recall <- legacy$cs_distinct_recall <- .5
legacy$stratum <- do.call(paste, c(legacy[c("family", "n_gwas", "n_ref", "annotation_type", "enrichment_fold", "S", "phi", "p_causal", "nominal_p")], sep = ":"))
tasks <- do.call(rbind, lapply(c("abf", "sbayesrc"), function(method) {
  z <- rbind(base, legacy); z$method <- method
  z$pip_event <- if (method == "abf") "method_model_sparse_inclusion" else "any_nonzero_mixture_component"
  z$probability_target_aligned <- method == "abf"
  z$brier_any_nonzero_total_effect <- z$log_loss_any_nonzero_total_effect <- if (method == "sbayesrc") .2 else NA_real_
  if (method == "sbayesrc") {
    failed <- z$dataset_id == "gaussian_fixture_5"
    z$success[failed] <- FALSE; z$ap[failed] <- z$brier[failed] <- z$log_loss[failed] <- NA_real_
    z$brier_any_nonzero_total_effect[failed] <- z$log_loss_any_nonzero_total_effect[failed] <- NA_real_
    z$recall[failed] <- z$cs_distinct_recall[failed] <- 0
  }
  z
}))
summary <- do.call(rbind, lapply(split(seq_len(nrow(tasks)), .paper_regional_key(tasks)), function(index) {
  z <- tasks[index, ]; r <- z[1, c("stage", "stratum", "method")]
  r$n_attempted <- nrow(z); r$n_successful <- sum(z$success)
  r$n_ap_eligible_attempts <- sum(z$n_observed_causal > 0)
  r$n_ap_eligible_successful <- sum(z$success & z$n_observed_causal > 0)
  r$success_rate <- mean(z$success); r$ap_successful <- mean(z$ap, na.rm = TRUE)
  r$brier_successful <- mean(z$brier, na.rm = TRUE); r$fdr_micro <- .2
  r$log_loss_successful <- mean(z$log_loss, na.rm = TRUE)
  r$n_probability_target_misaligned_successful <- sum(z$success & !z$probability_target_aligned)
  r
}))
rownames(summary) <- NULL
contrasts <- do.call(rbind, lapply(unique(tasks$stratum), function(stratum) {
  z <- tasks[tasks$stratum == stratum, ]; a <- z[z$method == "abf", ]; b <- z[z$method == "sbayesrc", ]
  r <- data.frame(stage = a$stage[1], stratum = stratum, method_a = "abf", method_b = "sbayesrc",
    metric = c("ap", "brier", "recall", "cs_distinct_recall"), difference = .12345678901234566,
    lower = -.1, upper = .2, n_attempted_pairs = nrow(a), n_clusters = 1L, probability_targets_aligned = FALSE)
  r$n_pairs <- vapply(r$metric, function(metric) sum(is.finite(a[[metric]] - b[[metric]])), 1L)
  r
}))
path <- file.path(workspace, "task-metrics.csv"); schema_path <- file.path(workspace, "schema.rds")
write_tasks <- function(x, target = path) {
  schema <- list(schema = 1L, columns = names(x), storage_classes = .campaign_metric_classes(x),
    missing_marker = "<FIXTURE_MISSING>", numeric_precision = "17 significant decimal digits")
  saveRDS(schema, schema_path)
  .campaign_csv_block(x, target, append = FALSE, missing = schema$missing_marker)
}
write_tasks(tasks)
typed <- NULL
invisible(.paper_regional_scan(path, .paper_regional_schema(schema_path), function(z) typed <<- rbind(typed, z), chunk_size = 3L))
rownames(typed) <- rownames(tasks) <- NULL
metadata <- setdiff(names(tasks), c("ap", "brier", "log_loss", "recall", "cs_distinct_recall",
  "brier_any_nonzero_total_effect", "log_loss_any_nonzero_total_effect"))
check(identical(typed[metadata], tasks[metadata]))
check(identical(is.finite(typed$brier), is.finite(tasks$brier)))
check(identical(typed$ap, tasks$ap))
view <- .paper_regional_publication(path, schema_path, summary, contrasts, chunk_size = 3L)
check(identical(summary, view$rows[names(summary)]))
check(identical(contrasts, view$contrasts[names(contrasts)]))
check(view$streaming$peak_chunk_rows == 3L && view$streaming$peak_stratum_rows == 10L)
census <- view$census
get <- function(stage, method, endpoint) census[census$stage == stage & census$method == method & census$endpoint == endpoint, ]
z <- get("corrected_refit", "abf", "all_attempted")
check(z$actual_p_min == 8 && z$actual_p_max == 10 && z$n_dropped_candidates == 3)
check(z$synthetic_n_gwas_min == 100 && z$synthetic_n_ref_min == 20)
check(z$locus_target_prevalence_min == .2 && z$locus_target_prevalence_max == .25)
check(z$candidate_gap_fraction_max == .2)
z <- get("challenge", "abf", "all_attempted")
check(z$n_generating_null_loci == 1 && z$n_no_observed_cause_loci == 2 && z$n_missing_causal == 1)
check(grepl("gap_not_attributed_to_monomorphism", z$size_basis, fixed = TRUE))
check(get("challenge", "abf", "ap_successful")$n_loci == 3)
check(get("challenge", "sbayesrc", "ap_successful")$n_loci == 2)
check(get("challenge", "sbayesrc", "ap_observed_positive_attempts_zero_for_failure")$n_loci == 3)
check(get("challenge", "sbayesrc", "recall_all_attempts")$n_loci == 4)
check(get("challenge", "abf", "brier_successful")$n_loci == 5)
check(get("challenge", "abf", "log_loss_successful")$n_loci == 5)
check(get("challenge", "abf", "log_loss_successful")$n_infinite_endpoint_scores == 1)
check(grepl("actual_simulated_gwas", get("corrected_refit", "abf", "all_attempted")$n_basis, fixed = TRUE))
check(grepl("gaussian_score_design_parameter", get("challenge", "abf", "all_attempted")$n_basis, fixed = TRUE))
check(get("challenge", "abf", "all_attempted")$variant_n_distribution_status == "unavailable_in_locus_metric_csv")
paired <- view$contrasts[view$contrasts$stage == "challenge", ]
check(paired$regional_n_loci_a[paired$metric == "ap"] == 2)
check(paired$regional_n_loci_a[paired$metric == "recall"] == 4)
check(paired$regional_n_missing_causal_a[paired$metric == "ap"] == 0)
check(paired$regional_n_missing_causal_a[paired$metric == "recall"] == 1)
check(!get("corrected_refit", "sbayesrc", "brier_any_nonzero_total_effect_successful")$all_positive_target)
check(get("challenge", "sbayesrc", "brier_any_nonzero_total_effect_successful")$n_all_positive_loci == 2)
check(is.na(get("challenge", "sbayesrc", "brier_any_nonzero_total_effect_successful")$probability_target_aligned))
check(view$rows$regional_brier_probability_target_aligned[view$rows$stage == "challenge" & view$rows$method == "abf"])
check(!view$rows$regional_brier_probability_target_aligned[view$rows$stage == "challenge" & view$rows$method == "sbayesrc"])

# Counts, not sparse_inf/p_causal labels, identify the all-positive exception.
one <- .paper_regional_normalize(tasks[tasks$dataset_id == "gaussian_fixture_4" & tasks$method == "sbayesrc", ])
dense <- .paper_regional_describe(one, "any_nonzero_total_effect")
check(dense$all_positive_target && dense$constant_one_ap == 1 && dense$constant_one_brier == 0)
boundary <- .paper_regional_describe(.paper_regional_normalize(tasks[tasks$family != "challenge" & tasks$method == "sbayesrc", ]), "any_nonzero_total_effect")
check(!boundary$all_positive_target && is.na(boundary$constant_one_ap))
variable <- tasks[tasks$family == "challenge" & tasks$method == "abf", ]
variable$scenario <- "meta_known_overlap"; variable$stratum <- paste("challenge", variable$scenario, variable$p, variable$n, sep = ":")
variable <- .paper_regional_describe(.paper_regional_normalize(variable))
check(variable$n_known_variable_n_track_loci == 5L && variable$synthetic_n_gwas_min == 50L)
check(variable$variant_n_distribution_status == "unavailable_in_locus_metric_csv" && grepl("public/inputs.rds", variable$n_evidence_route, fixed = TRUE))
check(grepl("no individual rows", .paper_regional_n_caption(view$rows[view$rows$stage == "challenge", ]), fixed = TRUE))
check(grepl("simulated GWAS individuals", .paper_regional_n_caption(view$rows[view$rows$stage == "corrected_refit", ]), fixed = TRUE))
unknown <- tasks
unknown$probability_target_aligned[unknown$method == "abf" & unknown$dataset_id == "gaussian_fixture_1"] <- NA
write_tasks(unknown)
unknown_view <- .paper_regional_publication(path, schema_path, summary, contrasts)
check(is.na(unknown_view$rows$regional_brier_probability_target_aligned[unknown_view$rows$stage == "challenge" & unknown_view$rows$method == "abf"]))
check(unknown_view$rows$regional_brier_n_alignment_unknown[unknown_view$rows$stage == "challenge" & unknown_view$rows$method == "abf"] == 1L)
write_tasks(tasks)
reject(.paper_regional_publication(path, schema_path, summary, max_groups = 2L), "group memory bound")
reject(.paper_regional_publication(path, schema_path, summary, max_stratum_rows = 5L), "stratum memory bound")
reject(.paper_regional_publication(path, schema_path, rbind(summary, summary[1, ])), "Repeated regional summary")
wrong <- summary; wrong$n_attempted[1] <- wrong$n_attempted[1] + 1
reject(.paper_regional_publication(path, schema_path, wrong), "summary denominators")
wrong <- contrasts; wrong$n_pairs[1] <- wrong$n_pairs[1] + 1
reject(.paper_regional_publication(path, schema_path, summary, wrong), "contributing-pair denominator")
write_tasks(rbind(tasks, tasks[1, ]))
reject(.paper_regional_publication(path, schema_path, summary, chunk_size = 3L), "Duplicate regional locus")
wrong <- tasks; wrong$actual_p[wrong$dataset_id == "gaussian_fixture_1" & wrong$method == "sbayesrc"] <- 3L
write_tasks(wrong)
reject(.paper_regional_publication(path, schema_path, summary, contrasts), "paired candidate or truth")
wrong <- tasks; wrong$p[1] <- 6L
write_tasks(wrong)
reject(.paper_regional_publication(path, schema_path, summary), "stratum disagrees")
wrong <- tasks; wrong$n_missing_causal[1] <- 1L; write_tasks(wrong)
reject(.paper_regional_publication(path, schema_path, summary), "truth or candidate universe")
wrong <- tasks; wrong$score_target[1] <- "unknown_event"; write_tasks(wrong)
reject(.paper_regional_publication(path, schema_path, summary), "primary probability target")
wrong <- tasks; wrong$probability_target_aligned <- as.integer(wrong$probability_target_aligned); write_tasks(wrong)
reject(.paper_regional_publication(path, schema_path, summary), "probability alignment")
write_tasks(tasks)
schema <- .paper_regional_schema(schema_path); schema$columns[1] <- "foreign_id"; saveRDS(schema, schema_path)
reject(.paper_regional_publication(path, schema_path, summary), "schema column identities")

# Lots of unselected rows do not increase the retained stratum/group bound.
foreign <- tasks[rep(1, 100003L), ]; foreign$method <- "unselected"
write_tasks(rbind(foreign, tasks))
large <- .paper_regional_publication(path, schema_path, summary, contrasts, chunk_size = 10000L)
check(identical(view$census, large$census) && large$streaming$retained_rows == nrow(tasks))
check(large$streaming$rows_scanned == 100003L + nrow(tasks) && large$streaming$peak_chunk_rows <= 10000L)
write_tasks(tasks)

# One bounded typed partition may span all selected stages, but the census
# and existing paired masks remain stage-specific. No excluded stage appears.
three_tasks <- tasks[tasks$family != "challenge", ]
three_tasks <- do.call(rbind, lapply(c("original_refit", "corrected_refit", "corrected_resimulation"), function(stage) {
  z <- three_tasks; z$stage <- stage; z
}))
three_summary <- do.call(rbind, lapply(unique(three_tasks$stage), function(stage) {
  z <- summary[summary$stage == "corrected_refit", ]; z$stage <- stage; z
}))
three_pairs <- do.call(rbind, lapply(unique(three_tasks$stage), function(stage) {
  z <- contrasts[contrasts$stage == "corrected_refit", ]; z$stage <- stage; z
}))
write_tasks(three_tasks)
calls <- list()
three <- .paper_regional_publication(path, schema_path, three_summary, three_pairs, chunk_size = 2L,
  on_stratum = function(d, stratum) calls[[stratum]] <<- d)
check(length(calls) == 1L && nrow(calls[[1L]]) == 18L && three$streaming$peak_stratum_rows == 18L)
check(identical(sort(unique(calls[[1L]]$stage)), sort(unique(three_tasks$stage))))
check(all(c(".actual", ".nominal", ".nonzero", "panel_cluster_id") %in% names(calls[[1L]])))
check(identical(three_summary, three$rows[names(three_summary)]))
check(identical(three_pairs, three$contrasts[names(three_pairs)]))
reject(.paper_regional_publication(path, schema_path, three_summary, max_stratum_rows = 17L), "stratum memory bound")
calls <- list()
selected <- .paper_regional_publication(path, schema_path, three_summary[three_summary$stage == "corrected_refit", ],
  on_stratum = function(d, stratum) calls[[stratum]] <<- d)
check(length(calls) == 1L && nrow(calls[[1L]]) == 6L && all(calls[[1L]]$stage == "corrected_refit"))
check(selected$streaming$retained_rows == 6L)
reject(.paper_regional_publication(path, schema_path, summary, on_stratum = 1), "typed-stratum callback")
write_tasks(tasks)
single <- .paper_regional_publication(path, schema_path, summary[summary$method == "abf", ], contrasts[FALSE, ])
check(!nrow(single$contrasts) && single$streaming$retained_rows == 8L)
check(!nrow(.paper_regional_unassessed(summary, contrasts[FALSE, ])$contrasts))

# Sealed actual plotting CLI gate, export manifests and preservation tests.
campaign <- file.path(workspace, "campaign"); dir.create(campaign)
analysis <- file.path(campaign, "analysis"); dir.create(analysis)
plan <- list(profile = "integration", source = "Synthetic software fixture; no method-performance evidence",
  runtime_fingerprint = list(test_fixture = TRUE))
seal <- function() {
  plain <- plan; plain$plan_sha256 <- NULL; plan$plan_sha256 <<- .campaign_object_sha(plain)
  saveRDS(plan, file.path(campaign, "plan.rds"))
  .campaign_analysis_manifest(analysis, plan, "complete")
}
write.csv(summary, file.path(analysis, "regime-summary.csv"), row.names = FALSE)
write.csv(contrasts, file.path(analysis, "paired-comparisons.csv"), row.names = FALSE)
stopifnot(file.copy(path, file.path(analysis, "task-metrics.csv")))
stopifnot(file.copy(schema_path, file.path(analysis, "task-metrics-schema.rds")))
seal()
export <- function(name, extra = "--allow-incomplete") {
  log <- file.path(workspace, paste0(name, ".log"))
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "paper/plot_results.R",
    shQuote(analysis), shQuote(file.path(workspace, name)), extra), stdout = log, stderr = log)
  list(status = status, log = paste(readLines(log, warn = FALSE), collapse = "\n"))
}
before <- .campaign_artifact_manifest(c(list.files(analysis, full.names = TRUE), file.path(campaign, "plan.rds")))
result <- export("figures")
check(result$status == 0)
check(export("single-method", extra = c("--allow-incomplete", "--methods", "abf"))$status == 0)
manifest <- jsonlite::fromJSON(file.path(workspace, "figures/figure-manifest.json"))
check(!manifest$scientific_full_profile && manifest$regional_truth_census$available)
check(grepl("SYNTHETIC SOFTWARE FIXTURE", manifest$label, fixed = TRUE))
check(all(normalizePath(file.path(analysis, c("task-metrics.csv", "task-metrics-schema.rds"))) %in% manifest$inputs$path))
check(file.exists(file.path(workspace, "figures/regional-target-census.csv")))
after <- .campaign_artifact_manifest(before$path)
check(identical(before, after))
plan$profile <- "full"; seal()
result <- export("explicit-fixture-full", extra = character())
check(result$status != 0 && grepl("synthetic fixtures", result$log, fixed = TRUE))
plan$runtime_fingerprint$test_fixture <- FALSE; unlink(file.path(analysis, "task-metrics-schema.rds")); seal()
result <- export("missing-full", extra = character())
check(result$status != 0 && grepl("require sealed locus metrics", result$log, fixed = TRUE))
plan$profile <- "integration"; plan$runtime_fingerprint$test_fixture <- TRUE; seal()
result <- export("missing-development")
check(result$status == 0)
manifest <- jsonlite::fromJSON(file.path(workspace, "missing-development/figure-manifest.json"))
check(!manifest$regional_truth_census$available && grepl("unassessed", manifest$label, fixed = TRUE))

# The CLI separately calculates specified stage differences only when all
# required stages survive the export selection; stored method contrasts stay.
write_tasks(three_tasks)
write.csv(three_summary, file.path(analysis, "regime-summary.csv"), row.names = FALSE)
write.csv(three_pairs, file.path(analysis, "paired-comparisons.csv"), row.names = FALSE)
stopifnot(file.copy(path, file.path(analysis, "task-metrics.csv"), overwrite = TRUE))
stopifnot(file.copy(schema_path, file.path(analysis, "task-metrics-schema.rds"), overwrite = TRUE))
seal()
check(export("three-stages")$status == 0)
stage_manifest <- jsonlite::fromJSON(file.path(workspace, "three-stages/figure-manifest.json"))
check(stage_manifest$regional_stage_analysis$n_contrast_metric_rows == 24L)
check(!stage_manifest$regional_stage_analysis$registered_full_plan && !stage_manifest$scientific_full_profile)
stage_table <- read.csv(file.path(workspace, "three-stages/regional-stage-attribution.csv"))
check(nrow(stage_table) == 24L && all(stage_table$n_attempted_pairs == 3L))
check(all(grepl("development_or_incomplete", stage_table$analysis_status, fixed = TRUE)))
stored <- read.csv(file.path(workspace, "three-stages/regional-paired-publication.csv"))
check(identical(stored$lower, three_pairs$lower) && identical(stored$upper, three_pairs$upper))
check(export("one-stage", extra = c("--allow-incomplete", "--setting-regex", "corrected_refit"))$status == 0)
one_manifest <- jsonlite::fromJSON(file.path(workspace, "one-stage/figure-manifest.json"))
check(one_manifest$regional_stage_analysis$n_contrast_metric_rows == 0L)
one_selection <- read.csv(file.path(workspace, "one-stage/regional-stage-selection.csv"))
check(all(one_selection$status == "all_three_stages_not_selected; unassessed"))
check(!file.exists(file.path(workspace, "one-stage/regional-stage-attribution.csv")))
cat("Passed", checks, "regional publication checks; synthetic software fixtures only, no outcomes recomputed.\n")
