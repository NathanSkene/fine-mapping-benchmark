#!/usr/bin/env Rscript
# Analytical and publication-export fixtures only; no method-performance evidence.
source("paper/read_selected_csv.R")
source("paper/native_publication.R")
checks <- 0L
check <- function(ok) { stopifnot(ok); checks <<- checks + 1L }
reject <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  check(inherits(error, "error") && grepl(pattern, conditionMessage(error), fixed = TRUE))
}
workspace <- Sys.getenv("FMB_PAPER_FIXTURE_DIR", unset = tempfile("paper-native-publication-"))
if (dir.exists(workspace)) stop("Use a new publication-test fixture directory")
dir.create(workspace)
if (dir.exists(".local-lib")) .libPaths(c(".local-lib", .libPaths()))
for (path in list.files("R", full.names = TRUE)) if (endsWith(path, ".R")) source(path)
source("scripts/native/campaign/helpers.R")

# The exact source scorer proves the no-data all-positive baseline and AP's
# invariance to rank. These are identities, not simulated benchmark findings.
for (p in list(rep(1, 4), rep(.5, 4), c(.9, .6, .2, .1), c(.1, .2, .6, .9)))
  check(identical(.native_event_metrics(p, rep(TRUE, 4))$ap, 1))
baseline <- .native_event_metrics(rep(1, 4), rep(TRUE, 4))
check(baseline$brier == 0 && baseline$log_loss == 0 && baseline$recall == 1)

x <- expand.grid(model = c("sparse", "sparse_inf"), draw = 1:2,
  method = c("abf", "sbayesrc_native"), event = c("sparse_component_membership", "any_nonzero_effect"),
  stringsAsFactors = FALSE)
x$draw_id <- paste(x$model, x$draw, sep = "_draw")
x$pi <- .001; x$h2 <- .3; x$exposure <- "none"
x$pip_event <- ifelse(x$method == "abf", "modeled_sparse_effect_membership", "any_nonzero_mixture_component")
pieces <- lapply(seq_len(nrow(x)), function(i) {
  d <- x[rep(i, 4), ]; d$input_id <- paste(d$draw_id, c(1, 1, 2, 2), sep = "_phenotype")
  d$region_id <- rep(1:2, 2)
  truth <- if (d$model[1] == "sparse_inf" && d$event[1] == "any_nonzero_effect") rep(4L, 4) else c(1L, 0L, 2L, 0L)
  probability <- if (d$method[1] == "abf") .5 else .8
  scores <- do.call(rbind, lapply(truth, function(n) .native_event_metrics(rep(probability, 4), seq_len(4) <= n)))
  d$success <- d$complete <- TRUE
  d$probability_target_aligned <- if (d$method[1] == "abf") d$event == "sparse_component_membership" else d$event == "any_nonzero_effect"
  cbind(d, scores)
})
tasks <- do.call(rbind, pieces)
key <- .paper_native_key(tasks)
summary <- do.call(rbind, lapply(split(seq_len(nrow(tasks)), key), function(index) {
  z <- tasks[index, ]; d <- z[1, c(.paper_native_keys, "draw")]
  d$n_planned <- d$n_complete <- d$n_successful <- nrow(z)
  d$ap_successful <- mean(z$ap, na.rm = TRUE); d$brier_successful <- mean(z$brier)
  d$recall_all_completed <- mean(z$recall, na.rm = TRUE)
  d$all_probability_targets_aligned <- all(z$probability_target_aligned)
  d$all_native_accuracy_eligible <- if (d$method == "sbayesrc_native") TRUE else NA
  d$posterior_conditioning <- if (d$method == "abf") "exactly_one_modeled_causal_variant" else "joint_whole_panel_mixture_model"
  d$native_null_model_in_support <- if (d$method == "abf") FALSE else NA
  d
}))
rownames(summary) <- NULL
path <- file.path(workspace, "tasks.csv"); write.csv(tasks, path, row.names = FALSE)
census <- .paper_native_truth_census(path, summary, chunk_size = 3L)
check(identical(summary, census[names(summary)]))
all_positive <- census$model == "sparse_inf" & census$event == "any_nonzero_effect"
check(all(census$target_support[all_positive] == "all_positive"))
check(all(census$census_n_negative[all_positive] == 0 & census$census_n_positive[all_positive] == 16))
check(all(census$ap_localisation_status[all_positive] == "no_localisation_discrimination"))
check(all(census$recall_localisation_status[all_positive] == "no_localisation_discrimination"))
check(all(census$constant_one_ap[all_positive] == 1 & census$constant_one_brier[all_positive] == 0))
check(all(census$ap_localisation_status[!all_positive] == "within_locus_discrimination_available"))
check(all(census$census_n_all_negative_loci[!all_positive] == 2 & census$census_n_mixed_loci[!all_positive] == 2))
check(identical(.paper_native_degenerate(census, "ap_successful"), all_positive))
check(identical(.paper_native_degenerate(census, "success_rate"), rep(FALSE, nrow(census))))
reject(.paper_native_truth_census(path, summary, max_groups = 2L), "group memory bound")
reject(.paper_native_truth_census(path, rbind(summary, summary[1, ])), "Repeated native")
wrong <- summary; wrong$n_planned[1] <- 5L
reject(.paper_native_truth_census(path, wrong), "denominators")
wrong <- summary; wrong$all_probability_targets_aligned[1] <- !wrong$all_probability_targets_aligned[1]
reject(.paper_native_truth_census(path, wrong), "alignment")
bad <- tasks; bad$n_positive[1] <- 5
bad_path <- file.path(workspace, "bad.csv"); write.csv(bad, bad_path, row.names = FALSE)
reject(.paper_native_truth_census(bad_path, summary), "truth counts")
bad$n_positive[1] <- 1.5; write.csv(bad, bad_path, row.names = FALSE)
reject(.paper_native_truth_census(bad_path, summary), "truth counts")
bad <- tasks; bad$success[1] <- NA; write.csv(bad, bad_path, row.names = FALSE)
reject(.paper_native_truth_census(bad_path, summary), "census flags")

# Both pooled classes can exist without any within-locus ranking problem.
one <- summary[1, ]; alternate <- tasks[key == .paper_native_key(one), ]
alternate$n_positive <- c(0L, 4L, 0L, 4L)
alternate$probability_target_aligned <- one$all_probability_targets_aligned
write.csv(alternate, bad_path, row.names = FALSE)
whole_loci <- .paper_native_truth_census(bad_path, one, chunk_size = 1L)
check(whole_loci$target_support == "positive_and_negative_loci_only")
check(whole_loci$ap_localisation_status == "no_localisation_discrimination" && !whole_loci$all_positive_target)
alternate$n_positive <- 0L; write.csv(alternate, bad_path, row.names = FALSE)
null <- .paper_native_truth_census(bad_path, one)
check(null$target_support == "all_negative" && null$ap_localisation_status == "no_positive_scored_loci")
alternate$n_positive <- c(4L, 1L, 4L, 0L); alternate$success <- c(TRUE, FALSE, TRUE, TRUE)
one$n_successful <- 3L; write.csv(alternate, bad_path, row.names = FALSE)
conditional <- .paper_native_truth_census(bad_path, one, chunk_size = 2L)
check(conditional$ap_localisation_status == "no_localisation_discrimination" &&
      conditional$recall_localisation_status == "mix_of_degenerate_and_discriminating_loci")

# Large unselected input is scanned, while only the selected bounded census is retained.
foreign <- tasks[rep(1L, 100003L), ]; foreign$draw_id <- "unselected_draw"
large_path <- file.path(workspace, "large.csv")
write.csv(rbind(foreign, tasks), large_path, row.names = FALSE)
streamed <- .paper_native_truth_census(large_path, summary, chunk_size = 10000L)
check(identical(census[names(streamed)], streamed))
unlink(large_path)

paired <- expand.grid(model = c("sparse", "sparse_inf"), event = c("sparse_component_membership", "any_nonzero_effect"),
  metric = c("ap_successful", "brier_successful", "recall_all_completed"), stringsAsFactors = FALSE)
paired$pi <- .001; paired$h2 <- .3; paired$method_a <- "sbayesrc_native"; paired$method_b <- "abf"
paired$exposure_a <- paired$exposure_b <- "none"
paired$difference <- vapply(seq_len(nrow(paired)), function(i) {
  z <- paired[i, ]; scope <- summary$model == z$model & summary$event == z$event
  mean(summary[[z$metric]][scope & summary$method == "sbayesrc_native"]) -
    mean(summary[[z$metric]][scope & summary$method == "abf"])
}, 1.0)
# Both artificial draws are identical, so every resampled mean is identical.
paired$lower <- paired$upper <- paired$difference
paired$n_paired_genotype_draws <- 2L
census$aligned <- census$all_probability_targets_aligned
view <- .paper_native_paired_view(paired, census)
check(all(vapply(names(paired), function(name) identical(paired[[name]], view[[name]]), TRUE)))
check(all(view$pip_event_a == "any_nonzero_mixture_component" & view$pip_event_b == "modeled_sparse_effect_membership"))
check(all(!view$probability_events_both_aligned))
check(all(view$probability_target_aligned_a == (view$event == "any_nonzero_effect")))
check(all(view$probability_target_aligned_b == (view$event == "sparse_component_membership")))
check(all(grepl("event_misaligned", view$probability_interpretation, fixed = TRUE)))
degenerate <- view$model == "sparse_inf" & view$event == "any_nonzero_effect"
check(all(view$localisation_interpretation[degenerate & view$metric != "brier_successful"] == "no_localisation_discrimination"))
check(all(view$census_n_positive_a[degenerate] == 32 & view$census_n_negative_b[degenerate] == 0))
check(all(grepl("Brier=0", view$all_positive_no_data_baseline[degenerate], fixed = TRUE)))
old <- .paper_native_unassessed(summary); old$aligned <- old$all_probability_targets_aligned
unknown <- .paper_native_paired_view(paired, old)
check(all(unknown$localisation_interpretation[unknown$metric != "brier_successful"] == "truth_census_unassessed"))
check(all(unknown$all_positive_no_data_baseline == "all_positive_status_unassessed"))

# Cross-draw support must use summed counts on precisely the source's usable
# pairs. Status concatenation and all-selected support can both hide degeneracy.
cross_tasks <- tasks[tasks$model == "sparse" & tasks$event == "sparse_component_membership", ]
cross_summary <- summary[summary$model == "sparse" & summary$event == "sparse_component_membership", ]
rebuild_cross <- function(raw) {
  for (i in seq_len(nrow(raw))) {
    probability <- if (raw$method[i] == "abf") .5 else .8
    score <- .native_event_metrics(rep(probability, 4), seq_len(4) <= raw$n_positive[i])
    if (!raw$success[i]) {
      score$ap <- score$brier <- score$log_loss <- NA_real_
      score$recall <- if (raw$n_positive[i] > 0) 0 else NA_real_
    }
    raw[i, names(score)] <- score
  }
  d <- cross_summary
  mean_or_na <- function(v) if (!length(v) || all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  for (i in seq_len(nrow(d))) {
    z <- raw[.paper_native_key(raw) == .paper_native_key(d[i, ]), ]
    d$n_successful[i] <- sum(z$success); d$n_complete[i] <- sum(z$complete)
    d$ap_successful[i] <- mean_or_na(z$ap[z$success])
    d$brier_successful[i] <- mean_or_na(z$brier[z$success])
    d$recall_all_completed[i] <- mean_or_na(z$recall[z$complete])
  }
  write.csv(raw, bad_path, row.names = FALSE)
  d <- .paper_native_truth_census(bad_path, d, chunk_size = 3L)
  d$aligned <- d$all_probability_targets_aligned
  d
}
cross_tasks$n_positive <- ifelse(cross_tasks$draw == 1, 0L, 4L)
cross <- rebuild_cross(cross_tasks)
pair <- paired[paired$model == "sparse" & paired$event == "sparse_component_membership" &
                 paired$metric %in% c("ap_successful", "recall_all_completed"), ]
pair$n_paired_genotype_draws <- 1L; pair$difference <- pair$lower <- pair$upper <- 0
cross_view <- .paper_native_paired_view(pair, cross)
check(all(cross_view$localisation_interpretation == "no_localisation_discrimination"))
check(all(cross_view$selected_scope_target_support_a == "positive_and_negative_loci_only"))
check(all(cross_view$n_publication_usable_draw_pairs == 1L & cross_view$census_n_mixed_loci_a == 0))
check(all(cross_view$all_positive_target_both_sides))
check(all(vapply(names(pair), function(name) identical(pair[[name]], cross_view[[name]]), TRUE)))

# A mixed second draw excluded for nonfinite AP cannot suppress the first
# draw's exact all-positive paired annotation. Completed failures remain in
# recall's declared source endpoint until incompleteness excludes the pair.
cross_tasks$n_positive <- ifelse(cross_tasks$draw == 1, 4L, 2L)
failed <- cross_tasks$draw == 2 & cross_tasks$method == "abf"
cross_tasks$success[failed] <- FALSE
cross <- rebuild_cross(cross_tasks)
ap_pair <- pair[pair$metric == "ap_successful", ]
excluded <- .paper_native_paired_view(ap_pair, cross)
check(excluded$localisation_interpretation == "no_localisation_discrimination")
check(excluded$census_n_mixed_loci_a == 0 && excluded$selected_scope_census_n_mixed_loci_a == 4)
check(excluded$ap_localisation_status_a == "no_localisation_discrimination" &&
      excluded$selected_scope_ap_localisation_status_a == "mix_of_degenerate_and_discriminating_loci")
cross_tasks$complete[failed] <- FALSE
cross <- rebuild_cross(cross_tasks)
excluded <- .paper_native_paired_view(pair, cross)
check(all(excluded$localisation_interpretation == "no_localisation_discrimination" &
            excluded$n_publication_usable_draw_pairs == 1))
wrong_pair <- pair; wrong_pair$n_paired_genotype_draws <- 2L
reject(.paper_native_paired_view(wrong_pair, cross), "stored paired-draw denominator")
reject(.paper_native_paired_view(pair, rbind(cross, cross[1, ])), "ambiguous draw identities")
no_pairs <- pair; no_pairs$n_paired_genotype_draws <- 0L
no_pairs$difference <- no_pairs$lower <- no_pairs$upper <- NA_real_
cross$n_complete <- 0L
empty_view <- .paper_native_paired_view(no_pairs, cross)
check(all(empty_view$localisation_interpretation == "no_usable_paired_draws" &
          empty_view$census_n_loci_a == 0 & is.na(empty_view$probability_target_aligned_a)))

# Actual sealed-analysis CLI export tests. Every successful fixture is labelled
# development-only. Full-profile fixtures below are used only to prove refusal.
campaign <- file.path(workspace, "campaign"); dir.create(campaign)
analysis <- file.path(campaign, "analysis"); dir.create(analysis)
plan <- list(profile = "integration", source = "Synthetic publication software fixture; no method-performance evidence",
  runtime_fingerprint = list(test_fixture = TRUE))
plan$plan_sha256 <- .campaign_object_sha(plan)
saveRDS(plan, file.path(campaign, "plan.rds"), version = 2)
write.csv(summary, file.path(analysis, "draw-metrics.csv"), row.names = FALSE)
write.csv(tasks, file.path(analysis, "task-metrics.csv"), row.names = FALSE)
write.csv(paired, file.path(analysis, "paired-draw-comparisons.csv"), row.names = FALSE)
.native_seal_analysis(analysis, plan, plan$runtime_fingerprint, "complete", list(test_fixture = TRUE), Sys.time())
healthy_inputs <- c(list.files(analysis, full.names = TRUE), file.path(campaign, "plan.rds"))
healthy_bytes <- lapply(healthy_inputs, function(path) readBin(path, "raw", n = file.info(path)$size))
export <- function(output, extra = character()) {
  log <- file.path(workspace, paste0(basename(output), ".log"))
  result <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "paper/plot_results.R", shQuote(analysis), shQuote(output), extra), stdout = log, stderr = log)
  list(status = result, log = paste(readLines(log, warn = FALSE), collapse = "\n"))
}
out <- file.path(workspace, "figures")
result <- export(out, "--allow-incomplete")
if (result$status != 0) stop(result$log)
check(result$status == 0 && all(file.exists(file.path(out, c("settings.pdf", "paired-contrasts.pdf",
  "native-draw-publication.csv", "native-target-census.csv", "native-paired-publication.csv")))))
out_rows <- read.csv(file.path(out, "native-draw-publication.csv"))
source_rows <- read.csv(file.path(analysis, "draw-metrics.csv"))
check(all(vapply(names(source_rows), function(name) identical(source_rows[[name]], out_rows[[name]]), TRUE)))
out_pairs <- read.csv(file.path(out, "native-paired-publication.csv"))
source_pairs <- read.csv(file.path(analysis, "paired-draw-comparisons.csv"))
check(all(vapply(names(source_pairs), function(name) identical(source_pairs[[name]], out_pairs[[name]]), TRUE)))
manifest <- jsonlite::fromJSON(file.path(out, "figure-manifest.json"))
check(!manifest$scientific_full_profile && grepl("DEVELOPMENT ONLY", manifest$label, fixed = TRUE))
check(grepl("SYNTHETIC SOFTWARE FIXTURE: no method-performance evidence", manifest$label, fixed = TRUE))
check(manifest$native_truth_census_available && nrow(manifest$plotting_helpers) == 2L)
check(normalizePath(file.path(analysis, "task-metrics.csv")) %in% manifest$inputs$path)
.campaign_validate_artifacts(manifest$outputs)
check(TRUE)
gated <- export(file.path(workspace, "gated"))
check(gated$status != 0 && grepl("complete full-profile", gated$log, fixed = TRUE))

# Scientific mode cannot silently proceed without the exact sealed truth census.
unlink(file.path(analysis, "task-metrics.csv"))
plan$profile <- "full"; plan$plan_sha256 <- NULL; plan$plan_sha256 <- .campaign_object_sha(plan)
saveRDS(plan, file.path(campaign, "plan.rds"), version = 2)
.native_seal_analysis(analysis, plan, plan$runtime_fingerprint, "complete", list(test_fixture = TRUE), Sys.time())
gated <- export(file.path(workspace, "fixture-full-profile"))
check(gated$status != 0 && grepl("explicit synthetic fixtures", gated$log, fixed = TRUE))
plan$runtime_fingerprint$test_fixture <- FALSE
plan$plan_sha256 <- NULL; plan$plan_sha256 <- .campaign_object_sha(plan)
saveRDS(plan, file.path(campaign, "plan.rds"), version = 2)
.native_seal_analysis(analysis, plan, plan$runtime_fingerprint, "complete", list(test_fixture = TRUE), Sys.time())
gated <- export(file.path(workspace, "missing-census"))
check(gated$status != 0 && grepl("sealed locus-level truth counts", gated$log, fixed = TRUE))

# A changed source table is refused before a figure or derivative is written.
write.csv(tasks, file.path(analysis, "task-metrics.csv"), row.names = FALSE)
plan$profile <- "integration"; plan$plan_sha256 <- NULL; plan$plan_sha256 <- .campaign_object_sha(plan)
saveRDS(plan, file.path(campaign, "plan.rds"), version = 2)
.native_seal_analysis(analysis, plan, plan$runtime_fingerprint, "complete", list(test_fixture = TRUE), Sys.time())
cat("\n", file = file.path(analysis, "task-metrics.csv"), append = TRUE)
gated <- export(file.path(workspace, "tampered"), "--allow-incomplete")
check(gated$status != 0 && !dir.exists(file.path(workspace, "tampered")))

# Restore the exact healthy input bytes, retaining all negative-case logs.
for (i in seq_along(healthy_inputs)) writeBin(healthy_bytes[[i]], healthy_inputs[i])
.campaign_validate_artifacts(manifest$inputs)
check(!dir.exists(file.path(campaign, "analysis.lock")))

cat("Passed", checks, "native publication checks. Diagnostic fixture output:", out, "\n")
