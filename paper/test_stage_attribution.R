#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))
source(file.path(repo, "paper", "stage_attribution.R"))
library(testthat)
fixture <- function(method = "abf") {
  d <- expand.grid(stage = c("original_refit", "corrected_refit", "corrected_resimulation"),
    dataset_id = c("scenario1", "scenario2"), region_id = 1:2, stringsAsFactors = FALSE)
  d$stratum <- "legacy_reconstruction:200:NA:binary:10:2:0.1:NA:50"
  d$family <- "legacy_reconstruction"; d$method <- method
  d$panel_cluster_id <- "sealed_original_simulation_hash"
  d$cluster_id <- paste0(d$panel_cluster_id, ":", d$region_id)
  d$success <- TRUE; d$actual_p <- ifelse(d$region_id == 1L, 48L, 50L)
  d$nominal_p <- 50L; d$n_gwas <- 200L; d$n_ref <- NA_integer_
  d$n_causal <- d$n_observed_causal <- 2L
  d$n_any_nonzero_total_effect <- 2L
  d$ap <- c(original_refit = .2, corrected_refit = .6, corrected_resimulation = .5)[d$stage]
  d$brier <- c(original_refit = .1, corrected_refit = .08, corrected_resimulation = .12)[d$stage]
  d$recall <- c(original_refit = .2, corrected_refit = .4, corrected_resimulation = .3)[d$stage]
  d$cs_distinct_recall <- d$recall
  d$probability_target_aligned <- TRUE
  d
}
pick <- function(x, metric, contrast = "implementation_and_acceptance")
  x[x$metric == metric & x$contrast == contrast, , drop = FALSE]

test_that("implementation and generating-experiment changes have separate directions and events", {
  z <- .paper_stage_contrasts(fixture(), n_boot = 100L)
  expect_equal(nrow(z), 12L)
  expect_equal(pick(z, "ap")$difference, .4)
  expect_equal(pick(z, "brier")$difference, -.02)
  expect_equal(pick(z, "recall")$difference, .2)
  expect_equal(pick(z, "ap", "generating_experiment")$difference, -.1)
  expect_equal(pick(z, "brier", "generating_experiment")$difference, .04)
  expect_true(all(z$n_attempted_pairs == 4L & z$n_pairs == 4L & z$n_clusters == 2L))
  expect_equal(pick(z, "ap")$lower, .4)
  expect_equal(pick(z, "ap")$upper, .4)
  expect_true(all(z$difference_direction == "stage_a_minus_stage_b"))
  expect_equal(pick(z, "ap")$fixed, "observed_inputs_and_truth")
  expect_equal(pick(z, "ap", "generating_experiment")$fixed, "genotypes_and_design")
  expect_false(any(z$interval_nominal_coverage_validated))
  expect_true(all(grepl("descriptive_percentile_envelope", z$interval_coverage_status)))
  expect_true(all(grepl("no_statistical_distinguishability_claim", z$few_cluster_warning)))
})

test_that("the attempted-recovery and complete-fit masks remain distinct", {
  d <- fixture(); failed <- d$stage == "original_refit" & d$dataset_id == "scenario1" & d$region_id == 1L
  d$success[failed] <- FALSE; d$ap[failed] <- d$brier[failed] <- NA_real_
  d$recall[failed] <- d$cs_distinct_recall[failed] <- 0
  z <- .paper_stage_contrasts(d, n_boot = 100L)
  expect_equal(pick(z, "ap")$n_pairs, 3L)
  expect_equal(pick(z, "brier")$n_pairs, 3L)
  expect_equal(pick(z, "recall")$n_pairs, 4L)
  expect_equal(pick(z, "ap_observed_positive_attempts_zero_for_failure")$n_pairs, 4L)
  expect_equal(pick(z, "ap_observed_positive_attempts_zero_for_failure")$difference, .45)
  expect_equal(pick(z, "availability")$difference, .25)
  expect_equal(pick(z, "ap")$n_recovered_stage_a, 1L)
  expect_equal(pick(z, "ap")$n_success_both, 3L)
  expect_equal(pick(z, "ap")$n_failed_both, 0L)
  null <- d$dataset_id == "scenario2" & d$region_id == 2L
  d$n_causal[null] <- d$n_observed_causal[null] <- 0L
  d$ap[null] <- d$recall[null] <- d$cs_distinct_recall[null] <- NA_real_
  z <- .paper_stage_contrasts(d, n_boot = 100L)
  expect_equal(pick(z, "ap")$n_pairs, 2L)
  expect_equal(pick(z, "brier")$n_pairs, 3L)
  expect_equal(pick(z, "ap_observed_positive_attempts_zero_for_failure")$n_pairs, 3L)
  expect_equal(pick(z, "recall")$n_pairs, 3L)
  expect_equal(pick(z, "availability")$n_pairs, 4L)
})

test_that("shared-prior methods require whole-panel independent clusters", {
  z <- .paper_stage_contrasts(fixture("fb_xregion"), n_boot = 100L)
  expect_true(all(z$n_clusters == 1L & z$cluster_field == "panel_cluster_id"))
  expect_true(all(is.na(z$lower) & is.na(z$upper)))
  expect_true(all(grepl("point_difference_only", z$inference, fixed = TRUE)))
  expect_equal(pick(z, "ap")$difference, .4)
  expect_true(all(z$interval_coverage_status == "no_interval"))
})

test_that("two-cluster bootstrap envelopes do not establish nominal coverage", {
  # Two independent mean-zero cluster effects +/-0.25 give equal signs with
  # probability 1/2. Their bootstrap envelopes exclude zero on those samples.
  # The point/interval arithmetic is retained, with no 95%-coverage claim.
  signs <- expand.grid(a = c(-1, 1), b = c(-1, 1))
  excludes_zero <- logical(nrow(signs))
  for (i in seq_len(nrow(signs))) {
    d <- fixture(); d$ap <- .5
    changed <- d$stage == "corrected_refit"
    d$ap[changed] <- .5 + .25 * c(signs$a[i], signs$b[i])[d$region_id[changed]]
    z <- pick(.paper_stage_contrasts(d, n_boot = 1000L), "ap")
    excludes_zero[i] <- z$lower > 0 || z$upper < 0
    expect_false(z$interval_nominal_coverage_validated)
    expect_match(z$few_cluster_warning, "no_statistical_distinguishability_claim", fixed = TRUE)
  }
  expect_equal(mean(excludes_zero), .5)
})

test_that("event alignment follows exactly the contributing probability-score pairs", {
  d <- fixture(); failed <- d$stage == "original_refit" & d$dataset_id == "scenario1" & d$region_id == 1L
  d$success[failed] <- FALSE; d$ap[failed] <- d$brier[failed] <- NA_real_
  d$recall[failed] <- d$cs_distinct_recall[failed] <- 0
  d$probability_target_aligned[failed] <- FALSE
  z <- .paper_stage_contrasts(d, n_boot = 100L)
  expect_true(pick(z, "brier")$probability_targets_aligned)
  d$probability_target_aligned[d$stage == "corrected_refit" & d$region_id == 2L] <- NA
  z <- .paper_stage_contrasts(d, n_boot = 100L)
  expect_false(pick(z, "brier")$probability_targets_aligned)
  expect_match(pick(z, "brier")$interpretation, "operational Brier diagnostic only", fixed = TRUE)
  expect_true(is.na(pick(z, "ap")$probability_targets_aligned))
})

test_that("fixed input/truth identities fail closed, while resimulation can change truth", {
  for (field in c("panel_cluster_id", "cluster_id", "actual_p", "nominal_p", "n_gwas", "n_ref",
                  "n_causal", "n_observed_causal", "n_any_nonzero_total_effect")) {
    d <- fixture(); i <- which(d$stage == "corrected_refit")[1L]
    if (is.character(d[[field]])) d[[field]][i] <- paste0(d[[field]][i], "_changed") else
      if (field == "n_ref") d[[field]][i] <- 100L else d[[field]][i] <- d[[field]][i] + 1L
    expect_error(.paper_stage_contrasts(d, n_boot = 100L), "fixed|counts", info = field)
  }
  d <- fixture(); changed <- d$stage == "corrected_resimulation"
  d$n_causal[changed] <- d$n_observed_causal[changed] <- 3L
  d$n_any_nonzero_total_effect[changed] <- 50L
  expect_equal(pick(.paper_stage_contrasts(d, n_boot = 100L), "ap", "generating_experiment")$n_pairs, 4L)
})

test_that("generating-experiment AP eligibility is paired separately from availability", {
  d <- fixture(); null <- d$stage == "corrected_resimulation" & d$region_id == 2L
  d$n_causal[null] <- d$n_observed_causal[null] <- d$n_any_nonzero_total_effect[null] <- 0L
  d$ap[null] <- d$recall[null] <- d$cs_distinct_recall[null] <- NA_real_
  z <- .paper_stage_contrasts(d, n_boot = 100L)
  expect_equal(pick(z, "ap", "generating_experiment")$n_pairs, 2L)
  expect_equal(pick(z, "ap_observed_positive_attempts_zero_for_failure", "generating_experiment")$n_pairs, 2L)
  expect_equal(pick(z, "availability", "generating_experiment")$n_pairs, 4L)
  expect_equal(pick(z, "brier", "generating_experiment")$n_pairs, 4L)
  expect_equal(pick(z, "recall", "generating_experiment")$n_pairs, 2L)
})

test_that("cluster resampling preserves mean-locus weighting and phenotype reuse", {
  d <- fixture(); d <- d[d$dataset_id == "scenario1", ]
  one <- d[d$region_id == 1L, ]; two <- d[d$region_id == 2L, ]
  d <- rbind(one, transform(one, dataset_id = "scenario2"),
             transform(one, dataset_id = "scenario3"), two)
  d$ap <- .5; changed <- d$stage == "corrected_refit"
  d$ap[changed] <- ifelse(d$region_id[changed] == 1L, .9, .3)
  z <- pick(.paper_stage_contrasts(d, n_boot = 1000L), "ap")
  expect_equal(z$difference, .25) # (3*0.4 - 0.2)/4, not equal-cluster 0.1.
  expect_equal(z$lower, -.2); expect_equal(z$upper, .4)
  expect_equal(z$n_attempted_pairs, 4L); expect_equal(z$n_clusters, 2L)
  repeated <- do.call(rbind, lapply(seq_len(10L), function(i) {
    x <- fixture(); x <- x[x$dataset_id == "scenario1", ]; x$dataset_id <- paste0("phenotype", i); x
  }))
  z <- .paper_stage_contrasts(repeated, n_boot = 100L)
  expect_true(all(z$n_attempted_pairs == 20L & z$n_clusters == 2L))
  expect_true(all(grepl("at_most_two", z$few_cluster_warning)))
  repeated$method <- "fb_pooled"
  z <- .paper_stage_contrasts(repeated, n_boot = 100L)
  expect_true(all(z$n_attempted_pairs == 20L & z$n_clusters == 1L))
  expect_true(all(is.na(z$lower) & is.na(z$upper)))
})

test_that("missing stages are explicit unmatched attempts and cannot become zeros", {
  d <- fixture(); d <- d[d$stage != "original_refit", ]
  z <- .paper_stage_contrasts(d, n_boot = 100L)
  x <- pick(z, "ap")
  expect_equal(x$n_stage_a, 4L); expect_equal(x$n_stage_b, 0L)
  expect_equal(x$n_attempted_pairs, 0L); expect_equal(x$n_unmatched_stage_a, 4L)
  expect_true(is.na(x$difference)); expect_equal(x$n_pairs, 0L)
  expect_equal(pick(z, "ap", "generating_experiment")$n_pairs, 4L)
})

test_that("input order and execution do not perturb the caller's random stream", {
  d <- fixture(); set.seed(7L); before <- .Random.seed
  z <- .paper_stage_contrasts(d, n_boot = 100L)
  expect_identical(.Random.seed, before)
  expect_identical(.paper_stage_contrasts(d[rev(seq_len(nrow(d))), ], n_boot = 100L), z)
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(.paper_stage_contrasts(d, n_boot = 100L))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("malformed or oversized strata cannot silently enter an attribution table", {
  d <- fixture()
  expect_error(.paper_stage_contrasts(rbind(d, d[1L, ]), n_boot = 100L), "Duplicate")
  expect_error(.paper_stage_contrasts(d, max_stratum_rows = 2L), "memory bound")
  expect_error(.paper_stage_contrasts(d, n_boot = 20.5), "specification")
  expect_error(.paper_stage_contrasts(d, seed = NA_real_), "specification")
  x <- d; x$stratum[1L] <- "other"; expect_error(.paper_stage_contrasts(x), "one regional")
  x <- d; x$family <- "challenge"; expect_error(.paper_stage_contrasts(x), "Gaussian")
  x <- d; x$stage[1L] <- "archived_rescore"; expect_error(.paper_stage_contrasts(x), "stage identity")
  x <- d; x$success[1L] <- NA; expect_error(.paper_stage_contrasts(x), "logical")
  x <- d; x$brier[1L] <- Inf; expect_error(.paper_stage_contrasts(x), "range")
  x <- d; x$ap[1L] <- NA; expect_error(.paper_stage_contrasts(x), "assessable AP")
  x <- d; x$success[1L] <- FALSE; expect_error(.paper_stage_contrasts(x), "successful-fit")
  x$ap[1L] <- x$brier[1L] <- NA_real_; expect_error(.paper_stage_contrasts(x), "zero attempted recovery")
  x <- d; x$n_ref[1L] <- 0L; expect_error(.paper_stage_contrasts(x), "reference")
  x <- d; x$cluster_id[1L] <- ""; expect_error(.paper_stage_contrasts(x), "identity")
  expect_equal(.paper_stage_contrasts(d[FALSE, ]), data.frame())
})
