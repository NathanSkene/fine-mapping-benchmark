#!/usr/bin/env Rscript
# Software fixtures for selection/provenance gates, not performance evidence.
source("paper/stage_attribution.R")
source("paper/stage_publication.R")
library(testthat)
specification <- "paper/stage-attribution-spec.json"
spec <- jsonlite::fromJSON(specification)
fixture <- function() {
  d <- expand.grid(stage = c("original_refit", "corrected_refit", "corrected_resimulation"),
    region_id = 1:2, method = c("abf", "fb_pooled"), stringsAsFactors = FALSE)
  d$stratum <- "legacy_reconstruction:200:NA:binary:10:2:0.1:NA:50"
  d$family <- "legacy_reconstruction"; d$dataset_id <- "scenario1"
  d$panel_cluster_id <- "sealed_original_simulation_hash"
  d$cluster_id <- paste0(d$panel_cluster_id, ":", d$region_id)
  d$success <- d$probability_target_aligned <- TRUE
  d$actual_p <- 48L; d$nominal_p <- 50L; d$n_gwas <- 200L; d$n_ref <- NA_integer_
  d$n_causal <- d$n_observed_causal <- 2L
  d$ap <- c(original_refit = .2, corrected_refit = .6, corrected_resimulation = .5)[d$stage]
  d$brier <- .1; d$recall <- d$cs_distinct_recall <- .5
  d
}
selected <- function(d) unique(d[c("stage", "stratum", "method")])
plan <- list(profile = "integration", plan_sha256 = paste(rep("a", 64L), collapse = ""))

test_that("development adapters preserve stage scope and label new analytical differences", {
  d <- fixture(); rows <- selected(d)
  adapter <- .paper_stage_publication(plan, rows, specification, collection_complete = TRUE)
  expect_null(adapter$on_stratum(d, d$stratum[1L]))
  view <- adapter$finish()
  expect_equal(nrow(view$contrasts), 24L)
  expect_equal(nrow(view$selection), 2L)
  expect_false(view$metadata$registered_full_plan)
  expect_false(view$metadata$scientific_full_profile)
  expect_true(view$metadata$development_results_already_seen)
  expect_true(all(grepl("development_or_incomplete", view$contrasts$analysis_status)))
  expect_true(all(view$contrasts$n_attempted_pairs == 2L))
  expect_true(all(is.na(view$contrasts$lower[view$contrasts$method == "fb_pooled"])))
  expect_error(adapter$on_stratum(d, d$stratum[1L]), "Repeated")
})

test_that("excluded stages are not transformed into missing or failed scientific attempts", {
  d <- fixture(); d <- d[d$stage == "corrected_refit", ]
  adapter <- .paper_stage_publication(plan, selected(d), specification)
  adapter$on_stratum(d, d$stratum[1L]); view <- adapter$finish()
  expect_equal(nrow(view$contrasts), 0L)
  expect_true(all(view$selection$status == "all_three_stages_not_selected; unassessed"))
  expect_false(any(grepl("failed", view$selection$status)))
  d <- fixture(); rows <- selected(d); rows <- rows[rows$method != "fb_pooled" | rows$stage == "corrected_refit", ]
  d <- d[d$method != "fb_pooled" | d$stage == "corrected_refit", ]
  adapter <- .paper_stage_publication(plan, rows, specification)
  adapter$on_stratum(d, d$stratum[1L]); view <- adapter$finish()
  expect_equal(unique(view$contrasts$method), "abf")
  expect_equal(view$selection$status[view$selection$method == "fb_pooled"], "all_three_stages_not_selected; unassessed")
})

test_that("full-profile registered analysis is bound to the exact frozen plan", {
  d <- fixture(); rows <- selected(d)
  registered_plan <- list(profile = "full", plan_sha256 = spec$regional_full_plan_sha256)
  adapter <- .paper_stage_publication(registered_plan, rows, specification,
    scientific = TRUE, collection_complete = TRUE)
  adapter$on_stratum(d, d$stratum[1L]); view <- adapter$finish()
  expect_true(view$metadata$registered_full_plan && view$metadata$scientific_full_profile)
  expect_true(all(view$contrasts$analysis_status == "registered_full_profile_stage_analysis"))
  other_plan <- plan; other_plan$profile <- "full"
  adapter <- .paper_stage_publication(other_plan, rows, specification, scientific = TRUE)
  adapter$on_stratum(d, d$stratum[1L]); view <- adapter$finish()
  expect_equal(nrow(view$contrasts), 0L)
  expect_true(all(view$selection$status == "unregistered_full_plan_not_analysed"))
  expect_error(.paper_stage_publication(plan, rows, specification, scientific = TRUE), "full profile")
})

test_that("unmatched locus attempts are only permitted in explicitly incomplete evidence", {
  d <- fixture(); rows <- selected(d)
  d <- d[!(d$stage == "original_refit" & d$region_id == 1L), ]
  adapter <- .paper_stage_publication(plan, rows, specification, collection_complete = TRUE)
  expect_error(adapter$on_stratum(d, d$stratum[1L]), "unmatched attempted")
  adapter <- .paper_stage_publication(plan, rows, specification)
  adapter$on_stratum(d, d$stratum[1L]); view <- adapter$finish()
  z <- view$contrasts[view$contrasts$contrast == "implementation_and_acceptance", ]
  expect_true(all(z$n_unmatched_stage_a == 1L & z$n_unmatched_stage_b == 0L))
  expect_true(all(z$n_attempted_pairs == 1L))
})

test_that("Gaussian tasks have no version-stage comparison", {
  d <- fixture(); d <- d[d$stage == "original_refit", ]; d$family <- d$stage <- "challenge"
  adapter <- .paper_stage_publication(plan, selected(d), specification)
  adapter$on_stratum(d, d$stratum[1L]); view <- adapter$finish()
  expect_equal(nrow(view$contrasts), 0L)
  expect_true(all(view$selection$status == "gaussian_has_no_implementation_stage_pair"))
})

test_that("changed specifications and unbounded output cannot silently enter exports", {
  d <- fixture(); rows <- selected(d)
  changed <- spec; changed$bootstrap$draws <- 1001L
  path <- tempfile(fileext = ".json"); on.exit(unlink(path))
  jsonlite::write_json(changed, path, auto_unbox = TRUE)
  expect_error(.paper_stage_publication(plan, rows, path), "specification differs")
  expect_error(.paper_stage_publication(plan, rows[c(1L, 1L), ], specification), "summary identities")
  expect_error(.paper_stage_publication(plan, rows, specification, scientific = NA), "Invalid")
  adapter <- .paper_stage_publication(plan, rows, specification, max_output_rows = 11L)
  expect_error(adapter$on_stratum(d, d$stratum[1L]), "output memory bound")
  expect_equal(nrow(.paper_stage_publication(plan, rows, specification)$finish()$contrasts), 0L)
})
