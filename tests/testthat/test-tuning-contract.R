source(testthat::test_path("..", "..", "scripts", "tuning", "tuning_contract.R"), local = TRUE)

testthat::test_that("tuning uses a common locus panel with authoritative region metadata", {
  genotypes <- lapply(c(10L, 20L, 40L), function(p) list(p = p))
  regions <- lapply(seq_along(genotypes), function(k) list(region_id = k, truth = list(causal_indices = k)))
  sim <- list(genotypes = genotypes,
              scenarios = list(list(scenario_id = 1L, regions = regions),
                               list(scenario_id = 2L, regions = regions)),
              params = list(p = c(12L, 22L, 42L), p_actual = c(10L, 20L, 40L),
                            n_regions = 3L, n_scenarios = 2L))
  a <- subset_tuning_scenario(sim, 1L, 2L)
  b <- subset_tuning_scenario(sim, 2L, 2L)
  testthat::expect_identical(a$params$p, c(12L, 22L))
  testthat::expect_identical(a$params$p_actual, c(10L, 20L))
  testthat::expect_identical(a$params$source_region_indices, b$params$source_region_indices)
  testthat::expect_identical(b$scenarios[[1]]$scenario_id, 1L)
  testthat::expect_identical(b$params$source_scenario_index, 2L)
  testthat::expect_identical(b$params$n_regions, 2L)
  testthat::expect_identical(length(b$scenarios[[1]]$regions), 2L)
  testthat::expect_identical(sim$params$n_regions, 3L)
  testthat::expect_error(subset_tuning_scenario(sim, 1L, rotate_regions = TRUE), "confounds")
  testthat::expect_error(subset_tuning_scenario(sim, 1L, 4L), "exceeds")
  testthat::expect_error(subset_tuning_scenario(sim, 1.5), "integer")
  testthat::expect_error(subset_tuning_scenario(sim, 3L), "outside")
  sim$params$p_actual <- c(10L, 20L)
  testthat::expect_error(subset_tuning_scenario(sim, 1L), "metadata")
})

testthat::test_that("tuning ranking cannot silently substitute pooled or success-only AP", {
  testthat::expect_identical(tuning_attempted_ap(list(ap = 1, ap_macro = .8,
                                                    ap_macro_attempted = .4)), .4)
  testthat::expect_error(tuning_attempted_ap(list(ap = 1)), "upgrade")
  testthat::expect_error(tuning_attempted_ap(list(ap_macro_attempted = 1.1)), "Invalid")
  testthat::expect_true(is.na(tuning_attempted_ap(list(ap_macro_attempted = NA_real_))))
})
