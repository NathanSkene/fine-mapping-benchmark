test_that("known overlap is a normalized study-inclusion Gram", {
  case <- .challenge_one("meta_known_overlap", p = 12L, n = 10000L,
    signal_z = 6, prior_sd = .05, reference_n = 30L,
    split = "final", serial = 1L, seed = 120361066L)
  input <- case$input
  included <- input$n_per_variant == input$n
  study_weights <- cbind(1, as.numeric(included)) / sqrt(1 + included)
  expected <- tcrossprod(study_weights)
  values <- eigen(input$sample_overlap, symmetric = TRUE, only.values = TRUE)$values
  expect_equal(input$sample_overlap, expected, tolerance = 1e-14)
  expect_equal(diag(input$sample_overlap), rep(1, 12L))
  expect_gte(min(values), -1e-12)
  expect_equal(sum(values > 1e-10), 2L)
  expect_true(validate_challenge_input(input))
  expect_true(all(case$truth$LD_true != 0))
  expect_equal(input$sample_overlap,
    case$truth$score_covariance_true / case$truth$LD_true, tolerance = 1e-14)

  # Disclosure of overlap must not alter the sampled scores or noise model.
  unknown <- .challenge_one("meta_missingness", p = 12L, n = 10000L,
    signal_z = 6, prior_sd = .05, reference_n = 30L,
    split = "final", serial = 1L, seed = 120361066L)
  expect_null(unknown$input$sample_overlap)
  expect_identical(input$z, unknown$input$z)
  expect_identical(case$truth$score_covariance_true, unknown$truth$score_covariance_true)

  # A genuine invalid covariance remains invalid; no PSD tolerance was relaxed.
  bad <- input
  bad$sample_overlap[1, 2] <- bad$sample_overlap[2, 1] <- 1.1
  expect_error(validate_challenge_input(bad), "Invalid LD matrix")
})

test_that("p1000 overlap survives zero and subnormal LD entries", {
  # Actual final generation-5 case: seed 20361016, serial 50, final offset 1e8.
  case <- .challenge_one("meta_known_overlap", p = 1000L, n = 10000L,
    signal_z = 6, prior_sd = .05, reference_n = 30L,
    split = "final", serial = 50L, seed = 120361066L)
  input <- case$input
  R <- case$truth$LD_true
  zeros <- R == 0
  subnormal <- R > 0 & R < .Machine$double.xmin
  expect_true(any(zeros))
  expect_true(any(subnormal))
  expect_true(all(is.finite(input$sample_overlap)))
  expect_equal(diag(input$sample_overlap), rep(1, 1000L))
  expect_true(validate_challenge_input(input))

  included <- input$n_per_variant == input$n
  study_weights <- cbind(1, as.numeric(included)) / sqrt(1 + included)
  expect_equal(input$sample_overlap, tcrossprod(study_weights), tolerance = 1e-14)
  expect_equal(case$truth$score_covariance_true, R * input$sample_overlap, tolerance = 1e-14)
  expect_true(all(case$truth$score_covariance_true[zeros] == 0))
  old <- case$truth$score_covariance_true / R
  expect_true(all(is.nan(old[zeros])))
  # Nonzero is insufficient: division also amplifies subnormal roundoff.
  expect_gt(max(abs(old[subnormal] - input$sample_overlap[subnormal])), .1)
  normal <- R >= .Machine$double.xmin
  expect_equal(old[normal], input$sample_overlap[normal], tolerance = 1e-14)
})
