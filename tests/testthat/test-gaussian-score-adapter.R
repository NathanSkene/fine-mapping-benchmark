test_that("Gaussian SuSiE uses the matched single-effect Gaussian likelihood", {
  skip_if_not_installed("susieR")
  z <- c(3, 1, 0)
  n <- 100
  v <- .1
  fit <- suppressMessages(run_susie(z, diag(3), n, L = 1,
    statistic_type = "normal_score", prior_variance = v,
    estimate_prior_variance = FALSE, min_abs_corr = 0))
  log_bf <- -.5 * log1p(n * v) + .5 * z^2 * n * v / (1 + n * v)
  weights <- exp(log_bf - max(log_bf)); weights <- weights / sum(weights)
  expect_null(fit$error)
  expect_equal(fit$pip, weights, tolerance = 1e-7)
  expect_equal(fit$additional$posterior_mean,
    weights * z * sqrt(n) * v / (1 + n * v), tolerance = 1e-7)
  expect_false(fit$params$effective_estimate_residual_variance)
  old <- suppressMessages(run_susie(z, diag(3), n, L = 1,
    prior_variance = v, estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE, min_abs_corr = 0))
  direct <- susieR::susie_rss(z = z, R = diag(3), n = n, L = 1,
    prior_variance = v, estimate_prior_variance = FALSE,
    estimate_residual_variance = FALSE, min_abs_corr = 0)
  expect_equal(old$pip, as.numeric(direct$pip))
})

test_that("public region adapters retain the declared likelihood without manual flags", {
  skip_if_not_installed("susieR")
  input <- generate_challenge(n_replicates = 1L, p = 12L, n = 100L, seed = 730L)$inputs[[1L]]
  input$annotations[,] <- 0
  args <- list(L = 1L, estimate_prior_variance = FALSE, min_abs_corr = 0)
  for (pair in list(list(wrapper = run_susie_region, runner = run_susie),
                   list(wrapper = run_polyfun_est_region, runner = run_polyfun_est),
                   list(wrapper = run_polyfun_ldsc_region, runner = run_polyfun_ldsc))) {
    direct <- suppressMessages(do.call(pair$runner, c(list(z = input$z, LD = input$LD,
      n = input$n, statistic_type = "normal_score"), args)))
    pred <- suppressMessages(challenge_region_adapter(pair$wrapper, args)(input))
    expect_equal(pred$pip, direct$pip, tolerance = 1e-7)
  }
  expect_error(challenge_region_adapter(run_susie_inf_region)(input), "unsupported_input")
  expect_error(challenge_region_adapter(run_finemap_inf_region)(input), "unsupported_input")
})

test_that("annotation-prior wrappers preserve the Gaussian likelihood type", {
  skip_if_not_installed("susieR")
  z <- c(3, 1, 0)
  args <- list(z = z, LD = diag(3), n = 100, L = 1,
    statistic_type = "normal_score", estimate_prior_variance = FALSE,
    min_abs_corr = 0)
  reference <- suppressMessages(do.call(.fm_susie_rss,
    c(list(z = z, R = diag(3), n = 100, statistic_type = "normal_score"),
      list(L = 1, estimate_prior_variance = FALSE, min_abs_corr = 0))))
  for (method in list(run_polyfun_est, run_polyfun_oracle, run_polyfun_ldsc)) {
    fit <- suppressMessages(do.call(method, args))
    expect_null(fit$error)
    expect_equal(fit$pip, as.numeric(reference$pip), tolerance = 1e-7)
  }
})
