# Dependency-free scoring fixtures: truth and expected values are hand computed.
.ev_fixture <- function(pips, causal = rep(list(1L), length(pips)),
                        sets = lapply(causal, function(x) if (length(x)) list(x) else list()),
                        iters = seq_along(pips), p = vapply(pips, length, integer(1))) {
  if (length(unique(p)) != 1L) stop("This fixture uses one fixed variant universe")
  scenarios <- lapply(seq_along(pips), function(i) {
    list(scenario_id = i, S = length(causal[[i]]), phi = 0.3, iter = iters[i],
         regions = list(list(truth = list(causal_indices = causal[[i]],
                                          beta_true = rep(0, p[i])))))
  })
  fits <- lapply(seq_along(pips), function(i) {
    list(scenario_id = i, region_id = 1L, iter = iters[i], S = length(causal[[i]]),
         phi = 0.3, pip = pips[[i]], credible_sets = sets[[i]], runtime_seconds = i)
  })
  list(simulation = list(scenarios = scenarios, genotypes = list(list(LD = diag(p[1L]),
       p = p[1L], maf = rep(0.2, p[1L]))), params = list(p = p[1L], p_actual = p[1L])),
       results = list(methods_run = "fixture", fixture = list(results = fits)))
}
.ev_score <- function(x, ...) {
  evaluate_methods(x$simulation, x$results, verbose = FALSE, ...)$fixture$global
}

test_that("exact AP handles endpoint perfection, compressed scores and ties", {
  expect_equal(.compute_ap_exact(c(1, 0), c(TRUE, FALSE)), 1)
  expect_equal(.compute_ap_exact(c(0.0033, 0.0032, 0.0031), c(FALSE, TRUE, FALSE)), 1 / 2)
  expect_equal(.compute_ap_exact(rep(0.5, 4), c(FALSE, TRUE, FALSE, TRUE)), 1 / 2)
  expect_equal(.compute_ap_exact(rep(0.5, 4), c(TRUE, TRUE, FALSE, FALSE)), 1 / 2)
  expect_equal(.compute_ap_exact(c(0.9, 0.8, 0.8, 0.1), c(FALSE, TRUE, FALSE, TRUE)),
               (1 / 3 + 1 / 2) / 2)
  expect_equal(.compute_ap_exact(0.3, TRUE), 1)
  expect_true(is.na(.compute_ap_exact(c(0.9, 0.1), c(FALSE, FALSE))))
  expect_error(.compute_ap_exact(c(0.5, NA), c(TRUE, FALSE)), "finite")
  expect_error(.compute_ap_exact(c(0.5, 0.1), TRUE), "matching")
})

test_that("public AP is exact and grid independent; pooled and macro AP differ", {
  x <- .ev_fixture(list(c(0.9, 0.5, 0.1), c(0.6, 0.4, 0.2)),
                   causal = list(1L, 3L))
  fine <- .ev_score(x)
  coarse <- .ev_score(x, pip_thresholds = c(0, 1))
  expect_equal(fine$ap, 0.7)
  expect_equal(fine$ap_macro, 2 / 3)
  expect_equal(fine$ap, coarse$ap)
  expect_equal(fine$ap_macro, coarse$ap_macro)
  expect_equal(fine$ap_n_fits, 2L)
})

test_that("threshold curves match brute force including ties and endpoints", {
  set.seed(71)
  p <- sample(c(0, 0.1, 0.5, 0.9, 1), 200, replace = TRUE)
  y <- seq_along(p) %% 11 == 0
  t <- c(1, 0.5, 0, 0.51, 0.9, 0.9)
  curve <- .fdr_power_curve(p, y, t)
  expect_equal(curve$tp, vapply(t, function(z) sum(y & p >= z), integer(1)))
  expect_equal(curve$fp, vapply(t, function(z) sum(!y & p >= z), integer(1)))
  expect_equal(curve$fn, sum(y) - curve$tp)
  expect_equal(curve$recall, curve$tp / sum(y))
})

test_that("null loci retain false positives and undefined recall", {
  x <- .ev_fixture(list(c(0.9, 0.1)), causal = list(integer()), sets = list(list(1L)))
  g <- .ev_score(x, pip_thresholds = c(0, 0.5, 1))
  expect_equal(g$fdr_power_curve$fp, c(2, 1, 0))
  expect_equal(g$fdr_power_curve$fdr, c(1, 1, 0))
  expect_equal(g$fdr_power_curve$precision, c(0, 0, 1))
  expect_true(all(is.na(g$fdr_power_curve$recall)))
  expect_true(is.na(g$ap))
  expect_equal(g$cs_coverage, 0)
  expect_equal(g$null_region_false_call_rate, 1)
  expect_equal(g$null_region_cs_rate, 1)
  expect_equal(g$n_null_fits, 1L)
})

test_that("failures and omitted outputs remain in attempted power denominators", {
  x <- .ev_fixture(list(c(0.9, 0.1), c(0.9, 0.1)))
  x$results$fixture$results[[2]]$error <- "timeout"
  g <- .ev_score(x, pip_thresholds = c(0.5, 1))
  expect_equal(g$fdr_power_curve$power[1], 1)
  expect_equal(g$fdr_power_curve_attempted$power[1], 1 / 2)
  expect_equal(g$fdr_power_curve_attempted$fn[1], 1)
  expect_equal(g$cs_power, 1)
  expect_equal(g$cs_power_attempted, 1 / 2)
  expect_equal(g$ap_macro_attempted, 1 / 2)
  expect_equal(g$success_rate, 1 / 2)
  expect_equal(g$n_failed, 1L)
  x$results$fixture$results <- x$results$fixture$results[1]
  omitted <- .ev_score(x)
  expect_equal(omitted$n_fits, 2L)
  expect_equal(omitted$n_missing, 1L)
  expect_equal(omitted$ap_macro_attempted, 1 / 2)
  x$results$fixture$results <- list()
  none <- .ev_score(x)
  expect_equal(none$n_failed, 2L)
  expect_equal(none$n_missing, 2L)
  expect_null(none$fdr_power_curve)
  expect_equal(none$fdr_power_curve_attempted$power, rep(0, 201))
  expect_equal(none$ap_macro_attempted, 0)
  expect_equal(none$cs_power_attempted, 0)
})

test_that("invalid submissions fail explicitly rather than changing the truth universe", {
  x <- .ev_fixture(list(c(0.9, 0.1)))
  for (p in list(0.9, c(0.9, NA), c(0.9, Inf), c(1.01, 0), c(-0.01, 0))) {
    bad <- x
    bad$results$fixture$results[[1]]$pip <- p
    g <- .ev_score(bad)
    expect_equal(g$n_invalid, 1L)
    expect_equal(g$n_failed, 1L)
    expect_equal(nrow(g$validation_errors), 1L)
    expect_equal(g$cs_power_attempted, 0)
  }
  x$results$fixture$results[[1]]$credible_sets <- list(c(1, 3))
  expect_equal(.ev_score(x)$n_invalid, 1L)
  x$results$fixture$results[[1]]$credible_sets <- list(integer())
  expect_equal(.ev_score(x)$n_invalid, 1L)
})

test_that("attempt identity and authoritative scenario labels cannot be gamed", {
  x <- .ev_fixture(list(c(0.9, 0.1)))
  x$results$fixture$results[[1]]$S <- 900L
  x$results$fixture$results[[1]]$iter <- 900L
  ev <- evaluate_methods(x$simulation, x$results, verbose = FALSE)
  expect_named(ev$fixture$by_S, "1")
  expect_equal(ev$fixture$global$n_replicate_clusters, 1L)
  x$results$fixture$results <- rep(x$results$fixture$results, 2)
  expect_error(.ev_score(x), "Duplicate")
  x$results$fixture$results <- x$results$fixture$results[1]
  x$results$fixture$results[[1]]$region_id <- 2L
  expect_error(.ev_score(x), "outside")
})

test_that("overlapping credible sets use maximum one-to-one signal matching", {
  x <- .ev_fixture(list(c(0.8, 0.7, 0.1)), causal = list(c(1L, 2L)),
                   sets = list(list(c(1L, 2L), 1L)))
  g <- .ev_score(x)
  expect_equal(g$cs_matched_precision, 1)
  expect_equal(g$cs_matched_power, 1)
  x$results$fixture$results[[1]]$credible_sets <- list(c(1L, 1L), 1L)
  duplicate <- .ev_score(x)
  expect_equal(duplicate$cs_coverage, 1)
  expect_equal(duplicate$cs_matched_precision, 1 / 2)
  expect_equal(duplicate$cs_matched_power, 1 / 2)
  expect_equal(duplicate$cs_duplicate_sets, 1L)
  expect_equal(duplicate$cs_size_mean, 1)
  x$results$fixture$results[[1]]$credible_sets <- list()
  empty <- .ev_score(x)
  expect_true(is.na(empty$cs_coverage))
  expect_equal(empty$cs_power, 0)
})

test_that("proper scores retain raw within-bin error and impossible-event loss", {
  x <- .ev_fixture(list(c(0.1, 0.4)))
  g <- .ev_score(x, n_pip_cal_bins = 1L)
  expect_equal(g$brier_score, (0.9^2 + 0.4^2) / 2)
  expect_equal(g$log_loss, (-log(0.1) - log(0.6)) / 2)
  b <- g$pip_calibration
  m <- .murphy_decomposition(b$n, b$n_causal, b$sum_pip, b$sum_pip_sq, b$sum_causal_pip)
  expect_equal(m$bs, 0.3125)
  expect_equal(m$bs_raw, 0.485)
  expect_equal(m$bs_raw, g$brier_score)
  expect_true(is.na(.murphy_decomposition(b$n, b$n_causal, b$sum_pip)$bs_raw))
  certainty <- .ev_fixture(list(c(0, 1)))
  expect_equal(.ev_score(certainty)$log_loss, Inf)
  perfect <- .ev_fixture(list(c(1, 0)))
  expect_equal(.ev_score(perfect)$log_loss, 0)
  expect_equal(.ev_score(perfect)$brier_score, 0)
})

test_that("ranking diagnostics do not break ties by genomic index", {
  a <- .rank_set_metrics(rep(0.5, 4), c(TRUE, FALSE, FALSE, TRUE))
  b <- .rank_set_metrics(rep(0.5, 4), c(FALSE, TRUE, TRUE, FALSE))
  expect_equal(a, b)
  expect_equal(a$prec_at_1, 1 / 2)
  expect_equal(a$prec_at_S, 1 / 2)
  expect_equal(a$set_size_5, 4L)
  expect_equal(a$set_size_95, 4L)
  expect_error(.rank_set_metrics(c(0.5, 1.2), c(TRUE, FALSE)), "finite")
  expect_equal(.honesty(1, 1, 1)$rho, Inf)
})

test_that("SE resamples replicate clusters and the actual pooled estimator", {
  x <- .ev_fixture(list(c(0.9, 0.1), c(0.1, 0.9), c(0.2, 0.8)),
                   iters = c(1L, 2L, 2L))
  # Cluster 1 Brier=.01; cluster 2 Brier=(.81+.64)/2=.725.
  # Delete-one jackknife K=2 gives |.725-.01|/2, independent of how
  # many LD-correlated variant observations are expanded inside each cluster.
  g <- .ev_score(x)
  expect_equal(g$n_replicate_clusters, 2L)
  expect_equal(g$brier_score, (0.01 + 0.81 + 0.64) / 3)
  expect_equal(g$brier_score_se, (0.725 - 0.01) / 2)
  expect_equal(g$se_method, "delete-one-iter cluster jackknife")
  x$simulation$scenarios[[3]]$iter <- NA_integer_
  expect_true(is.na(.ev_score(x)$brier_score_se))
})

test_that("evaluation rejects malformed scoring parameters", {
  x <- .ev_fixture(list(c(0.9, 0.1)))
  expect_error(.ev_score(x, pip_thresholds = c(0, NA)), "finite")
  expect_error(.ev_score(x, pip_thresholds = c(-0.1, 1)), "finite")
  expect_error(.ev_score(x, n_pip_cal_bins = 0), "positive integer")
  expect_error(.ev_score(x, n_pip_cal_bins = 1.5), "positive integer")
})


test_that("single-threshold and single-bin uncertainty preserve matrix dimensions", {
  x <- .ev_fixture(list(c(0.9, 0.1), c(0.1, 0.9)))
  g <- .ev_score(x, pip_thresholds = 0.5, n_pip_cal_bins = 1L)
  expect_equal(nrow(g$fdr_power_curve), 1L)
  expect_equal(g$fdr_power_curve$power_se, 1 / 2)
  expect_equal(g$pip_calibration$frac_causal_se, 0)
})
