test_that("challenge inputs never expose latent mechanisms or generator seeds", {
  panel <- generate_challenge(n_replicates = 1, p = 12)
  expect_length(panel$inputs, nrow(challenge_scenarios()))
  for (input in panel$inputs) {
    expect_true(validate_challenge_input(input))
    expect_false(any(c("truth", "seed", "scenario", "causal_indices") %in% names(input)))
  }
  expect_equal(generate_challenge(n_replicates = 1, p = 12), panel)
  final <- generate_challenge(n_replicates = 1, p = 12, split = "final")
  expect_false(any(names(final$inputs) %in% names(panel$inputs)))
  expect_false(identical(final$inputs[[1]]$z, panel$inputs[[1]]$z))
  set.seed(991); before <- .Random.seed
  generate_challenge("null", n_replicates = 1, p = 12)
  expect_identical(.Random.seed, before)
})

test_that("conditional Gaussian sampler preserves covariance including singular twins", {
  R <- .challenge_ar(5, .7)
  set.seed(70); x <- .challenge_mvn(R, 10000)
  expect_lt(max(abs(stats::cov(x) - R)), .06)
  R[, 2] <- R[, 1]; R[2, ] <- R[1, ]; R[2, 2] <- 1
  twins <- .challenge_mvn(R, 100)
  expect_lt(max(abs(twins[, 1] - twins[, 2])), 1e-7)
  expect_error(.challenge_mvn(matrix(c(1, 2, 2, 1), 2)), "semidefinite")
})

test_that("single-effect posterior matches a hand-computed null-aware Bayes factor", {
  input <- generate_challenge("null", n_replicates = 1, p = 12)$inputs[[1]]
  input$LD <- diag(12); input$z <- seq(-2, 3, length.out = 12)
  W <- .05^2; shrink <- input$n * W / (1 + input$n * W)
  bf <- exp(-.5 * log1p(input$n * W) + .5 * shrink * input$z^2)
  norm <- .5 + sum(.5 / 12 * bf)
  pred <- challenge_exact(input, prior_k = c(.5, .5, 0))
  expect_equal(pred$pip, .5 / 12 * bf / norm, tolerance = 1e-12)
  expect_equal(pred$null_probability, .5 / norm, tolerance = 1e-12)
  expect_equal(sum(pred$pip) + pred$null_probability, 1, tolerance = 1e-12)
})

test_that("exact reference respects allele recoding and arbitrary input permutation", {
  input <- generate_challenge("multi_sparse", n_replicates = 1, p = 12)$inputs[[1]]
  original <- challenge_exact(input)
  sign <- rep(c(-1, 1), 6)
  recoded <- input; recoded$z <- input$z * sign
  recoded$LD <- input$LD * outer(sign, sign)
  expect_equal(challenge_exact(recoded)$pip, original$pip, tolerance = 1e-10)
  perm <- 12:1; shuffled <- input
  for (field in c("variant_id", "z", "beta_hat", "se")) shuffled[[field]] <- input[[field]][perm]
  shuffled$LD <- input$LD[perm, perm]; shuffled$annotations <- input$annotations[perm, , drop = FALSE]
  expect_equal(challenge_exact(shuffled)$pip[perm], original$pip, tolerance = 1e-10)
})

test_that("perfect LD cannot identify which twin caused the association", {
  panel <- generate_challenge("perfect_ld", n_replicates = 1, p = 12)
  input <- panel$inputs[[1]]; pred <- challenge_exact(input)
  expect_equal(input$z[1], input$z[2], tolerance = 1e-6)
  expect_equal(pred$pip[1], pred$pip[2], tolerance = 1e-7)
  forced <- challenge_forced_one(input)
  expect_lte(max(forced$pip[1:2]), .5 + 1e-7)
})

test_that("matched prior posterior passes an independent Bayesian calibration check", {
  # This checks the model under its own random-S and normal-effect prior.
  # Fixed-S stress strata have no universal Bayesian calibration guarantee.
  panel <- generate_challenge("matched_prior", n_replicates = 240, p = 12, seed = 221)
  pred <- run_challenge(panel$inputs, list(exact = challenge_exact))$exact
  residual <- vapply(names(panel$inputs), function(id)
    length(panel$truth[[id]]$causal_variant_id) - sum(pred[[id]]$pip), 0)
  expect_lt(abs(mean(residual)), 4 * stats::sd(residual) / sqrt(length(residual)))
  hits <- vapply(names(panel$inputs), function(id) any(vapply(
    pred[[id]]$credible_configurations,
    function(s) setequal(s, panel$truth[[id]]$causal_variant_id), TRUE)), TRUE)
  expect_gte(mean(hits), .9)
})

test_that("missing causal variants limit recovery and never become negative truth", {
  panel <- generate_challenge("missing_causal", n_replicates = 1, p = 12)
  input <- panel$inputs[[1]]; truth <- panel$truth[[1]]
  expect_equal(length(input$z), 11)
  expect_equal(truth$n_missing_causal, 1)
  expect_false(any(truth$causal_variant_id %in% input$variant_id))
  result <- evaluate_challenge(panel, run_challenge(panel$inputs))
  expect_true(all(is.na(result$tasks$ap)))
  expect_equal(result$tasks$recall, c(0, 0))
  expect_equal(result$summary$missing_causal_fraction, c(1, 1))
})

test_that("meta-analysis effective sample sizes and overlap are internally coherent", {
  panel <- generate_challenge("meta_missingness", n_replicates = 1, p = 12)
  input <- panel$inputs[[1]]; truth <- panel$truth[[1]]
  expect_equal(input$se, 1 / sqrt(input$n_per_variant))
  expect_equal(input$beta_hat / input$se, input$z)
  expect_equal(diag(truth$score_covariance_true), rep(1, 12))
  expect_equal(truth$score_covariance_true[1, 2] / truth$LD_true[1, 2], 1 / sqrt(2))
  expect_error(challenge_exact(input), "sample-overlap")
  known <- generate_challenge("meta_known_overlap", n_replicates = 1, p = 12)$inputs[[1]]
  expect_true(all(is.finite(challenge_exact(known)$pip)))
})

test_that("failures, omitted tasks, invalid PIPs and duplicates remain visible", {
  panel <- generate_challenge("single_sparse", n_replicates = 2, p = 12)
  good <- run_challenge(panel$inputs, list(exact = challenge_exact))
  good$exact[[1]]$pip[1] <- 2
  result <- evaluate_challenge(panel, good)
  expect_equal(result$summary$success_rate, .5)
  expect_equal(result$tasks$recall[1], 0)
  expect_match(result$tasks$error[1], "PIPs")
  good$exact <- good$exact[0]
  result <- evaluate_challenge(panel, good)
  expect_equal(result$summary$success_rate, 0)
  expect_true(is.na(result$summary$brier_successful))
  expect_equal(result$summary$recall_all_attempts, 0)
  bad <- run_challenge(panel$inputs, list(exact = challenge_exact))
  bad$exact <- c(bad$exact, bad$exact[1])
  expect_error(evaluate_challenge(panel, bad), "duplicate")
})

test_that("configuration collections need declared levels and native sets need IDs", {
  panel <- generate_challenge("null", n_replicates = 1, p = 12)
  input <- panel$inputs[[1]]; pred <- challenge_exact(input)
  pred$configuration_coverage <- NULL
  expect_error(.challenge_validate_prediction(pred, input), "declared")
  pred <- challenge_exact(input); pred$credible_sets <- list(1L)
  expect_error(.challenge_validate_prediction(pred, input), "credible_sets")
  expect_equal(.challenge_distinct_hits(list(c("a", "b"), "a"), c("a", "b")), 2)
  expect_equal(.challenge_distinct_hits(list("a", "a"), c("a", "b")), 1)
  pred <- challenge_exact(input); pred$credible_configurations <- NULL
  pred$configuration_coverage <- c(.5, .95)
  result <- evaluate_challenge(panel, list(bad = stats::setNames(list(pred), input$task_id)))
  expect_equal(nrow(result$tasks), 1)
  expect_equal(result$summary$success_rate, 0)
})

test_that("nested truth leakage and inconsistent assessor labels are rejected", {
  panel <- generate_challenge("multi_ancestry_shared", n_replicates = 1, p = 12)
  leaked <- panel$inputs[[1]]; leaked$cohorts[[1]]$causal_variant_id <- "secret"
  expect_error(validate_challenge_input(leaked), "fields")
  bad <- panel; bad$truth[[1]]$n_missing_causal <- -1
  expect_error(evaluate_challenge(bad, run_challenge(panel$inputs)), "inconsistent")
  names(panel$inputs)[1] <- "wrong"
  expect_error(run_challenge(panel$inputs), "identifiers")
})

test_that("ranking is tie-invariant and proper scores distinguish confidence", {
  expect_equal(.challenge_ap(c(1, 0, 1, 0), rep(.5, 4)), .5)
  expect_equal(.challenge_ap(c(0, 1, 0, 1), rep(.5, 4)), .5)
  expect_true(is.na(.challenge_ap(c(0, 0), c(.9, .1))))
  expect_equal(.challenge_ap(c(1, 0), c(1, 0)), 1)
  panel <- generate_challenge("null", n_replicates = 2, p = 12)
  zero <- function(x) list(task_id = x$task_id, variant_id = x$variant_id,
    pip = rep(0, length(x$z)), status = "ok", credible_sets = list())
  high <- function(x) { p <- zero(x); p$pip[1] <- 1; p }
  result <- evaluate_challenge(panel, run_challenge(panel$inputs, list(zero = zero, high = high)))
  expect_equal(result$summary$null_false_call_rate_successful, c(1, 0))
  expect_equal(result$summary$brier_successful, c(1 / 12, 0))
  expect_true(all(is.finite(result$tasks$log_loss_clipped)))
})

test_that("paired bootstrap preserves locus pairs and honest availability", {
  panel <- generate_challenge("null", n_replicates = 4, p = 12)
  result <- evaluate_challenge(panel, run_challenge(panel$inputs))
  same <- challenge_compare(result, "exact", "exact", "brier", "null", 30)
  expect_equal(unlist(same[c("difference", "lower", "upper")]), c(difference = 0, lower = 0, upper = 0))
  expect_equal(same$n_clusters, 4)
  a <- challenge_compare(result, "exact", "forced_one", "brier", "null", 30)
  b <- challenge_compare(result, "exact", "forced_one", "brier", "null", 30)
  expect_equal(a, b)
})

test_that("release integrity rejects tampering and keeps truth outside public tree", {
  parent <- tempfile(); dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE))
  panel <- generate_challenge("null", n_replicates = 1, p = 12)
  pub <- file.path(parent, "public"); priv <- file.path(parent, "private")
  expect_error(write_challenge_release(panel, pub, file.path(pub, "truth")), "separate")
  write_challenge_release(panel, pub, priv)
  read <- read_challenge_release(pub)
  expect_identical(read$inputs, panel$inputs)
  expect_false(any(c("truth.rds", "generator.rds") %in% list.files(pub)))
  expect_error(write_challenge_release(panel, pub, priv), "overwrite")
  writeBin(as.raw(1:4), file.path(pub, "inputs.rds"))
  expect_error(read_challenge_release(pub), "hash mismatch")
})
