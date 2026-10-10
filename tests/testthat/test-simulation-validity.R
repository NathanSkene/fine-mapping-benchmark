# Deterministic scientific invariants, independent of external genotype files.
simulation_validity_genotypes <- function(n = 80L, p = 6L, regions = 2L) {
  set.seed(310)
  lapply(seq_len(regions), function(i) {
    X <- scale(matrix(rnorm(n * p), n, p))
    list(X = X, X_raw = X, LD = cor(X), n = n, p = p,
         region_id = i, variant_ids = paste0("r", i, "v", seq_len(p)))
  })
}

test_that("summary statistics equal intercept OLS for uncentred dosages", {
  set.seed(12)
  X <- matrix(rbinom(300, 2, .3), 100, 3)
  y <- 5 + as.numeric(X %*% c(.5, -.3, .1)) + rnorm(100)
  got <- compute_summary_statistics(X, y)
  expected <- vapply(seq_len(ncol(X)), function(j) {
    coef(summary(lm(y ~ X[, j])))[2, c(1, 2, 3)]
  }, numeric(3))
  expect_equal(got$beta_hat, unname(expected[1, ]), tolerance = 1e-12)
  expect_equal(got$se, unname(expected[2, ]), tolerance = 1e-12)
  expect_equal(got$z, unname(expected[3, ]), tolerance = 1e-12)
  shifted <- compute_summary_statistics(sweep(X, 2, c(10, -5, 3), "+"), y + 100)
  expect_equal(got$z, shifted$z, tolerance = 1e-12)
  expect_equal(got$statistic_type, "marginal_ols_t")
  expect_equal(got$degrees_of_freedom, 98L)
})

test_that("the t-to-correlation identity matches standardized OLS", {
  g <- simulation_validity_genotypes(regions = 1L)[[1]]
  set.seed(4)
  y <- as.numeric(scale(g$X[, 1] + rnorm(g$n)))
  ss <- compute_summary_statistics(g$X, y)
  expect_equal(ss$z / sqrt(g$n - 2 + ss$z^2), ss$beta_hat, tolerance = 1e-12)
  expect_error(compute_summary_statistics(matrix(1, 4, 2), 1:4), "nonzero variance")
  expect_error(compute_summary_statistics(matrix(1:4, 2), 1:2), "three individuals")
})

test_that("sample genetic variation and target PVE are recorded separately", {
  g <- simulation_validity_genotypes(regions = 1L)[[1]]
  set.seed(2)
  fit <- generate_phenotype_sparse(g$X, c(1, rep(0, g$p - 1)), .4, "normal", 1L)
  expect_equal(fit$pve_expected, .4, tolerance = 1e-12)
  expect_equal(fit$pve_realised, var(g$X[, 1]) / var(fit$y), tolerance = 1e-12)
  expect_gt(abs(fit$pve_realised - .4), 1e-6)
  v <- fit$variance_components
  expect_equal(v$phenotype, v$genetic + v$residual_realised +
                 2 * v$genetic_residual_covariance, tolerance = 1e-12)
})

test_that("background covariance and all effects survive phenotype simulation", {
  g <- simulation_validity_genotypes(regions = 1L)
  fit <- simulate_phenotypes(g, S = 2, phi = .4, model = "sparse_inf",
                            p_causal = .6, inf_model = "susie_inf", seed = 2,
                            verbose = FALSE)[[1]]
  t <- fit$truth
  expect_equal(t$beta_total, t$beta_sparse + t$beta_infinitesimal)
  expect_equal(t$estimand, "sparse_component_membership")
  expect_true(all(t$beta_infinitesimal != 0))
  v <- t$variance_components
  expect_equal(var(as.numeric(g[[1]]$X %*% t$beta_total)), v$genetic, tolerance = 1e-12)
  expect_equal(v$genetic, v$sparse + v$infinitesimal +
                 2 * v$sparse_infinitesimal_covariance, tolerance = 1e-12)
  expect_gt(abs(v$sparse_infinitesimal_covariance), 1e-6)
  expect_equal(t$pve_expected, .4, tolerance = 1e-12)
})

test_that("p_causal one does not demand an infinitesimal component", {
  g <- simulation_validity_genotypes(regions = 1L)
  fit <- simulate_phenotypes(g, S = g[[1]]$p, phi = .2, model = "sparse_inf",
                            p_causal = 1, inf_model = "beatrice", seed = 8,
                            verbose = FALSE)[[1]]
  expect_equal(fit$truth$beta_infinitesimal, rep(0, g[[1]]$p))
  expect_equal(fit$truth$pve_expected, .2, tolerance = 1e-12)
})

test_that("normal effect variance is a nuisance scale at fixed PVE", {
  g <- simulation_validity_genotypes(regions = 1L)
  a <- simulate_phenotypes(g, S = 2, phi = .2, effect_variance = .36,
                          seed = 10, verbose = FALSE)[[1]]
  b <- simulate_phenotypes(g, S = 2, phi = .2, effect_variance = 1,
                          seed = 10, verbose = FALSE)[[1]]
  expect_equal(a$z, b$z, tolerance = 1e-12)
  expect_equal(a$truth$beta_sparse / .6, b$truth$beta_sparse, tolerance = 1e-12)
})

test_that("random enrichment defines one shared mapping across regions", {
  g <- simulation_validity_genotypes(regions = 3L)
  fit <- simulate_phenotypes(g, S = 2, phi = .2, annotations = "binary",
                             n_annotations = 2, enrichment = NULL,
                             seed = 10, verbose = FALSE)
  expect_equal(fit[[1]]$truth$enrichment, fit[[2]]$truth$enrichment)
  expect_equal(fit[[1]]$truth$enrichment, fit[[3]]$truth$enrichment)
  expect_true(all(fit[[1]]$truth$enrichment >= 2 & fit[[1]]$truth$enrichment <= 10))
})

test_that("local parameter-grid scenarios keep their sampled enrichment architecture", {
  g <- simulation_validity_genotypes(regions = 2L)
  sandbox <- new.env(parent = environment(run_simulation))
  sandbox$simulate_genotypes <- function(...) g
  runner <- run_simulation
  environment(runner) <- sandbox
  sim <- runner(n_regions = 2L, n = 80L, p = 6L, n_iter = 2L, S = 1L,
                phi = c(.1, .2), annotations = "binary", n_annotations = 2L,
                seed = 7, genetic_map_dir = NULL, verbose = FALSE)
  for (scenario in sim$scenarios) {
    for (region in scenario$regions) {
      expect_equal(region$truth$enrichment, sim$params$enrichment_actual)
    }
  }
})

test_that("causal probability saturation preserves requested genome-wide pi", {
  probabilities <- gwfm_calibrate_causal_probabilities(c(0, 20, 0), .8)
  expect_equal(probabilities, c(.7, 1, .7), tolerance = 1e-10)
  expect_equal(mean(probabilities), .8, tolerance = 1e-10)
  expect_true(all(probabilities >= 0 & probabilities <= 1))
  g <- simulation_validity_genotypes()
  g[[1]]$annotations_matrix <- matrix(1, g[[1]]$p, 1)
  g[[2]]$annotations_matrix <- matrix(0, g[[2]]$p, 1)
  set.seed(9)
  a <- gwfm_assign_causal_variants(g, .1, "binary", 10, 1)
  expect_equal(mean(unlist(lapply(a, `[[`, "causal_probabilities"))), .1, tolerance = 1e-10)
  expect_gt(mean(a[[1]]$causal_probabilities), mean(a[[2]]$causal_probabilities))
})

test_that("genome-wide null sparse draws retain requested background effects", {
  g <- simulation_validity_genotypes()
  a <- lapply(g, function(x) list(causal_indices = integer(0)))
  set.seed(14)
  expect_warning(fit <- gwfm_simulate_joint_phenotype(g, a, .3, "sparse_inf", .5,
                        "beatrice", "normal", .36, sum(vapply(g, `[[`, integer(1), "p"))),
                 "No nonzero sparse")
  expect_equal(fit$beta_scaled, lapply(g, function(x) rep(0, x$p)))
  expect_true(all(unlist(fit$beta_infinitesimal) != 0))
  expect_equal(fit$h2_expected, .3, tolerance = 1e-12)
  genetic <- Reduce("+", lapply(seq_along(g), function(i) {
    as.numeric(g[[i]]$X %*% (fit$beta_scaled[[i]] + fit$beta_infinitesimal[[i]]))
  }))
  expect_equal(fit$h2_realised, var(genetic) / var(fit$y), tolerance = 1e-12)
})

test_that("constant heterozygote columns are excluded from LD", {
  G <- cbind(rep(0, 4), rep(1, 4), rep(2, 4), c(0, 1, 1, 2))
  expect_equal(genotype_variable_columns(G), c(FALSE, FALSE, FALSE, TRUE))
})

test_that("null chi-squared intercept is not genetic heritability", {
  A <- cbind(c(0, 1, 0, 1), c(1, 0, 0, 1))
  LD <- diag(4)
  tau <- .fit_ldsc_stratified(rep(1, 4), .ldscore_matrix(A, LD), rowSums(LD^2), n = 100)
  expect_equal(tau, rep(0, 3))
  expect_equal(.tau_annot_to_prior_weights(A, tau), rep(.25, 4))
  expected_tau <- c(.001, .002, .003)
  ell <- .ldscore_matrix(A, LD)
  chi2 <- 1 + 100 * as.numeric(cbind(rowSums(LD^2), ell) %*% expected_tau)
  expect_equal(.fit_ldsc_stratified(chi2, ell, rowSums(LD^2), n = 100),
               expected_tau, tolerance = 1e-7)
})

test_that("LD-score prior training holds out complete chromosomes", {
  g <- simulation_validity_genotypes(p = 8L, regions = 3L)
  set.seed(1)
  for (i in seq_along(g)) {
    g[[i]]$chrom <- c(1, 1, 2)[i]
    g[[i]]$annotations_matrix <- matrix(rep(c(0, 1), 4), 8, 1)
  }
  r <- lapply(g, function(x) list(z = seq_len(x$p) / 2))
  original <- run_polyfun_ldsc_scenario_setup(g, r, list())$pooled_tau
  r[[1]]$z <- r[[1]]$z * 10
  r[[2]]$z <- r[[2]]$z * 20
  changed <- run_polyfun_ldsc_scenario_setup(g, r, list())$pooled_tau
  expect_equal(original[["1"]], changed[["1"]])
  expect_equal(original[["2"]], changed[["2"]])
  expect_equal(attr(original, "held_out_unit"), "chromosome")
})

test_that("rank scores are not advertised as posterior credible sets", {
  fit <- run_marginal_z(c(1, 2, 3))
  expect_equal(fit$credible_sets, list())
  expect_false(fit$additional$is_posterior_probability)
  expect_true(length(fit$additional$pip_mass_set) > 0)
  expect_error(run_marginal_z(c(1, NA_real_)), "non-empty numeric")
})

test_that("sampler validates posterior sample count and labels surrogate identity", {
  skip_if_not_installed("nnet")
  expect_error(run_sbayesrc(c(1, 2), diag(2), 50, n_iter = 10, burn_in = 10),
               "leave posterior samples")
  fit <- run_sbayesrc(c(1, 2), diag(2), 50, n_iter = 30, burn_in = 10, seed = 1)
  expect_equal(fit$credible_sets, list())
  expect_false(fit$additional$upstream_equivalent)
  expect_equal(fit$additional$implementation, "sbayesrc_style_surrogate")
  expect_true(all(fit$additional$component_assignments %in% 0:4))
})

test_that("genome-wide metadata works with verbose region dispatch", {
  g <- simulation_validity_genotypes(regions = 1L)
  sim <- list(genotypes = g, params = list(model = "sparse", n = 80L),
              scenarios = list(list(scenario_id = 1L, pi = .1, h2 = .2, iter = 1L,
                                    regions = list(list(z = seq_len(g[[1]]$p))))))
  expect_message(out <- run_methods(sim, "marginal_z", verbose = TRUE), "pi=")
  expect_equal(out$marginal_z$n_failed, 0L)
  expect_equal(out$marginal_z$results[[1]]$pi, .1)
  expect_equal(out$marginal_z$results[[1]]$h2, .2)
})

test_that("saved experiment labels change when LD or method arguments change", {
  params <- list(model = "sparse", n = 80L, n_iter = 1L, seed = 1L, phi_values = .2)
  first <- .make_sim_label(params)
  expect_identical(first, .make_sim_label(params))
  params$n_ref <- 20L
  expect_false(identical(first, .make_sim_label(params)))
  expect_false(identical(.make_sim_label(params), .make_sim_label(params, list(susie = list(L = 1)))))
})

simulation_validity_spy_runner <- function(region_runner, setup = NULL) {
  sandbox <- new.env(parent = environment(run_methods))
  sandbox$.FM_REGISTRY <- list(spy = "run_spy_region")
  sandbox$run_spy_region <- region_runner
  if (!is.null(setup)) sandbox$run_spy_scenario_setup <- setup
  runner <- run_methods
  environment(runner) <- sandbox
  runner
}

test_that("observed-data adapters and setup hooks do not receive truth", {
  g <- simulation_validity_genotypes(regions = 1L)
  sim <- list(genotypes = g, params = list(), scenarios = list(list(
    scenario_id = 1L, S = 1L, phi = .2, iter = 1L,
    regions = list(list(z = seq_len(g[[1]]$p), truth = list(causal_indices = 1L)))
  )))
  runner <- simulation_validity_spy_runner(
    function(region_geno, region_pheno, ...) {
      if (!is.null(region_pheno$truth)) stop("Truth reached adapter")
      list(pip = rep(.2, region_geno$p), credible_sets = list())
    },
    function(genotypes, regions, user_args) {
      if (!is.null(regions[[1]]$truth)) stop("Truth reached setup")
      list()
    }
  )
  out <- runner(sim, "spy", verbose = FALSE)
  expect_equal(out$spy$n_failed, 0L)
  expect_equal(out$spy$results[[1]]$benchmark_access, "observed_data")
})

test_that("scenario setup errors fail fits instead of changing the estimator", {
  g <- simulation_validity_genotypes(regions = 1L)
  sim <- list(genotypes = g, params = list(), scenarios = list(list(
    scenario_id = 1L, S = 1L, phi = .2, iter = 1L,
    regions = list(list(z = seq_len(g[[1]]$p)))
  )))
  runner <- simulation_validity_spy_runner(
    function(...) stop("Adapter must not run after failed setup"),
    function(...) stop("Pilot convergence failure")
  )
  out <- runner(sim, "spy", verbose = FALSE)
  expect_equal(out$spy$n_failed, 1L)
  expect_match(out$spy$results[[1]]$error, "Scenario setup failed: Pilot convergence failure")
})

test_that("nonfinite or misaligned method outputs are failures", {
  g <- simulation_validity_genotypes(regions = 1L)
  sim <- list(genotypes = g, params = list(), scenarios = list(list(
    scenario_id = 1L, S = 1L, phi = .2, iter = 1L,
    regions = list(list(z = seq_len(g[[1]]$p)))
  )))
  runner <- simulation_validity_spy_runner(function(...) list(pip = c(1, NA_real_), credible_sets = list()))
  out <- runner(sim, "spy", verbose = FALSE)
  expect_equal(out$spy$n_failed, 1L)
  expect_match(out$spy$results[[1]]$error, "invalid or misaligned PIPs")
})

test_that("observed-data and oracle wrappers run on an annotated synthetic scenario", {
  skip_if_not_installed("susieR")
  skip_if_not_installed("nnet")
  g <- simulation_validity_genotypes(n = 150L, p = 20L, regions = 3L)
  for (i in seq_along(g)) {
    g[[i]]$annotations_matrix <- matrix(rep(c(0, 1), 10), 20, 1)
    g[[i]]$annotation_proportions <- .5
  }
  regions <- simulate_phenotypes(g, S = 2, phi = .2, annotations = "binary",
                                 n_annotations = 1, enrichment = 2,
                                 seed = 23, verbose = FALSE)
  sim <- list(genotypes = g, params = list(model = "sparse", n = 150L),
              scenarios = list(list(scenario_id = 1L, S = 2L, phi = .2,
                                    iter = 1L, regions = regions)))
  methods <- c("susie", "polyfun_est", "polyfun_ldsc", "polyfun_oracle", "sbayesrc")
  args <- setNames(rep(list(list(L = 2)), 4), methods[1:4])
  args$sbayesrc <- list(n_iter = 50, burn_in = 20,
                        scenario_pilot_iter = 40, scenario_pilot_burn = 15, seed = 3)
  suppressMessages(suppressWarnings(out <- run_methods(sim, methods, args, verbose = FALSE)))
  expect_true(all(vapply(out[methods], function(x) x$n_failed == 0L, logical(1))))
  expect_equal(out$polyfun_oracle$results[[1]]$benchmark_access, "oracle_enrichment")
  expect_equal(out$polyfun_ldsc$results[[1]]$additional$prior_training_unit, "region")
  expect_false(out$sbayesrc$results[[1]]$additional$upstream_equivalent)
})

test_that("SBayesRC surrogate pilots freeze the neutral prior before pooling labels", {
  skip_if_not_installed("nnet")
  g <- simulation_validity_genotypes(p = 6L)
  for (i in seq_along(g)) g[[i]]$annotations_matrix <- matrix(rep(c(0, 1), 3), 6, 1)
  observed <- lapply(g, function(r) list(z = rep(1, r$p)))
  calls <- list()
  sandbox <- new.env(parent = environment(run_sbayesrc_scenario_setup))
  sandbox$run_sbayesrc <- function(...) {
    args <- list(...)
    calls[[length(calls) + 1L]] <<- args
    list(additional = list(component_assignments = c(0L, 0L, 0L, 0L, 0L, 1L)))
  }
  setup <- run_sbayesrc_scenario_setup
  environment(setup) <- sandbox
  setup(g, observed, list(statistic_type = "normal_score", residual_variance = .7))
  expect_equal(length(calls), 2L)
  for (args in calls) {
    expect_identical(args$statistic_type, "normal_score")
    expect_equal(args$residual_variance, .7)
    expect_equal(args$pooled_gamma$alpha, rep(log((.03 / 4) / .97), 4))
    expect_equal(args$pooled_gamma$gamma, matrix(0, 1, 4))
  }
})

test_that("orthogonal surrogate Gibbs PIPs agree with exact mixture Bayes factors", {
  z <- c(1, 5)
  n <- 100
  sigma2 <- c(.05, .005)
  beta_hat <- z / sqrt(n - 2 + z^2)
  precision <- n - 1
  log_bf <- outer(beta_hat^2, .5 * precision^2 * sigma2 / (1 + precision * sigma2))
  log_bf <- sweep(log_bf, 2, -.5 * log(1 + precision * sigma2), "+")
  marginal_slab <- rowSums(exp(log_bf)) * (.03 / length(sigma2))
  exact_pip <- marginal_slab / (.97 + marginal_slab)
  fit <- run_sbayesrc(z, diag(2), n, sigma2_scale = sigma2,
                       n_iter = 12000L, burn_in = 1000L, seed = 418)
  expect_equal(fit$pip, exact_pip, tolerance = .015)
})

test_that("fixed-mixture SBayesRC runs when nnet is unavailable and learning fails explicitly", {
  sandbox <- new.env(parent = environment(run_sbayesrc))
  sandbox$requireNamespace <- function(package, ...) {
    if (identical(package, "nnet")) FALSE else base::requireNamespace(package, ...)
  }
  full_setup <- setup_sbayesrc
  environment(full_setup) <- sandbox
  sandbox$setup_sbayesrc <- full_setup
  runner <- run_sbayesrc
  environment(runner) <- sandbox
  expect_error(full_setup(), "nnet")

  neutral <- runner(c(1, 3), diag(2), 80, n_iter = 100, burn_in = 40, seed = 6)
  expect_true(all(is.finite(neutral$pip)))
  expect_equal(neutral$additional$prior_source, "no_annotations")
  annotation <- matrix(c(0, 1), 2, 1)
  fixed <- list(alpha = rep(log((.03 / 4) / .97), 4), gamma = matrix(0, 1, 4))
  annotated <- runner(c(1, 3), diag(2), 80, annotations = annotation,
                       pooled_gamma = fixed, n_iter = 100, burn_in = 40, seed = 6)
  expect_equal(annotated$pip, neutral$pip)
  expect_equal(annotated$additional$prior_source, "pooled_scenario_gamma")
  expect_error(runner(c(1, 3), diag(2), 80, annotations = annotation,
                       n_iter = 100, burn_in = 40), "nnet")

  setup <- run_sbayesrc_scenario_setup
  environment(setup) <- sandbox
  g <- simulation_validity_genotypes(p = 2L)
  observed <- lapply(g, function(r) list(z = c(1, 3)))
  expect_equal(setup(g, observed, list()), list())
  sandbox$run_sbayesrc <- runner
  sandbox$run_sbayesrc_scenario_setup <- setup
  region_runner <- run_sbayesrc_region
  environment(region_runner) <- sandbox
  sandbox$run_sbayesrc_region <- region_runner
  dispatch <- run_methods
  environment(dispatch) <- sandbox
  sim <- list(genotypes = g, params = list(model = "sparse", n = 80L),
              scenarios = list(list(scenario_id = 1L, S = 1L, phi = .2,
                                    iter = 1L, regions = observed)))
  controls <- list(n_iter = 30, burn_in = 10, seed = 6)
  neutral_dispatch <- dispatch(sim, "sbayesrc", list(sbayesrc = controls), verbose = FALSE)
  expect_equal(neutral_dispatch$sbayesrc$n_failed, 0L)
  for (i in seq_along(g)) g[[i]]$annotations_matrix <- annotation
  expect_equal(setup(g, observed, list(pooled_gamma = fixed)), list(pooled_gamma = fixed))
  expect_error(setup(g, observed, list()), "nnet")
  sim$genotypes <- g
  supplied_dispatch <- dispatch(sim, "sbayesrc", list(sbayesrc = c(controls, list(pooled_gamma = fixed))),
                                  verbose = FALSE)
  expect_equal(supplied_dispatch$sbayesrc$n_failed, 0L)
  learned_dispatch <- dispatch(sim, "sbayesrc", list(sbayesrc = controls), verbose = FALSE)
  expect_equal(learned_dispatch$sbayesrc$n_failed, 2L)
  expect_match(learned_dispatch$sbayesrc$results[[1]]$error, "Scenario setup failed:.*nnet")
})
