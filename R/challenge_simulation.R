# A dependency-light, conditional-Gaussian challenge. This is deliberately a
# second generator family, not a replacement for haplotype-based simulations.

.challenge_scalar <- function(x, name, lower, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower ||
      (integer && x != floor(x))) stop(name, " has an invalid value", call. = FALSE)
}

.challenge_rng <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
          else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  force(code)
}

.challenge_mvn <- function(R, draws = 1L) {
  ee <- eigen(R, symmetric = TRUE)
  if (min(ee$values) < -1e-8) stop("Covariance is not positive semidefinite")
  # Numerical eigensolvers can turn exact zero eigenvalues into tiny positives.
  # Remove only roundoff-scale eigenvalues, preserving the singular experiment.
  ee$values[ee$values < max(1, max(ee$values)) * 1e-12] <- 0
  B <- sweep(ee$vectors, 2L, sqrt(pmax(ee$values, 0)), "*")
  matrix(stats::rnorm(nrow(R) * draws), draws, nrow(R)) %*% t(B)
}

.challenge_ar <- function(p, rho) rho^abs(outer(seq_len(p), seq_len(p), "-"))

#' Built-in scientific challenge scenarios
#'
#' These names describe assessor-side mechanisms. They are never sent to
#' predictors because, for example, a public null label would reveal truth.
#' @return Data frame of generator mechanisms and assessed questions.
#' @export
challenge_scenarios <- function() {
  data.frame(
    scenario = c("matched_prior", "null", "single_sparse", "multi_sparse", "many_sparse", "opposing_effects",
                 "sparse_inf", "ld_reference_small", "ld_reference_shift",
                 "annotation_informative", "annotation_irrelevant",
                 "annotation_reversal", "missing_causal", "perfect_ld",
                 "allele_flip", "allele_recoding", "selected_null",
                 "selected_signal", "meta_missingness", "meta_known_overlap", "multi_ancestry_shared",
                 "multi_ancestry_heterogeneous"),
    question = c("Calibration under the exact baseline generative prior",
                 "Does the method avoid confident discoveries under no signal?",
                 "Single-signal localisation", "Multiple sparse signals",
                 "Architecture outside the reference model's two-effect support",
                 "Marginal signal cancellation between correlated opposite effects",
                 "Sparse targets in a nonzero polygenic background",
                 "Finite independent reference sampling", "Population LD shift",
                 "Useful annotations without oracle parameters",
                 "Fitting-only irrelevant annotation controls",
                 "Negative transfer under annotation reversal",
                 "Unobserved causal variants and attainable recall",
                 "Unidentifiable variant identity under identical LD",
                 "Inconsistent allele orientation", "Consistent allele recoding",
                 "Nulls conditional on discovery selection",
                 "Signals conditional on discovery selection",
                 "Variant-wise sample size and overlap in meta-analysis",
                 "Covariance-aware inference when study overlap is available",
                 "Complementary LD in independent cohorts with shared effects",
                 "Shared-effect misspecification across independent cohorts"),
    stringsAsFactors = FALSE
  )
}

.challenge_one <- function(scenario, p, n, signal_z, prior_sd, reference_n,
                           split, serial, seed) {
  .challenge_rng(seed, {
    R <- .challenge_ar(p, stats::runif(1, 0.35, 0.85))
    if (scenario == "opposing_effects") R <- .challenge_ar(p, .97)
    if (scenario == "perfect_ld") {
      # Both columns are the same random variable, not merely high r^2.
      R[, 2] <- R[, 1]; R[2, ] <- R[1, ]; R[2, 2] <- 1
    }
    annotations <- cbind(activity = stats::rnorm(p), decoy = stats::rnorm(p))
    k <- if (scenario %in% c("null", "selected_null")) 0L else
      if (scenario %in% c("multi_sparse", "sparse_inf")) 2L else 1L
    if (scenario == "matched_prior") k <- sample(0:2, 1, prob = c(.5, .3, .2))
    if (scenario == "many_sparse") k <- 5L
    weights <- rep(1, p)
    if (scenario == "annotation_informative")
      weights <- exp(annotations[, 1] - max(annotations[, 1]))
    if (scenario == "annotation_reversal")
      weights <- exp(-annotations[, 1] - max(-annotations[, 1]))
    causal <- if (k) sample.int(p, k, prob = weights) else integer(0)
    if (scenario == "perfect_ld") causal <- sample(1:2, 1)
    if (scenario == "opposing_effects") causal <- 1:2
    beta_sparse <- numeric(p)
    if (length(causal)) beta_sparse[causal] <- if (scenario == "matched_prior")
      stats::rnorm(length(causal), sd = prior_sd) else
      sample(c(-1, 1), length(causal), replace = TRUE) * signal_z / sqrt(n)
    if (scenario == "selected_signal") beta_sparse[causal] <- sign(beta_sparse[causal]) * 3 / sqrt(n)
    if (scenario == "opposing_effects") beta_sparse[causal] <- c(1, -1) * signal_z / sqrt(n)
    beta_background <- numeric(p)
    if (scenario == "sparse_inf")
      beta_background <- stats::rnorm(p, sd = 1.25 / sqrt(n))
    beta_total <- beta_sparse + beta_background
    mean_z <- as.numeric(sqrt(n) * R %*% beta_total)
    selected <- scenario %in% c("selected_null", "selected_signal")
    selection_cutoff <- if (scenario == "selected_null") 3.5 else 5
    selection_draws <- 0L
    repeat {
      z <- mean_z + as.numeric(.challenge_mvn(R))
      selection_draws <- selection_draws + 1L
      if (!selected || max(abs(z)) >= selection_cutoff) break
      if (selection_draws >= 100000L) stop("Selection rejection cap reached")
    }
    R_fit <- R
    n_per_variant <- NULL
    sample_overlap <- NULL
    noise_covariance <- R
    if (scenario %in% c("meta_missingness", "meta_known_overlap")) {
      # Two independent equal-size study scores; alternating variants are
      # observed in one versus both studies. Covariance is R * overlap.
      included <- rep(c(TRUE, FALSE), length.out = p)
      z1 <- as.numeric(sqrt(n / 2) * R %*% beta_total) + as.numeric(.challenge_mvn(R))
      z2 <- as.numeric(sqrt(n / 2) * R %*% beta_total) + as.numeric(.challenge_mvn(R))
      z <- ifelse(included, (z1 + z2) / sqrt(2), z1)
      n_per_variant <- ifelse(included, n, n / 2)
      overlap <- outer(as.numeric(included), as.numeric(included), "*") + 1
      noise_covariance <- R * overlap / sqrt(outer(1 + included, 1 + included))
      # This is the normalized Gram of the two study-inclusion vectors.
      # Recovering it by dividing by LD fails at zero or subnormal LD entries.
      if (scenario == "meta_known_overlap")
        sample_overlap <- overlap / sqrt(outer(1 + included, 1 + included))
    }
    if (scenario == "ld_reference_small") {
      ref <- .challenge_mvn(R, reference_n)
      R_fit <- stats::cor(ref)  # may be singular: a valid benchmark condition
    }
    if (scenario == "ld_reference_shift") R_fit <- .challenge_ar(p, 0.2)
    all_ids <- paste0("sim:", split, ":", serial, ":v", seq_len(p))
    observed <- seq_len(p)
    if (scenario == "missing_causal") observed <- setdiff(observed, causal[1])
    # Canonical effect orientation changes require simultaneous LD row/column
    # sign changes. The error arm deliberately changes z only.
    orientation <- rep(1, p)
    if (scenario %in% c("allele_flip", "allele_recoding")) {
      orientation[if (length(causal)) causal[1] else 1L] <- -1
      z <- z * orientation
      if (scenario == "allele_recoding") {
        R_fit <- R_fit * outer(orientation, orientation)
        R <- R * outer(orientation, orientation)
        noise_covariance <- noise_covariance * outer(orientation, orientation)
        beta_sparse <- beta_sparse * orientation
        beta_background <- beta_background * orientation
      }
    }
    cohorts <- NULL
    if (scenario %in% c("multi_ancestry_shared", "multi_ancestry_heterogeneous")) {
      R2 <- .challenge_ar(p, 0.15)
      beta2 <- beta_total
      if (scenario == "multi_ancestry_heterogeneous") beta2[causal] <- -beta2[causal]
      z2 <- as.numeric(sqrt(n) * R2 %*% beta2) + as.numeric(.challenge_mvn(R2))
      cohorts <- list(list(z = z, LD = R_fit, n = n, cohort_id = "A"),
                      list(z = z2, LD = R2, n = n, cohort_id = "B"))
    }
    task_id <- sprintf("%s_%05d", split, serial)
    summary_n <- if (is.null(n_per_variant)) rep(n, length(observed)) else n_per_variant[observed]
    input <- list(task_id = task_id, variant_id = all_ids[observed],
                  build = "synthetic-v1", z = z[observed], n = n,
                  beta_hat = z[observed] / sqrt(summary_n),
                  se = 1 / sqrt(summary_n),
                  LD = R_fit[observed, observed, drop = FALSE],
                  annotations = annotations[observed, , drop = FALSE],
                  cohorts = cohorts, n_per_variant = n_per_variant, sample_overlap = sample_overlap)
    truth <- list(task_id = task_id, scenario = scenario, split = split,
                  cluster_id = task_id, causal_variant_id = all_ids[causal],
                  all_variant_id = all_ids, observed_variant_id = all_ids[observed],
                  beta_sparse = stats::setNames(beta_sparse, all_ids),
                  beta_background = stats::setNames(beta_background, all_ids),
                  LD_true = R[observed, observed, drop = FALSE],
                  score_covariance_true = noise_covariance[observed, observed, drop = FALSE],
                  null = !length(causal), n_missing_causal = sum(!causal %in% observed),
                  selection_draws = selection_draws,
                  selection_rule = if (selected) paste0("max(abs(z)) >= ", selection_cutoff) else "none",
                  seed = seed)
    list(input = input, truth = truth)
  })
}

#' Generate a challenge panel with separate predictor inputs and assessor truth
#'
#' The summary-statistic experiment is z | beta ~ N(sqrt(n) R beta, R),
#' conditional on fixed LD. Effects are in residual-standard-deviation units.
#' No claim of realistic genotypes or exact marginal-OLS statistics is made.
#' @param scenarios Character vector from challenge_scenarios().
#' @param n_replicates Independent loci per mechanism.
#' @param p Number of variants, at least 12.
#' @param n Sample-size parameter in the score-statistic likelihood.
#' @param seed Assessor seed; do not release a final-panel seed before closing.
#' @param split One of development or final. Final generation does not confer
#'   blindness; independent custody and unavailable truth are also required.
#' @param signal_z Noncentral signal size for stress scenarios.
#' @param prior_sd Normal effect prior SD in the matched-prior scenario.
#' @param reference_n Independent reference observations in the small-LD arm.
#' @return List with inputs, truth and private generator configuration.
#' @export
generate_challenge <- function(scenarios = challenge_scenarios()$scenario,
                               n_replicates = 10L, p = 40L, n = 10000L,
                               seed = 20261010L, split = c("development", "final"),
                               signal_z = 6, prior_sd = 0.05, reference_n = 30L) {
  split <- match.arg(split)
  .challenge_scalar(n_replicates, "n_replicates", 1, TRUE)
  .challenge_scalar(p, "p", 12, TRUE); .challenge_scalar(n, "n", 3, TRUE)
  .challenge_scalar(seed, "seed", 0, TRUE)
  .challenge_scalar(signal_z, "signal_z", .Machine$double.eps)
  .challenge_scalar(prior_sd, "prior_sd", .Machine$double.eps)
  .challenge_scalar(reference_n, "reference_n", 3, TRUE)
  if (!length(scenarios) || anyDuplicated(scenarios) ||
      any(!scenarios %in% challenge_scenarios()$scenario)) stop("Invalid scenarios")
  count <- length(scenarios) * n_replicates
  if (seed + count + 100000000 > .Machine$integer.max) stop("Seed is too large")
  # Randomise order so task identifiers do not encode a causal mechanism.
  order_scenarios <- .challenge_rng(seed, sample(rep(scenarios, each = n_replicates)))
  offset <- if (split == "final") 100000000L else 0L
  cases <- lapply(seq_len(count), function(i) .challenge_one(
    order_scenarios[i], p, n, signal_z, prior_sd, reference_n, split, i, seed + offset + i))
  inputs <- lapply(cases, `[[`, "input"); truth <- lapply(cases, `[[`, "truth")
  names(inputs) <- names(truth) <- vapply(inputs, `[[`, "", "task_id")
  list(inputs = inputs, truth = truth, config = list(
    specification = "conditional-gaussian-v1", scenarios = scenarios,
    n_replicates = n_replicates, p = p, n = n, seed = seed, split = split,
    signal_z = signal_z, prior_sd = prior_sd, reference_n = reference_n))
}

#' Exact small-locus Bayesian reference with an explicit null model
#'
#' Enumerates zero, one and two-effect configurations with Gaussian effects.
#' This is exact for its declared model, not an oracle or universal ceiling.
#' Native output is a credible collection of configurations, not a set made by
#' accumulating marginal PIPs. Independent multi-cohort evidence is combined
#' under the declared shared-effect model.
#' @param input A predictor input from generate_challenge().
#' @param prior_sd Effect prior standard deviation.
#' @param prior_k Probabilities for 0, 1, 2 causal variants; sums to one.
#' @param coverage Posterior mass for the credible configuration collection.
#' @param annotation_coefficient Prespecified activity coefficient (zero default).
#' @return Prediction compatible with evaluate_challenge().
#' @export
challenge_exact <- function(input, prior_sd = 0.05, prior_k = c(.5, .3, .2),
                            coverage = .95, annotation_coefficient = 0) {
  validate_challenge_input(input)
  if (!is.null(input$n_per_variant) && is.null(input$sample_overlap))
    stop("Variant-wise n requires sample-overlap metadata for this adapter")
  .challenge_scalar(prior_sd, "prior_sd", .Machine$double.eps)
  if (length(prior_k) != 3L || any(!is.finite(prior_k)) || any(prior_k < 0) ||
      abs(sum(prior_k) - 1) > 1e-10) stop("prior_k must be probabilities for 0:2")
  if (length(coverage) != 1L || !is.finite(coverage) || coverage <= 0 || coverage >= 1)
    stop("coverage must lie in (0,1)")
  if (length(annotation_coefficient) != 1L || !is.finite(annotation_coefficient))
    stop("Invalid annotation coefficient")
  p <- length(input$z)
  if (p > 300L) stop("Exact enumeration is limited to 300 variants; use a scalable adapter")
  sets <- c(list(integer(0)), lapply(seq_len(p), function(j) j),
            utils::combn(seq_len(p), 2L, simplify = FALSE))
  k <- lengths(sets)
  weights <- if (annotation_coefficient == 0) rep(0, p) else
    annotation_coefficient * input$annotations[, "activity"]
  log_subset <- vapply(sets, function(s) sum(weights[s]), 0)
  log_prior <- rep(-Inf, length(sets))
  for (kk in 0:2) {
    ii <- which(k == kk)
    aa <- log_subset[ii]; norm <- max(aa) + log(sum(exp(aa - max(aa))))
    if (prior_k[kk + 1] > 0) log_prior[ii] <- log(prior_k[kk + 1]) + aa - norm
  }
  cohorts <- if (is.null(input$cohorts)) list(list(z = input$z, LD = input$LD, n = input$n))
             else input$cohorts
  W <- prior_sd^2
  information_full <- matrix(0, p, p); score_full <- numeric(p)
  if (!is.null(input$n_per_variant)) {
    Sigma <- input$LD * input$sample_overlap
    ee <- eigen(Sigma, symmetric = TRUE)
    inverse_values <- ifelse(ee$values > max(ee$values) * 1e-10, 1 / ee$values, 0)
    inverse <- sweep(ee$vectors, 2, inverse_values, "*") %*% t(ee$vectors)
    H <- sweep(input$LD, 1, sqrt(input$n_per_variant), "*")
    information_full <- crossprod(H, inverse %*% H)
    score_full <- as.numeric(crossprod(H, inverse %*% input$z))
  } else for (cc in cohorts) {
    information_full <- information_full + cc$n * cc$LD
    score_full <- score_full + sqrt(cc$n) * cc$z
  }
  log_bf <- vapply(sets, function(s) {
    if (!length(s)) return(0)
    information <- information_full[s, s, drop = FALSE]; score <- score_full[s]
    A <- diag(length(s)) + W * information
    as.numeric(-.5 * determinant(A, logarithm = TRUE)$modulus +
               .5 * W * crossprod(score, solve(A, score)))
  }, 0)
  lp <- log_prior + log_bf
  posterior <- exp(lp - max(lp)); posterior <- posterior / sum(posterior)
  pip <- numeric(p)
  for (i in seq_along(sets)) pip[sets[[i]]] <- pip[sets[[i]]] + posterior[i]
  # A sum of probabilities can exceed one by floating-point roundoff alone.
  pip <- pmin(1, pmax(0, pip))
  oo <- order(posterior, decreasing = TRUE)
  last <- which(cumsum(posterior[oo]) >= coverage)[1]
  config_ids <- lapply(sets[oo[seq_len(last)]], function(s) input$variant_id[s])
  list(task_id = input$task_id, variant_id = input$variant_id, pip = pip,
       status = "ok", credible_sets = list(),
       credible_configurations = config_ids, configuration_coverage = coverage,
       null_probability = posterior[1],
       implementation = "exact-0:2-shared-gaussian", oracle = FALSE)
}

#' Force-one-cause rank baseline for exposing null overconfidence
#' @param input A challenge input.
#' @return Prediction under a deliberately restrictive one-effect model.
#' @export
challenge_forced_one <- function(input) {
  challenge_exact(input, prior_k = c(0, 1, 0))
}
