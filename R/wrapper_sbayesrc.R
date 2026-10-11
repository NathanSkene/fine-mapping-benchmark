# =============================================================================
# wrapper_sbayesrc.R
#
# SBayesRC-style Bayesian regression with annotation-informed mixture prior,
# implemented in R for the fine-mapping benchmark's per-locus regime.
#
# This is a local experimental surrogate, not the upstream SBayesRC package.
# All outputs identify implementation = sbayesrc_style_surrogate and
# upstream_equivalent = FALSE. A per-region pilot, one pooled multinomial
# regression, and final per-region chains do not reproduce joint genome-wide
# SBayesRC posterior inference. See docs/methods.md for scientific limits.
#
# Model
# -----
# For each region i, effect sizes beta_{i,j} follow a spike + K-normal-slabs
# mixture:
#
#     beta_{i,j} | comp_{i,j} = k ~ Normal(0, sigma2_k),  k = 1, ..., K
#     beta_{i,j} | comp_{i,j} = 0 = 0                     (spike)
#
# where sigma2_1 > sigma2_2 > ... > sigma2_K is a fixed grid of variance
# scales (default: c(0.05, 0.005, 5e-4, 5e-5), corresponding roughly to
# large/medium/small/tiny per-SNP heritability).
#
# The per-SNP mixture proportions are annotation-modulated via a multinomial
# logit with the spike (k=0) as the reference class:
#
#     log(pi_{i,j,k} / pi_{i,j,0}) = alpha_k + A_{i,j}^T gamma_k,
#     pi_{i,j,k} = softmax_k(alpha_k + A_{i,j}^T gamma_k)
#
# alpha and gamma are SHARED across regions - this is the "pool annotations
# across regions to estimate the shared annotation-to-prior mapping" step
# from the plan (Sec 0.1 of the auto-research doc).
#
# Summary-stat likelihood
# -----------------------
# For standardized X and y, OLS t statistics imply the exact marginal slope
# beta_hat_j = t_j / sqrt(n - 2 + t_j^2). The likelihood uses fixed residual
# variance (default 1) with precision (n - 1) / residual_variance. Genuine
# approximate normal scores may use z/sqrt(n) via statistic_type=normal_score.
# The Gibbs conditional precision is based on X'X, not the marginal variance
# of beta_hat. A Gibbs sweep updates sampled mixture class and effect size.
# A scenario hook runs pilots with a fixed sparse annotation-neutral prior,
# then pools sampled component assignments to fit a fixed final prior.
# Single-region calls can still use the documented unstable local refit.
# Without annotations the mixture proportions stay at their initialization.
# =============================================================================


# =============================================================================
# setup_sbayesrc()
# =============================================================================

#' Check full annotation-learning support for the SBayesRC surrogate
#'
#' Checks \code{nnet}, required for local or pooled annotation regression.
#' Fixed-mixture sampling without annotations, or with supplied
#' \code{pooled_gamma}, uses base R and does not call this full-support check.
#'
#' @return Invisible TRUE.
#' @export
setup_sbayesrc <- function() {
  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("The `nnet` package is required for sbayesrc (used for the ",
         "shared multinomial annotation regression). It ships with R by ",
         "default; install with install.packages('nnet') if missing.",
         call. = FALSE)
  }
  invisible(TRUE)
}


# =============================================================================
# run_sbayesrc(): single-region interface
# =============================================================================

#' Run SBayesRC-style Bayesian regression on a single region
#'
#' Fits a K-component mixture-of-normals prior (with the point mass at zero
#' as the spike) using a Gibbs sampler on the summary-stat likelihood
#' r_j | beta_j ~ Normal(beta_j, 1/n). When \code{pooled_gamma} is supplied
#' (as it is under \code{run_methods()} via
#' \code{run_sbayesrc_scenario_setup()}) the per-SNP mixture weights use the
#' cross-region annotation regression fit; otherwise they are refit on this
#' region alone.
#'
#' When \code{annotations} is NULL, the prior collapses to region-common
#' mixture proportions (fixed at their initial values), and the method runs.
#' This fixed-prior path and supplied \code{pooled_gamma} need only base R;
#' learning local or scenario-pooled annotation coefficients requires nnet.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD correlation matrix (p x p).
#' @param n Integer. Sample size.
#' @param annotations Matrix or NULL. Functional annotations (p x m).
#' @param sigma2_scale Numeric vector of length K. Per-component prior
#'   variance grid. Default \code{c(0.05, 0.005, 5e-4, 5e-5)}.
#' @param n_iter Integer. Total Gibbs iterations. Default 300.
#' @param burn_in Integer. Iterations to discard before averaging. Default 150.
#' @param gamma_update_every Integer. Refit the annotation regression
#'   (alpha, gamma) every this many Gibbs iterations. Default 10.
#' @param pooled_gamma List with \code{alpha}, \code{gamma} matrices from
#'   \code{run_sbayesrc_scenario_setup()}, or NULL. When supplied, no
#'   region-local annotation refit is performed - the shared coefficients
#'   are used as-is.
#' @param seed Integer or NULL. Sampler seed.
#' @param variant_ids Character or NULL. Passed through unused.
#' @param statistic_type Character. \code{marginal_ols_t} (default) uses the
#'   exact t-to-standardized-slope transform; \code{normal_score} uses z/sqrt(n).
#' @param residual_variance Positive scalar. Fixed residual variance on the
#'   standardized phenotype scale. Default 1 is an approximation.
#' @param ... Ignored (wrapper compatibility).
#'
#' @return List with pip, credible_sets, method = "sbayesrc", params,
#'   runtime_seconds, additional (posterior_mean_beta, prior_source, etc.).
#' @export
run_sbayesrc <- function(z, LD, n,
                         annotations = NULL,
                         sigma2_scale = c(0.05, 0.005, 5e-4, 5e-5),
                         n_iter = 300L,
                         burn_in = 150L,
                         gamma_update_every = 10L,
                         pooled_gamma = NULL,
                         seed = NULL,
                         variant_ids = NULL,
                         statistic_type = c("marginal_ols_t", "normal_score"),
                         residual_variance = 1,
                         ...) {
  if (!is.null(annotations) && is.null(pooled_gamma)) setup_sbayesrc()
  t0 <- Sys.time()
  if (!is.null(seed)) set.seed(seed)

  p <- length(z)
  stopifnot("LD must be p x p" = is.matrix(LD) && all(dim(LD) == c(p, p)))
  K <- length(sigma2_scale)
  statistic_type <- match.arg(statistic_type)
  stopifnot(
    "z and LD must be finite and nonempty" = p >= 1L && all(is.finite(z)) && all(is.finite(LD)),
    "LD must be a symmetric correlation matrix" =
      max(abs(LD - t(LD))) < 1e-8 && max(abs(diag(LD) - 1)) < 1e-8,
    "n must be an integer >= 3" = length(n) == 1L && is.finite(n) && n >= 3 && n == floor(n),
    "slab variances must be finite and positive" = K >= 1L && all(is.finite(sigma2_scale)) && all(sigma2_scale > 0),
    "burn_in must leave posterior samples" =
      length(n_iter) == 1L && length(burn_in) == 1L && is.finite(n_iter) &&
      is.finite(burn_in) && n_iter == floor(n_iter) && burn_in == floor(burn_in) &&
      burn_in >= 0 && n_iter > burn_in,
    "gamma_update_every must be a positive integer" =
      length(gamma_update_every) == 1L && is.finite(gamma_update_every) &&
      gamma_update_every >= 1 && gamma_update_every == floor(gamma_update_every),
    "residual_variance must be finite and positive" =
      length(residual_variance) == 1L && is.finite(residual_variance) && residual_variance > 0
  )

  prior_source <- if (is.null(annotations)) {
    "no_annotations"
  } else if (!is.null(pooled_gamma)) {
    "pooled_scenario_gamma"
  } else {
    "single_region_gamma"
  }

  # --- Marginal beta_hat + per-SNP residual state ---------------------------
  # Marginal OLS t statistics imply exactly cor(X, y) = t/sqrt(df+t^2)
  # for a standardized phenotype. z/sqrt(n) is only a small-effect approximation.
  beta_hat <- if (statistic_type == "marginal_ols_t") z / sqrt(n - 2 + z^2) else z / sqrt(n)
  likelihood_precision <- (if (statistic_type == "marginal_ols_t") n - 1 else n) / residual_variance
  beta     <- rep(0, p)                # current effect estimates
  comp     <- rep(0L, p)               # current mixture assignments (0..K)
  # Running R %*% beta so per-SNP residuals are cheap: init 0 (beta = 0)
  R_beta   <- rep(0, p)

  # --- Annotation setup ------------------------------------------------------
  has_annot <- !is.null(annotations)
  if (has_annot) {
    A <- as.matrix(annotations)
    m_ann <- ncol(A)
    if (nrow(A) != p) {
      stop("annotations must have p rows to match z", call. = FALSE)
    }
    if (!is.null(pooled_gamma)) {
      alpha <- as.numeric(pooled_gamma$alpha)
      gamma <- as.matrix(pooled_gamma$gamma)
      stopifnot(length(alpha) == K, nrow(gamma) == m_ann, ncol(gamma) == K)
    } else {
      # Sparse initial prior: ~97% mass on the spike, ~3% distributed
      # roughly evenly across the K slabs. alpha_k = log((0.03/K) / 0.97).
      # This matches the sparse-model regime of the benchmark (S causal
      # variants out of p, S << p) and gives the sampler a reasonable
      # starting point before the multinomial refit learns the true
      # sparsity from the data.
      alpha <- rep(log((0.03 / K) / 0.97), K)
      gamma <- matrix(0, nrow = m_ann, ncol = K)  # per-annotation slopes
    }
  } else {
    A <- NULL; m_ann <- 0L
    alpha <- rep(log((0.03 / K) / 0.97), K)   # same sparse initial prior
    gamma <- matrix(0, nrow = 0L, ncol = K)
  }

  # Silent-degradation guard: with annotations but no pooled gamma we fall back
  # to the unstable in-loop refit (see the WARNING in the Gibbs loop below).
  if (has_annot && is.null(pooled_gamma)) .sbayesrc_warn_no_pooled()

  # Prior mixture proportions per SNP: pi[j, ] over (0, 1, ..., K)
  pi_mat <- .sbayesrc_priors_from_gamma(A, alpha, gamma, K, p = p)

  # --- Storage for post-burn PIP estimation ---------------------------------
  pip_running   <- rep(0, p)
  beta_running  <- rep(0, p)
  n_kept        <- 0L

  # --- Main Gibbs loop ------------------------------------------------------
  for (it in seq_len(n_iter)) {
    # Sweep over SNPs
    swp <- .sbayesrc_gibbs_sweep(
      beta_hat = beta_hat, LD = LD, n = likelihood_precision,
      beta = beta, comp = comp, R_beta = R_beta,
      pi_mat = pi_mat, sigma2_scale = sigma2_scale
    )
    beta   <- swp$beta
    comp   <- swp$comp
    R_beta <- swp$R_beta

    # Update annotation regression periodically (only when pooled_gamma is NULL
    # and we have annotations). Under run_methods()/scenario_setup the pooled
    # coefficients are supplied and this branch is skipped.
    #
    # WARNING (verified empirically 2026-07-29): this in-loop refit path is
    # UNSTABLE. On a benchmark-like locus (p=300, n=1000, phi=0.05, S=3) it
    # inflates total posterior mass ~15x relative to the pooled-gamma path
    # (mass ratio ~55 vs ~3.7) EVEN WITH THE TRUE LD - so the inflation is
    # caused by the annotation refit, not by LD mis-specification. It is only
    # reached when run_sbayesrc_scenario_setup() returned nothing (fewer than 2
    # regions, absent/ragged annotations, or a failed pilot), which previously
    # degraded results silently. It now warns; see .sbayesrc_warn_no_pooled().
    if (has_annot && is.null(pooled_gamma) &&
        (it %% gamma_update_every == 0L)) {
      fit <- tryCatch(.sbayesrc_fit_gamma(comp, A, K), error = function(e) NULL)
      if (!is.null(fit)) {
        alpha  <- fit$alpha
        gamma  <- fit$gamma
        pi_mat <- .sbayesrc_priors_from_gamma(A, alpha, gamma, K, p = p)
      }
    }

    # Accumulate post-burn samples
    if (it > burn_in) {
      pip_running  <- pip_running  + (comp > 0L)
      beta_running <- beta_running + beta
      n_kept       <- n_kept + 1L
    }
  }

  pip  <- if (n_kept > 0L) pip_running  / n_kept else rep(0, p)
  bhat <- if (n_kept > 0L) beta_running / n_kept else rep(0, p)

  # A set containing 95% of SUMMED marginal PIPs is a posterior mass summary,
  # not a 95% credible set for an effect or the complete causal configuration.
  ord <- order(pip, decreasing = TRUE)
  mass_set <- integer(0)
  if (sum(pip) > 0) {
    cutoff <- which(cumsum(pip[ord]) >= 0.95 * sum(pip))[1L]
    mass_set <- as.integer(ord[seq_len(cutoff)])
  }
  cs_list <- list()

  list(
    pip             = pip,
    credible_sets   = cs_list,
    method          = "sbayesrc",
    input_type      = "summary",
    params          = list(
      sigma2_scale = sigma2_scale, n_iter = n_iter, burn_in = burn_in,
      gamma_update_every = gamma_update_every, K = K,
      statistic_type = statistic_type, residual_variance = residual_variance,
      implementation = "sbayesrc_style_surrogate", upstream_equivalent = FALSE
    ),
    runtime_seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    additional      = list(
      posterior_mean_beta = bhat,
      component_assignments = comp,
      pip_event = "any_nonzero_mixture_component",
      pip_mass_set = mass_set,
      credible_set_semantics = "not_provided",
      implementation = "sbayesrc_style_surrogate",
      upstream_equivalent = FALSE,
      prior_source        = prior_source,
      alpha               = alpha,
      gamma               = gamma,
      variant_ids         = variant_ids
    )
  )
}


# =============================================================================
# run_sbayesrc_region(): adapter for run_methods()
# =============================================================================

#' Adapter for run_methods() - forwards to run_sbayesrc
#'
#' @param region_geno One element of \code{simulation$genotypes}.
#' @param region_pheno One element of a scenario's \code{regions}.
#' @param pooled_gamma List or NULL. Shared (alpha, gamma) from
#'   \code{run_sbayesrc_scenario_setup()}. When supplied, no region-local
#'   annotation refit happens.
#' @param ... Forwarded to \code{run_sbayesrc()}.
#' @param statistic_type Optional summary type inherited from the region;
#'   defaults to \code{marginal_ols_t} when no type is declared.
#'
#' @return Output of run_sbayesrc.
#' @export
run_sbayesrc_region <- function(region_geno, region_pheno,
                                pooled_gamma = NULL, ..., statistic_type = NULL) {
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  run_sbayesrc(
    z            = region_pheno$z,
    LD           = region_geno$LD,
    n            = region_geno$n,
    annotations  = A,
    pooled_gamma = pooled_gamma,
    statistic_type = statistic_type %||% region_pheno$statistic_type %||% "marginal_ols_t",
    variant_ids  = region_geno$variant_ids,
    ...
  )
}


# =============================================================================
# Scenario-level setup: pool across regions, estimate shared (alpha, gamma)
# =============================================================================

#' Fit a pooled annotation regression for sbayesrc across the scenario
#'
#' Runs a short-burn Gibbs on every region under a fixed sparse prior,
#' pools the resulting per-SNP component assignments across all
#' regions, and does one final multinomial regression on the pooled data to
#' obtain the shared \code{(alpha, gamma)} - i.e., the annotation-to-prior
#' mapping learned from the whole scenario.
#'
#' Downstream calls (per region) use those shared coefficients as
#' \code{pooled_gamma}, so each region's Gibbs skips the per-region
#' annotation refit and its priors are set from the cross-region model. This
#' pools information across regions, but does not implement upstream SBayesRC
#' joint posterior inference.
#'
#' Returns an empty list (no forwarding) when there are fewer than 2 regions,
#' when annotation matrices are missing/inconsistent, or when the pilot Gibbs
#' fails - the per-region wrapper then does its own single-region refit.
#'
#' @param genotypes List of region_geno.
#' @param regions List of region_pheno for this scenario.
#' @param user_args User's method_args for sbayesrc.
#'
#' @return List with a single element \code{pooled_gamma} = list(alpha, gamma),
#'   or an empty list.
#' @export
run_sbayesrc_scenario_setup <- function(genotypes, regions, user_args) {
  if (!is.null(user_args$pooled_gamma)) return(list(pooled_gamma = user_args$pooled_gamma))
  n_regions <- length(regions)
  if (n_regions < 2L) return(list())

  # Resolve tuning knobs (allow user override of sigma2_scale, K, etc.)
  sigma2_scale       <- user_args$sigma2_scale %||% c(0.05, 0.005, 5e-4, 5e-5)
  pilot_iter         <- user_args$scenario_pilot_iter %||% 100L
  pilot_burn         <- user_args$scenario_pilot_burn %||% 40L
  K <- length(sigma2_scale)

  # Validate that every region has an annotation matrix of consistent width
  A_list <- vector("list", n_regions)
  comp_list <- vector("list", n_regions)
  for (i in seq_len(n_regions)) {
    A_i <- genotypes[[i]]$annotations_matrix
    if (is.null(A_i)) A_i <- regions[[i]]$annotations_matrix
    if (is.null(A_i) || !is.matrix(A_i)) return(list())
    A_list[[i]] <- A_i
  }
  ncols <- vapply(A_list, ncol, integer(1))
  if (length(unique(ncols)) != 1L) return(list())
  # Learning the shared coefficients requires nnet. Do not let a missing
  # optional dependency become a silently degenerate pooled prior through
  # the numerical-fit fallback in .sbayesrc_fit_gamma(). Supplied coefficients
  # return above, and unannotated scenarios return without this requirement.
  setup_sbayesrc()

  # A pilot must not fit annotations repeatedly to its own sampled labels:
  # that empirical-Bayes feedback can inflate inclusion and contaminate the
  # pooled fit before it sees another region. Freeze the same sparse neutral
  # prior for all pilots, then estimate annotations once on the pooled draws.
  pilot_prior <- list(alpha = rep(log((0.03 / K) / 0.97), K),
                      gamma = matrix(0, nrow = ncols[1L], ncol = K))
  for (i in seq_len(n_regions)) {
    A_i <- A_list[[i]]
    z_i <- regions[[i]]$z
    if (is.null(z_i) || length(z_i) != nrow(A_i)) return(list())
    LD_i <- genotypes[[i]]$LD
    n_i  <- genotypes[[i]]$n
    if (is.null(LD_i) || is.null(n_i)) return(list())

    pilot <- tryCatch(
      run_sbayesrc(
        z = z_i, LD = LD_i, n = n_i, annotations = A_i,
        sigma2_scale = sigma2_scale, n_iter = pilot_iter, burn_in = pilot_burn,
        statistic_type = user_args$statistic_type %||% regions[[i]]$statistic_type %||% "marginal_ols_t",
        residual_variance = user_args$residual_variance %||% 1,
        pooled_gamma = pilot_prior
      ),
      error = function(e) NULL
    )
    if (is.null(pilot)) stop("SBayesRC surrogate pilot fit failed.", call. = FALSE)
    # Use a real posterior Gibbs assignment. Slab membership is a sampled
    # latent variable, not a rank bin of thresholded posterior-mean effects.
    comp_i <- pilot$additional$component_assignments
    comp_list[[i]] <- comp_i
  }

  # Pool and fit the shared (alpha, gamma)
  A_pooled    <- do.call(rbind, A_list)
  comp_pooled <- unlist(comp_list, use.names = FALSE)
  fit <- tryCatch(
    .sbayesrc_fit_gamma(comp_pooled, A_pooled, K),
    error = function(e) NULL
  )
  if (is.null(fit)) stop("SBayesRC surrogate pooled prior fit failed.", call. = FALSE)

  list(pooled_gamma = list(alpha = fit$alpha, gamma = fit$gamma))
}


# =============================================================================
# Internals
# =============================================================================

# One Gibbs sweep over SNPs in a single region.
# Updates in place (returns new copies of) beta, comp, and R_beta = LD %*% beta
# so per-SNP residuals stay cheap to compute.
.sbayesrc_gibbs_sweep <- function(beta_hat, LD, n, beta, comp, R_beta,
                                  pi_mat, sigma2_scale) {
  p <- length(beta_hat)
  K <- length(sigma2_scale)
  # Precompute per-component posterior variance and log-BF scale factor
  # (which do not depend on r_j).
  v_k    <- 1 / (1 / sigma2_scale + n)           # posterior variance
  # log p(r|k) - log p(r|k=0) has two pieces:
  #   log-sd ratio  = 0.5 * log((1/n) / (sigma2_k + 1/n))
  #                 = -0.5 * log(1 + n * sigma2_k)
  #   quadratic in r = 0.5 * (n - 1/(sigma2_k + 1/n)) * r^2
  #                  = 0.5 * n^2 * sigma2_k / (1 + n * sigma2_k) * r^2
  log_sd_ratio <- -0.5 * log(1 + n * sigma2_scale)
  quad_scale   <- 0.5 * n^2 * sigma2_scale / (1 + n * sigma2_scale)

  for (j in seq_len(p)) {
    # Residual: beta_hat_j - sum_{k != j} R_{j,k} beta_k
    #         = beta_hat_j - (R_beta[j] - R[j,j] * beta[j])
    #         = beta_hat_j - R_beta[j] + beta[j]     (R[j,j] = 1)
    r_j <- beta_hat[j] - R_beta[j] + beta[j]

    # log-likelihood contribution per component (relative to k=0)
    log_lik <- c(0, log_sd_ratio + quad_scale * r_j^2)   # length K+1

    # log-priors from pi_mat (already normalized rows over 0..K)
    log_pri <- log(pmax(pi_mat[j, ], 1e-30))

    log_post <- log_pri + log_lik
    log_post <- log_post - max(log_post)
    w_post   <- exp(log_post); w_post <- w_post / sum(w_post)

    new_k <- sample.int(K + 1L, size = 1L, prob = w_post) - 1L
    if (new_k == 0L) {
      new_beta <- 0
    } else {
      m_k      <- v_k[new_k] * n * r_j
      new_beta <- stats::rnorm(1L, mean = m_k, sd = sqrt(v_k[new_k]))
    }

    # Update R_beta with the change in beta_j
    delta_j <- new_beta - beta[j]
    if (delta_j != 0) R_beta <- R_beta + LD[, j] * delta_j
    beta[j] <- new_beta
    comp[j] <- new_k
  }

  list(beta = beta, comp = comp, R_beta = R_beta)
}


# Compute per-SNP prior mixture proportions pi[j, ] over (0, 1, ..., K)
# from alpha and gamma (with k=0 spike as reference class).
# When annotations are NULL, the row-common softmax over (0, alpha) is broadcast
# to all p rows.
.sbayesrc_priors_from_gamma <- function(A, alpha, gamma, K, p = NULL) {
  if (is.null(A) || nrow(gamma) == 0L) {
    logits <- matrix(alpha, nrow = 1L, ncol = K)
  } else {
    logits <- sweep(as.matrix(A) %*% gamma, 2, alpha, "+")
    if (is.null(p)) p <- nrow(logits)
  }
  ext <- cbind(0, logits)
  ext <- ext - apply(ext, 1L, max)
  ex  <- exp(ext)
  pi_row <- ex / rowSums(ex)
  if (is.null(A) || nrow(gamma) == 0L) {
    if (is.null(p)) stop(".sbayesrc_priors_from_gamma needs p when A is NULL")
    matrix(rep(as.numeric(pi_row), each = p), nrow = p, ncol = K + 1L)
  } else {
    pi_row
  }
}


# Fit a multinomial logistic regression of comp (0..K) on A. Returns
# alpha (K-vector) and gamma (m_ann x K matrix) with k=0 (spike) as the
# reference class.
#
# Two failure modes have to be guarded for at benchmark scale:
#   (a) EMPTY CLASSES: with small p some slab labels never appear. The
#       correct behaviour is "give this class ~0 prior mass" (very negative
#       alpha), NOT the padded 0 that softmax would then convert into equal
#       prior with the spike — the latter puts non-causal PIPs at ~0.8.
#   (b) DIVERGENT FITS: with few non-spike observations, unregularized
#       multinom overshoots. Heavy L2 (decay = 1) plus a coefficient cap
#       keeps priors sane. This is a Gibbs-sampler auxiliary update, not a
#       final scientific estimate, so aggressive regularization is fine.
.sbayesrc_fit_gamma <- function(comp, A, K,
                                empty_alpha = -12,
                                coef_cap    =   6) {
  m_ann <- ncol(A)
  y <- factor(comp, levels = as.character(0:K))
  # Empirical class counts — used both to detect empty classes and to
  # sanity-check the fitted intercepts.
  counts <- table(y)

  fit <- tryCatch(
    suppressMessages(suppressWarnings(nnet::multinom(
      y ~ ., data = data.frame(y = y, as.data.frame(A)),
      trace = FALSE, MaxNWts = 1e5, decay = 1.0, maxit = 200L
    ))),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(alpha = rep(empty_alpha, K),
                gamma = matrix(0, nrow = m_ann, ncol = K)))
  }
  coefs <- tryCatch(stats::coef(fit), error = function(e) NULL)
  if (is.null(coefs)) {
    return(list(alpha = rep(empty_alpha, K),
                gamma = matrix(0, nrow = m_ann, ncol = K)))
  }
  # coefs shape:
  #   - (K x (1 + m_ann)) matrix with rownames "1".."K" when >=2 non-reference
  #     classes are present in the data;
  #   - plain (1 + m_ann) numeric vector when only ONE non-reference class was
  #     fit (nnet::multinom collapses to a binary logistic). The class it
  #     belongs to has NO rowname on the vector, so infer it from `comp`.
  if (is.null(dim(coefs))) {
    present <- setdiff(sort(unique(as.integer(as.character(y)))), 0L)
    lone_class <- if (length(present) == 1L) present else 1L
    coefs <- matrix(coefs, nrow = 1L,
                    dimnames = list(as.character(lone_class),
                                    names(coefs)))
  }

  # Reindex to always have rows 1..K in that order, padding empty classes
  # with `empty_alpha` (intercept) and zero (annotation slopes) so that
  # they get essentially zero prior mass — not the softmax-uniform prior
  # a raw zero-padding would produce.
  full <- matrix(0, nrow = K, ncol = ncol(coefs),
                 dimnames = list(as.character(1:K), colnames(coefs)))
  seen <- rownames(coefs)
  if (!is.null(seen)) full[seen, ] <- coefs
  # Empty classes: fill intercept column with a strongly negative value
  unseen <- setdiff(as.character(1:K), seen)
  if (length(unseen) > 0L) full[unseen, 1L] <- empty_alpha

  # Cap magnitudes to prevent divergence when the data is scarce. Use
  # `full[] <-` so the matrix dim/dimnames survive (pmax on a matrix
  # returns a plain vector, which breaks the [, 1L] extraction below).
  full[] <- pmax(-coef_cap, pmin(coef_cap, full))

  alpha <- as.numeric(full[, 1L])
  gamma <- t(full[, -1L, drop = FALSE])
  list(alpha = alpha, gamma = gamma)
}


# Warn (once per session) that sbayesrc is running WITHOUT pooled annotation
# coefficients and has therefore fallen back to the unstable in-loop refit.
# Verified 2026-07-29: that path inflates posterior mass ~15x relative to the
# pooled path even under the true LD, so a silent fallback would corrupt any
# annotated-arm result without producing an error.
.sbayesrc_warned <- new.env(parent = emptyenv())
.sbayesrc_warn_no_pooled <- function() {
  if (isTRUE(.sbayesrc_warned$done)) return(invisible(NULL))
  .sbayesrc_warned$done <- TRUE
  warning("sbayesrc: annotations supplied but pooled_gamma is NULL, so the ",
          "unstable in-loop annotation refit is being used (empirically ~15x ",
          "posterior-mass inflation vs the pooled path). This happens when ",
          "run_sbayesrc_scenario_setup() returned nothing - check that the ",
          "scenario has >=2 regions with consistent annotation matrices.",
          call. = FALSE)
  invisible(NULL)
}
