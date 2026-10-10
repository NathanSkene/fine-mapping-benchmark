# Typed Gaussian score likelihood: z | beta ~ N(sqrt(n) R beta, sigma^2 R).
# With theta=sqrt(n)*beta, the fixed-noise SuSiE RSS likelihood uses
# XtX=R and Xty=z. Passing n to susie_rss instead invokes its OLS-t
# transformation, which is appropriate only for marginal_ols_t inputs.
.fm_susie_rss <- function(z, R, n, statistic_type = "marginal_ols_t",
                          residual_variance = 1, ...) {
  statistic_type <- match.arg(statistic_type, c("marginal_ols_t", "normal_score"))
  args <- list(...)
  if (statistic_type == "marginal_ols_t") {
    return(do.call(susieR::susie_rss, c(list(z = z, R = R, n = n), args)))
  }
  if (!is.null(args$unmappable_effects) &&
      match.arg(args$unmappable_effects, c("none", "inf", "ash")) == "inf") {
    stop(paste0("unsupported_input: normal_score SuSiE-inf needs a supported fixed-noise ",
      "infinitesimal variance estimator. susieR estimates tau2 jointly with phenotype ",
      "residual variance and requires observed phenotype sufficient statistics; ",
      "disabling this update silently removes the infinitesimal component."), call. = FALSE)
  }
  if (length(n) != 1L || !is.finite(n) || n <= 0 ||
      length(residual_variance) != 1L || !is.finite(residual_variance) || residual_variance <= 0)
    stop("Gaussian scores need positive n and residual_variance.", call. = FALSE)
  # This Gaussian experiment has a specified noise covariance. y'y is not
  # observed, so estimating a phenotype residual variance would be fictitious.
  args$estimate_residual_variance <- FALSE
  args$residual_variance <- residual_variance
  args$prior_variance <- n * (args$prior_variance %||% .2)
  fit <- withCallingHandlers(
    do.call(susieR::susie_rss, c(list(z = z, R = R, n = NULL), args)),
    warning = function(w) {
      if (grepl("Providing the sample size", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  attr(fit, "score_effect_divisor") <- sqrt(n)
  attr(fit, "score_likelihood") <- "z_given_beta_normal_sqrt_n_R_beta_sigma2_R"
  fit
}

.fm_susie_posterior_mean <- function(fit) {
  colSums(fit$alpha * fit$mu) / (attr(fit, "score_effect_divisor") %||% 1)
}

# OLS t determines the sample correlation exactly. The author Python
# infinitesimal implementation uses XtX=n*R, Xty=sqrt(n)*z, yty=n*meansq.
# Population-SD standardizing X, and centering/scaling y to the chosen meansq,
# therefore gives z=sqrt(n*meansq)*r. No unobserved phenotype norm is invented.
.fm_ols_t_to_score <- function(z, n, meansq = 1) {
  if (!is.numeric(z) || any(!is.finite(z)) || length(n) != 1L ||
      !is.finite(n) || n <= 2 || n != floor(n) || length(meansq) != 1L ||
      !is.finite(meansq) || meansq <= 0)
    stop("OLS-t reconstruction needs finite statistics, integer n > 2, and positive meansq.", call. = FALSE)
  # This form avoids overflowing t^2 for large finite statistics.
  correlation <- sign(z) / sqrt(1 + (sqrt(n - 2) / abs(z))^2)
  sqrt(n) * sqrt(meansq) * correlation
}
