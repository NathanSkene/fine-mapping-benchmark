# =============================================================================
# wrapper_polyfun_ldsc.R
#
# Corrected LD-score PolyFun-style prior (§0.6 Error 1, Option b in the
# autoresearch plan). Companion to (not replacement of) polyfun_est.
#
# What polyfun_est gets wrong
# ---------------------------
# polyfun_est regresses chi^2_j on each variant's OWN annotation vector:
#
#     chi2_j = tau_0 + sum_k tau_k * A_{j,k} + noise
#
# Under LD this is statistically incorrect. chi2_j is inflated by tagging
# of nearby causal variants, so the fit attributes signal to annotations
# carried by tagging variants rather than causal ones.
#
# What polyfun_ldsc does instead
# ------------------------------
# Follows the canonical S-LDSC / PolyFun model:
#
#     E[chi2_j] = 1 + N * sum_c tau_c * l_{j,c},
#     l_{j,c} = sum_k r_{j,k}^2 * A_{k,c}
#
# where l_{j,c} is the "annotation-c LD score" of variant j - the sum of
# squared LD to every variant carrying annotation c. tau is estimated by
# weighted non-negative least squares. Per-SNP prior variance is then
#
#     sigma^2(j) proportional to sum_c tau_c * A_{j,c}
#
# (using annotation membership, not LD scores, at the prior-conversion
# step - matching canonical PolyFun).
#
# Held-out chromosome or region prior training
# ---------------------------
# With complete chromosome metadata, the scenario hook excludes every region
# on the target chromosome. Otherwise it excludes the target region only.
# The output records which grouping was possible; regions are not silently
# called chromosomes. A single available group cannot provide a held-out fit.
#
# When run standalone or on a single-region input, tau is fitted from
# that region's data alone (no LOCO possible) and a warning is recorded
# in prior_source.
#
# When annotations are absent, the wrapper falls back to a uniform prior
# and reports prior_source = "uniform_fallback".
# =============================================================================


# =============================================================================
# setup_polyfun_ldsc()
# =============================================================================

#' Check that susieR is available for polyfun_ldsc
#'
#' @return Invisible TRUE if susieR is available.
#' @export
setup_polyfun_ldsc <- function() {
  if (!requireNamespace("susieR", quietly = TRUE)) {
    stop(
      "susieR is required for polyfun_ldsc but is not installed.\n",
      "Install it with:\n",
      "  remotes::install_github('stephenslab/susieR')",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# =============================================================================
# run_polyfun_ldsc(): single-region interface
# =============================================================================

#' Run corrected LD-score PolyFun on a single region
#'
#' Estimates per-annotation contributions by regressing chi-squared on
#' annotation-weighted LD scores (S-LDSC), converts them into per-SNP prior
#' variances, and runs \code{susieR::susie_rss} with those priors.
#'
#' Typically called via \code{run_polyfun_ldsc_region()} from
#' \code{run_methods()}, which uses the scenario-setup hook to compute
#' LOCO-fitted \code{pooled_tau} values per region.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD correlation matrix (p x p).
#' @param n Integer. Sample size.
#' @param annotations Matrix or NULL. Functional annotation matrix (p x m).
#'   NULL runs with a uniform prior (reports prior_source =
#'   "uniform_fallback").
#' @param pooled_tau Numeric vector, named list, or NULL. If numeric, used
#'   directly as the (length m + 1) coefficient vector. If a named list
#'   (as returned by \code{run_polyfun_ldsc_scenario_setup()}), the entry
#'   keyed by \code{region_id} is looked up. Default: NULL.
#' @param region_id Character/integer or NULL. Identifier used to look up
#'   the correct per-region tau vector inside a list-valued
#'   \code{pooled_tau}. Ignored when \code{pooled_tau} is numeric or NULL.
#'   Automatically forwarded by \code{run_polyfun_ldsc_region()}.
#' @param L Integer. Max single-effect components. Default: 10.
#' @param coverage Numeric. Credible-set coverage. Default: 0.95.
#' @param min_abs_corr Numeric. Purity threshold. Default: 0.5.
#' @param max_iter Integer. IBSS iterations. Default: 100.
#' @param estimate_residual_variance Logical. Default: TRUE.
#' @param estimate_prior_variance Logical. Default: TRUE.
#' @param variant_ids Character or NULL. Optional variant labels.
#' @param statistic_type \code{marginal_ols_t} (default) or \code{normal_score}.
#' @param residual_variance Specified Gaussian score noise variance; default 1.
#' @param ... Ignored (for wrapper compatibility).
#'
#' @return List with pip, credible_sets, method = "polyfun_ldsc", params,
#'   runtime_seconds, additional (prior_weights, tau, prior_source, ...).
#' @export
run_polyfun_ldsc <- function(z, LD, n,
                             annotations = NULL,
                             pooled_tau  = NULL,
                             region_id   = NULL,
                             L           = 10,
                             coverage    = 0.95,
                             min_abs_corr = 0.5,
                             max_iter    = 100,
                             estimate_residual_variance = TRUE,
                             estimate_prior_variance    = TRUE,
                             variant_ids = NULL,
                             statistic_type = "marginal_ols_t", residual_variance = 1,
                             ...) {
  setup_polyfun_ldsc()
  t0 <- Sys.time()

  p <- length(z)
  stopifnot("LD must be p x p" = is.matrix(LD) && all(dim(LD) == c(p, p)))

  prior_source <- if (is.null(annotations)) {
    "uniform_fallback"
  } else if (!is.null(pooled_tau)) {
    "loco_scenario_setup"
  } else {
    "single_region_ldsc"
  }

  # `pooled_tau` may arrive as either:
  #   (a) a numeric vector - used directly (single-region / no-LOCO caller), or
  #   (b) a named list of per-region tau vectors, keyed by region_id
  #       (LOCO caller from run_polyfun_ldsc_scenario_setup). In that case
  #       we look up this region's own tau via the `region_id` arg.
  tau <- NULL
  if (!is.null(annotations)) {
    if (is.list(pooled_tau) && !is.data.frame(pooled_tau)) {
      key <- if (!is.null(region_id)) as.character(region_id) else NA_character_
      if (!is.na(key) && !is.null(pooled_tau[[key]])) {
        tau <- as.numeric(pooled_tau[[key]])
      }
    } else if (!is.null(pooled_tau)) {
      tau <- as.numeric(pooled_tau)
    }
    if (is.null(tau)) {
      ldsc_mat <- .ldscore_matrix(annotations, LD)
      tau      <- .fit_ldsc_stratified(z^2, ldsc_mat,
                                       baseline_ldscore = rowSums(LD * LD), n = n)
      if (is.list(pooled_tau)) prior_source <- "loco_lookup_missed"
    }
  }

  prior_weights <- if (is.null(tau)) {
    rep(1 / p, p)
  } else {
    .tau_annot_to_prior_weights(annotations, tau)
  }

  fit <- .fm_susie_rss(
    z            = z, R = LD, n = n, L = L, coverage = coverage,
    statistic_type = statistic_type, residual_variance = residual_variance,
    min_abs_corr = min_abs_corr, max_iter = max_iter,
    prior_weights = prior_weights,
    estimate_residual_variance = estimate_residual_variance,
    estimate_prior_variance    = estimate_prior_variance
  )

  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  cs_list <- if (is.null(fit$sets$cs)) list() else fit$sets$cs
  cs_list <- lapply(cs_list, as.integer)

  list(
    pip             = as.numeric(fit$pip),
    credible_sets   = cs_list,
    method          = "polyfun_ldsc",
    input_type      = "summary",
    params          = list(
      L = L, coverage = coverage, min_abs_corr = min_abs_corr,
      max_iter = max_iter,
      estimate_residual_variance = estimate_residual_variance,
      statistic_type = statistic_type, score_residual_variance = residual_variance,
      effective_estimate_residual_variance = if (statistic_type == "normal_score") FALSE else estimate_residual_variance,
      estimate_prior_variance    = estimate_prior_variance
    ),
    runtime_seconds = runtime,
    additional      = list(
      tau           = tau,
      prior_weights = prior_weights,
      prior_source  = prior_source,
      prior_training_unit = attr(pooled_tau, "held_out_unit") %||% "not_cross_fitted",
      variant_ids   = variant_ids
    )
  )
}


# =============================================================================
# run_polyfun_ldsc_region(): the adapter called by run_methods()
# =============================================================================

#' Adapter for run_methods() - forwards to run_polyfun_ldsc
#'
#' @param region_geno One element of \code{simulation$genotypes}.
#' @param region_pheno One element of a scenario's regions.
#' @param pooled_tau Numeric or NULL. LOCO-fitted per-region coefficients
#'   supplied by the scenario_setup hook.
#' @param ... Forwarded to run_polyfun_ldsc().
#' @param statistic_type Optional summary type inherited from the region;
#'   defaults to \code{marginal_ols_t} when no type is declared.
#'
#' @return Output of run_polyfun_ldsc.
#' @export
run_polyfun_ldsc_region <- function(region_geno, region_pheno,
                                    pooled_tau = NULL, ..., statistic_type = NULL) {
  A <- region_geno$annotations_matrix
  if (is.null(A)) A <- region_pheno$annotations_matrix
  run_polyfun_ldsc(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    annotations = A,
    pooled_tau  = pooled_tau,
    statistic_type = statistic_type %||% region_pheno$statistic_type %||% "marginal_ols_t",
    region_id   = region_geno$region_id,
    variant_ids = region_geno$variant_ids,
    ...
  )
}


# =============================================================================
# Scenario-level setup hook: LOCO tau per region
# =============================================================================

#' Fit LOCO tau values across the scenario for polyfun_ldsc
#'
#' For each region i, fit S-LDSC tau on the pooled (chi^2, LD-score) rows
#' from all regions except i, and return the per-region tau vectors so
#' each region's prior is built from data that never touched it. When
#' fewer than two regions carry an annotation matrix, returns an empty
#' list and the per-region wrapper falls back to a single-region fit.
#'
#' @param genotypes List of region_geno objects.
#' @param regions List of region_pheno objects for the current scenario.
#' @param user_args User's method_args for polyfun_ldsc (unused here).
#'
#' @return A list with one element - \code{pooled_tau} - itself a list of
#'   per-region tau vectors, or an empty list to signal "not applicable".
#' @export
run_polyfun_ldsc_scenario_setup <- function(genotypes, regions, user_args) {
  if (!is.null(user_args$pooled_tau)) return(list(pooled_tau = user_args$pooled_tau))
  n_regions <- length(regions)
  if (n_regions < 2L) return(list())

  chi2_list <- vector("list", n_regions)
  ell_list  <- vector("list", n_regions)
  ncols     <- integer(n_regions)
  baseline_list <- vector("list", n_regions)
  sample_sizes <- vector("list", n_regions)

  for (i in seq_len(n_regions)) {
    A_i <- genotypes[[i]]$annotations_matrix
    if (is.null(A_i)) A_i <- regions[[i]]$annotations_matrix
    if (is.null(A_i) || !is.matrix(A_i)) return(list())

    z_i <- regions[[i]]$z
    if (is.null(z_i) || length(z_i) != nrow(A_i)) return(list())

    LD_i <- genotypes[[i]]$LD
    if (is.null(LD_i) || !all(dim(LD_i) == c(nrow(A_i), nrow(A_i)))) return(list())

    chi2_list[[i]] <- z_i^2
    ell_list[[i]]  <- .ldscore_matrix(A_i, LD_i)
    baseline_list[[i]] <- rowSums(LD_i * LD_i)
    sample_sizes[[i]] <- rep(genotypes[[i]]$n, length(z_i))
    ncols[i]       <- ncol(A_i)
  }

  if (length(unique(ncols)) != 1L) return(list())

  chromosomes <- vapply(genotypes, function(g) as.character(g$chrom %||% NA_character_), character(1))
  have_chromosomes <- all(!is.na(chromosomes) & nzchar(chromosomes))
  groups <- if (have_chromosomes) sub("^chr", "", tolower(chromosomes)) else as.character(seq_len(n_regions))
  if (length(unique(groups)) < 2L) return(list())

  # Hold out the entire chromosome when metadata exists; otherwise a region.
  # Key the result by region_id so the region wrapper can look up its own tau
  # (run_methods()'s scenario_setup merge is scenario-wide, not per-region).
  per_region_tau <- list()
  for (i in seq_len(n_regions)) {
    idx <- which(groups != groups[i])
    chi2_pooled <- unlist(chi2_list[idx], use.names = FALSE)
    ell_pooled  <- do.call(rbind, ell_list[idx])
    tau_i <- tryCatch(.fit_ldsc_stratified(chi2_pooled, ell_pooled,
                         baseline_ldscore = unlist(baseline_list[idx], use.names = FALSE),
                         n = unlist(sample_sizes[idx], use.names = FALSE)),
                      error = function(e) NULL)
    if (is.null(tau_i)) stop("Held-out LD-score prior estimation failed.", call. = FALSE)
    key_i <- as.character(genotypes[[i]]$region_id %||% i)
    per_region_tau[[key_i]] <- tau_i
  }

  attr(per_region_tau, "held_out_unit") <- if (have_chromosomes) "chromosome" else "region"
  list(pooled_tau = per_region_tau)
}


# =============================================================================
# Internal helpers
# =============================================================================

# Annotation-c LD score of variant j:
#   l_{j,c} = sum_k r_{j,k}^2 * A_{k,c}
# Computed as R^2 %*% A where R^2 is the element-wise squared LD matrix.
.ldscore_matrix <- function(A, LD) {
  as.matrix((LD * LD) %*% A)
}

# Restricted S-LDSC-style fit with a fixed null/confounding intercept of 1.
# The first genetic coefficient multiplies the ALL-SNP LD score. It must not
# be the regression intercept: the latter includes null sampling variation and
# population confounding and has no per-variant genetic interpretation.
# This small-sample NNLS comparator is not the upstream PolyFun implementation.
.fit_ldsc_stratified <- function(chi2, ldscore_matrix, baseline_ldscore = NULL, n = 1) {
  ell <- as.matrix(ldscore_matrix)
  m <- ncol(ell)
  if (is.null(baseline_ldscore)) baseline_ldscore <- rep(1, nrow(ell))
  if (length(n) == 1L) n <- rep(n, nrow(ell))
  if (length(chi2) != nrow(ell) || length(baseline_ldscore) != nrow(ell) ||
      length(n) != nrow(ell) || any(!is.finite(n)) || any(n <= 0)) {
    stop("LD-score fit inputs have inconsistent lengths or invalid sample sizes.", call. = FALSE)
  }
  y <- (as.numeric(chi2) - 1) / n
  X <- cbind(baseline_ldscore, ell)
  total_ldsc <- as.numeric(baseline_ldscore)
  w <- 1 / pmax(total_ldsc, 1e-8)
  bad <- !is.finite(y) | !is.finite(w) | apply(!is.finite(X), 1, any)
  if (any(bad)) { y <- y[!bad]; X <- X[!bad, , drop = FALSE]; w <- w[!bad] }
  if (!length(y)) return(rep(0, m + 1L))
  .nnls_active_set(X, y, w)
}

# Convex NNLS by coordinate minimization. Unlike irreversible sign-dropping,
# coordinates may re-enter after another coefficient changes. Validate the
# projected gradient (the KKT conditions) before calling a fit converged.
.nnls_active_set <- function(X, y, w, max_iter = 10000L, tolerance = 1e-10) {
  Xw <- X * sqrt(w)
  yw <- y * sqrt(w)
  q <- ncol(Xw)
  beta <- rep(0, q)
  residual <- yw
  norms <- colSums(Xw * Xw)
  scale <- max(.Machine$double.eps, max(abs(as.numeric(crossprod(Xw, yw)))))
  for (iter in seq_len(max_iter)) {
    for (j in which(norms > 0)) {
      old <- beta[j]
      beta[j] <- max(0, old + sum(Xw[, j] * residual) / norms[j])
      residual <- residual - Xw[, j] * (beta[j] - old)
    }
    gradient <- as.numeric(crossprod(Xw, residual))
    violation <- ifelse(beta > 0, abs(gradient), pmax(gradient, 0))
    if (max(violation) <= tolerance * scale) return(beta)
  }
  stop("LD-score NNLS did not satisfy its KKT convergence criterion.", call. = FALSE)
}

# Turn tau into per-SNP prior weights via sigma^2_j = sum_c tau_c * A_{j,c}
# (with tau_0 as the all-variant genetic baseline). Floored to keep priors strictly
# positive and normalised to sum to 1 so susieR treats it as a probability
# distribution over variants.
.tau_annot_to_prior_weights <- function(A, tau) {
  p <- nrow(A)
  m <- ncol(A)
  stopifnot(length(tau) == m + 1L)
  intercept <- tau[1L]
  slope     <- tau[-1L]
  sigma2 <- intercept + as.numeric(A %*% slope)
  floor_val <- max(.Machine$double.eps, 1e-8 * mean(abs(sigma2)))
  sigma2 <- pmax(sigma2, floor_val)
  sigma2 / sum(sigma2)
}
