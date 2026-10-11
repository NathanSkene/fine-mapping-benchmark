# =============================================================================
# marginal_z.R
#
# Wrapper for the marginal-z baseline.
#
# This is a model-free baseline that ranks variants by the magnitude of their
# marginal z-score and converts the ranking into normalised "PIPs":
#
#   pip_j = |z_j| / sum_k |z_k|
#
# These normalized scores sum to one but are not posterior probabilities.
# They support a ranking baseline; calibration, Bayesian credible-set coverage
# and PIP-threshold discovery claims do not apply to them. A score-mass set is
# provided separately for descriptive use, never as a posterior credible set.
# A fitted method can legitimately rank worse in a finite sample or under
# model mismatch; there is no universal dominance guarantee.
#
# No external packages or binaries are required.
#
# This file provides:
#   - run_marginal_z()        : runs the baseline on a single region (explicit inputs)
#   - run_marginal_z_region() : adapter called by run_methods()
# =============================================================================


# =============================================================================
# Run marginal_z on a single region
# =============================================================================

#' Run the marginal-z baseline on a single region
#'
#' Computes normalized absolute z ranking scores in the legacy \code{pip}
#' field. A descriptive score-mass set appears in \code{additional}; posterior
#' credible sets are not provided.
#'
#' This is a model-free baseline — no LD matrix, sample size, or prior
#' assumption is used. Its purpose is to anchor reported ranking performance.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param coverage Numeric. Descriptive score-mass threshold. Default: 0.95.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Normalised |z|-based PIPs that sum
#'     to 1.}
#'   \item{credible_sets}{Empty list: this ranking model has no posterior credible sets.}
#'   \item{method}{Character. Always \code{"marginal_z"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used (\code{coverage}).}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds (effectively
#'     zero — included for interface consistency with the other methods).}
#'   \item{additional}{Score semantics and a descriptive score-mass set.}
#' }
#'
#' @export
run_marginal_z <- function(z, coverage = 0.95) {

  # --- Validate ---------------------------------------------------------------

  p <- length(z)

  stopifnot(
    "z must be a non-empty numeric vector" =
      is.numeric(z) && p >= 1 && all(is.finite(z)),
    "coverage must be a single number in (0, 1)" =
      is.numeric(coverage) && length(coverage) == 1 &&
      coverage > 0 && coverage < 1
  )

  # --- Compute PIPs -----------------------------------------------------------

  start_time <- proc.time()

  abs_z <- abs(z)

  if (sum(abs_z) == 0) {
    # Degenerate case: all z-scores are zero. Spread PIP uniformly.
    pip <- rep(1 / p, p)
  } else {
    pip <- abs_z / sum(abs_z)
  }

  # --- Derive credible set ----------------------------------------------------
  # Greedy: sort variants by PIP descending, accumulate until >= coverage.

  ord        <- order(pip, decreasing = TRUE)
  cumulative <- cumsum(pip[ord])
  n_in_cs    <- which(cumulative >= coverage)[1]
  if (is.na(n_in_cs)) n_in_cs <- p   # safety: include all if PIPs don't sum to 1
  cs_indices <- sort(ord[seq_len(n_in_cs)])

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = list(),
    method          = "marginal_z",
    input_type      = "summary",
    params          = list(coverage = coverage),
    runtime_seconds = elapsed,
    additional      = list(
      score_semantics = "normalized_absolute_z_rank_score",
      is_posterior_probability = FALSE,
      credible_set_semantics = "not_provided",
      pip_mass_set = as.integer(cs_indices)
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run the marginal-z baseline on a single region from simulation data
#'
#' Thin adapter that extracts z-scores from the simulation's
#' \code{region_pheno} object and calls \code{\link{run_marginal_z}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes}. Not
#'   used by marginal_z but required for a consistent adapter signature.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}.
#' @param ... Additional arguments passed to \code{\link{run_marginal_z}}
#'   (e.g. \code{coverage}).
#'
#' @return The output of \code{\link{run_marginal_z}}.
#' @export
run_marginal_z_region <- function(region_geno, region_pheno, ...) {
  run_marginal_z(z = region_pheno$z, ...)
}
