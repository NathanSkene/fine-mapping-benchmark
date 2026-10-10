# =============================================================================
# simulate_phenotypes.R
#
# Simulate phenotypes, compute summary statistics, and record ground truth
# for fine-mapping benchmarking.
#
# Takes the output of simulate_genotypes() and adds:
#   - Phenotype vector y
#   - Summary statistics (z-scores, beta-hat, se)
#   - LD matrix
#   - Ground truth (causal indices, effect sizes, etc.)
#   - Annotations (if applicable)
# =============================================================================


#' Simulate phenotypes and summary statistics for fine-mapping benchmarking
#'
#' Takes the output of \code{\link{simulate_genotypes}} and, for each region,
#' selects causal variants, generates a phenotype vector, computes marginal
#' summary statistics, and records the ground truth.
#'
#' @param genotypes List. Output from \code{\link{simulate_genotypes}}.
#' @param S Integer or integer vector. Number of causal variants per region.
#'   If scalar, the same value is used for all regions. If a vector, must
#'   have length equal to the number of regions. Default: 1.
#' @param phi Numeric or numeric vector. Proportion of variance explained
#'   (PVE) by the genetic component. If scalar, same for all regions.
#'   Must be in (0, 1). Default: 0.1.
#' @param model Character. Genetic architecture model. Either \code{"sparse"}
#'   (standard sparse model, y = Xb + e) or \code{"sparse_inf"} (sparse +
#'   infinitesimal effects). Default: "sparse".
#' @param p_causal Numeric. Proportion of genetic variance attributable to the
#'   sparse (causal) component. Only used when \code{model = "sparse_inf"}.
#'   Must be in (0, 1]. Default: 0.5.
#' @param inf_model Character. Which infinitesimal formulation to use. Either
#'   \code{"beatrice"} (noncausal variants only) or \code{"susie_inf"} (all
#'   variants). Only used when \code{model = "sparse_inf"}. Default: "beatrice".
#' @param effect_distribution Character. Distribution for causal effect sizes.
#'   \code{"normal"} draws from N(0, effect_variance).
#'   \code{"equal"} assigns equal effect magnitudes with random signs.
#'   LD means equal magnitudes do not imply equal variance contributions.
#'   Default: "normal".
#' @param effect_variance Numeric. Variance of the normal effect size
#'   distribution. Only used when \code{effect_distribution = "normal"}.
#'   Default: 0.36 (i.e. sd = 0.6, following SuSiE).
#' @param annotations Character or matrix. Controls functional annotation
#'   simulation. \code{"none"} for no annotations (causal variants selected
#'   uniformly at random). \code{"binary"} for synthetic binary annotations.
#'   \code{"continuous"} for synthetic continuous annotations (drawn from
#'   N(0,1)). Or a user-supplied matrix (p x m) of annotation values.
#'   Default: "none".
#' @param n_annotations Integer. Number of annotation categories. Only used
#'   when \code{annotations} is \code{"binary"} or \code{"continuous"}.
#'   Default: 3.
#' @param annotation_proportions Numeric, scalar, vector, or NULL. Controls
#'   the proportion of variants with value 1 in each binary annotation.
#'   If NULL, proportions are drawn randomly from Uniform(0.01, 0.30) for
#'   each annotation. If scalar, that proportion is used for all annotations.
#'   If a vector, must have length \code{n_annotations}. Only used when
#'   \code{annotations = "binary"}. Default: NULL.
#' @param enrichment Numeric, scalar, vector, or NULL. Fold-enrichment of
#'   each annotation for causal variant selection. If NULL, enrichments are
#'   drawn from Uniform(2, 10) for each annotation. If scalar, that
#'   enrichment is used for all annotations. If a vector, must have length
#'   \code{n_annotations}. Values must be > 0 (values < 1 indicate
#'   depletion). Default: NULL.
#' @param annotation_correlation Numeric scalar in \code{[0, 1]}. Approximate
#'   pairwise correlation between binary annotations that share the same
#'   enrichment fold. \code{0} (default) recovers the original independent
#'   generation. When positive, annotations at the same enrichment level are
#'   generated via a shared-factor latent Gaussian construction and
#'   thresholded to preserve each column's marginal frequency; annotations
#'   with different enrichment folds remain uncorrelated. Only affects
#'   \code{annotations = "binary"}; ignored when \code{enrichment} is NULL
#'   or has an unexpected length. The realised binary correlation is
#'   attenuated relative to the latent target by thresholding.
#' @param seed Integer or NULL. Random seed. Default: NULL.
#' @param save Logical. If TRUE, save the returned list (genotypes augmented
#'   with phenotypes) as an \code{.rds} file inside \code{output_dir}. The
#'   filename encodes the key simulation parameters and the seed. Default: FALSE.
#' @param output_dir Character. Directory in which to save the result when
#'   \code{save = TRUE}. Created automatically if it does not exist.
#'   Default: \code{"results"}.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return The input \code{genotypes} list, with additional fields appended
#'   to each region:
#'   \describe{
#'     \item{y}{Phenotype vector (n x 1).}
#'     \item{z}{Marginal z-scores (p x 1).}
#'     \item{beta_hat}{Marginal effect size estimates (p x 1).}
#'     \item{se}{Standard errors of marginal effects (p x 1).}
#'     \item{LD}{LD (correlation) matrix (p x p).}
#'     \item{annotations_matrix}{Annotation matrix (p x m), or NULL if none.}
#'     \item{truth}{List containing ground truth:
#'       \describe{
#'         \item{causal_indices}{Integer vector of causal variant indices.}
#'         \item{causal_effects}{Numeric vector of true effect sizes.}
#'         \item{beta_true}{Legacy p-length sparse-component effect vector.}
#'         \item{beta_sparse, beta_infinitesimal, beta_total}{Explicit effect decomposition.}
#'         \item{estimand}{Sparse-component membership, not all nonzero effects.}
#'         \item{pve_expected, variance_components}{Expected PVE and sample variance ledger.}
#'         \item{pve}{Realised PVE.}
#'         \item{S}{Number of causal variants.}
#'         \item{phi}{Target PVE.}
#'         \item{model}{Genetic architecture model used.}
#'         \item{effect_distribution}{Effect size distribution used.}
#'         \item{annotation_type}{Type of annotations used.}
#'         \item{enrichment}{Enrichment values used (or NULL).}
#'         \item{annotation_proportions}{Proportions used for binary annotations (or NULL).}
#'       }
#'     }
#'   }
#'
#' @examples
#' \dontrun{
#' geno <- simulate_genotypes(n_regions = 2, n = 200, p = 100, seed = 1)
#'
#' # Sparse model, no annotations
#' sim <- simulate_phenotypes(geno, S = 3, phi = 0.2, seed = 1)
#'
#' # Sparse + infinitesimal, with binary annotations
#' sim <- simulate_phenotypes(
#'   geno, S = 2, phi = 0.3,
#'   model = "sparse_inf", p_causal = 0.5,
#'   annotations = "binary", n_annotations = 3,
#'   seed = 1
#' )
#' }
#'
#' @export
simulate_phenotypes <- function(genotypes,
                                S = 1,
                                phi = 0.1,
                                model = "sparse",
                                p_causal = 0.5,
                                inf_model = "beatrice",
                                effect_distribution = "normal",
                                effect_variance = 0.36,
                                annotations = "none",
                                n_annotations = 3,
                                annotation_proportions = NULL,
                                enrichment = NULL,
                                annotation_correlation = 0,
                                seed = NULL,
                                save = FALSE,
                                output_dir = "results",
                                verbose = TRUE) {

  stopifnot(
    "annotation_correlation must be a scalar in [0, 1]" =
      is.numeric(annotation_correlation) &&
      length(annotation_correlation) == 1L &&
      !is.na(annotation_correlation) &&
      annotation_correlation >= 0 && annotation_correlation <= 1
  )

  # --- Input validation -------------------------------------------------------

  n_regions <- length(genotypes)
  if (!is.list(genotypes) || n_regions < 1L) {
    stop("genotypes must contain at least one region.", call. = FALSE)
  }

  # S: scalar or vector
  S <- validate_per_region_param(S, n_regions, "S", integer_valued = TRUE)

  if (any(S < 1L)) stop("S must be positive; use genome-wide simulation for null regions.", call. = FALSE)

  # phi: scalar or vector
  phi <- validate_per_region_param(phi, n_regions, "phi")
  stopifnot(
    "phi must be in (0, 1)" = all(phi > 0 & phi < 1)
  )

  # model
  model <- match.arg(model, choices = c("sparse", "sparse_inf"))

  # p_causal
  if (model == "sparse_inf") {
    stopifnot(
      "p_causal must be a single number in (0, 1]" =
        is.numeric(p_causal) && length(p_causal) == 1 &&
        p_causal > 0 && p_causal <= 1
    )
  }

  # inf_model
  inf_model <- match.arg(inf_model, choices = c("beatrice", "susie_inf"))

  # effect_distribution
  effect_distribution <- match.arg(effect_distribution, choices = c("normal", "equal"))

  # effect_variance
  stopifnot(
    "effect_variance must be a positive number" =
      is.numeric(effect_variance) && length(effect_variance) == 1 &&
      effect_variance > 0
  )

  # annotations
  user_annotation_matrix <- NULL
  if (is.matrix(annotations)) {
    if (!is.numeric(annotations) || ncol(annotations) < 1L || any(!is.finite(annotations))) {
      stop("annotations must be a finite numeric matrix with at least one column.", call. = FALSE)
    }
    user_annotation_matrix <- annotations
    n_annotations <- ncol(annotations)
    annotation_type <- "user_supplied"
  } else if (is.character(annotations)) {
    annotation_type <- match.arg(annotations, choices = c("none", "binary", "continuous"))
  } else {
    stop("annotations must be 'none', 'binary', 'continuous', or a matrix.", call. = FALSE)
  }

  # n_annotations
  if (annotation_type %in% c("binary", "continuous")) {
    stopifnot(
      "n_annotations must be a positive integer" =
        is.numeric(n_annotations) && length(n_annotations) == 1 &&
        n_annotations == floor(n_annotations) && n_annotations >= 1
    )
    n_annotations <- as.integer(n_annotations)
  }

  # annotation_proportions (binary only)
  if (annotation_type == "binary" && !is.null(annotation_proportions)) {
    annotation_proportions <- validate_annotation_param(
      annotation_proportions, n_annotations, "annotation_proportions"
    )
    stopifnot(
      "annotation_proportions must be in (0, 1)" =
        all(annotation_proportions > 0 & annotation_proportions < 1)
    )
  }

  # enrichment
  if (annotation_type != "none" && !is.null(enrichment)) {
    enrichment <- validate_annotation_param(
      enrichment, n_annotations, "enrichment"
    )
    stopifnot(
      "enrichment must be positive" = all(enrichment > 0)
    )
  }

  # Check S <= p for each region
  for (i in seq_len(n_regions)) {
    if (S[i] > genotypes[[i]]$p) {
      stop(
        sprintf(
          "Region %d: S = %d but only %d SNPs available. S must be <= p.",
          i, S[i], genotypes[[i]]$p
        ),
        call. = FALSE
      )
    }
  }

  # --- Set seed ---------------------------------------------------------------

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # An annotation-to-causal mapping is one scenario-level architecture.
  # Drawing folds separately in each region confounded pooled-prior methods.
  if (annotation_type != "none" && is.null(enrichment)) {
    enrichment <- runif(n_annotations, min = 2, max = 10)
  }

  # --- Simulate each region ---------------------------------------------------

  for (i in seq_len(n_regions)) {
    if (verbose) {
      message(sprintf(
        "Region %d/%d: S=%d, phi=%.3f, model=%s",
        i, n_regions, S[i], phi[i], model
      ))
    }

    p_i <- genotypes[[i]]$p
    n_i <- genotypes[[i]]$n
    X_i <- genotypes[[i]]$X  # standardised genotype matrix

    # --- Simulate annotations for this region ---------------------------------
    # If a pre-computed annotation matrix is already stored in this region
    # (e.g. generated once by run_simulation and shared across scenarios),
    # reuse it directly instead of generating a new one.

    if (!is.null(genotypes[[i]]$annotations_matrix) &&
        annotation_type %in% c("binary", "continuous")) {
      A_i     <- genotypes[[i]]$annotations_matrix
      props_i <- genotypes[[i]]$annotation_proportions
    } else {
      annot_result <- simulate_annotations_for_region(
        p = p_i,
        annotation_type = annotation_type,
        n_annotations = if (annotation_type %in% c("binary", "continuous")) n_annotations else 0,
        annotation_proportions = annotation_proportions,
        user_annotation_matrix = user_annotation_matrix,
        annotation_correlation = annotation_correlation,
        enrichment = enrichment
      )
      A_i     <- annot_result$matrix
      props_i <- annot_result$proportions
    }

    if (n_i < 3L || !is.matrix(X_i) || any(dim(X_i) != c(n_i, p_i)) ||
        any(!is.finite(X_i)) || any(apply(X_i, 2, var) <= 0)) {
      stop("Each region needs at least three individuals and finite, variable genotype columns.",
           call. = FALSE)
    }

    # --- Select causal variants -----------------------------------------------

    causal_result <- select_causal_variants(
      p = p_i,
      S = S[i],
      annotation_matrix = A_i,
      enrichment = enrichment,
      n_annotations = if (annotation_type %in% c("binary", "continuous")) n_annotations else 0,
      annotation_type = annotation_type
    )

    causal_indices <- causal_result$causal_indices
    enrichment_used <- causal_result$enrichment

    # --- Draw effect sizes ----------------------------------------------------

    beta_true <- rep(0, p_i)

    if (effect_distribution == "normal") {
      beta_true[causal_indices] <- rnorm(S[i], mean = 0, sd = sqrt(effect_variance))
    } else if (effect_distribution == "equal") {
      # Effect sizes set later during variance calibration
      beta_true[causal_indices] <- sample(c(-1, 1), S[i], replace = TRUE)
    }

    # --- Generate phenotype ---------------------------------------------------

    if (model == "sparse") {
      pheno_result <- generate_phenotype_sparse(
        X = X_i,
        beta_true = beta_true,
        phi = phi[i],
        effect_distribution = effect_distribution,
        causal_indices = causal_indices
      )
    } else if (model == "sparse_inf") {
      pheno_result <- generate_phenotype_sparse_inf(
        X = X_i,
        beta_true = beta_true,
        phi = phi[i],
        p_causal = p_causal,
        inf_model = inf_model,
        effect_distribution = effect_distribution,
        causal_indices = causal_indices
      )
    }

    y_i <- pheno_result$y
    beta_true <- pheno_result$beta_true  # may have been rescaled for "equal"

    # --- Compute summary statistics -------------------------------------------

    sumstats <- compute_summary_statistics(X = X_i, y = y_i)

    # --- Store results --------------------------------------------------------

    genotypes[[i]]$y <- y_i
    genotypes[[i]]$z <- sumstats$z
    genotypes[[i]]$beta_hat <- sumstats$beta_hat
    genotypes[[i]]$se <- sumstats$se
    genotypes[[i]]$statistic_type <- sumstats$statistic_type
    genotypes[[i]]$degrees_of_freedom <- sumstats$degrees_of_freedom
    # Only compute LD if not already present (e.g. pre-computed by run_simulation)
    if (is.null(genotypes[[i]]$LD)) {
      genotypes[[i]]$LD <- cor(X_i)
    }
    genotypes[[i]]$annotations_matrix <- A_i

    genotypes[[i]]$truth <- list(
      causal_indices = causal_indices,
      causal_effects = beta_true[causal_indices],
      beta_true = beta_true,  # legacy field: sparse-component effects
      beta_sparse = beta_true,
      beta_infinitesimal = pheno_result$beta_infinitesimal,
      beta_total = beta_true + pheno_result$beta_infinitesimal,
      estimand = "sparse_component_membership",
      pve = pheno_result$pve_realised,
      pve_expected = pheno_result$pve_expected,
      variance_components = pheno_result$variance_components,
      S = S[i],
      phi = phi[i],
      model = model,
      p_causal = if (model == "sparse_inf") p_causal else NULL,
      inf_model = if (model == "sparse_inf") inf_model else NULL,
      effect_distribution = effect_distribution,
      effect_variance = effect_variance,
      annotation_type = annotation_type,
      enrichment = enrichment_used,
      annotation_proportions = props_i
    )

    if (verbose) {
      message(sprintf(
        "  Causal indices: {%s}, realised PVE: %.4f",
        paste(causal_indices, collapse = ", "),
        pheno_result$pve_realised
      ))
    }
  }

  if (verbose) {
    message("Phenotype simulation complete.")
  }

  # --- Save to disk (optional) ------------------------------------------------

  if (save) {
    stopifnot(
      "output_dir must be a single character string" =
        is.character(output_dir) && length(output_dir) == 1L
    )
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    seed_tag <- if (!is.null(seed)) paste0("seed", seed) else "noseed"
    S_tag    <- paste(sort(unique(S)),   collapse = "-")
    phi_tag  <- paste(sort(unique(phi)), collapse = "-")
    fname    <- sprintf("phenotypes_%s_S%s_phi%s_%s.rds",
                        model, S_tag, phi_tag, seed_tag)
    fpath    <- file.path(output_dir, fname)
    saveRDS(genotypes, file = fpath)
    if (verbose) message(sprintf("Phenotypes saved to: %s", fpath))
  }

  # --- Return -----------------------------------------------------------------

  genotypes
}


# =============================================================================
# Internal: validate a parameter that can be scalar or per-region vector
# =============================================================================

validate_per_region_param <- function(x, n_regions, name, integer_valued = FALSE) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    stop(sprintf("%s must be nonempty and finite numeric.", name), call. = FALSE)
  }
  if (integer_valued && any(x != floor(x))) {
    stop(sprintf("%s must be integer-valued.", name), call. = FALSE)
  }
  if (length(x) == 1) {
    x <- rep(x, n_regions)
  } else if (length(x) != n_regions) {
    stop(
      sprintf(
        "If %s is a vector, it must have length n_regions (%d). Got length %d.",
        name, n_regions, length(x)
      ),
      call. = FALSE
    )
  }
  if (integer_valued) as.integer(x) else as.numeric(x)
}


# =============================================================================
# Internal: validate annotation parameter (scalar or vector of length m)
# =============================================================================

validate_annotation_param <- function(x, n_annotations, name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    stop(sprintf("%s must be nonempty and finite numeric.", name), call. = FALSE)
  }
  if (length(x) == 1) {
    x <- rep(x, n_annotations)
  } else if (length(x) != n_annotations) {
    stop(
      sprintf(
        "%s must be a scalar or vector of length n_annotations (%d). Got length %d.",
        name, n_annotations, length(x)
      ),
      call. = FALSE
    )
  }
  as.numeric(x)
}


# =============================================================================
# Internal: simulate annotations for a single region
# =============================================================================

simulate_annotations_for_region <- function(p,
                                            annotation_type,
                                            n_annotations,
                                            annotation_proportions,
                                            user_annotation_matrix,
                                            annotation_correlation = 0,
                                            enrichment = NULL) {

  if (annotation_type == "none") {
    return(list(matrix = NULL, proportions = NULL))
  }

  if (annotation_type == "user_supplied") {
    if (nrow(user_annotation_matrix) != p) {
      stop(
        sprintf(
          "User-supplied annotation matrix has %d rows but region has %d SNPs.",
          nrow(user_annotation_matrix), p
        ),
        call. = FALSE
      )
    }
    return(list(matrix = user_annotation_matrix, proportions = NULL))
  }

  if (annotation_type == "binary") {
    # Determine proportions
    if (is.null(annotation_proportions)) {
      props <- runif(n_annotations, min = 0.01, max = 0.30)
    } else {
      props <- annotation_proportions
    }

    A <- .binary_annotations(p, n_annotations, props,
                             annotation_correlation, enrichment)
    colnames(A) <- paste0("annot_", seq_len(n_annotations))
    return(list(matrix = A, proportions = props))
  }

  if (annotation_type == "continuous") {
    # Draw from N(0, 1) following Funmap
    A <- matrix(rnorm(p * n_annotations), nrow = p, ncol = n_annotations)
    colnames(A) <- paste0("annot_", seq_len(n_annotations))

    return(list(matrix = A, proportions = NULL))
  }
}


# =============================================================================
# Internal: binary annotation matrix with optional within-enrichment-group
# correlation.
#
# When annotation_correlation == 0 (default) OR the enrichment grouping is
# unknown (NULL / wrong length), columns are drawn independently via rbinom
# - recovers the original behaviour exactly.
#
# When annotation_correlation > 0 and enrichment is supplied, columns are
# grouped by identical enrichment fold. Within each group of size >= 2 we
# generate a latent Gaussian field with compound-symmetric correlation rho
# via a shared-factor construction:
#
#     Z_k = sqrt(rho) * F + sqrt(1 - rho) * eps_k,  F, eps_k ~ N(0, 1) iid
#
# giving marginal N(0,1) and pairwise cor(Z_k, Z_l) = rho for k != l in the
# same group; independent across groups. Each column is then thresholded at
# qnorm(1 - prop_k) so its marginal P(A_k = 1) is exactly prop_k.
#
# NB: thresholding attenuates the induced Bernoulli correlation relative to
# the latent rho - this is expected, documented in the tests, and swept as
# an axis by the autoresearch simulation grid rather than solved for
# analytically.
.binary_annotations <- function(p, n_annotations, props,
                                annotation_correlation, enrichment) {
  A <- matrix(0L, nrow = p, ncol = n_annotations)

  # Enrichment-based grouping: identical enrichment fold -> same group.
  # Fall back to "each column is its own group" (i.e. independent) when
  # enrichment isn't supplied or when the requested correlation is 0.
  if (annotation_correlation == 0 || is.null(enrichment) ||
      length(enrichment) != n_annotations) {
    groups <- as.list(seq_len(n_annotations))
  } else {
    groups <- split(seq_len(n_annotations), enrichment)
  }

  rho <- annotation_correlation
  for (grp in groups) {
    if (length(grp) == 1L || rho == 0) {
      for (k in grp) {
        A[, k] <- rbinom(p, size = 1L, prob = props[k])
      }
    } else {
      F <- rnorm(p)   # shared factor for this enrichment group
      for (k in grp) {
        eps_k    <- rnorm(p)
        latent_k <- sqrt(rho) * F + sqrt(1 - rho) * eps_k
        # Threshold at qnorm(1 - prop_k) => P(latent > cut) = prop_k
        A[, k]   <- as.integer(latent_k > stats::qnorm(1 - props[k]))
      }
    }
  }
  A
}


# =============================================================================
# Internal: select causal variants
# =============================================================================

select_causal_variants <- function(p,
                                   S,
                                   annotation_matrix,
                                   enrichment,
                                   n_annotations,
                                   annotation_type) {

  # --- No annotations: uniform random selection -------------------------------

  if (annotation_type == "none" || is.null(annotation_matrix)) {
    causal_indices <- sort(sample(p, S))
    return(list(causal_indices = causal_indices, enrichment = NULL))
  }

  # --- With annotations: weighted selection -----------------------------------

  m <- ncol(annotation_matrix)

  # Determine enrichment values
  if (is.null(enrichment)) {
    enrichment_vals <- runif(m, min = 2, max = 10)
  } else {
    enrichment_vals <- enrichment
  }

  # Compute unnormalised weights: w_j = exp(sum_k A_jk * log(enrichment_k))
  log_enrichment <- log(enrichment_vals)
  log_weights <- as.numeric(annotation_matrix %*% log_enrichment)
  weights <- exp(log_weights - max(log_weights))  # subtract max for numerical stability

  # Normalise to probabilities
  probs <- weights / sum(weights)

  # Sample S causal indices without replacement
  causal_indices <- sort(sample(p, S, replace = FALSE, prob = probs))

  return(list(causal_indices = causal_indices, enrichment = enrichment_vals))
}


# =============================================================================
# Internal: generate phenotype under sparse model
# =============================================================================

generate_phenotype_sparse <- function(X, beta_true, phi, effect_distribution,
                                      causal_indices) {

  n <- nrow(X)

  # Compute genetic signal
  g_sparse <- as.numeric(X %*% beta_true)
  var_g <- var(g_sparse)

  # Handle "equal" effect distribution: effects have random signs and unit magnitude
  # magnitudes); the PVE target is achieved by calibrating sigma^2, not by
  # rescaling effects. A new sign draw is tried only if Var(Xb) == 0.
  if (effect_distribution == "equal") {
    S <- length(causal_indices)
    # Set effects so that Var(Xb) gives a reasonable signal,
    # then calibrate via sigma^2
    # Start with unit effects, compute variance, then proceed
    if (var_g == 0) {
      # All effects are 1 (placeholder), but Var(Xb) = 0 is unlikely
      # with standardised genotypes. If it happens, redraw.
      warning("Var(Xb) = 0 with equal effects. Trying random signs.", call. = FALSE)
      beta_true[causal_indices] <- sample(c(-1, 1), S, replace = TRUE)
      g_sparse <- as.numeric(X %*% beta_true)
      var_g <- var(g_sparse)
    }
  }

  # Check for degenerate case
  if (var_g < .Machine$double.eps) {
    stop(
      "Var(Xb) is effectively zero. The drawn effect sizes produced no signal. ",
      "This can happen with very small S or very small effect_variance. ",
      "Try increasing effect_variance or S.",
      call. = FALSE
    )
  }

  # Calibrate residual variance to achieve target PVE
  # phi = Var(Xb) / (Var(Xb) + sigma^2)
  # => sigma^2 = Var(Xb) * (1 - phi) / phi
  sigma2 <- var_g * (1 - phi) / phi

  # Draw phenotype
  e <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  y <- g_sparse + e

  ledger <- phenotype_variance_ledger(g_sparse, rep(0, n), e, sigma2)

  list(
    y = y,
    beta_true = beta_true,
    beta_infinitesimal = rep(0, ncol(X)),
    pve_realised = ledger$pve_realised,
    pve_expected = ledger$pve_expected,
    variance_components = ledger$components
  )
}


# =============================================================================
# Internal: generate phenotype under sparse + infinitesimal model
# =============================================================================

generate_phenotype_sparse_inf <- function(X, beta_true, phi, p_causal,
                                          inf_model, effect_distribution,
                                          causal_indices) {

  n <- nrow(X)
  p <- ncol(X)
  S <- length(causal_indices)

  # --- Sparse component -------------------------------------------------------

  g_sparse <- as.numeric(X %*% beta_true)
  var_g_sparse <- var(g_sparse)

  # Only retry with random signs when the user picked "equal" effects — for
  # normal effects, zero variance means the draws cancelled exactly (vanishing
  # probability) or S = 0, neither of which is fixed by reshuffling signs.
  if (effect_distribution == "equal" && var_g_sparse < .Machine$double.eps) {
    warning("Var(Xb) = 0 with equal effects. Trying random signs.", call. = FALSE)
    beta_true[causal_indices] <- sample(c(-1, 1), S, replace = TRUE)
    g_sparse <- as.numeric(X %*% beta_true)
    var_g_sparse <- var(g_sparse)
  }

  if (var_g_sparse < .Machine$double.eps) {
    stop(
      "Var(Xb) is effectively zero for the sparse component. ",
      "Try increasing effect_variance or S.",
      call. = FALSE
    )
  }

  # Retain all coefficients so sparse membership and nonzero-effect truth
  # are distinguishable. No background is required at the p_causal = 1 boundary.
  alpha <- rep(0, p)
  if (p_causal < 1) {
    inf_indices <- if (inf_model == "beatrice") {
      setdiff(seq_len(p), causal_indices)
    } else seq_len(p)
    if (length(inf_indices) == 0L) {
      stop("No variants are available for the requested infinitesimal component.", call. = FALSE)
    }
    alpha[inf_indices] <- rnorm(length(inf_indices), sd = 1 / sqrt(length(inf_indices)))
  }
  g_inf <- as.numeric(X %*% alpha)
  var_g_inf <- var(g_inf)
  if (p_causal < 1 && var_g_inf < .Machine$double.eps) {
    stop("Var(g_inf) is effectively zero for the infinitesimal component.", call. = FALSE)
  }

  sparse_scale <- sqrt(p_causal * phi / var_g_sparse)
  inf_scale <- if (p_causal < 1) sqrt((1 - p_causal) * phi / var_g_inf) else 0
  beta_true <- beta_true * sparse_scale
  alpha <- alpha * inf_scale
  g_sparse_scaled <- g_sparse * sparse_scale
  g_inf_scaled <- g_inf * inf_scale
  g_total <- g_sparse_scaled + g_inf_scaled
  var_g_total <- var(g_total)
  if (var_g_total < .Machine$double.eps) {
    stop("The sparse and infinitesimal signals cancel to zero variance.", call. = FALSE)
  }

  # LD induces covariance between the two signals. Calibrate using Var(g_total),
  # preserving this covariance rather than assuming the variances add.
  sigma2 <- var_g_total * (1 - phi) / phi
  e <- rnorm(n, sd = sqrt(sigma2))
  y <- g_total + e
  ledger <- phenotype_variance_ledger(g_sparse_scaled, g_inf_scaled, e, sigma2)
  list(
    y = y,
    beta_true = beta_true,
    beta_infinitesimal = alpha,
    pve_realised = ledger$pve_realised,
    pve_expected = ledger$pve_expected,
    variance_components = ledger$components
  )
}

# Finite-sample signal/noise covariance means Var(g)/Var(y) need not be <= 1.
# The expected PVE conditions on the drawn genetic signal and uses the noise
# distribution's variance; neither quantity is the sample R-squared.
phenotype_variance_ledger <- function(g_sparse, g_infinitesimal, residual,
                                      residual_variance) {
  g <- g_sparse + g_infinitesimal
  vg <- var(g)
  vy <- var(g + residual)
  list(
    pve_expected = if (vg == 0) 0 else vg / (vg + residual_variance),
    pve_realised = if (vg == 0) 0 else vg / vy,
    components = list(
      sparse = var(g_sparse), infinitesimal = var(g_infinitesimal),
      sparse_infinitesimal_covariance = stats::cov(g_sparse, g_infinitesimal),
      genetic = vg, residual_expected = residual_variance,
      residual_realised = var(residual), genetic_residual_covariance = stats::cov(g, residual),
      phenotype = vy
    )
  )
}


# =============================================================================
# Internal: compute marginal summary statistics
# =============================================================================

compute_summary_statistics <- function(X, y) {
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) < 3L || ncol(X) < 1L ||
      length(y) != nrow(X) || any(!is.finite(X)) || any(!is.finite(y))) {
    stop("Summary statistics require finite X and y with at least three individuals.", call. = FALSE)
  }
  n <- nrow(X)
  # Explicit centering is essential when standardise = FALSE: the model
  # includes an intercept whether genotypes are dosages or standardised.
  X_centered <- sweep(X, 2, colMeans(X), "-")
  y_centered <- y - mean(y)
  xty <- as.numeric(crossprod(X_centered, y_centered))
  xtx <- colSums(X_centered * X_centered)
  if (any(xtx <= 0)) stop("Genotype columns must have nonzero variance.", call. = FALSE)
  beta_hat <- xty / xtx
  rss <- sum(y_centered * y_centered) - beta_hat * xty
  # Numerical cancellation near a perfect fit can produce tiny negative RSS.
  rss <- pmax(rss, 0)
  se <- sqrt(rss / (n - 2) / xtx)
  z <- beta_hat / se
  list(beta_hat = beta_hat, se = se, z = z,
       statistic_type = "marginal_ols_t", degrees_of_freedom = n - 2L)
}
