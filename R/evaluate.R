# =============================================================================
# evaluate.R
#
# Evaluation metrics for the fine-mapping benchmark.
#
# Usage:
#   eval_out <- evaluate_methods(simulation, results)
#   eval_out <- evaluate_methods(simulation, results, save = TRUE, output_dir = "results/my_run")
#
# Output structure:
#   eval_out$<method>$global          — metrics pooled across all scenarios
#   eval_out$<method>$by_S            — named list, one entry per S value
#   eval_out$<method>$by_phi          — named list, one entry per phi value
#   eval_out$<method>$by_p_causal     — named list (sparse_inf only); NULL otherwise
#   eval_out$<method>$by_region_size  — named list, one entry per nominal region
#                                       length (e.g. "100","200","400","500",
#                                       "1000"). Region size is a real scenario
#                                       axis — fine-mapping difficulty scales with
#                                       it — so metrics must never be pooled
#                                       across lengths.
#   eval_out$<method>$by_causal_maf   — named list, keys "rare" / "low" / "common"
#                                       (bins on the minimum causal-variant MAF
#                                       per region; NULL if no fits have causal
#                                       MAFs available)
#   eval_out$<method>$by_true_annotation_type — named list, keys "none" /
#                                       "binary" / "continuous" / "user_supplied"
#                                       (true annotation regime used in
#                                       simulation; scenario-level, so typically
#                                       single-bin per evaluate_methods() call
#                                       and becomes meaningful only when eval
#                                       objects from differently-annotated sims
#                                       are merged)
#
# Each stratum contains:
#   fdr_power_curve  — data.frame: threshold, tp, fp, fn, fdr, power, precision, recall
#   auprc            — scalar (area under precision-recall curve)
#   pip_calibration  — data.frame: bin, bin_lower, bin_upper, bin_mid, n, n_causal,
#                                  mean_pip, frac_causal
#   cs_coverage      — proportion of reported CSs containing ≥1 true causal variant
#   cs_power         — proportion of true causal variants captured by any CS
#   cs_size_median   — median variants per CS
#   cs_size_mean     — mean variants per CS
#   n_cs_reported    — total CSs reported
#   runtime_mean     — mean runtime in seconds (successful fits only)
#   runtime_sd       — SD of runtime
#   n_fits           — total fits attempted
#   n_failed         — fits that errored
# =============================================================================


#' Evaluate fine-mapping methods against ground truth
#'
#' @param simulation  Output of \code{run_simulation()}.
#' @param results     Output of \code{run_methods()}.
#' @param pip_thresholds Numeric vector of PIP thresholds for the power/FDR
#'   curve. Default: seq(0, 1, by = 0.005).
#' @param n_pip_cal_bins Integer. Number of equal-width bins for PIP
#'   calibration. Default: 10.
#' @param save Logical. If TRUE, save the full evaluation object as
#'   \code{evaluation.rds} inside \code{output_dir}, and write a flat
#'   per-method summary table as \code{evaluation_summary.csv}.
#'   Default: FALSE.
#' @param output_dir Character. Directory in which to save results when
#'   \code{save = TRUE}. Created automatically if it does not exist.
#'   Default: \code{"results"}.
#' @param verbose Logical. Print progress. Default: TRUE.
#'
#' @return A named list, one element per method, each containing
#'   \code{global}, \code{by_S}, \code{by_phi}, and \code{by_p_causal}
#'   sub-lists.  See file header for the fields in each stratum.
#'
#' @export
evaluate_methods <- function(simulation,
                             results,
                             pip_thresholds = seq(0, 1, by = 0.005),
                             n_pip_cal_bins = 10L,
                             save           = FALSE,
                             output_dir     = "results",
                             verbose        = TRUE) {

  stopifnot(
    "simulation must contain 'scenarios'" = !is.null(simulation$scenarios),
    "results must contain 'methods_run'"  = !is.null(results$methods_run)
  )

  if (!is.numeric(pip_thresholds) || !length(pip_thresholds) ||
      any(!is.finite(pip_thresholds)) || any(pip_thresholds < 0 | pip_thresholds > 1))
    stop("pip_thresholds must be finite probabilities in [0, 1].", call. = FALSE)
  if (!is.numeric(n_pip_cal_bins) || length(n_pip_cal_bins) != 1L ||
      !is.finite(n_pip_cal_bins) || n_pip_cal_bins < 1 ||
      n_pip_cal_bins != as.integer(n_pip_cal_bins))
    stop("n_pip_cal_bins must be a positive integer.", call. = FALSE)

  methods <- intersect(results$methods_run, names(results))
  if (length(methods) == 0) {
    stop("No recognised methods found in results.", call. = FALSE)
  }

  output <- list()

  for (method in methods) {
    if (verbose) message(sprintf("Evaluating %s...", method))

    fits_raw <- .complete_fit_grid(results[[method]]$results, simulation, method)
    fits     <- .annotate_fits_with_truth(fits_raw, simulation)

    output[[method]] <- list(
      global      = .metrics_with_se(fits, pip_thresholds, n_pip_cal_bins),
      by_S                    = .stratify_metrics(fits, "S",           pip_thresholds, n_pip_cal_bins),
      by_phi                  = .stratify_metrics(fits, "phi",         pip_thresholds, n_pip_cal_bins),
      by_p_causal             = .stratify_metrics(fits, "p_causal",    pip_thresholds, n_pip_cal_bins),
      by_region_size          = .stratify_metrics(fits, "region_size", pip_thresholds, n_pip_cal_bins),
      by_causal_maf           = .stratify_metrics_by_maf(fits, pip_thresholds, n_pip_cal_bins),
      by_true_annotation_type = .stratify_metrics_by_annotation_type(fits, pip_thresholds, n_pip_cal_bins)
    )
  }

  output$methods_evaluated   <- methods
  output$simulation_params   <- simulation$params
  output$pip_thresholds_used <- pip_thresholds
  output$evaluation_version <- "2.0"
  output$uncertainty <- "delete-one-iter cluster jackknife, conditional on simulated loci"

  # --- Save to disk -----------------------------------------------------------

  if (save) {
    stopifnot(
      "output_dir must be a single character string" =
        is.character(output_dir) && length(output_dir) == 1
    )
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    # Full evaluation object
    rds_path <- file.path(output_dir, "evaluation.rds")
    saveRDS(output, file = rds_path)
    if (verbose) message(sprintf("  Saved evaluation object: %s", rds_path))

    # Flat summary table — one row per method, global metrics only
    summary_rows <- lapply(methods, function(m) {
      g <- output[[m]]$global
      data.frame(
        method              = m,
        auprc               = g$auprc,
        auprc_se            = g$auprc_se,
        ap                  = g$ap,
        ap_se               = g$ap_se,
        ap_macro            = g$ap_macro,
        ap_macro_se         = g$ap_macro_se,
        ap_macro_attempted  = g$ap_macro_attempted,
        ap_macro_attempted_se = g$ap_macro_attempted_se,
        success_rate        = g$success_rate,
        n_invalid           = g$n_invalid,
        n_missing           = g$n_missing,
        brier_score         = g$brier_score,
        log_loss            = g$log_loss,
        cs_power_attempted  = g$cs_power_attempted,
        cs_matched_precision = g$cs_matched_precision,
        cs_matched_power    = g$cs_matched_power,
        cs_duplicate_sets   = g$cs_duplicate_sets,
        null_region_false_call_rate = g$null_region_false_call_rate,
        null_region_cs_rate = g$null_region_cs_rate,
        cs_coverage         = g$cs_coverage,
        cs_coverage_se      = g$cs_coverage_se,
        cs_power            = g$cs_power,
        cs_power_se         = g$cs_power_se,
        cs_size_median      = g$cs_size_median,
        cs_size_median_se   = g$cs_size_median_se,
        cs_size_mean        = g$cs_size_mean,
        cs_size_mean_se     = g$cs_size_mean_se,
        n_cs_reported       = g$n_cs_reported,
        runtime_mean        = g$runtime_mean,
        runtime_mean_se     = g$runtime_mean_se,
        runtime_sd          = g$runtime_sd,
        n_fits              = g$n_fits,
        n_failed            = g$n_failed,
        stringsAsFactors    = FALSE
      )
    })
    summary_df  <- do.call(rbind, summary_rows)
    csv_path    <- file.path(output_dir, "evaluation_summary.csv")
    write.csv(summary_df, file = csv_path, row.names = FALSE)
    if (verbose) message(sprintf("  Saved summary table:     %s", csv_path))
  }

  output
}


# =============================================================================
# Internal: annotate fits with ground truth
# =============================================================================

# Complete the simulation design before scoring: a missing output is an
# unsuccessful attempt, not permission to drop a hard locus from the denominator.
.complete_fit_grid <- function(fits, simulation, method) {
  expected <- unlist(lapply(seq_along(simulation$scenarios), function(sc) {
    paste(sc, seq_along(simulation$scenarios[[sc]]$regions), sep = ":")
  }), use.names = FALSE)
  keys <- vapply(fits, function(f) {
    ids <- c(f$scenario_id, f$region_id)
    if (length(ids) != 2L || !is.numeric(ids) || any(!is.finite(ids)) ||
        any(ids != floor(ids)) || any(ids < 1L))
      stop("Each result needs valid scenario_id and region_id indices.", call. = FALSE)
    paste(ids, collapse = ":")
  }, character(1))
  if (anyDuplicated(keys))
    stop("Duplicate scenario/region result: attempts cannot be counted twice.", call. = FALSE)
  if (any(!keys %in% expected))
    stop("Result references a scenario/region outside the simulation.", call. = FALSE)
  lapply(expected, function(key) {
    j <- match(key, keys)
    if (!is.na(j)) return(fits[[j]])
    ids <- as.integer(strsplit(key, ":", fixed = TRUE)[[1L]])
    sc <- simulation$scenarios[[ids[1L]]]
    list(scenario_id = ids[1L], region_id = ids[2L], method = method,
         S = sc$S, phi = sc$phi, iter = sc$iter, p_causal = sc$p_causal,
         pip = numeric(), credible_sets = list(), runtime_seconds = NA_real_,
         error = "Missing method result", missing_result = TRUE)
  })
}

# Invalid method outputs become explicit failures. Never silently clip PIPs,
# remove nonfinite variants or infer the variant count from a truncated output.
.validate_evaluation_fit <- function(f) {
  if (!is.null(f$error)) return(f)
  issue <- NULL
  if (!is.numeric(f$pip) || !is.null(dim(f$pip)) ||
      length(f$pip) != f$n_variants || any(!is.finite(f$pip)) ||
      any(f$pip < 0 | f$pip > 1)) {
    issue <- "PIPs must be a finite [0,1] vector matching the simulated variant count"
  } else if (!is.null(f$credible_sets) && !is.list(f$credible_sets)) {
    issue <- "credible_sets must be a list of variant-index vectors"
  } else if (any(vapply(f$credible_sets, function(cs) {
    !is.numeric(cs) || !length(cs) || any(!is.finite(cs)) ||
      any(cs != floor(cs) | cs < 1 | cs > f$n_variants)
  }, logical(1)))) {
    issue <- "Every reported credible set needs nonempty in-range integer indices"
  }
  if (!is.null(issue)) {
    f$validation_error <- issue
    f$error <- paste("Invalid method output:", issue)
  }
  f
}

.annotate_fits_with_truth <- function(fits, simulation) {
  lapply(fits, function(f) {
    sc    <- simulation$scenarios[[f$scenario_id]]
    truth <- sc$regions[[f$region_id]]$truth
    f$causal_indices <- truth$causal_indices
    # The variant universe comes from the simulation, not the submitted PIPs.
    geno <- simulation$genotypes[[f$region_id]]
    f$n_variants <- if (!is.null(geno$LD)) nrow(geno$LD) else if (!is.null(geno$X)) {
      ncol(geno$X)
    } else if (!is.null(geno$p)) as.integer(geno$p) else if (!is.null(truth$beta_true)) {
      length(truth$beta_true)
    } else if (!is.null(simulation$params$p_actual)) {
      as.integer(simulation$params$p_actual[f$region_id])
    } else length(f$pip)
    if (length(f$n_variants) != 1L || !is.finite(f$n_variants) || f$n_variants < 1L ||
        any(!is.finite(f$causal_indices)) || any(f$causal_indices != floor(f$causal_indices)) ||
        any(f$causal_indices < 1L | f$causal_indices > f$n_variants) ||
        anyDuplicated(f$causal_indices))
      stop("Simulation has an invalid variant universe or causal indices.", call. = FALSE)
    # Scenario labels are authoritative; submitted labels cannot move a hard fit
    # into a more favourable stratum or change the uncertainty clusters.
    for (nm in c("S", "phi", "p_causal", "iter"))
      if (!is.null(sc[[nm]])) f[[nm]] <- sc[[nm]]

    # Nominal region length (the requested per-region p, NOT the realised
    # post-MAF-filter count in params$p_actual). Fine-mapping difficulty scales
    # steeply with region size, so this is a genuine scenario axis: metrics must
    # be stratified by it, never pooled across region lengths. Indexed by
    # region_id; falls back to p_actual then the realised variant count when the
    # nominal vector is unavailable (e.g. legacy simulations).
    p_nominal <- simulation$params$p
    rid       <- f$region_id
    f$region_size <- if (!is.null(p_nominal) && length(p_nominal) >= rid &&
                         !is.na(p_nominal[rid])) {
      as.integer(p_nominal[rid])
    } else {
      p_actual <- simulation$params$p_actual
      if (!is.null(p_actual) && length(p_actual) >= rid && !is.na(p_actual[rid]))
        as.integer(p_actual[rid]) else length(f$pip)
    }

    # Attach causal-variant MAFs and a coarse stratification bin so the
    # evaluator can stratify performance by how rare the rarest causal
    # variant is in the region (the rarest causal is the bottleneck for
    # fine-mapping). MAFs live on genotypes (shared across scenarios
    # for the same region).
    geno_maf <- simulation$genotypes[[f$region_id]]$maf
    if (!is.null(geno_maf) && length(truth$causal_indices) > 0L) {
      f$causal_maf     <- geno_maf[truth$causal_indices]
      f$min_causal_maf <- min(f$causal_maf)
    } else {
      f$causal_maf     <- numeric(0)
      f$min_causal_maf <- NA_real_
    }
    f$causal_maf_bin <- .maf_bin(f$min_causal_maf)

    # True annotation type used by the simulation. Scenario-level (the
    # whole sim was generated with one annotation regime), but attached
    # to each fit so the evaluator can stratify across simulations that
    # have been merged into one results object.
    f$true_annotation_type <- if (!is.null(simulation$params$annotation_type)) {
      as.character(simulation$params$annotation_type)
    } else {
      NA_character_
    }
    .validate_evaluation_fit(f)
  })
}


# =============================================================================
# Internal: bin a MAF value into "rare" / "low" / "common"
# =============================================================================

# Bins follow the standard human-genetics partition. `NA_real_` input
# (e.g. no causal variants in a region) maps to NA_character_ so the
# fit is excluded from MAF stratification rather than counted in any bin.
.maf_bin <- function(maf) {
  if (is.null(maf) || length(maf) == 0L || is.na(maf)) return(NA_character_)
  if (maf <= 0.01)  return("rare")
  if (maf <= 0.05)  return("low")
  "common"
}


# =============================================================================
# Internal: stratify fits by a scenario-level variable and compute metrics
# =============================================================================

.stratify_metrics <- function(fits, by_var, pip_thresholds, n_pip_cal_bins) {
  values <- sort(unique(sapply(fits, function(f) {
    v <- f[[by_var]]
    if (is.null(v)) NA_real_ else v
  })))
  values <- values[!is.na(values)]

  if (length(values) == 0) return(NULL)

  result <- lapply(values, function(v) {
    subset <- Filter(function(f) {
      fv <- f[[by_var]]
      !is.null(fv) && !is.na(fv) && fv == v
    }, fits)
    .metrics_with_se(subset, pip_thresholds, n_pip_cal_bins)
  })
  names(result) <- as.character(values)
  result
}


# =============================================================================
# Internal: stratify fits by causal-variant MAF bin
#
# Groups fits by `causal_maf_bin` (attached in .annotate_fits_with_truth) and
# computes metrics per bin. Bins are kept in the canonical order
# "rare" -> "low" -> "common" rather than alphabetical, so downstream plots
# read left-to-right by increasing MAF.
#
# Returns NULL when no fit has a valid bin assignment (e.g. simulation with
# no causal variants assigned to any region, or genotypes without an attached
# `maf` field).
# =============================================================================

.stratify_metrics_by_maf <- function(fits, pip_thresholds, n_pip_cal_bins) {
  bin_order <- c("rare", "low", "common")

  present <- unique(vapply(fits, function(f) {
    b <- f$causal_maf_bin
    if (is.null(b) || length(b) == 0L) NA_character_ else as.character(b[[1L]])
  }, character(1L)))
  present <- present[!is.na(present)]
  present <- bin_order[bin_order %in% present]

  if (length(present) == 0L) return(NULL)

  result <- lapply(present, function(b) {
    subset <- Filter(function(f) {
      bv <- f$causal_maf_bin
      !is.null(bv) && !is.na(bv) && bv == b
    }, fits)
    .metrics_with_se(subset, pip_thresholds, n_pip_cal_bins)
  })
  names(result) <- present
  result
}


# =============================================================================
# Internal: stratify fits by the true annotation type used in simulation
#
# Groups fits by `true_annotation_type` (attached in
# .annotate_fits_with_truth from simulation$params$annotation_type) and
# computes metrics per type. Kept in canonical order
# "none" -> "binary" -> "continuous" -> "user_supplied" so plots read from
# the misspecified annotation-null case through to richer annotation
# regimes.
#
# Within a single evaluate_methods() call this typically produces a
# single bin, because annotation_type is scenario-level. The axis
# becomes meaningful when evaluation objects from differently-annotated
# simulations are merged.
# =============================================================================

.stratify_metrics_by_annotation_type <- function(fits, pip_thresholds,
                                                  n_pip_cal_bins) {
  type_order <- c("none", "binary", "continuous", "user_supplied")

  present <- unique(vapply(fits, function(f) {
    v <- f$true_annotation_type
    if (is.null(v) || length(v) == 0L) NA_character_ else as.character(v[[1L]])
  }, character(1L)))
  present <- present[!is.na(present)]
  present <- c(type_order[type_order %in% present],
               setdiff(present, type_order))  # tolerate unknown labels

  if (length(present) == 0L) return(NULL)

  result <- lapply(present, function(t) {
    subset <- Filter(function(f) {
      tv <- f$true_annotation_type
      !is.null(tv) && !is.na(tv) && tv == t
    }, fits)
    .metrics_with_se(subset, pip_thresholds, n_pip_cal_bins)
  })
  names(result) <- present
  result
}


# =============================================================================
# Internal: compute all metrics for a collection of annotated fits
# =============================================================================

.metrics_for_fits <- function(fits, pip_thresholds, n_pip_cal_bins) {
  n_fits <- length(fits)
  valid <- Filter(function(f) is.null(f$error), fits)
  total_causal <- sum(vapply(fits, function(f) length(f$causal_indices), integer(1)))
  n_failed <- n_fits - length(valid)
  labels <- function(f) {
    y <- logical(f$n_variants); y[f$causal_indices] <- TRUE; y
  }
  pip_vals <- as.numeric(unlist(lapply(valid, function(f) f$pip), use.names = FALSE))
  is_causal <- as.logical(unlist(lapply(valid, labels), use.names = FALSE))
  curve <- if (length(valid)) .fdr_power_curve(pip_vals, is_causal, pip_thresholds) else NULL
  attempted <- .fdr_power_curve(pip_vals, is_causal, pip_thresholds,
                                total_causal = total_causal)
  ap_fit <- vapply(valid, function(f) .compute_ap_exact(f$pip, labels(f)), numeric(1))
  signal_fits <- sum(vapply(fits, function(f) length(f$causal_indices) > 0L, logical(1)))
  ap_macro <- if (any(!is.na(ap_fit))) mean(ap_fit, na.rm = TRUE) else NA_real_
  cal <- if (length(valid)) .pip_calibration(pip_vals, is_causal, n_pip_cal_bins) else NULL
  cs <- .cs_metrics(valid)
  cs_attempted <- .cs_metrics(lapply(fits, function(f) {
    if (!is.null(f$error)) f$credible_sets <- list()
    f
  }))
  rts <- vapply(valid, function(f) {
    rt <- f$runtime_seconds
    if (is.numeric(rt) && length(rt) == 1L && is.finite(rt) && rt >= 0) rt else NA_real_
  }, numeric(1))
  rts <- rts[!is.na(rts)]
  null_fits <- Filter(function(f) !length(f$causal_indices), valid)
  # The null operating point is predeclared and reported alongside its denominator.
  null_any_call <- vapply(null_fits, function(f) any(f$pip >= 0.5), logical(1))
  invalid <- Filter(function(f) !is.null(f$validation_error), fits)
  validation_errors <- if (length(invalid)) do.call(rbind, lapply(invalid, function(f) {
    data.frame(scenario_id = f$scenario_id, region_id = f$region_id,
               error = f$validation_error, stringsAsFactors = FALSE)
  })) else data.frame(scenario_id = integer(), region_id = integer(), error = character())
  y <- as.integer(is_causal)
  brier <- if (length(y)) mean((pip_vals - y)^2) else NA_real_
  prevalence <- if (length(y)) mean(y) else NA_real_
  log_loss <- if (length(y)) {
    # No clipping: a confidently impossible observed event has infinite loss.
    loss <- numeric(length(y))
    loss[y == 1L] <- -log(pip_vals[y == 1L])
    loss[y == 0L] <- -log1p(-pip_vals[y == 0L])
    mean(loss)
  } else NA_real_
  list(
    fdr_power_curve = curve,
    fdr_power_curve_attempted = attempted,
    auprc = if (is.null(curve)) NA_real_ else .compute_auprc(curve$precision, curve$recall),
    ap = .compute_ap_exact(pip_vals, is_causal),
    ap_macro = ap_macro,
    ap_macro_attempted = if (signal_fits) sum(ap_fit, na.rm = TRUE) / signal_fits else NA_real_,
    ap_n_fits = sum(!is.na(ap_fit)),
    pip_calibration = cal,
    brier_score = brier,
    brier_skill_score = if (!is.na(prevalence) && prevalence > 0 && prevalence < 1)
      1 - brier / (prevalence * (1 - prevalence)) else NA_real_,
    log_loss = log_loss,
    ece = if (is.null(cal)) NA_real_ else
      sum(cal$n * abs(cal$mean_pip - cal$frac_causal), na.rm = TRUE) / sum(cal$n),
    cs_coverage = cs$coverage,
    cs_power = cs$power,
    cs_power_attempted = cs_attempted$power,
    cs_matched_precision = cs$matched_precision,
    cs_matched_power = cs$matched_power,
    cs_duplicate_sets = cs$duplicate_sets,
    cs_size_median = cs$size_median,
    cs_size_mean = cs$size_mean,
    n_cs_reported = cs$n_cs,
    runtime_mean = if (length(rts)) mean(rts) else NA_real_,
    runtime_sd = if (length(rts) > 1L) stats::sd(rts) else NA_real_,
    n_fits = n_fits,
    n_failed = n_failed,
    n_invalid = length(invalid),
    n_missing = sum(vapply(fits, function(f) isTRUE(f$missing_result), logical(1))),
    success_rate = if (n_fits) length(valid) / n_fits else NA_real_,
    validation_errors = validation_errors,
    n_null_fits = length(null_fits),
    null_region_false_call_rate = if (length(null_fits)) mean(null_any_call) else NA_real_,
    null_region_cs_rate = if (length(null_fits))
      mean(vapply(null_fits, function(f) length(f$credible_sets) > 0L, logical(1))) else NA_real_
  )
}

# =============================================================================
# Internal: delete-one-replicate cluster jackknife for the actual aggregate
# estimator. Variants in LD and fits sharing one iter are never independent
# resampling units. Undefined delete-one statistics yield an undefined SE.
# =============================================================================

.metrics_with_se <- function(fits, pip_thresholds, n_pip_cal_bins) {
  main <- .metrics_for_fits(fits, pip_thresholds, n_pip_cal_bins)
  scalar_fields <- c("auprc", "ap", "ap_macro", "ap_macro_attempted",
                     "cs_coverage", "cs_power", "cs_power_attempted",
                     "cs_matched_precision", "cs_matched_power",
                     "cs_size_median", "cs_size_mean", "runtime_mean",
                     "brier_score", "ece", "success_rate",
                     "null_region_false_call_rate", "null_region_cs_rate")
  cluster <- vapply(fits, function(f) {
    if (is.null(f$iter) || length(f$iter) != 1L || is.na(f$iter)) NA_character_ else as.character(f$iter)
  }, character(1))
  iter_vals <- sort(unique(cluster[!is.na(cluster)]))
  K <- length(iter_vals)
  main$n_replicate_clusters <- K
  main$se_method <- "delete-one-iter cluster jackknife"
  main$se_conditioning <- "simulated loci/genotypes and scenario grid fixed"
  # Missing cluster identifiers make the assumed resampling unit ambiguous.
  can_estimate <- K >= 2L && !anyNA(cluster)
  deleted <- if (can_estimate) lapply(iter_vals, function(it) {
    .metrics_for_fits(fits[cluster != it], pip_thresholds, n_pip_cal_bins)
  }) else list()
  .jackknife_se <- function(x) {
    # Dropping undefined leave-one-out statistics would change the estimand.
    if (length(x) != K || K < 2L || any(!is.finite(x))) return(NA_real_)
    sqrt((K - 1) / K * sum((x - mean(x))^2))
  }
  for (fld in scalar_fields) {
    vals <- vapply(deleted, function(m) as.numeric(m[[fld]]), numeric(1))
    main[[paste0(fld, "_se")]] <- .jackknife_se(vals)
  }
  for (name in c("fdr_power_curve", "fdr_power_curve_attempted")) {
    if (is.null(main[[name]])) next
    for (fld in c("power", "precision")) {
      n <- nrow(main[[name]])
      if (can_estimate) {
        mat <- matrix(vapply(deleted, function(m) {
          if (is.null(m[[name]])) rep(NA_real_, n) else m[[name]][[fld]]
        }, numeric(n)), nrow = n)
        main[[name]][[paste0(fld, "_se")]] <- apply(mat, 1L, .jackknife_se)
      } else main[[name]][[paste0(fld, "_se")]] <- rep(NA_real_, n)
    }
  }
  if (!is.null(main$pip_calibration)) {
    n <- nrow(main$pip_calibration)
    if (can_estimate) {
      mat <- matrix(vapply(deleted, function(m) {
        if (is.null(m$pip_calibration)) rep(NA_real_, n) else m$pip_calibration$frac_causal
      }, numeric(n)), nrow = n)
      main$pip_calibration$frac_causal_se <- apply(mat, 1L, .jackknife_se)
    } else main$pip_calibration$frac_causal_se <- rep(NA_real_, n)
  }
  main
}

# =============================================================================
# Internal: power / FDR curve
#
# Uses a sort + cumsum approach: O(N log N + T) rather than O(N * T).
# Pools all variants across all fits in the supplied collection.
# =============================================================================

.fdr_power_curve <- function(pip_vals, is_causal, thresholds,
                             total_causal = sum(is_causal)) {

  N            <- length(pip_vals)

  if (N == 0) {
    return(data.frame(
      threshold = thresholds,
      tp = 0L, fp = 0L, fn = total_causal,
      fdr = 0, power = if (total_causal > 0) 0 else NA_real_,
      precision = 1, recall = if (total_causal > 0) 0 else NA_real_
    ))
  }

  # Sort by decreasing PIP; build cumulative TP and selection counts
  ord          <- order(pip_vals, decreasing = TRUE)
  sorted_pips  <- pip_vals[ord]
  cum_tp       <- c(0L, cumsum(as.integer(is_causal[ord])))  # length N+1

  # For threshold t, n_selected = #{pip >= t}
  # sorted_pips is non-increasing => -sorted_pips is non-decreasing
  # findInterval(-t, -sorted_pips) = largest k such that -sorted_pips[k] <= -t
  #                                = largest k such that sorted_pips[k] >= t
  #                                = n_selected(t)
  neg_sorted <- -sorted_pips
  n_sel      <- findInterval(-thresholds, neg_sorted)
  n_sel      <- pmin(pmax(n_sel, 0L), N)

  tp <- cum_tp[n_sel + 1L]
  fp <- n_sel - tp
  fn <- total_causal - tp

  data.frame(
    threshold = thresholds,
    tp        = tp,
    fp        = fp,
    fn        = fn,
    fdr       = ifelse(n_sel > 0L, fp / n_sel,          0),
    power     = if (total_causal > 0) tp / total_causal else NA_real_,
    precision = ifelse(n_sel > 0L, tp / n_sel,          1),
    recall    = if (total_causal > 0) tp / total_causal else NA_real_
  )
}


# =============================================================================
# Internal: AUPRC via the trapezoidal rule
# =============================================================================

.compute_auprc <- function(precision, recall) {
  if (length(recall) < 2 || any(!is.finite(recall)) ||
      any(!is.finite(precision))) return(NA_real_)

  # Sort by increasing recall; at recall ties keep the higher precision
  ord   <- order(recall, -precision)
  prec  <- precision[ord]
  rec   <- recall[ord]

  # Trapezoidal rule: sum of trapezoids.
  # Legacy threshold-grid trapezoid retained for numerical comparison only.
  # Its interpolation and missing endpoints can bias in either direction;
  # use exact AP for new ranking comparisons.
  sum(diff(rec) * (head(prec, -1) + tail(prec, -1)) / 2)
}


#' Average precision on an already supplied precision-recall grid
#'
#' This compatibility helper is exact only when the grid contains every
#' distinct score operating point. evaluate_methods() uses .compute_ap_exact().
#'
#' @param precision,recall Numeric vectors from the FDR/power curve.
#' @return Scalar average precision, or NA if the curve is degenerate.
#' @keywords internal
.compute_ap <- function(precision, recall) {
  if (length(recall) < 2 || any(!is.finite(recall)) ||
      any(!is.finite(precision))) return(NA_real_)
  ord  <- order(recall, -precision)
  prec <- precision[ord]
  rec  <- recall[ord]
  # sum over steps: (R_i - R_{i-1}) * P_i, with R_0 = 0
  sum(diff(c(0, rec)) * prec)
}


# =============================================================================
# Internal: PIP calibration
#
# Bins all (pip, is_causal) pairs into n_bins equal-width intervals on [0,1].
# Returns mean_pip (expected) vs frac_causal (observed) per bin.
# =============================================================================

.pip_calibration <- function(pip_vals, is_causal, n_bins = 10L) {

  breaks  <- seq(0, 1, length.out = n_bins + 1L)
  bin_ids <- findInterval(pip_vals, breaks, rightmost.closed = TRUE)
  bin_ids <- pmin(pmax(bin_ids, 1L), n_bins)

  ns          <- tabulate(bin_ids, nbins = n_bins)
  n_causal    <- tapply(as.integer(is_causal), bin_ids, sum)
  n_causal    <- as.integer(n_causal[as.character(seq_len(n_bins))])
  n_causal[is.na(n_causal)] <- 0L

  # Mean PIP per bin (only over non-empty bins)
  mean_pip <- vapply(seq_len(n_bins), function(b) {
    idx <- bin_ids == b
    if (!any(idx)) NA_real_ else mean(pip_vals[idx])
  }, numeric(1))

  frac_causal <- ifelse(ns > 0L, n_causal / ns, NA_real_)

  data.frame(
    bin         = seq_len(n_bins),
    bin_lower   = breaks[-length(breaks)],
    bin_upper   = breaks[-1],
    bin_mid     = (breaks[-length(breaks)] + breaks[-1]) / 2,
    n           = ns,
    n_causal    = n_causal,
    mean_pip    = mean_pip,
    frac_causal = frac_causal,
    sum_pip = vapply(seq_len(n_bins), function(b) sum(pip_vals[bin_ids == b]), numeric(1)),
    sum_pip_sq = vapply(seq_len(n_bins), function(b) sum(pip_vals[bin_ids == b]^2), numeric(1)),
    sum_causal_pip = vapply(seq_len(n_bins), function(b)
      sum(pip_vals[bin_ids == b] * is_causal[bin_ids == b]), numeric(1))
  )
}


# =============================================================================
# Internal: credible set metrics
#
# Coverage : proportion of reported CSs containing ≥1 true causal variant
# Power    : proportion of true causal variants captured by any CS
# Size     : median and mean number of variants per CS
# =============================================================================

.cs_metrics <- function(valid_fits) {
  cs_hits <- logical()
  cs_sizes <- integer()
  causal_captured <- logical()
  n_matched <- 0L
  duplicate_sets <- 0L
  for (f in valid_fits) {
    # Count distinct variants inside each set; retain repeated reported sets in
    # the denominator so duplicate reporting cannot improve matched precision.
    sets <- lapply(f$credible_sets, function(cs) sort(unique(as.integer(cs))))
    keys <- vapply(sets, function(cs) paste(cs, collapse = ","), character(1))
    duplicate_sets <- duplicate_sets + sum(duplicated(keys))
    all_cs_variants <- unique(unlist(sets))
    cs_hits <- c(cs_hits, vapply(sets, function(cs) any(cs %in% f$causal_indices), logical(1)))
    cs_sizes <- c(cs_sizes, lengths(sets))
    causal_captured <- c(causal_captured, f$causal_indices %in% all_cs_variants)
    n_matched <- n_matched + .maximum_cs_matching(sets, f$causal_indices)
  }
  n_cs <- length(cs_hits)
  n_causal <- length(causal_captured)
  list(coverage = if (n_cs) mean(cs_hits) else NA_real_,
       power = if (n_causal) mean(causal_captured) else NA_real_,
       matched_precision = if (n_cs) n_matched / n_cs else NA_real_,
       matched_power = if (n_causal) n_matched / n_causal else NA_real_,
       duplicate_sets = duplicate_sets,
       size_median = if (n_cs) stats::median(cs_sizes) else NA_real_,
       size_mean = if (n_cs) mean(cs_sizes) else NA_real_, n_cs = n_cs)
}

# Maximum bipartite matching: each reported set may recover one distinct causal
# and each causal may be credited once. A greedy matching is not sufficient when
# overlapping sets contain multiple signals.
.maximum_cs_matching <- function(sets, causal_indices) {
  if (!length(sets) || !length(causal_indices)) return(0L)
  edges <- lapply(sets, function(cs) which(causal_indices %in% cs))
  owner <- integer(length(causal_indices))
  augment <- function(i) {
    for (j in edges[[i]]) {
      if (seen[j]) next
      seen[j] <<- TRUE
      if (owner[j] == 0L || augment(owner[j])) {
        owner[j] <<- i
        return(TRUE)
      }
    }
    FALSE
  }
  matched <- 0L
  for (i in seq_along(sets)) {
    seen <- rep(FALSE, length(causal_indices))
    if (augment(i)) matched <- matched + 1L
  }
  matched
}
