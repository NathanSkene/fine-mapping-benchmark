# Compare versions of the experiment within one frozen regional design stratum.
# Callers must validate the sealed analysis before supplying its metric rows.
# This helper does not refit, rescore, tune or manufacture absent attempts.
.paper_stage_contrasts <- function(d, n_boot = 1000L, seed = 20261011L,
                                   max_stratum_rows = 100000L) {
  scalar_integer <- function(x, lower) length(x) == 1L && is.numeric(x) &&
    is.finite(x) && x == floor(x) && x >= lower && x <= .Machine$integer.max
  if (!scalar_integer(n_boot, 20L) || !scalar_integer(seed, 0L) ||
      !scalar_integer(max_stratum_rows, 1L)) stop("Invalid stage-contrast resource/seed specification")
  required <- c("stage", "stratum", "family", "method", "dataset_id", "region_id",
    "panel_cluster_id", "cluster_id", "success", "actual_p", "nominal_p", "n_gwas", "n_ref",
    "n_causal", "n_observed_causal", "ap", "brier", "recall", "cs_distinct_recall",
    "probability_target_aligned")
  if (!is.data.frame(d) || !all(required %in% names(d)))
    stop("Stage attribution requires the typed locus-level regional metric schema")
  if (nrow(d) > max_stratum_rows) stop("Stage attribution exceeds its explicit stratum memory bound")
  if (!nrow(d)) return(data.frame())
  if (anyNA(d$stratum) || length(unique(d$stratum)) != 1L ||
      anyNA(d$family) || length(unique(d$family)) != 1L || any(d$family == "challenge"))
    stop("Supply exactly one regional design stratum; Gaussian tasks have no implementation-stage pair")
  stages <- c("original_refit", "corrected_refit", "corrected_resimulation")
  if (anyNA(d$stage) || any(!d$stage %in% stages)) stop("Unknown implementation-stage identity")
  ids <- c("method", "dataset_id", "panel_cluster_id", "cluster_id")
  if (any(vapply(d[ids], function(x) !is.character(x) || anyNA(x) || any(!nzchar(x)), TRUE)) ||
      !is.numeric(d$region_id) || anyNA(d$region_id) || any(d$region_id < 1 | d$region_id != floor(d$region_id)))
    stop("Missing or malformed locus/cluster identity")
  if (!is.logical(d$success) || anyNA(d$success) || !is.logical(d$probability_target_aligned))
    stop("Stage attribution requires logical fit/event status")
  counts <- c("actual_p", "nominal_p", "n_gwas", "n_causal", "n_observed_causal")
  if (any(vapply(d[counts], function(x) !is.numeric(x) || any(!is.finite(x)) ||
      any(x < 0 | x != floor(x)), TRUE)) || any(d$actual_p < 1 | d$n_gwas < 1 |
      d$actual_p > d$nominal_p | d$n_observed_causal > d$actual_p |
      d$n_observed_causal > d$n_causal)) stop("Inconsistent retained candidate/truth counts")
  if (!is.numeric(d$n_ref) || any(!is.na(d$n_ref) &
      (!is.finite(d$n_ref) | d$n_ref < 1 | d$n_ref != floor(d$n_ref))))
    stop("Invalid reference sample count")
  if (any(vapply(d[c("ap", "brier", "recall", "cs_distinct_recall")],
      function(x) !is.numeric(x) || any(is.infinite(x)) || any(is.finite(x) & (x < 0 | x > 1)), TRUE)))
    stop("Stage metrics must retain their original probability-score range")
  if (any(d$success & d$n_observed_causal > 0 & !is.finite(d$ap)))
    stop("A successful observed-positive locus lacks assessable AP")
  if (any(!d$success & (is.finite(d$ap) | is.finite(d$brier))))
    stop("A failed locus contains successful-fit AP/Brier")
  if (any(!d$success & d$n_observed_causal > 0 &
      (!is.finite(d$recall) | d$recall != 0 | !is.finite(d$cs_distinct_recall) | d$cs_distinct_recall != 0)))
    stop("An observed-positive failed locus must retain zero attempted recovery")
  keys <- c("dataset_id", "region_id", "method")
  if (anyDuplicated(d[c("stage", keys)])) stop("Duplicate locus identity within an implementation stage")
  pooled <- c("fb_pooled", "fb_xregion", "polyfun_est", "polyfun_ldsc", "sbayesrc", "paintor")
  definitions <- list(
    list(a = "corrected_refit", b = "original_refit", type = "implementation_and_acceptance",
         fixed = "observed_inputs_and_truth", interpretation = paste(
           "Whole implementation/adapter/acceptance bundle on identical supplied inputs;",
           "not attribution to an individual repair without a separate ablation")),
    list(a = "corrected_resimulation", b = "corrected_refit", type = "generating_experiment",
         fixed = "genotypes_and_design", interpretation = paste(
           "Changed phenotype/truth experiment on retained genotypes;",
           "not an inference-only or historical-ranking correction")))
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (old_seed_exists) assign(".Random.seed", old_seed, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  output <- list()
  same <- function(x, y) all((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y))
  for (definition in definitions) for (method in sort(unique(d$method), method = "radix")) {
    a <- d[d$stage == definition$a & d$method == method, , drop = FALSE]
    b <- d[d$stage == definition$b & d$method == method, , drop = FALSE]
    # An absent stage stays an explicit unmatched census, rather than a zero score.
    paired <- merge(a, b, by = keys, suffixes = c("_a", "_b"), sort = TRUE)
    fixed <- c("panel_cluster_id", "cluster_id", "actual_p", "nominal_p", "n_gwas", "n_ref")
    if (definition$type == "implementation_and_acceptance")
      fixed <- c(fixed, "n_causal", "n_observed_causal", intersect("n_any_nonzero_total_effect", names(d)))
    for (field in fixed) if (!same(paired[[paste0(field, "_a")]], paired[[paste0(field, "_b")]]))
      stop("Stage pair changes a declared fixed input/truth field: ", field)
    cluster_field <- if (method %in% pooled) "panel_cluster_id" else "cluster_id"
    metrics <- c("availability", "ap", "ap_observed_positive_attempts_zero_for_failure",
                 "brier", "recall", "cs_distinct_recall")
    for (metric in metrics) {
      if (metric == "availability") {
        va <- as.numeric(paired$success_a); vb <- as.numeric(paired$success_b)
      } else if (metric == "ap_observed_positive_attempts_zero_for_failure") {
        va <- ifelse(paired$n_observed_causal_a > 0, ifelse(paired$success_a, paired$ap_a, 0), NA_real_)
        vb <- ifelse(paired$n_observed_causal_b > 0, ifelse(paired$success_b, paired$ap_b, 0), NA_real_)
      } else { va <- paired[[paste0(metric, "_a")]]; vb <- paired[[paste0(metric, "_b")]] }
      good <- is.finite(va - vb)
      delta <- (va - vb)[good]
      cluster <- paired[[paste0(cluster_field, "_a")]][good]
      blocks <- split(delta, cluster)
      interval <- rep(NA_real_, 2L)
      if (length(blocks) >= 2L) {
        # Fixed seed and sorted merge make results independent of input row order.
        set.seed(as.integer(seed))
        samples <- replicate(as.integer(n_boot), mean(unlist(blocks[
          sample.int(length(blocks), length(blocks), replace = TRUE)], use.names = FALSE)))
        interval <- stats::quantile(samples, c(.025, .975), names = FALSE)
      }
      aligned <- if (metric == "brier" && any(good)) all(
        c(paired$probability_target_aligned_a[good], paired$probability_target_aligned_b[good]) %in% TRUE) else NA
      output[[length(output) + 1L]] <- data.frame(
        stratum = d$stratum[1L], method = method, contrast = definition$type,
        stage_a = definition$a, stage_b = definition$b, difference_direction = "stage_a_minus_stage_b",
        fixed = definition$fixed, metric = metric,
        difference = if (length(delta)) mean(delta) else NA_real_, lower = interval[1L], upper = interval[2L],
        n_stage_a = nrow(a), n_stage_b = nrow(b), n_attempted_pairs = nrow(paired), n_pairs = sum(good),
        n_unmatched_stage_a = nrow(a) - nrow(paired), n_unmatched_stage_b = nrow(b) - nrow(paired),
        n_success_both = sum(paired$success_a & paired$success_b),
        n_failed_both = sum(!paired$success_a & !paired$success_b),
        n_recovered_stage_a = sum(paired$success_a & !paired$success_b),
        n_lost_stage_a = sum(!paired$success_a & paired$success_b),
        n_clusters = length(blocks), cluster_field = cluster_field,
        n_boot = as.integer(n_boot), bootstrap_seed = as.integer(seed),
        interval_nominal_coverage_validated = FALSE,
        interval_coverage_status = if (length(blocks) >= 2L)
          "descriptive_percentile_envelope; nominal_coverage_not_validated" else "no_interval",
        few_cluster_warning = if (length(blocks) <= 2L)
          "at_most_two_independent_declared_clusters; no_statistical_distinguishability_claim" else NA_character_,
        probability_targets_aligned = aligned,
        inference = if (length(blocks) >= 2L) "descriptive_paired_cluster_percentile_envelope; conditional_on_supplied_genotype_design" else
          "point_difference_only; fewer_than_two_independent_declared_clusters",
        interpretation = paste(definition$interpretation, if (metric == "brier" && !isTRUE(aligned))
          "Event/conditioning alignment unqualified: operational Brier diagnostic only." else ""),
        stringsAsFactors = FALSE)
    }
  }
  result <- do.call(rbind, output); rownames(result) <- NULL
  result
}
