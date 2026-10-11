# Publication adapter for the separately specified stage analysis. Callers
# establish analysis custody; this adapter never reads predictions or truth.
.paper_stage_publication <- function(plan, rows, specification_path, scientific = FALSE,
                                     collection_complete = FALSE, max_output_rows = 100000L) {
  if (!is.logical(scientific) || length(scientific) != 1L || is.na(scientific) ||
      !is.logical(collection_complete) || length(collection_complete) != 1L || is.na(collection_complete) ||
      length(max_output_rows) != 1L || !is.numeric(max_output_rows) ||
      !is.finite(max_output_rows) || max_output_rows < 1L || max_output_rows != floor(max_output_rows))
    stop("Invalid stage-publication status or output bound")
  spec <- jsonlite::fromJSON(specification_path)
  stages <- c("original_refit", "corrected_refit", "corrected_resimulation")
  metrics <- c("availability", "ap", "ap_observed_positive_attempts_zero_for_failure",
               "brier", "recall", "cs_distinct_recall")
  pooled <- c("fb_pooled", "fb_xregion", "polyfun_est", "polyfun_ldsc", "sbayesrc", "paintor")
  definitions <- data.frame(stage_a = stages[c(2L, 3L)], stage_b = stages[c(1L, 2L)],
    type = c("implementation_and_acceptance", "generating_experiment"),
    fixed = c("observed_inputs_and_truth", "genotypes_and_design"), stringsAsFactors = FALSE)
  valid_sha <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && grepl("^[a-f0-9]{64}$", x)
  if (!is.list(spec) || !identical(spec$schema, 1L) ||
      !identical(spec$analysis, "regional_stage_attribution_v1") ||
      !identical(spec$registered_before_full_campaign_outcomes, TRUE) ||
      !identical(spec$development_results_already_seen, TRUE) ||
      !valid_sha(spec$regional_full_plan_sha256) || !valid_sha(plan$plan_sha256) ||
      !identical(spec$metrics, metrics) || !identical(spec$pooled_methods, pooled) ||
      !identical(spec$pair_keys, c("dataset_id", "region_id", "method")) ||
      !identical(spec$direction, "stage_a_minus_stage_b") ||
      !is.data.frame(spec$contrasts) || !all(names(definitions) %in% names(spec$contrasts)) ||
      !identical(spec$contrasts[names(definitions)], definitions) ||
      !identical(spec$bootstrap$draws, 1000L) || !identical(spec$bootstrap$seed, 20261011L) ||
      !identical(spec$bootstrap$interval, "percentile 0.025/0.975") ||
      !identical(spec$limits$max_stratum_rows, 100000L) ||
      !identical(spec$limits$no_all_memory_full_campaign, TRUE))
    stop("Stage-publication specification differs from the implemented analysis contract")
  if (!is.data.frame(rows) || !all(c("stage", "stratum", "method") %in% names(rows)) ||
      anyNA(rows[c("stage", "stratum", "method")]) ||
      anyDuplicated(rows[c("stage", "stratum", "method")])) stop("Invalid selected stage-summary identities")
  registered <- identical(plan$plan_sha256, spec$regional_full_plan_sha256)
  if (scientific && !identical(plan$profile, "full")) stop("Scientific stage analysis requires the full profile")
  allowed <- registered || !scientific
  tables <- decisions <- list(); output_rows <- 0L; seen <- character()
  callback <- function(d, stratum) {
    if (length(stratum) != 1L || is.na(stratum) || stratum %in% seen ||
        !nrow(d) || anyNA(d$stratum) || any(d$stratum != stratum)) stop("Repeated or inconsistent stage-publication partition")
    seen <<- c(seen, stratum)
    gaussian <- all(d$family == "challenge")
    for (method in sort(unique(d$method), method = "radix")) {
      selected <- rows$stage[rows$stratum == stratum & rows$method == method]
      status <- if (gaussian) "gaussian_has_no_implementation_stage_pair" else if (!allowed)
        "unregistered_full_plan_not_analysed" else if (!all(stages %in% selected))
        "all_three_stages_not_selected; unassessed" else "analysed_selected_three_stage_pairs"
      decisions[[length(decisions) + 1L]] <<- data.frame(stratum = stratum, method = method,
        selected_stages = paste(sort(selected, method = "radix"), collapse = ";"), status = status,
        registered_full_plan = registered, scientific_full_profile = scientific, stringsAsFactors = FALSE)
      if (status != "analysed_selected_three_stage_pairs") next
      z <- d[d$method == method & d$stage %in% stages, , drop = FALSE]
      x <- .paper_stage_contrasts(z, n_boot = spec$bootstrap$draws,
        seed = spec$bootstrap$seed, max_stratum_rows = spec$limits$max_stratum_rows)
      if (collection_complete && any(x$n_unmatched_stage_a != 0L | x$n_unmatched_stage_b != 0L))
        stop("Complete stage collection has unmatched attempted locus identities")
      output_rows <<- output_rows + nrow(x)
      if (output_rows > max_output_rows) stop("Stage-publication output memory bound exceeded")
      x$registered_full_plan <- registered; x$scientific_full_profile <- scientific
      x$analysis_status <- if (scientific) "registered_full_profile_stage_analysis" else
        "development_or_incomplete_stage_analysis; no_full_campaign_claim"
      tables[[length(tables) + 1L]] <<- x
    }
    invisible(NULL)
  }
  finish <- function() list(
    contrasts = if (length(tables)) do.call(rbind, tables) else data.frame(),
    selection = if (length(decisions)) do.call(rbind, decisions) else data.frame(),
    metadata = list(available = TRUE, analysis = spec$analysis,
      registered_full_plan = registered, scientific_full_profile = scientific,
      development_results_already_seen = spec$development_results_already_seen,
      n_strata_visited = length(seen), n_contrast_metric_rows = output_rows,
      max_output_rows = max_output_rows, bootstrap = spec$bootstrap,
      interval_nominal_coverage_validated = FALSE,
      interval_interpretation = "Descriptive percentile envelopes; few independent genotype clusters cannot establish statistical distinguishability",
      interpretation = "Separate stage analysis of stored scores; original method contrasts and intervals are unchanged"))
  list(on_stratum = callback, finish = finish)
}
