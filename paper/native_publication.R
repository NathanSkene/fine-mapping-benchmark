# Publication annotations only. Collected scores and intervals are never changed.
.paper_native_keys <- c("draw_id", "model", "pi", "h2", "exposure", "method", "event", "pip_event")

.paper_native_key <- function(x, keys = .paper_native_keys) {
  if (!all(keys %in% names(x))) stop("Native publication view lacks identity fields")
  if (anyNA(x[keys]) || any(vapply(x[keys], function(v) any(grepl("\037", v, fixed = TRUE)), TRUE)))
    stop("Invalid native publication identity")
  do.call(paste, c(x[keys], sep = "\037"))
}

.paper_native_support <- function(n_loci, n_all_positive, n_all_negative, n_mixed) {
  ifelse(is.na(n_loci) | is.na(n_all_positive) | is.na(n_all_negative) | is.na(n_mixed) |
    n_loci == 0, "unassessed", ifelse(n_all_positive == n_loci, "all_positive",
    ifelse(n_all_negative == n_loci, "all_negative", ifelse(n_mixed > 0,
      "contains_within_locus_positive_and_negative_targets", "positive_and_negative_loci_only"))))
}

.paper_native_localisation <- function(n_positive_loci, n_mixed_loci) {
  ifelse(is.na(n_positive_loci) | is.na(n_mixed_loci), "unassessed",
    ifelse(n_positive_loci == 0, "no_positive_scored_loci", ifelse(n_mixed_loci == 0,
      "no_localisation_discrimination", ifelse(n_mixed_loci < n_positive_loci,
        "mix_of_degenerate_and_discriminating_loci", "within_locus_discrimination_available"))))
}

.paper_native_unassessed <- function(rows) {
  counts <- c("loci", "variants", "positive", "negative", "all_positive_loci", "all_negative_loci",
    "mixed_loci", "successful_loci", "complete_loci", "successful_positive_loci",
    "successful_mixed_loci", "complete_positive_loci", "complete_mixed_loci", "misaligned_loci")
  for (name in counts) rows[[paste0("census_n_", name)]] <- NA_real_
  rows$target_support <- rows$ap_localisation_status <- rows$recall_localisation_status <- "unassessed"
  rows$all_positive_target <- NA
  rows$constant_one_ap <- rows$constant_one_recall <- rows$constant_one_brier <- rows$constant_one_log_loss <- NA_real_
  rows
}

.paper_native_truth_census <- function(path, rows, chunk_size = 10000L, max_groups = 100000L) {
  if (!nrow(rows) || nrow(rows) > max_groups) stop("Native truth census exceeds its group memory bound")
  wanted <- .paper_native_key(rows)
  if (anyDuplicated(wanted)) stop("Repeated native draw-summary identity")
  fields <- c("loci", "variants", "positive", "negative", "all_positive_loci", "all_negative_loci",
    "mixed_loci", "successful_loci", "complete_loci", "successful_positive_loci",
    "successful_mixed_loci", "complete_positive_loci", "complete_mixed_loci", "misaligned_loci")
  totals <- matrix(0, nrow(rows), length(fields), dimnames = list(NULL, fields))
  required <- c(.paper_native_keys, "n_variants", "n_positive", "success", "complete", "probability_target_aligned")
  .paper_read_selected_csv(path, function(chunk) {
    if (!all(required %in% names(chunk))) stop("Sealed native task metrics lack truth-census fields")
    group <- match(.paper_native_key(chunk), wanted)
    x <- chunk[!is.na(group), , drop = FALSE]; group <- group[!is.na(group)]
    if (nrow(x)) {
      if (any(!is.finite(x$n_variants) | x$n_variants <= 0 | x$n_variants != floor(x$n_variants)) ||
          any(!is.finite(x$n_positive) | x$n_positive < 0 | x$n_positive > x$n_variants |
              x$n_positive != floor(x$n_positive))) stop("Invalid sealed native truth counts")
      if (anyNA(x[c("success", "complete", "probability_target_aligned")]) ||
          any(vapply(x[c("success", "complete", "probability_target_aligned")],
                     function(v) any(!v %in% c(TRUE, FALSE)), TRUE))) stop("Invalid native census flags")
      positive <- x$n_positive > 0
      mixed <- positive & x$n_positive < x$n_variants
      values <- cbind(loci = 1, variants = x$n_variants, positive = x$n_positive,
        negative = x$n_variants - x$n_positive,
        all_positive_loci = x$n_positive == x$n_variants, all_negative_loci = !positive,
        mixed_loci = mixed, successful_loci = x$success, complete_loci = x$complete,
        successful_positive_loci = x$success & positive, successful_mixed_loci = x$success & mixed,
        complete_positive_loci = x$complete & positive, complete_mixed_loci = x$complete & mixed,
        misaligned_loci = !x$probability_target_aligned)
      partial <- rowsum(values, group, reorder = FALSE)
      ids <- as.integer(rownames(partial)); totals[ids, ] <<- totals[ids, , drop = FALSE] + partial
    }
    # Only bounded group totals survive a chunk; no locus rows are retained.
    chunk[FALSE, , drop = FALSE]
  }, chunk_size = chunk_size, max_rows = 1L)
  if (any(totals[, "loci"] != rows$n_planned) ||
      any(totals[, "successful_loci"] != rows$n_successful))
    stop("Native truth census disagrees with selected draw-summary denominators")
  if ("n_complete" %in% names(rows) && any(totals[, "complete_loci"] != rows$n_complete))
    stop("Native truth census disagrees with completed denominators")
  if (any((totals[, "misaligned_loci"] == 0) != rows$all_probability_targets_aligned))
    stop("Native truth census disagrees with probability-event alignment")
  for (name in fields) rows[[paste0("census_n_", name)]] <- totals[, name]
  rows$target_support <- .paper_native_support(totals[, "loci"], totals[, "all_positive_loci"],
    totals[, "all_negative_loci"], totals[, "mixed_loci"])
  rows$all_positive_target <- rows$target_support == "all_positive"
  rows$ap_localisation_status <- .paper_native_localisation(totals[, "successful_positive_loci"], totals[, "successful_mixed_loci"])
  rows$recall_localisation_status <- .paper_native_localisation(totals[, "complete_positive_loci"], totals[, "complete_mixed_loci"])
  # Analytical no-data baseline only for the exact all-positive target census.
  rows$constant_one_ap <- rows$constant_one_recall <- ifelse(rows$all_positive_target, 1, NA_real_)
  rows$constant_one_brier <- rows$constant_one_log_loss <- ifelse(rows$all_positive_target, 0, NA_real_)
  rows
}

.paper_native_degenerate <- function(rows, metric) {
  if (metric == "ap_successful") return(rows$ap_localisation_status == "no_localisation_discrimination")
  if (metric == "recall_all_completed") return(rows$recall_localisation_status == "no_localisation_discrimination")
  if (metric == "brier_successful") return(rows$all_positive_target %in% TRUE)
  rep(FALSE, nrow(rows))
}

.paper_native_paired_view <- function(contrasts, rows) {
  if (!nrow(contrasts)) return(contrasts)
  scope <- c("model", "pi", "h2", "event", "method", "exposure")
  group <- .paper_native_key(rows, scope)
  groups <- split(seq_len(nrow(rows)), group)
  count_fields <- grep("^census_n_", names(rows), value = TRUE)
  collapse <- function(x) paste(sort(unique(as.character(x)), method = "radix"), collapse = ";")
  side <- function(x) {
    counts <- lapply(x[count_fields], function(v) if (anyNA(v)) NA_real_ else sum(v))
    support <- .paper_native_support(counts$census_n_loci, counts$census_n_all_positive_loci,
      counts$census_n_all_negative_loci, counts$census_n_mixed_loci)
    ap_status <- .paper_native_localisation(counts$census_n_successful_positive_loci, counts$census_n_successful_mixed_loci)
    recall_status <- .paper_native_localisation(counts$census_n_complete_positive_loci, counts$census_n_complete_mixed_loci)
    if (!nrow(x)) support <- ap_status <- recall_status <- "no_usable_paired_draws"
    result <- c(list(pip_event = if (!nrow(x)) "no_usable_paired_draws" else
      if ("pip_event" %in% names(x)) collapse(x$pip_event) else "not_declared",
      probability_target_aligned = if (!nrow(x) || anyNA(x$aligned)) NA else all(x$aligned),
      posterior_conditioning = if (!nrow(x)) "no_usable_paired_draws" else
        if ("posterior_conditioning" %in% names(x)) collapse(x$posterior_conditioning) else "not_declared",
      native_null_model_in_support = if (!"native_null_model_in_support" %in% names(x) ||
        !nrow(x) || anyNA(x$native_null_model_in_support) || length(unique(x$native_null_model_in_support)) != 1L) NA else x$native_null_model_in_support[1L],
      target_support = support,
      all_positive_target = if (support %in% c("unassessed", "no_usable_paired_draws")) NA else support == "all_positive",
      ap_localisation_status = ap_status, recall_localisation_status = recall_status), counts)
    as.data.frame(result, stringsAsFactors = FALSE)
  }
  metadata <- do.call(rbind, lapply(groups, function(index) side(rows[index, , drop = FALSE])))
  indices <- list()
  for (suffix in c("a", "b")) {
    x <- contrasts[c("model", "pi", "h2", "event", paste0("method_", suffix), paste0("exposure_", suffix))]
    names(x) <- scope
    index <- match(.paper_native_key(x, scope), names(groups))
    if (anyNA(index)) stop("Paired publication view has an unmatched comparison side")
    indices[[suffix]] <- index
    for (name in names(metadata)) contrasts[[paste0("selected_scope_", name, "_", suffix)]] <- metadata[[name]][index]
  }
  paired_metadata <- list(a = list(), b = list())
  known <- all(c("draw_id", "n_complete", "n_planned") %in% names(rows))
  contrasts$paired_support_scope <- if (known) "exact_source_metric_usable_draw_pairs" else "unassessed_legacy_pair_support"
  contrasts$n_publication_usable_draw_pairs <- NA_integer_
  for (i in seq_len(nrow(contrasts))) {
    a <- rows[groups[[indices$a[i]]], , drop = FALSE]
    b <- rows[groups[[indices$b[i]]], , drop = FALSE]
    if (known) {
      if (anyDuplicated(a$draw_id) || anyDuplicated(b$draw_id))
        stop("Paired publication support has ambiguous draw identities")
      metric <- contrasts$metric[i]
      if (!metric %in% names(rows)) stop("Paired publication support lacks the source metric")
      ids <- intersect(a$draw_id, b$draw_id)
      a <- a[match(ids, a$draw_id), , drop = FALSE]; b <- b[match(ids, b$draw_id), , drop = FALSE]
      # Match .native_paired_draws' existing inclusion rule; no output score,
      # difference or interval is recomputed or replaced.
      usable <- is.finite(a[[metric]] - b[[metric]]) & a$n_complete == a$n_planned & b$n_complete == b$n_planned
      if (anyNA(usable)) stop("Paired publication support has undefined completion flags")
      if (sum(usable) != contrasts$n_paired_genotype_draws[i])
        stop("Paired publication support disagrees with the stored paired-draw denominator")
      contrasts$n_publication_usable_draw_pairs[i] <- sum(usable)
      a <- a[usable, , drop = FALSE]; b <- b[usable, , drop = FALSE]
    }
    paired_metadata$a[[i]] <- side(a); paired_metadata$b[[i]] <- side(b)
  }
  for (suffix in c("a", "b")) {
    paired_side <- do.call(rbind, paired_metadata[[suffix]])
    for (name in names(paired_side)) contrasts[[paste0(name, "_", suffix)]] <- paired_side[[name]]
  }
  contrasts$probability_events_both_aligned <- contrasts$probability_target_aligned_a & contrasts$probability_target_aligned_b
  contrasts$all_positive_target_both_sides <- contrasts$all_positive_target_a & contrasts$all_positive_target_b
  # Classify the endpoint, never its measured value or apparent winning method.
  contrasts$localisation_interpretation <- vapply(seq_len(nrow(contrasts)), function(i) {
    metric <- contrasts$metric[i]
    if (!metric %in% c("ap_successful", "recall_all_completed")) return("probability_or_availability_endpoint")
    field <- if (metric == "ap_successful") "ap_localisation_status" else "recall_localisation_status"
    statuses <- c(contrasts[[paste0(field, "_a")]][i], contrasts[[paste0(field, "_b")]][i])
    if (any(statuses == "no_usable_paired_draws")) "no_usable_paired_draws" else
      if (all(statuses == "no_localisation_discrimination")) "no_localisation_discrimination" else
      if (any(grepl("unassessed", statuses, fixed = TRUE))) "truth_census_unassessed" else
      if (any(grepl("degenerate|no_localisation", statuses))) "contains_degenerate_scored_loci" else
        if (any(statuses == "no_positive_scored_loci")) "no_positive_scored_loci_on_at_least_one_side" else
          "within_locus_discrimination_available"
  }, "")
  contrasts$probability_interpretation <- ifelse(contrasts$probability_events_both_aligned %in% TRUE,
    "declared_target_scores; conditioning_retained; no_calibration_rank_inferred",
    "event_misaligned_or_undeclared_operational_scores; no_calibration_rank")
  contrasts$all_positive_no_data_baseline <- ifelse(is.na(contrasts$all_positive_target_both_sides),
    "all_positive_status_unassessed", ifelse(contrasts$all_positive_target_both_sides,
      "constant_PIP_1: AP=1; recall=1; Brier=0; log_loss=0; no_localisation",
      "not_exactly_all_positive_on_both_sides"))
  contrasts
}
