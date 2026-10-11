# Publication metadata only: never replace a collected score or interval.
.paper_regional_key <- function(x, fields = c("stage", "stratum", "method")) {
  if (!all(fields %in% names(x)) || anyNA(x[fields]) ||
      any(vapply(x[fields], function(v) any(!nzchar(v) | grepl("\037", v, fixed = TRUE)), TRUE)))
    stop("Invalid regional publication identity")
  do.call(paste, c(x[fields], sep = "\037"))
}

.paper_regional_schema <- function(path) {
  s <- readRDS(path)
  if (!is.list(s) || !identical(s$schema, 1L) || !is.character(s$columns) ||
      anyNA(s$columns) || anyDuplicated(s$columns) || !length(s$columns) ||
      !is.character(s$storage_classes) || length(s$storage_classes) != length(s$columns) ||
      any(!s$storage_classes %in% c("character", "logical", "integer", "double")) ||
      !is.character(s$missing_marker) || length(s$missing_marker) != 1L ||
      is.na(s$missing_marker) || !nzchar(s$missing_marker)) stop("Invalid regional metric CSV schema")
  if (!is.null(names(s$storage_classes)) && !identical(names(s$storage_classes), s$columns))
    stop("Regional metric schema column identities disagree")
  names(s$storage_classes) <- s$columns
  s
}

.paper_regional_scan <- function(path, schema, callback, chunk_size = 10000L) {
  if (length(chunk_size) != 1L || !is.finite(chunk_size) || chunk_size < 1L || chunk_size != floor(chunk_size))
    stop("Invalid regional CSV chunk memory bound")
  con <- file(path, "rt", encoding = "UTF-8"); on.exit(close(con), add = TRUE)
  first <- TRUE; n_chunks <- n_rows <- peak <- 0L
  repeat {
    options <- list(file = con, header = first, colClasses = "character", nrows = chunk_size,
      stringsAsFactors = FALSE, check.names = FALSE, na.strings = schema$missing_marker)
    if (!first) options$col.names <- schema$columns
    x <- do.call(utils::read.csv, options)
    if (!identical(names(x), schema$columns)) stop("Regional metric CSV header differs from sealed schema")
    first <- FALSE
    if (!nrow(x)) break
    for (name in names(x)) {
      token <- x[[name]]; kind <- schema$storage_classes[[name]]
      if (kind == "character") next
      if (kind == "logical") {
        if (any(!is.na(token) & !token %in% c("TRUE", "FALSE"))) stop("Invalid regional logical token")
        x[[name]] <- as.logical(token)
      } else {
        value <- suppressWarnings(as.numeric(token))
        if (any(is.na(value) & !is.na(token) & token != "NaN")) stop("Invalid regional numeric token")
        if (kind == "integer") {
          if (any(!is.na(value) & (!is.finite(value) | value != floor(value) | abs(value) > .Machine$integer.max)))
            stop("Invalid regional integer token")
          value <- as.integer(value)
        }
        x[[name]] <- value
      }
    }
    callback(x); n_chunks <- n_chunks + 1L; n_rows <- n_rows + nrow(x); peak <- max(peak, nrow(x))
    if (nrow(x) < chunk_size) break
  }
  list(chunks = n_chunks, rows_scanned = n_rows, peak_chunk_rows = peak)
}

.paper_regional_normalize <- function(x) {
  required <- c("dataset_id", "family", "stage", "stratum", "method", "success", "n_causal",
    "n_observed_causal", "ap", "brier", "recall", "cs_distinct_recall", "log_loss", "score_target",
    "pip_event", "probability_target_aligned", "posterior_conditioning")
  if (!all(required %in% names(x))) stop("Regional metric census lacks required fields")
  take <- function(a, b = NULL) {
    z <- if (a %in% names(x)) x[[a]] else rep(NA_real_, nrow(x))
    if (!is.null(b) && b %in% names(x)) {
      other <- x[[b]]
      if (any(!is.na(z) & !is.na(other) & z != other)) stop("Contradictory regional candidate counts")
      z[is.na(z)] <- other[is.na(z)]
    }
    z
  }
  challenge <- x$family == "challenge"
  x$.actual <- take("actual_p", "n_variants")
  x$.nominal <- take("nominal_p")
  if ("p" %in% names(x)) {
    if (any(challenge & !is.na(x$.nominal) & !is.na(x$p) & x$.nominal != x$p))
      stop("Contradictory Gaussian requested variant counts")
    x$.nominal[challenge] <- x$p[challenge]
  }
  x$.n <- take("n_gwas")
  if ("n" %in% names(x)) {
    if (any(challenge & !is.na(x$.n) & !is.na(x$n) & x$.n != x$n))
      stop("Contradictory Gaussian synthetic sample counts")
    x$.n[challenge] <- x$n[challenge]
  }
  x$.n_ref <- take("n_ref")
  x$.n_basis <- ifelse(challenge, "gaussian_score_design_parameter; no_individual_rows_simulated",
    "actual_simulated_gwas_individual_count")
  # The task-metric CSV records panel$config$n, not its variant-wise vector.
  # These are exact generator track names; no vector or distribution is
  # reconstructed from them for publication.
  scenario <- if ("scenario" %in% names(x)) x$scenario else rep(NA_character_, nrow(x))
  x$.variable_n_track <- challenge & scenario %in% c("meta_missingness", "meta_known_overlap")
  x$.variant_n_status <- ifelse(challenge, "unavailable_in_locus_metric_csv",
    "not_applicable; actual_simulated_sample_count")
  x$.n_evidence_route <- ifelse(challenge,
    "plan dataset-to-generation mapping -> datasets/generation_<generation_id>/public/inputs.rds; input$n and input$n_per_variant; verify public manifest",
    "sealed simulation genotype$n and params$n_ref; generation dimension audit")
  x$.missing <- if ("n_missing_causal" %in% names(x)) x$n_missing_causal else rep(NA_real_, nrow(x))
  x$.missing[!challenge & is.na(x$.missing)] <- 0L
  x$.nonzero <- take("n_any_nonzero_total_effect")
  x$.null <- x$n_causal == 0
  if ("null" %in% names(x) && any(!is.na(x$null) & x$null != x$.null)) stop("Contradictory generating-null labels")
  integer_count <- function(v, minimum = 0, allow_na = FALSE) {
    if ((!allow_na && anyNA(v)) || any(!is.na(v) & (!is.finite(v) | v < minimum | v != floor(v))))
      stop("Invalid regional truth, size or sample count")
  }
  for (v in list(x$.actual, x$.nominal)) integer_count(v, 1)
  integer_count(x$.n, 2); integer_count(x$.n_ref, 2, TRUE)
  for (v in list(x$n_causal, x$n_observed_causal, x$.missing)) integer_count(v)
  integer_count(x$.nonzero, 0, TRUE)
  if (any(x$.actual > x$.nominal | x$n_observed_causal > x$.actual |
          x$n_causal != x$n_observed_causal + x$.missing) ||
      any(!is.na(x$.nonzero) & x$.nonzero > x$.actual)) stop("Contradictory regional truth or candidate universe")
  if (!is.logical(x$success) || anyNA(x$success)) stop("Invalid regional success flags")
  if (!is.logical(x$probability_target_aligned))
    stop("Invalid regional probability alignment")
  if (anyNA(x$score_target) || any(x$score_target != "sparse_component_membership"))
    stop("Unknown regional primary probability target")
  for (name in c("ap", "brier", "recall", "cs_distinct_recall"))
    if (any(!is.na(x[[name]]) & (!is.finite(x[[name]]) | x[[name]] < 0 | x[[name]] > 1)))
      stop("Invalid regional bounded endpoint")
  if (any(x$success & x$n_observed_causal > 0 & !is.finite(x$ap)) ||
      any(x$n_observed_causal == 0 & !is.na(x$ap)) || any(!x$success & (!is.na(x$ap) | !is.na(x$brier))) ||
      any(x$n_causal == 0 & !is.na(x$recall)) ||
      any(!x$success & x$n_causal > 0 & (is.na(x$recall) | x$recall != 0)))
    stop("Regional endpoint availability contradicts scoring semantics")
  # Validate stratum labels against the collector's exact design rule. Size
  # remains nominal here; realized size is separate publication metadata.
  legacy_fields <- c("n_gwas", "n_ref", "annotation_type", "enrichment_fold", "S", "phi", "p_causal", "nominal_p")
  if (any(!challenge) && !all(legacy_fields %in% names(x))) stop("Regional legacy design fields missing")
  expected <- rep(NA_character_, nrow(x))
  if (any(!challenge)) expected[!challenge] <- do.call(paste, c(x[!challenge, c("family", legacy_fields)], sep = ":"))
  if (any(challenge)) {
    if (!all(c("scenario", "p", "n") %in% names(x))) stop("Gaussian design fields missing")
    expected[challenge] <- paste("challenge", x$scenario[challenge], x$p[challenge], x$n[challenge], sep = ":")
  }
  if (anyNA(expected) || any(expected != x$stratum)) stop("Regional stratum disagrees with design metadata")
  region <- if ("region_id" %in% names(x)) x$region_id else rep(NA_integer_, nrow(x))
  task <- if ("task_id" %in% names(x)) x$task_id else rep(NA_character_, nrow(x))
  integer_count(region[!challenge], 1)
  if (anyNA(task[challenge]) || any(!nzchar(task[challenge]))) stop("Gaussian task identity missing")
  x$.unit <- paste(x$dataset_id, ifelse(challenge, paste0("task:", task), paste0("region:", region)), sep = "\037")
  .paper_regional_key(x, c("dataset_id", "stage", "stratum", "method"))
  x$.size_basis <- ifelse(challenge, "public_gaussian_candidate_subset; gap_not_attributed_to_monomorphism",
    "retained_legacy_candidates; nominal_requested_design")
  x
}

.paper_regional_describe <- function(x, event = "sparse_component_membership") {
  positive <- if (event == "any_nonzero_total_effect") x$.nonzero else x$n_observed_causal
  if (anyNA(positive)) stop("Dense-event endpoint lacks aligned total-effect truth counts")
  distribution <- function(v, prefix) {
    v <- v[!is.na(v)]
    q <- if (length(v)) stats::quantile(v, c(0, .25, .5, .75, 1), names = FALSE) else rep(NA_real_, 5)
    stats::setNames(as.list(q), paste0(prefix, c("_min", "_q25", "_median", "_q75", "_max")))
  }
  n <- nrow(x); all_positive <- n > 0 && all(positive == x$.actual)
  mixed <- positive > 0 & positive < x$.actual
  collapse <- function(v) if (!length(v)) "no_contributing_loci" else paste(sort(unique(as.character(v)), method = "radix"), collapse = ";")
  answer <- c(list(n_loci = n, n_successful = sum(x$success), n_variants = sum(x$.actual),
    n_positive = sum(positive), n_generating_causal = sum(x$n_causal), n_missing_causal = sum(x$.missing),
    n_generating_null_loci = sum(x$.null), n_no_observed_cause_loci = sum(x$n_observed_causal == 0),
    n_all_positive_loci = sum(positive == x$.actual), n_mixed_loci = sum(mixed),
    n_dropped_candidates = sum(x$.nominal - x$.actual),
    variant_weighted_target_prevalence = if (n) sum(positive) / sum(x$.actual) else NA_real_,
    target_event = event, size_basis = collapse(x$.size_basis), pip_event = collapse(x$pip_event),
    n_basis = collapse(x$.n_basis), variant_n_distribution_status = collapse(x$.variant_n_status),
    n_known_variable_n_track_loci = sum(x$.variable_n_track), n_evidence_route = collapse(x$.n_evidence_route),
    posterior_conditioning = collapse(x$posterior_conditioning),
    probability_target_aligned = if (event != "sparse_component_membership" || !n ||
      anyNA(x$probability_target_aligned)) NA else all(x$probability_target_aligned),
    alignment_scope = "source_flag_for_sparse_component_membership; dense_event_alignment_not_inferred",
    n_alignment_unknown = sum(is.na(x$probability_target_aligned)), all_positive_target = all_positive,
    localisation_status = if (!n) "no_contributing_loci" else if (!any(positive > 0)) "no_positive_scored_loci" else
      if (!any(mixed)) "no_localisation_discrimination" else if (any(positive == x$.actual))
        "mix_of_degenerate_and_discriminating_loci" else "within_locus_discrimination_available",
    constant_one_ap = if (all_positive) 1 else NA_real_, constant_one_brier = if (all_positive) 0 else NA_real_,
    constant_one_log_loss = if (all_positive) 0 else NA_real_),
    distribution(x$.nominal, "nominal_p"), distribution(x$.actual, "actual_p"),
    distribution((x$.nominal - x$.actual) / x$.nominal, "candidate_gap_fraction"),
    distribution(x$.n, "synthetic_n_gwas"), distribution(x$.n_ref, "synthetic_n_ref"),
    distribution(positive / x$.actual, "locus_target_prevalence"))
  as.data.frame(answer, stringsAsFactors = FALSE)
}

.paper_regional_n_caption <- function(x) {
  if (!"regional_n_basis" %in% names(x))
    return("N labels synthetic experiment parameters; sample-size basis is unassessed. ")
  if (any(grepl("gaussian_score_design_parameter", x$regional_n_basis, fixed = TRUE)))
    return(paste0("N is a Gaussian score design parameter; no individual rows are simulated. ",
      "Variant-specific N distributions are unavailable in this locus CSV; see committed public inputs. "))
  "N counts simulated GWAS individuals; reference N counts simulated reference individuals when supplied. "
}

.paper_regional_masks <- function(x) {
  result <- list(all_attempted = rep(TRUE, nrow(x)), ap_successful = is.finite(x$ap),
    brier_successful = is.finite(x$brier), log_loss_successful = !is.na(x$log_loss),
    recall_all_attempts = is.finite(x$recall), cs_distinct_recall_all_attempts = is.finite(x$cs_distinct_recall),
    ap_observed_positive_attempts_zero_for_failure = x$n_observed_causal > 0,
    fdr_micro = rep(TRUE, nrow(x)))
  for (name in c("brier_any_nonzero_total_effect", "log_loss_any_nonzero_total_effect"))
    if (name %in% names(x)) result[[paste0(name, "_successful")]] <- !is.na(x[[name]])
  result
}

.paper_regional_unassessed <- function(rows, contrasts = NULL, reason = "sealed_metric_csv_or_schema_unavailable") {
  rows$regional_census_status <- "unassessed"; rows$regional_census_reason <- reason
  if (!is.null(contrasts)) contrasts$regional_census_status <- rep("unassessed", nrow(contrasts))
  list(rows = rows, contrasts = contrasts,
    census = data.frame(rows[c("stage", "stratum", "method")], endpoint = "unassessed", census_status = "unassessed", reason = reason),
    streaming = list(available = FALSE, reason = reason))
}

.paper_regional_write_csv <- function(x, path) {
  # Preserve the source point values in derived CSVs, including infinities.
  values <- lapply(x, function(v) {
    z <- if (is.double(v)) sprintf("%.17g", v) else as.character(v)
    z[is.na(v) & !is.nan(v)] <- "NA"
    z
  })
  values <- as.data.frame(values, stringsAsFactors = FALSE, optional = TRUE)
  names(values) <- names(x)
  utils::write.table(values, path, sep = ",", row.names = FALSE,
    quote = which(vapply(x, is.character, TRUE)), qmethod = "double", fileEncoding = "UTF-8")
}

.paper_regional_publication <- function(path, schema_path, rows, contrasts = NULL,
    chunk_size = 10000L, max_groups = 100000L, max_stratum_rows = 100000L,
    on_stratum = NULL) {
  bounds <- c(max_groups, max_stratum_rows)
  if (length(bounds) != 2L || any(!is.finite(bounds) | bounds < 1L | bounds != floor(bounds)) ||
      !nrow(rows) || nrow(rows) > max_groups) stop("Regional census group memory bound exceeded")
  wanted <- .paper_regional_key(rows)
  if (anyDuplicated(wanted)) stop("Repeated regional summary identity")
  if (!is.null(on_stratum) && !is.function(on_stratum)) stop("Invalid regional typed-stratum callback")
  # A partition contains exactly the selected stage/method rows in one design
  # stratum. The callback may reuse those typed rows; excluded stages are not
  # missing attempts and are never silently added here.
  strata <- unique(rows$stratum)
  schema <- .paper_regional_schema(schema_path)
  workspace <- tempfile("paper-regional-partitions-"); dir.create(workspace)
  on.exit(unlink(workspace, recursive = TRUE), add = TRUE)
  sizes <- integer(length(strata)); paths <- file.path(workspace, paste0(seq_along(strata), ".blocks"))
  scan <- .paper_regional_scan(path, schema, function(chunk) {
    selected <- .paper_regional_key(chunk) %in% wanted
    if (!any(selected)) return(invisible(NULL))
    x <- .paper_regional_normalize(chunk[selected, , drop = FALSE])
    group <- match(x$stratum, strata)
    for (g in unique(group)) {
      z <- x[group == g, , drop = FALSE]; sizes[g] <<- sizes[g] + nrow(z)
      if (sizes[g] > max_stratum_rows) stop("Regional census stratum memory bound exceeded")
      raw <- serialize(z, NULL, version = 2L); con <- file(paths[g], "ab")
      tryCatch({ writeBin(as.double(length(raw)), con, size = 8L, endian = "little"); writeBin(raw, con) },
        finally = close(con))
    }
  }, chunk_size)
  if (any(sizes == 0)) stop("Regional census has a missing selected design stratum")
  census <- list(); annotations <- list(); pairs <- if (is.null(contrasts)) NULL else contrasts
  for (g in seq_along(strata)) {
    con <- file(paths[g], "rb"); blocks <- list()
    tryCatch(repeat {
      size <- readBin(con, "double", n = 1L, size = 8L, endian = "little")
      if (!length(size)) break
      blocks[[length(blocks) + 1L]] <- unserialize(readBin(con, "raw", n = size))
    }, finally = close(con))
    x <- do.call(rbind, blocks)
    if (anyDuplicated(paste(x$stage, x$method, x$.unit, sep = "\037"))) stop("Duplicate regional locus identity across CSV chunks")
    for (idx in which(rows$stratum == strata[g])) {
      method <- rows$method[idx]
      z <- x[x$stage == rows$stage[idx] & x$method == method, , drop = FALSE]
      if (!nrow(z) || nrow(z) != rows$n_attempted[idx] || sum(z$success) != rows$n_successful[idx])
        stop("Regional census disagrees with summary denominators")
      if (sum(z$n_observed_causal > 0) != rows$n_ap_eligible_attempts[idx] ||
          sum(z$success & z$n_observed_causal > 0) != rows$n_ap_eligible_successful[idx])
        stop("Regional census disagrees with AP eligibility denominators")
      for (endpoint in names(.paper_regional_masks(z))) {
        mask <- .paper_regional_masks(z)[[endpoint]]
        event <- if (grepl("any_nonzero", endpoint, fixed = TRUE)) "any_nonzero_total_effect" else "sparse_component_membership"
        description <- .paper_regional_describe(z[mask, , drop = FALSE], event)
        description$stage <- rows$stage[idx]; description$stratum <- rows$stratum[idx]; description$method <- method
        description$endpoint <- endpoint; description$census_status <- "verified_selected_locus_counts"
        description$n_infinite_endpoint_scores <- if (endpoint %in% names(rows) && grepl("log_loss", endpoint, fixed = TRUE)) {
          name <- sub("_successful$", "", endpoint); sum(is.infinite(z[[name]][mask]))
        } else 0L
        census[[length(census) + 1L]] <- description
        if (endpoint == "all_attempted") annotations[[as.character(idx)]] <- description
        if (endpoint == "brier_successful") {
          annotations[[as.character(idx)]]$brier_probability_target_aligned <- description$probability_target_aligned
          annotations[[as.character(idx)]]$brier_n_alignment_unknown <- description$n_alignment_unknown
        }
      }
    }
    if (!is.null(pairs)) for (i in which(pairs$stratum == strata[g])) {
      a <- x[x$stage == pairs$stage[i] & x$method == pairs$method_a[i], , drop = FALSE]
      b <- x[x$stage == pairs$stage[i] & x$method == pairs$method_b[i], , drop = FALSE]
      if (!nrow(a) || !nrow(b) || !pairs$metric[i] %in% c("ap", "brier", "recall", "cs_distinct_recall"))
        stop("Regional paired census side or source metric missing")
      j <- match(a$.unit, b$.unit); a <- a[!is.na(j), , drop = FALSE]; b <- b[j[!is.na(j)], , drop = FALSE]
      for (field in c(".actual", ".nominal", ".n", ".n_ref", "n_causal", "n_observed_causal", ".missing"))
        if (!identical(unname(a[[field]]), unname(b[[field]]))) stop("Regional paired candidate or truth identities disagree")
      metric <- pairs$metric[i]; good <- is.finite(a[[metric]] - b[[metric]])
      if (nrow(a) != pairs$n_attempted_pairs[i] || sum(good) != pairs$n_pairs[i])
        stop("Regional paired census disagrees with stored contributing-pair denominator")
      for (side in c("a", "b")) {
        d <- .paper_regional_describe(if (side == "a") a[good, , drop = FALSE] else b[good, , drop = FALSE])
        for (name in names(d)) {
          field <- paste0("regional_", name, "_", side)
          if (!field %in% names(pairs)) pairs[[field]] <- d[[name]][rep(NA_integer_, nrow(pairs))]
          pairs[[field]][i] <- d[[name]]
        }
      }
    }
    if (!is.null(on_stratum)) on_stratum(x, stratum = strata[g])
  }
  all <- do.call(rbind, annotations[as.character(seq_len(nrow(rows)))])
  for (name in setdiff(names(all), c("stage", "stratum", "method", "endpoint", "census_status")))
    rows[[paste0("regional_", name)]] <- all[[name]]
  rows$regional_census_status <- "verified_selected_locus_counts"
  if (!is.null(pairs)) pairs$regional_census_status <- rep("verified_exact_finite_difference_pairs", nrow(pairs))
  list(rows = rows, contrasts = pairs, census = do.call(rbind, census),
    streaming = c(list(available = TRUE, retained_rows = sum(sizes), peak_stratum_rows = max(sizes),
      max_stratum_rows = max_stratum_rows, max_groups = max_groups), scan))
}
