# Predictor validation and locus-weighted scoring. Failed outputs are not PIPs.

#' Validate a public challenge input
#' @param input Public predictor input; no assessor fields are permitted.
#' @return Invisibly TRUE, or an error identifying the invalid contract.
#' @export
validate_challenge_input <- function(input) {
  required <- c("task_id", "variant_id", "build", "z", "n", "beta_hat", "se", "LD", "annotations")
  allowed <- c(required, "cohorts", "n_per_variant", "sample_overlap")
  if (!is.list(input) || !all(required %in% names(input)) ||
      any(!names(input) %in% allowed)) stop("Input fields violate the public contract")
  p <- length(input$variant_id)
  if (!p || !is.character(input$task_id) || length(input$task_id) != 1L || anyNA(input$task_id) ||
      !nzchar(input$task_id) || !is.character(input$variant_id) ||
      anyNA(input$variant_id) || any(!nzchar(input$variant_id)) || anyDuplicated(input$variant_id))
    stop("Invalid task/variant identifiers")
  if (!is.character(input$build) || length(input$build) != 1L || anyNA(input$build) ||
      !nzchar(input$build)) stop("Invalid build identity")
  .challenge_scalar(input$n, "n", 3)
  for (v in c("z", "beta_hat", "se"))
    if (!is.numeric(input[[v]]) || length(input[[v]]) != p || any(!is.finite(input[[v]])))
      stop("Invalid ", v)
  if (any(input$se <= 0)) stop("Nonpositive standard errors")
  check_ld <- function(R) {
    if (!is.matrix(R) || !is.numeric(R) || !identical(dim(R), c(p, p)) ||
        any(!is.finite(R)) || max(abs(R - t(R))) > 1e-8 ||
        max(abs(diag(R) - 1)) > 1e-8 || min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) < -1e-7)
      stop("Invalid LD matrix")
  }
  check_ld(input$LD)
  if (!is.matrix(input$annotations) || nrow(input$annotations) != p ||
      any(!is.finite(input$annotations))) stop("Invalid annotations")
  if (!is.null(input$cohorts)) {
    if (!is.list(input$cohorts) || !length(input$cohorts)) stop("Invalid cohorts")
    for (cc in input$cohorts) {
    if (!is.list(cc) || !all(c("z", "LD", "n", "cohort_id") %in% names(cc)) ||
        any(!names(cc) %in% c("z", "LD", "n", "cohort_id")) ||
        !is.character(cc$cohort_id) || length(cc$cohort_id) != 1L || anyNA(cc$cohort_id) ||
        !is.numeric(cc$z) || length(cc$z) != p || any(!is.finite(cc$z))) stop("Invalid cohort scores/fields")
    .challenge_scalar(cc$n, "cohort n", 3); check_ld(cc$LD)
    }
  }
  if (!is.null(input$n_per_variant) &&
      (length(input$n_per_variant) != p || any(!is.finite(input$n_per_variant)) ||
       any(input$n_per_variant < 3))) stop("Invalid variant-wise n")
  if (!is.null(input$sample_overlap)) {
    check_ld(input$sample_overlap)
    if (is.null(input$n_per_variant) || !is.null(input$cohorts) || any(input$sample_overlap < 0))
      stop("Invalid sample-overlap metadata")
  }
  invisible(TRUE)
}

.challenge_validate_inputs <- function(inputs) {
  if (!is.list(inputs) || !length(inputs) || is.null(names(inputs)) || anyDuplicated(names(inputs)) ||
      !identical(names(inputs), unname(vapply(inputs, function(x) x$task_id %||% "", ""))))
    stop("Input names must be unique task identifiers")
  for (x in inputs) validate_challenge_input(x)
}

.challenge_validate_panel <- function(panel) {
  .challenge_validate_inputs(panel$inputs)
  if (!is.list(panel$truth) || !identical(names(panel$inputs), names(panel$truth)))
    stop("Assessor truth must align with public inputs")
  for (id in names(panel$inputs)) {
    input <- panel$inputs[[id]]; truth <- panel$truth[[id]]
    if (!is.list(truth) || !identical(truth$task_id, id)) stop("Invalid assessor task ID")
    for (field in c("all_variant_id", "observed_variant_id", "causal_variant_id")) {
      ids <- truth[[field]]
      if (!is.character(ids) || anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids)))
        stop("Invalid assessor variant IDs")
    }
    if (!identical(truth$observed_variant_id, input$variant_id) ||
        any(!input$variant_id %in% truth$all_variant_id) ||
        any(!truth$causal_variant_id %in% truth$all_variant_id) ||
        !identical(truth$null, length(truth$causal_variant_id) == 0L) ||
        length(truth$n_missing_causal) != 1L || !is.finite(truth$n_missing_causal) ||
        truth$n_missing_causal != sum(!truth$causal_variant_id %in% input$variant_id))
      stop("Assessor truth is inconsistent with observed/causal variants")
    for (field in c("scenario", "cluster_id"))
      if (!is.character(truth[[field]]) || length(truth[[field]]) != 1L ||
          anyNA(truth[[field]]) || !nzchar(truth[[field]])) stop("Invalid assessor stratum/cluster")
    checked <- input; checked$LD <- truth$LD_true; validate_challenge_input(checked)
  }
  invisible(TRUE)
}

.challenge_validate_prediction <- function(pred, input) {
  if (!is.list(pred) || !identical(pred$task_id, input$task_id)) stop("Task ID mismatch")
  if (!identical(pred$status, "ok")) stop(pred$error %||% "Predictor did not succeed")
  if (!is.character(pred$variant_id) || anyDuplicated(pred$variant_id) ||
      !setequal(pred$variant_id, input$variant_id)) stop("Variant IDs must match exactly")
  if (!is.numeric(pred$pip) || length(pred$pip) != length(input$z) ||
      any(!is.finite(pred$pip)) || any(pred$pip < 0 | pred$pip > 1)) stop("Invalid PIPs")
  pred$pip <- pred$pip[match(input$variant_id, pred$variant_id)]
  pred$variant_id <- input$variant_id
  for (field in c("credible_sets", "credible_configurations")) if (!is.null(pred[[field]])) {
    if (!is.list(pred[[field]])) stop(field, " must be a list of variant ID sets")
    for (s in pred[[field]]) if (!is.character(s) || anyNA(s) ||
      anyDuplicated(s) || any(!s %in% input$variant_id) ||
      (field == "credible_sets" && !length(s))) stop("Invalid ", field)
  }
  if ((!is.null(pred$credible_configurations) || !is.null(pred$configuration_coverage)) &&
      (is.null(pred$credible_configurations) || !is.numeric(pred$configuration_coverage) ||
       length(pred$configuration_coverage) != 1L || !is.finite(pred$configuration_coverage) ||
       pred$configuration_coverage <= 0 || pred$configuration_coverage >= 1))
    stop("Configuration collections require a declared nominal coverage")
  if (!is.null(pred$null_probability) && (length(pred$null_probability) != 1L ||
      !is.finite(pred$null_probability) || pred$null_probability < 0 || pred$null_probability > 1))
    stop("Invalid null probability")
  if (!is.null(pred$null_probability) && max(pred$pip) + pred$null_probability > 1 + 1e-8)
    stop("Null probability contradicts marginal PIPs")
  if (!is.null(pred$null_probability) && 1 - pred$null_probability > sum(pred$pip) + 1e-8)
    stop("Null probability violates the PIP union bound")
  pred
}

.challenge_ap <- function(y, score) {
  if (!sum(y)) return(NA_real_)
  oo <- order(score, decreasing = TRUE); yy <- y[oo]; ss <- score[oo]
  ends <- c(which(diff(ss) != 0), length(ss))
  tp <- cumsum(yy)[ends]
  sum(diff(c(0, tp / sum(y))) * tp / ends)
}

.challenge_distinct_hits <- function(sets, causal) {
  if (!length(sets) || !length(causal)) return(0L)
  owner <- rep(0L, length(causal))
  visit <- function(i) {
    for (j in which(causal %in% sets[[i]])) {
      if (seen[j]) next
      seen[j] <<- TRUE
      if (owner[j] == 0L || visit(owner[j])) { owner[j] <<- i; return(TRUE) }
    }
    FALSE
  }
  hits <- 0L
  for (i in seq_along(sets)) { seen <- rep(FALSE, length(causal)); hits <- hits + visit(i) }
  hits
}

#' Run predictor functions using public input only
#' @param inputs Named list of public inputs.
#' @param methods Named list of functions taking one input and returning PIPs.
#' @return Named method predictions including failures and elapsed seconds.
#' @export
run_challenge <- function(inputs, methods = list(exact = challenge_exact,
                                                forced_one = challenge_forced_one)) {
  .challenge_validate_inputs(inputs)
  if (!length(inputs) || anyDuplicated(names(inputs)) || !length(methods) ||
      is.null(names(methods)) || any(!nzchar(names(methods))) || anyDuplicated(names(methods)) ||
      any(!vapply(methods, is.function, TRUE))) stop("Invalid methods/inputs")
  for (x in inputs) validate_challenge_input(x)
  out <- lapply(methods, function(fun) lapply(inputs, function(input) {
    start <- proc.time()["elapsed"]
    warnings <- character(0)
    pred <- tryCatch(withCallingHandlers(.challenge_validate_prediction(fun(input), input),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) list(task_id = input$task_id, status = "failed", error = conditionMessage(e)))
    pred$runtime_seconds <- unname(proc.time()["elapsed"] - start)
    pred$warnings <- warnings
    pred
  }))
  out
}

.challenge_safe_mean <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

.challenge_score <- function(pred, input, truth, threshold) {
  err <- NULL
  runtime <- if (is.list(pred)) pred$runtime_seconds %||% NA_real_ else NA_real_
  warning_count <- if (is.list(pred)) length(pred$warnings %||% character()) else 0L
  converged <- if (is.list(pred)) pred$converged %||% NA else NA
  if (!is.logical(converged) || length(converged) != 1L) converged <- NA
  if (length(runtime) != 1L || !is.numeric(runtime) || !is.finite(runtime) || runtime < 0)
    runtime <- NA_real_
  pred <- tryCatch(.challenge_validate_prediction(pred, input),
                   error = function(e) { err <<- conditionMessage(e); NULL })
  base <- list(task_id = input$task_id, scenario = truth$scenario,
               cluster_id = truth$cluster_id, null = truth$null,
               n_causal = length(truth$causal_variant_id),
               n_missing_causal = truth$n_missing_causal,
               success = is.null(err), error = err %||% "",
               converged = converged, warning_count = warning_count,
               runtime_seconds = runtime)
  if (!is.null(err)) return(c(base, list(ap = NA_real_, brier = NA_real_,
    log_loss_clipped = NA_real_, tp = 0, fp = 0, recall = if (base$n_causal) 0 else NA_real_,
    reported = FALSE, cs_hit_coverage = NA_real_, cs_distinct_recall = if (base$n_causal) 0 else NA_real_,
    configuration_hit = NA_real_, configuration_nominal = NA_real_, configuration_count = NA_real_,
    null_brier = NA_real_, ld_class_top_hit = NA_real_)))
  y <- as.numeric(input$variant_id %in% truth$causal_variant_id)
  pip <- pred$pip; calls <- pip >= threshold
  tp <- sum(y[calls]); fp <- sum(1 - y[calls]); eps <- 1e-15
  clipped <- pmin(1 - eps, pmax(eps, pip))
  cs <- pred$credible_sets %||% list()
  # Repeated sets do not create independent signal discoveries.
  if (length(cs)) cs <- cs[!duplicated(vapply(cs, function(s) paste(sort(s), collapse = "|"), ""))]
  cs_hits <- vapply(cs, function(s) any(s %in% truth$causal_variant_id), TRUE)
  config_hit <- if (is.null(pred$credible_configurations)) NA_real_ else
    as.numeric(any(vapply(pred$credible_configurations,
                         function(s) setequal(s, truth$causal_variant_id), TRUE)))
  top <- which(pip == max(pip))
  observed_causal <- which(y == 1)
  # A resolution diagnostic, never called a class posterior probability.
  equiv <- if (length(observed_causal)) which(apply(
    abs(truth$LD_true[, observed_causal, drop = FALSE]) >= 1 - 1e-8, 1, any)) else integer(0)
  c(base, list(ap = .challenge_ap(y, pip), brier = mean((pip - y)^2),
    log_loss_clipped = -mean(y * log(clipped) + (1 - y) * log1p(-clipped)),
    tp = tp, fp = fp, recall = if (base$n_causal) tp / base$n_causal else NA_real_,
    reported = any(calls), cs_hit_coverage = .challenge_safe_mean(cs_hits),
    cs_distinct_recall = if (base$n_causal) .challenge_distinct_hits(cs, truth$causal_variant_id) / base$n_causal else NA_real_,
    configuration_hit = config_hit,
    configuration_nominal = pred$configuration_coverage %||% NA_real_,
    configuration_count = if (is.null(pred$credible_configurations)) NA_real_ else length(pred$credible_configurations),
    null_brier = if (is.null(pred$null_probability)) NA_real_ else (pred$null_probability - truth$null)^2,
    ld_class_top_hit = if (length(equiv)) mean(top %in% equiv) else NA_real_))
}

#' Score all attempted challenge tasks and report separate mechanism strata
#' @param panel Output of generate_challenge, held by the assessor.
#' @param predictions Output of run_challenge or the same named contract.
#' @param threshold Prespecified PIP discovery threshold.
#' @return Task metrics, mechanism summary and variant-level reliability table.
#' @export
evaluate_challenge <- function(panel, predictions, threshold = .95) {
  .challenge_validate_panel(panel)
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0 || threshold >= 1)
    stop("threshold must be in (0,1)")
  if (!length(predictions) || is.null(names(predictions)) || anyDuplicated(names(predictions)))
    stop("Predictions must have unique method names")
  rows <- list(); reliability <- list()
  for (method in names(predictions)) {
    pp <- predictions[[method]]
    if (length(pp) && (is.null(names(pp)) || anyDuplicated(names(pp)) || any(!names(pp) %in% names(panel$inputs))))
      stop("Unexpected/duplicate task identifiers in predictions")
    for (id in names(panel$inputs)) {
      input <- panel$inputs[[id]]; truth <- panel$truth[[id]]
      row <- .challenge_score(pp[[id]], input, truth, threshold)
      rows[[length(rows) + 1L]] <- as.data.frame(c(list(method = method), row), stringsAsFactors = FALSE)
      if (row$success) {
        pred <- .challenge_validate_prediction(pp[[id]], input)
        reliability[[length(reliability) + 1L]] <- data.frame(
          method = method, scenario = truth$scenario, task_id = id,
          pip = pred$pip, causal = input$variant_id %in% truth$causal_variant_id)
      }
    }
  }
  tasks <- do.call(rbind, rows)
  groups <- split(tasks, interaction(tasks$method, tasks$scenario, drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(d) data.frame(
    method = d$method[1], scenario = d$scenario[1], n_attempted = nrow(d),
    success_rate = mean(d$success), ap_macro_successful = .challenge_safe_mean(d$ap),
    nonconverged_successful = sum(!d$converged[d$success], na.rm = TRUE),
    warning_count_attempted = sum(d$warning_count),
    brier_successful = .challenge_safe_mean(d$brier),
    log_loss_clipped_successful = .challenge_safe_mean(d$log_loss_clipped),
    recall_all_attempts = .challenge_safe_mean(d$recall),
    fdr_micro = if (sum(d$tp + d$fp)) sum(d$fp) / sum(d$tp + d$fp) else NA_real_,
    n_discoveries = sum(d$tp + d$fp),
    null_false_call_rate_successful = .challenge_safe_mean(d$reported[d$null & d$success]),
    cs_hit_coverage_successful = .challenge_safe_mean(d$cs_hit_coverage),
    cs_distinct_recall_all_attempts = .challenge_safe_mean(d$cs_distinct_recall),
    configuration_hit_rate_successful = .challenge_safe_mean(d$configuration_hit),
    configuration_nominal_successful = if (length(unique(stats::na.omit(d$configuration_nominal))) == 1L)
      unique(stats::na.omit(d$configuration_nominal)) else NA_real_,
    configuration_count_mean = .challenge_safe_mean(d$configuration_count),
    null_brier_successful = .challenge_safe_mean(d$null_brier),
    ld_class_top_hit_successful = .challenge_safe_mean(d$ld_class_top_hit),
    missing_causal_fraction = if (sum(d$n_causal)) sum(d$n_missing_causal) / sum(d$n_causal) else NA_real_,
    runtime_median_successful = stats::median(d$runtime_seconds[d$success], na.rm = TRUE),
    runtime_total_attempted = if (all(is.finite(d$runtime_seconds))) sum(d$runtime_seconds) else NA_real_,
    stringsAsFactors = FALSE)))
  rownames(tasks) <- rownames(summary) <- NULL
  rel <- if (length(reliability)) do.call(rbind, reliability) else data.frame()
  list(tasks = tasks, summary = summary, reliability = rel,
       specification = panel$config$specification, threshold = threshold,
       warning = "Conditional scores require the availability column; no cross-mechanism overall rank is defined.")
}

#' Paired cluster bootstrap of a prespecified task metric difference
#' @param evaluation Output of evaluate_challenge().
#' @param method_a,method_b Compared method names.
#' @param metric Numeric task metric, for example recall or brier.
#' @param scenario Mechanism stratum, required to avoid averaging incompatible tasks.
#' @param n_boot Number of bootstrap replicates.
#' @param seed Resampling seed.
#' @return Difference A minus B, percentile interval and paired availability.
#' @export
challenge_compare <- function(evaluation, method_a, method_b, metric, scenario,
                              n_boot = 1000L, seed = 1L) {
  .challenge_scalar(n_boot, "n_boot", 20, TRUE)
  d <- evaluation$tasks
  if (!method_a %in% d$method || !method_b %in% d$method || !scenario %in% d$scenario)
    stop("Unknown method or scenario")
  if (!metric %in% names(d) || !is.numeric(d[[metric]])) stop("Invalid metric")
  a <- d[d$method == method_a & d$scenario == scenario, ]
  b <- d[d$method == method_b & d$scenario == scenario, ]
  m <- merge(a[c("task_id", "cluster_id", metric)], b[c("task_id", metric)],
             by = "task_id", suffixes = c("_a", "_b"))
  good <- is.finite(m[[paste0(metric, "_a")]]) & is.finite(m[[paste0(metric, "_b")]])
  n_available <- sum(good); m <- m[good, ]
  diff <- m[[paste0(metric, "_a")]] - m[[paste0(metric, "_b")]]
  clusters <- split(seq_len(nrow(m)), m$cluster_id)
  if (length(clusters) < 2L) return(list(difference = .challenge_safe_mean(diff),
    lower = NA_real_, upper = NA_real_, n_pairs = n_available,
    n_attempted_pairs = length(good), n_clusters = length(clusters)))
  values <- .challenge_rng(seed, replicate(n_boot, {
    selected <- sample(seq_along(clusters), length(clusters), replace = TRUE)
    mean(diff[unlist(clusters[selected], use.names = FALSE)])
  }))
  ci <- stats::quantile(values, c(.025, .975), names = FALSE)
  list(difference = mean(diff), lower = ci[1], upper = ci[2],
       n_pairs = n_available, n_attempted_pairs = length(good), n_clusters = length(clusters))
}
