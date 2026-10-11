#!/usr/bin/env Rscript
# Plot frozen collector evidence; separately specified stage contrasts use
# stored locus scores and a fixed cluster bootstrap. No refitting or tuning.
main <- function() {
args <- commandArgs(TRUE)
if (length(args) < 2L)
  stop("Usage: Rscript paper/plot_results.R ANALYSIS OUTPUT [--allow-incomplete] [--setting-regex PATTERN] [--methods ID,ID] [--max-settings N]")
selection <- list(allow_incomplete = FALSE, setting_regex = "", methods = "", max_settings = 100L)
i <- 3L
while (i <= length(args)) {
  if (args[i] == "--allow-incomplete") { selection$allow_incomplete <- TRUE; i <- i + 1L; next }
  key <- sub("^--", "", args[i]); key <- gsub("-", "_", key, fixed = TRUE)
  if (!key %in% c("setting_regex", "methods", "max_settings") || i == length(args)) stop("Invalid plot option")
  selection[[key]] <- args[i + 1L]; i <- i + 2L
}
selection$max_settings <- as.integer(selection$max_settings)
if (is.na(selection$max_settings) || selection$max_settings < 1L) stop("Invalid setting limit")
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), ".."))
reader_source <- file.path(repo, "paper", "read_selected_csv.R")
source(reader_source)
native_publication_source <- file.path(repo, "paper", "native_publication.R")
source(native_publication_source)
regional_publication_source <- file.path(repo, "paper", "regional_publication.R")
source(regional_publication_source)
if (nzchar(Sys.getenv("FMB_R_LIB"))) .libPaths(c(Sys.getenv("FMB_R_LIB"), .libPaths()))
if (dir.exists(file.path(repo, ".local-lib"))) .libPaths(c(file.path(repo, ".local-lib"), .libPaths()))
for (path in list.files(file.path(repo, "R"), full.names = TRUE)) if (endsWith(path, ".R")) source(path)
analysis <- normalizePath(args[1L], mustWork = TRUE)
lock <- file.path(dirname(analysis), "analysis.lock")
.campaign_claim(lock)
on.exit(unlink(lock, recursive = TRUE), add = TRUE)
output <- args[2L]
if (dir.exists(output)) stop("Use a new figure output directory to preserve previous exports")
manifest_path <- file.path(analysis, "analysis-manifest.rds")
manifest <- readRDS(manifest_path)
if (!identical(manifest$schema, 1L)) stop("Unknown analysis manifest schema")
.campaign_validate_artifacts(manifest$artifacts)
manifest_json <- jsonlite::fromJSON(file.path(analysis, "analysis-manifest.json"))
if (!is.null(manifest_json$manifest_sha256)) {
  if (!identical(manifest_json$manifest_sha256, .challenge_sha256(manifest_path)))
    stop("Analysis manifest bytes differ from the readable commitment")
} else if (!.campaign_check_receipt(analysis, manifest_path, manifest$plan_sha256, "analysis")) {
  stop("Native analysis manifest lacks its terminal receipt")
}
plan <- readRDS(file.path(dirname(analysis), "plan.rds"))
if (!identical(plan$plan_sha256, manifest$plan_sha256)) stop("Analysis and campaign plan differ")
x <- plan; x$plan_sha256 <- NULL
if (!identical(plan$plan_sha256, .campaign_object_sha(x))) stop("Campaign plan bytes are incompatible with its commitment")
status <- manifest$collection_status %||% manifest$status
explicit_fixture <- is.list(plan$runtime_fingerprint) && isTRUE(plan$runtime_fingerprint$test_fixture)
scientific <- identical(status, "complete") && identical(plan$profile, "full") && !explicit_fixture
if (!scientific && !selection$allow_incomplete)
  stop("Publication exports require a complete full-profile collection; explicit synthetic fixtures are not scientific evidence; use --allow-incomplete for clearly marked diagnostics")
label <- if (scientific) "Full-profile collected evidence" else if (plan$profile != "full" || explicit_fixture)
  paste("DEVELOPMENT ONLY: profile", plan$profile, "collection", status) else
  paste("INCOMPLETE CAMPAIGN: collection", status)
synthetic_fixture <- explicit_fixture ||
  (is.character(plan$source) && length(plan$source) == 1L &&
   grepl("synthetic.*fixture|fixture.*synthetic", plan$source, ignore.case = TRUE))
if (synthetic_fixture) label <- paste(label, "| SYNTHETIC SOFTWARE FIXTURE: no method-performance evidence")
read_table <- function(name, select = identity, max_rows = 100000L) {
  path <- file.path(analysis, name)
  if (!file.exists(path) || !path %in% manifest$artifacts$path)
    stop("Required table is absent from the sealed analysis: ", name)
  .paper_read_selected_csv(path, select = select, max_rows = max_rows)
}
native <- file.exists(file.path(analysis, "draw-metrics.csv"))
validate_native_chunk <- function(chunk) {
  has_eligibility <- "all_native_accuracy_eligible" %in% names(chunk)
  if (!has_eligibility && scientific)
    stop("Publication native figures require explicit accuracy eligibility; recollect under the current contract")
  if (has_eligibility) {
    ineligible <- chunk$method == "sbayesrc_native" & !(chunk$all_native_accuracy_eligible %in% TRUE)
    fields <- intersect(c("ap_successful", "brier_successful", "recall_all_completed"), names(chunk))
    if (any(vapply(chunk[ineligible, fields, drop = FALSE], function(x) any(is.finite(x)), TRUE)))
      stop("Unqualified native accuracy appeared in ordinary comparison columns")
  }
}
requested_methods <- if (nzchar(selection$methods)) strsplit(selection$methods, ",", fixed = TRUE)[[1L]] else character()
if (anyDuplicated(requested_methods)) stop("Unknown or repeated method filter")
seen_methods <- seen_settings <- character()
rows <- read_table(if (native) "draw-metrics.csv" else "regime-summary.csv", function(chunk) {
  seen_methods <<- union(seen_methods, chunk$method)
  if (native) {
    # Validate all streamed native rows before filters can hide an invalid row.
    validate_native_chunk(chunk)
    keys <- c("model", "pi", "h2", "exposure", "event")
    chunk$setting <- apply(chunk[keys], 1L, function(v) paste(paste(keys, v, sep = "="), collapse = "; "))
  } else chunk$setting <- paste(chunk$stage, chunk$stratum, sep = "; ")
  if (length(requested_methods)) chunk <- chunk[chunk$method %in% requested_methods, , drop = FALSE]
  if (nzchar(selection$setting_regex)) chunk <- chunk[grepl(selection$setting_regex, chunk$setting, perl = TRUE), , drop = FALSE]
  seen_settings <<- union(seen_settings, chunk$setting)
  if (length(seen_settings) > selection$max_settings)
    stop("Export exceeds --max-settings; select prespecified settings using --setting-regex")
  chunk
})
if (any(!requested_methods %in% seen_methods)) stop("Unknown or repeated method filter")
if (!nrow(rows)) stop("No measured outcomes are available to plot")
if (native) {
  has_native_eligibility <- "all_native_accuracy_eligible" %in% names(rows)
  if (!has_native_eligibility && scientific)
    stop("Publication native figures require explicit accuracy eligibility; recollect under the current contract")
  if (has_native_eligibility) {
    ineligible <- rows$method == "sbayesrc_native" & !(rows$all_native_accuracy_eligible %in% TRUE)
    quality_fields <- intersect(c("ap_successful", "brier_successful", "recall_all_completed"), names(rows))
    if (any(vapply(rows[ineligible, quality_fields, drop = FALSE], function(x) any(is.finite(x)), TRUE)))
      stop("Unqualified native accuracy appeared in ordinary comparison columns")
  } else label <- paste(label, "| native accuracy eligibility unavailable")
  rows$aligned <- rows$all_probability_targets_aligned
  rows$success_rate <- rows$n_successful / rows$n_planned
  rows$availability_label <- paste0(rows$n_successful, "/", rows$n_planned)
  metrics <- c(success_rate = "Imported / planned", ap_successful = "AP: eligible fits",
               brier_successful = "Brier: eligible fits", recall_all_completed = "Recall: eligible fits")
  if (!has_native_eligibility) metrics <- c(success_rate = "Imported / planned",
    ap_successful = "AP: imported fits", brier_successful = "Brier: imported fits", recall_all_completed = "Recall: imported fits")
} else {
  rows$aligned <- rows$n_probability_target_misaligned_successful == 0L
  rows$availability_label <- paste0(rows$n_successful, "/", rows$n_attempted)
  metrics <- c(success_rate = "Successful / attempted", ap_successful = "AP: successful positive loci",
               brier_successful = "Brier: successful fits", fdr_micro = "FDR: pooled calls at PIP 0.95")
}
required <- c("method", "setting", names(metrics), "aligned")
if (!all(required %in% names(rows))) stop("Collector schema lacks plot fields: ", paste(setdiff(required, names(rows)), collapse = ", "))
if (!nrow(rows)) stop("No sealed outcomes match the requested setting/method filters")
if (length(unique(rows$setting)) > selection$max_settings)
  stop("Export exceeds --max-settings; select prespecified settings using --setting-regex")
native_census_path <- file.path(analysis, "task-metrics.csv")
native_census_available <- FALSE
if (native) {
  native_census_available <- file.exists(native_census_path) &&
    all(.paper_native_keys %in% names(rows))
  if (native_census_available) {
    if (!native_census_path %in% manifest$artifacts$path)
      stop("Native truth census is absent from the sealed analysis commitment")
    rows <- .paper_native_truth_census(native_census_path, rows)
  } else {
    if (scientific) stop("Publication native figures require sealed locus-level truth counts")
    rows <- .paper_native_unassessed(rows)
    label <- paste(label, "| truth census unavailable; localisation unassessed")
  }
}
contrast_name <- if (native) "paired-draw-comparisons.csv" else "paired-comparisons.csv"
contrast_path <- file.path(analysis, contrast_name)
contrasts <- NULL
if (file.exists(contrast_path)) {
  if (native) {
    key <- function(x) do.call(paste, c(x, sep = "|"))
    available <- unique(key(rows[c("model", "pi", "h2", "event", "method", "exposure")]))
    contrasts <- read_table(contrast_name, function(chunk) {
      a <- chunk[c("model", "pi", "h2", "event", "method_a", "exposure_a")]
      b <- chunk[c("model", "pi", "h2", "event", "method_b", "exposure_b")]
      chunk[key(a) %in% available & key(b) %in% available, , drop = FALSE]
    }, max_rows = 2400L)
  } else {
    contrasts <- read_table(contrast_name, function(chunk) {
      chunk[paste(chunk$stage, chunk$stratum, sep = "; ") %in% rows$setting &
        chunk$method_a %in% rows$method & chunk$method_b %in% rows$method, , drop = FALSE]
    }, max_rows = 2400L)
  }
  if (nrow(contrasts) > 2400L)
    stop("Export exceeds 2,400 paired contrasts; select prespecified settings/methods")
}
regional <- NULL
regional_census_inputs <- character()
stage_publication <- NULL
stage_publication_sources <- character()
if (!native) {
  regional_census_inputs <- file.path(analysis, c("task-metrics.csv", "task-metrics-schema.rds"))
  available <- all(file.exists(regional_census_inputs)) &&
    all(regional_census_inputs %in% manifest$artifacts$path)
  if (!available && scientific)
    stop("Publication regional figures require sealed locus metrics and their CSV schema for realized-size/truth census")
  stage_callback <- NULL
  if (available) {
    stage_publication_sources <- file.path(repo, "paper", c("stage_attribution.R", "stage_publication.R", "stage-attribution-spec.json"))
    source(stage_publication_sources[1L]); source(stage_publication_sources[2L])
    stage_adapter <- .paper_stage_publication(plan, rows, stage_publication_sources[3L],
      scientific = scientific, collection_complete = identical(status, "complete"))
    stage_callback <- stage_adapter$on_stratum
  }
  regional <- if (available) .paper_regional_publication(regional_census_inputs[1L],
    regional_census_inputs[2L], rows, contrasts, on_stratum = stage_callback) else .paper_regional_unassessed(rows, contrasts)
  stage_publication <- if (available) stage_adapter$finish() else list(metadata = list(
    available = FALSE, reason = "sealed_regional_metric_csv_or_schema_unavailable"))
  rows <- regional$rows; contrasts <- regional$contrasts
  if (available) rows$aligned <- rows$regional_brier_probability_target_aligned
  if (!available) {
    regional_census_inputs <- character()
    label <- paste(label, "| regional realized-size/truth census unavailable; unassessed")
  }
}
dir.create(output, recursive = TRUE)
inputs <- c(manifest_path, file.path(analysis, "analysis-manifest.json"),
            file.path(dirname(analysis), "plan.rds"),
            file.path(analysis, if (native) "draw-metrics.csv" else "regime-summary.csv"))
if (!native) {
  inputs <- c(inputs, regional_census_inputs)
  .paper_regional_write_csv(regional$census, file.path(output, "regional-target-census.csv"))
  .paper_regional_write_csv(rows, file.path(output, "regional-summary-publication.csv"))
  if (!is.null(contrasts)) .paper_regional_write_csv(contrasts, file.path(output, "regional-paired-publication.csv"))
  if (isTRUE(stage_publication$metadata$available)) {
    .paper_regional_write_csv(stage_publication$selection, file.path(output, "regional-stage-selection.csv"))
    if (nrow(stage_publication$contrasts))
      .paper_regional_write_csv(stage_publication$contrasts, file.path(output, "regional-stage-attribution.csv"))
  }
}
if (native) {
  if (native_census_available) inputs <- c(inputs, native_census_path)
  write.csv(rows, file.path(output, "native-draw-publication.csv"), row.names = FALSE)
  census_columns <- c(intersect(c(.paper_native_keys, "draw"), names(rows)),
    grep("^census_n_", names(rows), value = TRUE), "target_support", "all_positive_target",
    "ap_localisation_status", "recall_localisation_status",
    "constant_one_ap", "constant_one_recall", "constant_one_brier", "constant_one_log_loss")
  write.csv(rows[census_columns], file.path(output, "native-target-census.csv"), row.names = FALSE)
}
settings_path <- file.path(output, "settings.pdf")
grDevices::pdf(settings_path, width = 11.7, height = 8.3, onefile = TRUE, family = "Helvetica")
tryCatch({
  for (setting in sort(unique(rows$setting), method = "radix")) {
    d <- rows[rows$setting == setting, ]; methods <- sort(unique(d$method), method = "radix")
    layout(matrix(1:4, 1L)); par(mar = c(4, 1.5, 4, 1), oma = c(5, 13, 5, 1), cex = .85)
    for (i in seq_along(metrics)) {
      metric <- names(metrics)[i]
      y <- match(d$method, rev(methods))
      # Each native point is one genotype draw; no unregistered locus-level CI.
      if (native) y <- y + (d$draw - mean(range(d$draw))) * .09
      values <- d[[metric]]; valid <- is.finite(values)
      degenerate <- if (native) .paper_native_degenerate(d, metric) else rep(FALSE, nrow(d))
      heading <- metrics[i]
      if (native && all(degenerate)) heading <- switch(metric,
        ap_successful = "AP: no localisation discrimination",
        recall_all_completed = "Recall: all-positive call fraction",
        brier_successful = "Brier: all-positive target", heading)
      plot(NA, xlim = c(0, 1), ylim = c(.3, length(methods) + .7), xlab = "", ylab = "", yaxt = "n",
           main = paste(strwrap(heading, width = 22), collapse = "\n"), cex.main = .9, bty = "n")
      abline(v = seq(0, 1, .25), col = "#eeeeee")
      abline(h = seq_along(methods), col = "#f3f3f3")
      if (i == 1L) axis(2L, at = seq_along(methods), labels = rev(methods), las = 1L, tick = FALSE)
      misaligned <- metric == "brier_successful" & !(d$aligned %in% TRUE)
      ordinary <- valid & !degenerate
      points(values[ordinary & !misaligned], y[ordinary & !misaligned], pch = 16L, col = "#176b80", cex = .8)
      points(values[ordinary & misaligned], y[ordinary & misaligned], pch = 1L, col = "#777777", cex = .9)
      points(values[valid & degenerate & !misaligned], y[valid & degenerate & !misaligned], pch = 17L, col = "#a35b12", cex = .9)
      points(values[valid & degenerate & misaligned], y[valid & degenerate & misaligned], pch = 2L, col = "#777777", cex = .9)
      # Explicit absent outcomes prevent an unsupported method from vanishing.
      if (any(!valid)) points(rep(.025, sum(!valid)), y[!valid], pch = 4L, col = "#9b382a", cex = .7)
      if (!native && i == 1L) text(.98, y, d$availability_label, pos = 2L, cex = .65)
    }
    mtext(paste(strwrap(setting, width = 140), collapse = "\n"), outer = TRUE, side = 3L, line = 2.2, cex = .9)
    mtext(paste(strwrap(label, width = 155), collapse = "\n"), outer = TRUE, side = 3L, line = .1, cex = .72,
          col = if (scientific) "#444444" else "#9b382a")
    note <- if (native) "AP uses positive loci; points are genotype draws. Four-chain success does not establish MCMC convergence. " else
      paste0("Counts show successful / attempted locus fits. ", .paper_regional_n_caption(d), "Settings use requested p; ",
        if (isTRUE(regional$streaming$available)) paste0("retained p ", min(d$regional_actual_p_min), "–",
          max(d$regional_actual_p_max), "; endpoint-specific prevalence/gap distributions are in regional-target-census.csv. ") else
          "retained-size/truth census unassessed. ",
        "Reference-panel regimes are not LD-only paired inputs. ")
    if (native && any(d$ap_localisation_status == "no_localisation_discrimination" |
                      d$recall_localisation_status == "no_localisation_discrimination"))
      note <- paste0(note, "Triangles: all-positive scored loci; AP=1 at any ranking; recall is a call fraction, with no localisation discrimination. ")
    if (native && any(d$all_positive_target %in% TRUE))
      note <- paste0(note, "For an all-positive target, constant PIP=1 gives AP=1, recall=1, Brier=0 and log loss=0 without data. ")
    if (native && any(grepl("degenerate", d$ap_localisation_status, fixed = TRUE)))
      note <- paste0(note, "AP includes both all-positive and discriminating loci; see the exported truth census. ")
    legend <- if (native) "Open Brier points: event-misaligned operational scores. Cross: undefined/no successful score." else
      "Open Brier points: event alignment unqualified; operational scores. Cross: undefined/no successful score."
    mtext(paste(strwrap(paste0(note, legend),
                        width = 170), collapse = "\n"), outer = TRUE, side = 1L, line = .5, cex = .68)
  }
}, finally = grDevices::dev.off())

if (!is.null(contrasts) && nrow(contrasts)) {
  inputs <- c(inputs, contrast_path)
  if (native) {
    contrasts <- .paper_native_paired_view(contrasts, rows)
    write.csv(contrasts, file.path(output, "native-paired-publication.csv"), row.names = FALSE)
    keys <- c("model", "pi", "h2", "event", "metric")
    contrasts$setting <- apply(contrasts[keys], 1L, function(v) paste(paste(keys, v, sep = "="), collapse = "; "))
    contrasts$comparison <- paste0(contrasts$method_a, " [", contrasts$exposure_a, "] minus ",
                                  contrasts$method_b, " [", contrasts$exposure_b, "]")
    contrasts$n_clusters <- contrasts$n_paired_genotype_draws
    contrasts$aligned <- contrasts$metric != "brier_successful" |
      contrasts$probability_events_both_aligned %in% TRUE
  } else {
    contrasts$setting <- paste(contrasts$stage, contrasts$stratum, contrasts$metric, sep = "; ")
    contrasts$comparison <- paste(contrasts$method_a, "minus", contrasts$method_b)
    contrasts$aligned <- contrasts$metric != "brier" |
      (!is.na(contrasts$probability_targets_aligned) & contrasts$probability_targets_aligned)
  }
  needed <- c("difference", "lower", "upper", "n_clusters", "comparison", "setting")
  if (!all(needed %in% names(contrasts))) stop("Paired collector schema lacks required fields")
  grDevices::pdf(file.path(output, "paired-contrasts.pdf"), width = 11.7, height = 8.3, onefile = TRUE, family = "Helvetica")
  tryCatch({
    for (setting in sort(unique(contrasts$setting), method = "radix")) {
      d <- contrasts[contrasts$setting == setting, ]; d <- d[order(d$comparison, method = "radix"), ]
      # Split long pair tables into pages, retaining every declared contrast.
      for (ix in split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / 24L))) {
        z <- d[ix, ]; y <- rev(seq_len(nrow(z))); valid <- is.finite(z$difference)
        span <- range(c(0, z$difference[valid], z$lower[is.finite(z$lower)], z$upper[is.finite(z$upper)]))
        pad <- max(diff(span) * .08, .01); span <- span + c(-pad, pad)
        par(mar = c(5, 28, 2, 2), oma = c(if (native) 4 else 6, 0, 5, 0), cex = .95)
        plot(NA, xlim = span, ylim = c(.4, nrow(z) + .6), xlab = "",
             ylab = "", yaxt = "n", bty = "n")
        abline(v = 0, lty = 2, col = "#888888")
        axis(2L, at = y, labels = paste0(z$comparison, " (clusters=", z$n_clusters, ")"), las = 1L, tick = FALSE, cex.axis = .85)
        bounded <- valid & is.finite(z$lower) & is.finite(z$upper)
        degenerate <- if (native) z$localisation_interpretation == "no_localisation_discrimination" |
          (z$metric == "brier_successful" & z$all_positive_target_both_sides %in% TRUE) else rep(FALSE, nrow(z))
        col <- ifelse(!z$aligned, "#777777", ifelse(degenerate, "#a35b12", "#176b80"))
        segments(z$lower[bounded], y[bounded], z$upper[bounded], y[bounded], col = col[bounded])
        shapes <- ifelse(degenerate, ifelse(z$aligned, 17L, 2L), ifelse(z$aligned, 16L, 1L))
        points(z$difference[valid], y[valid], pch = shapes[valid], col = col[valid])
        if (any(!valid)) text(rep(span[1L] + pad / 2, sum(!valid)), y[!valid], "undefined", pos = 4L, cex = .7, col = "#9b382a")
        heading <- setting
        if (native && any(degenerate)) heading <- paste0(heading, "; all-positive scored loci: no localisation discrimination")
        mtext(paste(strwrap(heading, width = 115), collapse = "\n"), outer = TRUE,
              side = 3L, line = 2.2, cex = .85)
        mtext(paste(strwrap(label, width = 170), collapse = "\n"), outer = TRUE, side = 3L, line = .1, cex = .65,
              col = if (scientific) "#444444" else "#9b382a")
        mtext("Method A minus method B on the named endpoint", side = 1L, line = 3, cex = .85)
        note <- if (native) "Paired cluster intervals; direction depends on endpoint. Open grey points: event-misaligned operational Brier contrasts, not calibration ranks." else
          paste0("Descriptive cluster percentile envelopes; nominal coverage unvalidated with few genotype clusters; no multiplicity control. Direction depends on endpoint. ",
            "Open grey points: event alignment unqualified; operational Brier contrasts, not calibration ranks.")
        if (!native) note <- paste0(note,
          " Regional census uses exactly the finite-difference pairs; recall includes recorded failure zeros. Requested p differs from retained candidates where filtering applies. N is a synthetic experiment parameter; Gaussian per-variant N distributions require the committed public inputs. Reference-panel regimes are not LD-only pairs.")
        if (native && any(degenerate)) note <- paste0(note,
          " Triangles: no localisation resolution. For an all-positive target, constant PIP=1 attains AP=1, recall=1, Brier=0 and log loss=0 without data.")
        mtext(paste(strwrap(note, width = 175), collapse = "\n"),
              outer = TRUE, side = 1L, line = if (native) .5 else 2.2, cex = .65)
      }
    }
  }, finally = grDevices::dev.off())
}
outputs <- list.files(output, full.names = TRUE)
figure_manifest <- list(schema = 1L, plan_sha256 = manifest$plan_sha256, collection_status = status,
  profile = plan$profile, scientific_full_profile = scientific, label = label,
  selection = selection, n_settings_plotted = length(unique(rows$setting)),
  plotting_source_sha256 = .challenge_sha256(sub("^--file=", "", file_arg[1L])),
  plotting_helpers = .campaign_artifact_manifest(c(reader_source, native_publication_source,
    if (!native) c(regional_publication_source, stage_publication_sources))),
  native_truth_census_available = native && native_census_available,
  regional_truth_census = if (!native) regional$streaming else NULL,
  regional_stage_analysis = if (!native) stage_publication$metadata else NULL,
  inputs = .campaign_artifact_manifest(inputs), outputs = .campaign_artifact_manifest(outputs),
  interpretation = paste("Setting-specific endpoint views, alphabetical method order, no overall rank; existing method-paired intervals unchanged.",
    if (!native && isTRUE(stage_publication$metadata$available))
      "A separately specified regional stage table uses stored locus scores and a fixed paired cluster bootstrap; filters excluding a required stage leave that comparison unassessed. Percentile envelopes are descriptive, with nominal coverage unvalidated and no multiplicity control.",
    "Native truth census and paired event metadata are publication annotations; collected scores and intervals are unchanged.",
    "All-positive AP/recall do not demonstrate localisation; constant-one probabilities give zero all-positive proper-score loss."))
jsonlite::write_json(figure_manifest, file.path(output, "figure-manifest.json"), pretty = TRUE,
                    auto_unbox = TRUE, null = "null", na = "null")
cat("Exported frozen-evidence figures to", normalizePath(output), "\n")
}
main()
