#!/usr/bin/env Rscript
# Plot the collector's frozen evidence. No simulation, tuning or new bootstrap.
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
scientific <- identical(status, "complete") && identical(plan$profile, "full")
if (!scientific && !selection$allow_incomplete)
  stop("Publication exports require a complete full-profile collection; use --allow-incomplete for clearly marked diagnostics")
label <- if (scientific) "Full-profile collected evidence" else if (plan$profile != "full")
  paste("DEVELOPMENT ONLY: profile", plan$profile, "collection", status) else
  paste("INCOMPLETE CAMPAIGN: collection", status)
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
dir.create(output, recursive = TRUE)
inputs <- c(manifest_path, file.path(analysis, "analysis-manifest.json"),
            file.path(dirname(analysis), "plan.rds"),
            file.path(analysis, if (native) "draw-metrics.csv" else "regime-summary.csv"))
settings_path <- file.path(output, "settings.pdf")
grDevices::pdf(settings_path, width = 11.7, height = 8.3, onefile = TRUE, family = "Helvetica")
tryCatch({
  for (setting in sort(unique(rows$setting), method = "radix")) {
    d <- rows[rows$setting == setting, ]; methods <- sort(unique(d$method), method = "radix")
    layout(matrix(1:4, 1L)); par(mar = c(4, 1.5, 3, 1), oma = c(3, 13, 4, 1), cex = .85)
    for (i in seq_along(metrics)) {
      metric <- names(metrics)[i]
      y <- match(d$method, rev(methods))
      # Each native point is one genotype draw; no unregistered locus-level CI.
      if (native) y <- y + (d$draw - mean(range(d$draw))) * .09
      values <- d[[metric]]; valid <- is.finite(values)
      plot(NA, xlim = c(0, 1), ylim = c(.3, length(methods) + .7), xlab = "", ylab = "", yaxt = "n",
           main = paste(strwrap(metrics[i], width = 22), collapse = "\n"), cex.main = .95, bty = "n")
      abline(v = seq(0, 1, .25), col = "#eeeeee")
      abline(h = seq_along(methods), col = "#f3f3f3")
      if (i == 1L) axis(2L, at = seq_along(methods), labels = rev(methods), las = 1L, tick = FALSE)
      misaligned <- metric == "brier_successful" & !d$aligned
      points(values[valid & !misaligned], y[valid & !misaligned], pch = 16L, col = "#176b80", cex = .8)
      points(values[valid & misaligned], y[valid & misaligned], pch = 1L, col = "#777777", cex = .9)
      # Explicit absent outcomes prevent an unsupported method from vanishing.
      if (any(!valid)) points(rep(.025, sum(!valid)), y[!valid], pch = 4L, col = "#9b382a", cex = .7)
      if (!native && i == 1L) text(.98, y, d$availability_label, pos = 2L, cex = .65)
    }
    mtext(paste(strwrap(setting, width = 140), collapse = "\n"), outer = TRUE, side = 3L, line = 1.6, cex = .9)
    mtext(label, outer = TRUE, side = 3L, line = .1, cex = .72, col = if (scientific) "#444444" else "#9b382a")
    note <- if (native) "AP uses positive loci; points are genotype draws. Four-chain success does not establish MCMC convergence. " else "Counts show successful / attempted locus fits. "
    mtext(paste(strwrap(paste0(note, "Open Brier points: event-misaligned operational scores. Cross: undefined/no successful score."),
                        width = 170), collapse = "\n"), outer = TRUE, side = 1L, line = .5, cex = .68)
  }
}, finally = grDevices::dev.off())

if (!is.null(contrasts) && nrow(contrasts)) {
  inputs <- c(inputs, contrast_path)
  if (native) {
    keys <- c("model", "pi", "h2", "event", "metric")
    contrasts$setting <- apply(contrasts[keys], 1L, function(v) paste(paste(keys, v, sep = "="), collapse = "; "))
    contrasts$comparison <- paste0(contrasts$method_a, " [", contrasts$exposure_a, "] minus ",
                                  contrasts$method_b, " [", contrasts$exposure_b, "]")
    contrasts$n_clusters <- contrasts$n_paired_genotype_draws
    contrasts$aligned <- vapply(seq_len(nrow(contrasts)), function(i) {
      z <- contrasts[i, ]
      if (z$metric != "brier_successful") return(TRUE)
      scope <- rows$model == z$model & rows$pi == z$pi & rows$h2 == z$h2 & rows$event == z$event
      a <- scope & rows$method == z$method_a & rows$exposure == z$exposure_a
      b <- scope & rows$method == z$method_b & rows$exposure == z$exposure_b
      any(a) && any(b) && all(rows$aligned[a | b])
    }, TRUE)
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
        par(mar = c(5, 28, 5, 2), oma = c(3, 0, 0, 0), cex = .95)
        plot(NA, xlim = span, ylim = c(.4, nrow(z) + .6), xlab = "",
             ylab = "", yaxt = "n", bty = "n")
        abline(v = 0, lty = 2, col = "#888888")
        axis(2L, at = y, labels = paste0(z$comparison, " (clusters=", z$n_clusters, ")"), las = 1L, tick = FALSE, cex.axis = .85)
        bounded <- valid & is.finite(z$lower) & is.finite(z$upper)
        col <- ifelse(z$aligned, "#176b80", "#777777")
        segments(z$lower[bounded], y[bounded], z$upper[bounded], y[bounded], col = col[bounded])
        points(z$difference[valid], y[valid], pch = ifelse(z$aligned[valid], 16L, 1L), col = col[valid])
        if (any(!valid)) text(rep(span[1L] + pad / 2, sum(!valid)), y[!valid], "undefined", pos = 4L, cex = .7, col = "#9b382a")
        title(main = paste(strwrap(setting, width = 115), collapse = "\n"), cex.main = .9)
        mtext(label, side = 3L, line = .1, cex = .65, col = if (scientific) "#444444" else "#9b382a")
        mtext("Method A minus method B on the named endpoint", side = 1L, line = 3, cex = .85)
        mtext("Paired cluster intervals; direction depends on endpoint. Open grey points: event-misaligned operational Brier contrasts, not calibration ranks.",
              outer = TRUE, side = 1L, line = 1, cex = .65)
      }
    }
  }, finally = grDevices::dev.off())
}
outputs <- list.files(output, full.names = TRUE)
figure_manifest <- list(schema = 1L, plan_sha256 = manifest$plan_sha256, collection_status = status,
  profile = plan$profile, scientific_full_profile = scientific, label = label,
  selection = selection, n_settings_plotted = length(unique(rows$setting)),
  plotting_source_sha256 = .challenge_sha256(sub("^--file=", "", file_arg[1L])),
  plotting_helpers = .campaign_artifact_manifest(reader_source),
  inputs = .campaign_artifact_manifest(inputs), outputs = .campaign_artifact_manifest(outputs),
  interpretation = "Setting-specific endpoint views, alphabetical method order, no overall rank; existing paired intervals only")
jsonlite::write_json(figure_manifest, file.path(output, "figure-manifest.json"), pretty = TRUE,
                    auto_unbox = TRUE, null = "null", na = "null")
cat("Exported frozen-evidence figures to", normalizePath(output), "\n")
}
main()
