# Native campaign helpers; sources only public-input code, never dispatches jobs.
.native_profile <- function(profile) {
  switch(profile,
    full = list(n = 5000L, p = 200L, draws = 5L, architectures = c("sparse", "sparse_inf"),
      pi = c(.001, .01), h2 = c(.1, .3, .5), replicates = 3L, regions = 128L,
      n_iter = 25000L, burn_in = 5000L, thin = 10L),
    canary = list(n = 5000L, p = 200L, draws = 1L, architectures = c("sparse", "sparse_inf"),
      pi = .001, h2 = .3, replicates = 1L, regions = 128L,
      n_iter = 500L, burn_in = 100L, thin = 10L),
    integration = list(n = 200L, p = 24L, draws = 1L, architectures = "sparse",
      pi = .03, h2 = .2, replicates = 1L, regions = 2L,
      n_iter = 150L, burn_in = 50L, thin = 5L), stop("Unknown profile"))
}

.native_design <- function(profile, seed) {
  config <- .native_profile(profile)
  draws <- expand.grid(draw = seq_len(config$draws), model = config$architectures,
                       stringsAsFactors = FALSE)
  draws$draw_id <- sprintf("draw_%03d", seq_len(nrow(draws)))
  draws$seed <- seed + seq_len(nrow(draws)) * 10000L
  inputs <- list(); k <- 0L
  for (i in seq_len(nrow(draws))) for (pi in config$pi) for (h2 in config$h2)
    for (replicate in seq_len(config$replicates)) for (exposure in c("none", "informative", "permuted")) {
      k <- k + 1L
      inputs[[k]] <- data.frame(input_id = sprintf("input_%05d", k), draw_id = draws$draw_id[i],
        draw = draws$draw[i], model = draws$model[i], pi = pi, h2 = h2,
        replicate = replicate, exposure = exposure, stringsAsFactors = FALSE)
    }
  inputs <- do.call(rbind, inputs)
  chains <- inputs[rep(seq_len(nrow(inputs)), each = 4L), ]
  chains$chain <- rep(1:4, nrow(inputs))
  chains$chain_id <- sprintf("chain_%06d", seq_len(nrow(chains)))
  chains$seed <- seed + 1000000L + seq_len(nrow(chains)) * 37L
  list(config = config, draws = draws, inputs = inputs, chains = chains)
}

.native_source_snapshot <- function(repo) {
  paths <- c(list.files(file.path(repo, "R"), "\\.R$", full.names = TRUE),
             list.files(file.path(repo, "BEATRICE_annot_sparse"), "\\.py$", recursive = TRUE, full.names = TRUE),
             list.files(file.path(repo, "scripts/native/campaign"), full.names = TRUE))
  paths <- sort(paths[!dir.exists(paths)], method = "radix")
  data.frame(path = unname(substring(paths, nchar(repo) + 2L)),
             sha256 = unname(vapply(paths, .challenge_sha256, "")), row.names = NULL,
             stringsAsFactors = FALSE)
}

.native_assert_plan <- function(plan, repo, verify_source_data = FALSE) {
  x <- plan; x$plan_sha256 <- NULL
  if (!identical(plan$plan_sha256, .campaign_object_sha(x))) stop("Native plan changed")
  actual_source <- .native_source_snapshot(repo)
  if (!identical(plan$source, actual_source)) {
    changed <- merge(plan$source, actual_source, by = "path", all = TRUE,
                     suffixes = c("_expected", "_actual"))
    changed <- changed[is.na(changed$sha256_expected) | is.na(changed$sha256_actual) |
                       changed$sha256_expected != changed$sha256_actual, , drop = FALSE]
    stop("Native source changed in ", repo, ": ",
         if (nrow(changed)) paste(changed$path, collapse = ",") else "snapshot path attributes differ")
  }
  if (verify_source_data && (any(!file.exists(plan$data_snapshot$path)) ||
      !identical(unname(plan$data_snapshot$sha256), unname(vapply(plan$data_snapshot$path, .challenge_sha256, ""))))
    )
    stop("Native source VCF/map/region bytes changed")
  invisible(TRUE)
}

.native_write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, quote = FALSE, row.names = FALSE, sep = " ")
}

.native_variant_table <- function(genotypes) {
  xs <- lapply(seq_along(genotypes), function(i) {
    g <- genotypes[[i]]; m <- g$variant_metadata
    if (!is.data.frame(m) || !all(c("chrom", "position", "source_id", "ref", "alt") %in% names(m)) ||
        nrow(m) != ncol(g$X_raw) || !identical(g$counted_allele, "ALT"))
      stop("Verified retained-VCF ALT allele metadata is required")
    if (anyNA(m) || any(!m$ref %in% c("A", "C", "G", "T")) ||
        any(!m$alt %in% c("A", "C", "G", "T")) || any(m$ref == m$alt))
      stop("Native campaign requires biallelic A/C/G/T variants")
    chromosome <- as.integer(sub("^chr", "", as.character(m$chrom), ignore.case = TRUE))
    if (anyNA(chromosome) || any(!chromosome %in% 1:22) ||
        any(!is.finite(m$position)) || any(m$position < 1 | m$position != floor(m$position)))
      stop("Invalid autosomal variant coordinates")
    sd <- apply(g$X_raw, 2L, stats::sd); freq <- colMeans(g$X_raw) / 2
    if (any(!is.finite(sd) | sd <= 0) || any(!is.finite(freq) | freq <= 0 | freq >= 1))
      stop("Degenerate native allele frequency or genotype SD")
    if (max(abs(scale(g$X_raw) - g$X)) > 1e-7) stop("Genotype standardisation cannot be reversed")
    data.frame(SNP = paste(chromosome, m$position, m$ref, m$alt, sep = ":"),
       chrom = chromosome, position = m$position, source_id = m$source_id,
       A1 = m$alt, A2 = m$ref, freq = freq, genotype_sd = sd,
       region_id = i, region_label = g$region_label %||% as.character(i),
       local_index = seq_len(nrow(m)), source_vcf = g$vcf_source, stringsAsFactors = FALSE)
  })
  x <- do.call(rbind, xs)
  if (anyDuplicated(x$SNP)) stop("Overlapping source blocks duplicate modeled SNPs")
  x
}

.native_write_eigen <- function(R, path) {
  if (!is.matrix(R) || any(!is.finite(R)) || max(abs(R - t(R))) > 1e-8 ||
      max(abs(diag(R) - 1)) > 1e-8) stop("Invalid empirical correlation")
  e <- eigen(R, symmetric = TRUE)
  if (min(e$values) < -1e-7) stop("Empirical LD is not positive semidefinite")
  keep <- e$values > 1e-8
  if (!any(keep)) stop("No positive empirical LD eigenvalues")
  values <- e$values[keep]; vectors <- e$vectors[, keep, drop = FALSE]
  con <- file(path, "wb")
  tryCatch({
    writeBin(as.integer(c(nrow(R), length(values))), con, size = 4L, endian = "little")
    writeBin(c(sum(pmax(e$values, 0)), sum(values) / sum(pmax(e$values, 0))), con,
             size = 4L, endian = "little")
    writeBin(values, con, size = 4L, endian = "little")
    writeBin(as.numeric(vectors), con, size = 4L, endian = "little")
  }, finally = close(con))
  decoded <- .native_read_eigen(path)
  error <- max(abs(decoded - R))
  if (error > 2e-6) stop("Float32 empirical eigen LD reconstruction failed")
  data.frame(variants = nrow(R), retained_rank = length(values),
             variance_retained = sum(values) / sum(pmax(e$values, 0)), max_error = error)
}

.native_read_eigen <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  dimensions <- readBin(con, integer(), 2L, size = 4L, endian = "little")
  header <- readBin(con, numeric(), 2L, size = 4L, endian = "little")
  values <- readBin(con, numeric(), dimensions[2L], size = 4L, endian = "little")
  vectors <- matrix(readBin(con, numeric(), prod(dimensions), size = 4L, endian = "little"),
                    dimensions[1L], dimensions[2L])
  if (length(values) != dimensions[2L] || any(!is.finite(header))) stop("Truncated native eigen file")
  tcrossprod(sweep(vectors, 2L, sqrt(values), "*"))
}

.native_annotation_controls <- function(sim, seed) {
  useful <- lapply(sim$genotypes, function(g) g$annotations_matrix)
  if (any(vapply(useful, function(a) !is.matrix(a) || ncol(a) != 3L, TRUE)))
    stop("Native truth requires three informative binary annotations")
  set.seed(seed)
  permuted <- lapply(useful, function(a) a[sample.int(nrow(a)), , drop = FALSE])
  list(none = lapply(useful, function(a) matrix(numeric(), nrow(a), 0L)),
       informative = useful, permuted = permuted)
}

.native_prepare_draw <- function(sim, draw, plan, directory) {
  if (length(sim$genotypes) != plan$config$regions ||
      any(vapply(sim$genotypes, function(g)
        nrow(g$X_raw) != plan$config$n || ncol(g$X_raw) != plan$config$p, TRUE)))
    stop("Native genotype dimensions differ from the frozen profile")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  variants <- .native_variant_table(sim$genotypes)
  ld <- file.path(directory, "ld"); dir.create(ld, showWarnings = FALSE)
  block <- variants$region_id
  info <- data.frame(Chrom = variants$chrom, ID = variants$SNP,
    Index = seq_len(nrow(variants)) - 1L, GenPos = 0, PhysPos = variants$position,
    A1 = variants$A1, A2 = variants$A2, A1Freq = variants$freq,
    N = plan$config$n, Block = block)
  blocks <- do.call(rbind, lapply(unique(block), function(b) {
    idx <- which(block == b)
    data.frame(Block = b, Chrom = unique(variants$chrom[idx]), StartSnpIdx = min(idx) - 1L,
      StartSnpID = variants$SNP[min(idx)], EndSnpIdx = max(idx) - 1L,
      EndSnpID = variants$SNP[max(idx)], NumSnps = length(idx))
  }))
  .native_write_table(info, file.path(ld, "snp.info"))
  .native_write_table(blocks, file.path(ld, "ldm.info"))
  quality <- do.call(rbind, lapply(seq_along(sim$genotypes), function(i) {
    z <- .native_write_eigen(sim$genotypes[[i]]$LD, file.path(ld, paste0("block", i, ".eigen.bin")))
    z$block <- i; z
  }))
  write.csv(quality, file.path(directory, "eigen-quality.csv"), row.names = FALSE)
  write.csv(variants, file.path(directory, "variant-map.csv"), row.names = FALSE)
  controls <- .native_annotation_controls(sim, draw$seed + 91L)
  for (name in names(controls)) {
    a <- do.call(rbind, controls[[name]])
    ann <- data.frame(SNP = variants$SNP, Intercept = 1L)
    if (ncol(a)) { colnames(a) <- paste0("annotation", 1:3); ann <- cbind(ann, a) }
    .native_write_table(ann, file.path(directory, paste0("annotations-", name, ".txt")))
  }
  paths <- list()
  for (i in seq_along(sim$scenarios)) {
    sc <- sim$scenarios[[i]]
    if (any(vapply(sc$regions, function(r) !identical(r$y, sc$regions[[1L]]$y), TRUE)))
      stop("Native input does not have one shared phenotype")
    bhat <- unlist(lapply(sc$regions, `[[`, "beta_hat"), use.names = FALSE)
    se <- unlist(lapply(sc$regions, `[[`, "se"), use.names = FALSE)
    z <- bhat / se
    gwas <- data.frame(SNP = variants$SNP, A1 = variants$A1, A2 = variants$A2,
      freq = variants$freq, b = bhat / variants$genotype_sd,
      se = se / variants$genotype_sd, p = 2 * stats::pt(-abs(z), df = plan$config$n - 2),
      N = plan$config$n)
    file <- file.path(directory, sprintf("scenario_%03d.ma", i)); .native_write_table(gwas, file)
    matching <- plan$inputs[plan$inputs$draw_id == draw$draw_id & plan$inputs$pi == sc$pi &
      plan$inputs$h2 == sc$h2 & plan$inputs$replicate == sc$iter, ]
    if (nrow(matching) != 3L) stop("Native design and simulation scenarios do not align")
    paths[[i]] <- data.frame(input_id = matching$input_id, scenario = i, gwas_file = file,
      ld_dir = ld, annotation_file = file.path(directory, paste0("annotations-", matching$exposure, ".txt")))
  }
  .campaign_atomic_rds(sim, file.path(directory, "simulation.rds"))
  .campaign_atomic_rds(controls, file.path(directory, "annotation-controls.rds"))
  paths <- do.call(rbind, paths); write.csv(paths, file.path(directory, "inputs.csv"), row.names = FALSE)
  files <- sort(list.files(directory, recursive = TRUE, full.names = TRUE), method = "radix")
  list(terminal = TRUE, task = draw, plan_sha256 = plan$plan_sha256, inputs = paths,
    artifacts = data.frame(path = unname(files), sha256 = unname(vapply(files, .challenge_sha256, "")),
                           row.names = NULL, stringsAsFactors = FALSE),
    completed = Sys.time(), session = utils::sessionInfo())
}

.native_read_prepared <- function(plan, root, draw_id, include_simulation = TRUE,
                                  require_receipt = TRUE) {
  path <- file.path(root, "prepared", draw_id, "completed.rds")
  draw <- plan$draws[match(draw_id, plan$draws$draw_id), ]
  receipted <- .campaign_check_receipt(dirname(path), path, plan$plan_sha256, draw_id)
  if (require_receipt && !receipted) stop("Native preparation has no verified terminal receipt")
  if (!.campaign_validate_checkpoint(path, plan, draw, validate_artifacts = FALSE))
    stop("Native preparation is not complete")
  z <- readRDS(path)
  artifacts <- z$artifacts
  if (!include_simulation) artifacts <- artifacts[!grepl("\\.rds$", artifacts$path), ]
  .campaign_validate_artifacts(artifacts)
  z
}

.native_completed_checkpoint <- function(path, plan, task, id) {
  if (!.campaign_check_receipt(dirname(path), path, plan$plan_sha256, id)) return(FALSE)
  .campaign_validate_checkpoint(path, plan, task)
}

.native_archive_analysis <- function(root) {
  analysis <- file.path(root, "analysis")
  if (!dir.exists(analysis)) return(invisible(NULL))
  history <- file.path(root, "analysis-history")
  dir.create(history, showWarnings = FALSE)
  destination <- file.path(history, paste0(format(Sys.time(), "%Y%m%dT%H%M%OS6", tz = "UTC"),
                                         "-", Sys.getpid()))
  if (!file.rename(analysis, destination)) stop("Cannot archive previous native analysis")
  invisible(destination)
}

.native_seal_analysis <- function(analysis, plan, runtime, status, coverage, started,
                                  error = NULL) {
  files <- list.files(analysis, recursive = TRUE, full.names = TRUE)
  files <- files[!dir.exists(files) & !basename(files) %in%
                    c("analysis-manifest.rds", "analysis-manifest.json", "terminal.json")]
  manifest <- list(schema = 1L, plan_sha256 = plan$plan_sha256,
    collection_status = status, campaign_complete = identical(status, "complete"),
    error = error, started = started, completed = Sys.time(), coverage = coverage,
    source = plan$source, runtime_fingerprint = runtime,
    frozen_runtime_fingerprint = plan$runtime_fingerprint,
    posterior_version = as.character(utils::packageVersion("posterior")),
    artifacts = .campaign_artifact_manifest(files))
  checkpoint <- file.path(analysis, "analysis-manifest.rds")
  json_path <- file.path(analysis, "analysis-manifest.json")
  jsonlite::write_json(manifest, json_path,
                        auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  manifest$artifacts <- .campaign_artifact_manifest(c(files, json_path))
  .campaign_atomic_rds(manifest, checkpoint)
  .campaign_receipt(analysis, checkpoint, plan$plan_sha256, "analysis")
  invisible(manifest)
}

.native_read_core <- function(path, expected_draws) {
  if (!file.exists(path)) stop("Missing upstream CoreParameters text trace")
  x <- utils::read.table(path, header = TRUE, check.names = FALSE, comment.char = "", quote = "")
  if (nrow(x) != expected_draws || !ncol(x) || anyDuplicated(names(x)) ||
      any(vapply(x, function(y) !is.numeric(y) || any(!is.finite(y)), TRUE)))
    stop("CoreParameters trace has unexpected schema, draw count, or nonfinite values")
  # GCTB writes retained post-burn-in draws. Do not discard burn-in a second time.
  x
}

.native_diagnostics <- function(traces, expected_draws, thresholds = list(rhat = 1.01, ess = 400)) {
  if (length(traces) != 4L) stop("Four independent native chains are required")
  xs <- lapply(traces, .native_read_core, expected_draws = expected_draws)
  if (!all(vapply(xs, function(x) identical(names(x), names(xs[[1L]])), TRUE)))
    stop("CoreParameters differ across chains")
  if (!requireNamespace("posterior", quietly = TRUE)) stop("Install official posterior package for diagnostics")
  do.call(rbind, lapply(names(xs[[1L]]), function(name) {
    x <- do.call(cbind, lapply(xs, `[[`, name))
    rhat <- posterior::rhat(x); bulk <- posterior::ess_bulk(x); tail <- posterior::ess_tail(x)
    variances <- apply(x, 2L, stats::var)
    constant <- any(variances == 0)
    passed <- is.finite(rhat) && is.finite(bulk) && is.finite(tail) &&
      rhat <= thresholds$rhat && bulk >= thresholds$ess && tail >= thresholds$ess
    data.frame(parameter = name, n_draws_per_chain = nrow(x),
      rhat = rhat, bulk_ess = bulk, tail_ess = tail,
      minimum_within_chain_variance = min(variances),
      maximum_within_chain_variance = max(variances), constant_in_chain = constant,
      diagnostic_status = if (constant) "constant_unassessable" else
        if (!all(is.finite(c(rhat, bulk, tail)))) "nonfinite_unassessable" else
        if (passed) "thresholds_passed" else "thresholds_failed",
      assessed_pass = passed)
  }))
}

.native_event_metrics <- function(pip, y) {
  called <- pip >= .95
  data.frame(n_variants = length(y), n_positive = sum(y), tp = sum(called & y),
    fp = sum(called & !y), any_call = any(called),
    recall = if (any(y)) sum(called & y) / sum(y) else NA_real_,
    ap = .challenge_ap(y, pip), brier = mean((pip - y)^2),
    log_loss = -mean(ifelse(y, log(pip), log1p(-pip))),
    log_loss_clipped = -mean(ifelse(y, log(pmax(pip, 1e-15)), log(pmax(1 - pip, 1e-15)))),
    sum_pip = sum(pip))
}

.native_expose_annotations <- function(sim, scenario, exposure, controls) {
  sim <- .campaign_one_scenario(sim, scenario)
  for (i in seq_along(sim$genotypes)) {
    a <- if (exposure == "none") NULL else controls[[exposure]][[i]]
    sim$genotypes[[i]]$annotations_matrix <- a
    sim$scenarios[[1L]]$regions[[i]]$annotations_matrix <- a
  }
  sim$params$annotation_type <- if (exposure == "none") "none" else "binary"
  sim$params$campaign_annotation_exposure <- exposure
  sim
}

.native_metric_rows <- function(sim, input, method, pip = NULL, success = FALSE,
                                attempted = FALSE, complete = FALSE, error = "",
                                runtime = NA_real_, pip_event = "modeled_sparse_effect_membership") {
  if (success && (length(pip) != sum(vapply(sim$genotypes, `[[`, 1L, "p")) ||
      any(!is.finite(pip)) || any(pip < 0 | pip > 1))) stop("Invalid whole-panel PIPs")
  sizes <- vapply(sim$genotypes, `[[`, 1L, "p"); offsets <- c(0L, cumsum(sizes))
  sc <- sim$scenarios[[1L]]
  rows <- list(); k <- 0L
  for (i in seq_along(sim$genotypes)) {
    truth <- sc$regions[[i]]$truth
    sparse <- seq_len(sizes[i]) %in% truth$causal_indices
    nonzero <- (truth$beta_total %||% (truth$beta_true + truth$beta_infinitesimal)) != 0
    for (event in c("sparse_component_membership", "any_nonzero_effect")) {
      y <- if (event == "sparse_component_membership") sparse else nonzero
      if (success) {
        p <- pip[seq.int(offsets[i] + 1L, offsets[i + 1L])]
        score <- .native_event_metrics(p, y)
      } else score <- data.frame(n_variants = sizes[i], n_positive = sum(y), tp = 0L,
        fp = 0L, any_call = FALSE, recall = if (any(y)) 0 else NA_real_,
        ap = NA_real_, brier = NA_real_, log_loss = NA_real_, log_loss_clipped = NA_real_, sum_pip = NA_real_)
      k <- k + 1L
      rows[[k]] <- cbind(input, data.frame(method = method, region_id = i,
        source_locus = sim$genotypes[[i]]$vcf_source,
        genotype_cluster = input$draw_id, event = event, pip_event = pip_event,
        probability_target_aligned = if (pip_event %in% c("sparse_component_membership", "modeled_sparse_effect_membership"))
          event == "sparse_component_membership" else if (pip_event %in%
            c("any_nonzero_effect", "any_nonzero_mixture_component")) event == "any_nonzero_effect" else FALSE,
        success = success,
        attempted = attempted, complete = complete, error = error, runtime_seconds = runtime), score)
    }
  }
  do.call(rbind, rows)
}

.native_regional_rows <- function(sim, input, method, result = NULL, attempted = FALSE) {
  if (is.null(result)) return(.native_metric_rows(sim, input, method, attempted = attempted,
                                                error = "missing_checkpoint"))
  fits <- .complete_fit_grid(result[[method]]$results, sim, method)
  rows <- lapply(seq_along(fits), function(i) {
    mini <- sim; mini$genotypes <- sim$genotypes[i]
    mini$scenarios[[1L]]$regions <- sim$scenarios[[1L]]$regions[i]
    f <- fits[[i]]
    event <- f$additional$pip_event %||% if (method == "sbayesrc") "any_nonzero_mixture_component" else
      "modeled_sparse_effect_membership"
    x <- .native_metric_rows(mini, input, method, pip = f$pip, success = is.null(f$error),
      attempted = TRUE, complete = TRUE, error = f$error %||% "", runtime = f$runtime_seconds %||% NA_real_,
      pip_event = event)
    x$region_id <- i
    if (method == "marginal_z") {
      x$brier <- x$log_loss <- x$log_loss_clipped <- NA_real_; x$pip_event <- "normalized_absolute_z_score"
      x$probability_target_aligned <- FALSE
    }
    x
  })
  do.call(rbind, rows)
}

.native_aggregate_draw <- function(rows) {
  mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  keys <- c("draw_id", "draw", "model", "pi", "h2", "exposure", "method", "event", "pip_event")
  groups <- split(seq_len(nrow(rows)), do.call(paste, c(rows[keys], sep = "|")))
  do.call(rbind, lapply(groups, function(indices) {
    x <- rows[indices, ]
    cbind(x[1L, keys], data.frame(n_planned = nrow(x), n_attempted = sum(x$attempted),
      n_complete = sum(x$complete), n_successful = sum(x$success),
      completion_rate = mean(x$complete), success_rate_planned = mean(x$success),
      n_ap_eligible_completed = sum(x$complete & x$n_positive > 0),
      n_ap_eligible_successful = sum(x$success & x$n_positive > 0),
      ap_successful = mean_or_na(x$ap[x$success]), brier_successful = mean_or_na(x$brier[x$success]),
      log_loss_successful = mean_or_na(x$log_loss[x$success]),
      n_infinite_log_loss = sum(is.infinite(x$log_loss)),
      recall_all_completed = mean_or_na(x$recall[x$complete]),
      recall_planned_to_date = mean_or_na(x$recall),
      ap_failure_zero_completed = mean_or_na(ifelse(x$success[x$complete & x$n_positive > 0],
        x$ap[x$complete & x$n_positive > 0], 0)),
      tp = sum(x$tp), fp = sum(x$fp), n_successful_null = sum(x$success & x$n_positive == 0),
      n_called_successful_null = sum(x$success & x$n_positive == 0 & x$any_call),
      sum_pip_successful = mean_or_na(x$sum_pip[x$success]),
      all_probability_targets_aligned = all(x$probability_target_aligned)))
  }))
}

.native_paired_draws <- function(a, b, metric, n_boot = 1000L, seed = 1L) {
  paired <- merge(a[c("draw_id", metric, "n_complete", "n_planned")],
    b[c("draw_id", metric, "n_complete", "n_planned")], by = "draw_id", suffixes = c("_a", "_b"))
  d <- paired[[paste0(metric, "_a")]] - paired[[paste0(metric, "_b")]]
  usable <- is.finite(d) & paired$n_complete_a == paired$n_planned_a &
    paired$n_complete_b == paired$n_planned_b
  d <- d[usable]; interval <- rep(NA_real_, 2L)
  if (length(d) >= 2L) {
    set.seed(seed)
    interval <- quantile(replicate(n_boot, mean(sample(d, length(d), replace = TRUE))),
                         c(.025, .975), names = FALSE)
  }
  data.frame(metric = metric, difference = if (length(d)) mean(d) else NA_real_,
    lower = interval[1L], upper = interval[2L], n_paired_genotype_draws = length(d),
    n_planned_draw_pairs = nrow(paired), conditional_success_metric = grepl("successful", metric))
}
