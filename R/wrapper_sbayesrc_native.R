# Official GCTB SBayesRC: one external-file call for the complete supplied panel.
# Source audited at GCTB cc7fa7d765c83a89c6375946cf77fe50ba1a317e;
# the release executable is fingerprinted separately, not attributed to that commit.
# This adapter deliberately does not convert per-locus simulated matrices into
# purported genome-wide inputs or call the experimental R surrogate.

#' Run official SBayesRC jointly on a prepared panel
#'
#' Calls the official GCTB program once using its eigen-LD, COJO summary and
#' annotation inputs. Inputs must already be harmonized and frozen. This is a
#' separate interface from the region dispatcher in \code{run_methods()}.
#' A declared genome-wide run must contain all 22 autosomes; that declaration
#' alone does not establish high-density coverage or inferential convergence.
#'
#' @param gwas_file Path to prepared COJO-format summary statistics.
#' @param ld_dir Path containing \code{snp.info}, \code{ldm.info}, and every
#'   referenced \code{block<BLOCK>.eigen.bin} file.
#' @param annotation_file Path to annotations: SNP, Intercept, then categories.
#' @param output_dir New or empty directory for native outputs and provenance.
#' @param scale Required declaration: \code{genome_wide}, \code{chromosome},
#'   or \code{subset}. The entire supplied panel is jointly modeled in all cases.
#' @param provenance Named list with nonempty \code{gwas_source},
#'   \code{ld_source}, \code{annotation_source}, \code{genome_build},
#'   \code{ancestry}, \code{panel_id}, and \code{preprocessing} descriptions.
#' @param gctb_path Official GCTB executable, or its name on PATH.
#' @param expected_version Required exact upstream version; default 2.5.5.
#' @param variant_ids Optional unique output subset in the requested order.
#'   This changes only imported output, never the full panel fitted upstream.
#' @param n_iter Total native MCMC iterations. Default 3000.
#' @param burn_in Native burn-in iterations. Default 1000.
#' @param seed Positive integer native random seed. Default 22.
#' @param threads Positive integer thread count. Default 1.
#' @param save_mcmc Save upstream text parameter draws for convergence diagnostics.
#'   Default FALSE. SNP-effect binary draws are not requested by this option.
#' @param thin Positive integer native posterior sampling interval. Default 10.
#'   Always passed to GCTB, including when text output is disabled.
#' @param mcmc_scope Text retention policy, \code{all} (default) or \code{core}.
#'   With \code{save_mcmc = TRUE}, \code{core} retains core, annotation and block
#'   text traces while redirecting seven full-panel text arrays to \code{/dev/null}.
#'   This Linux-only policy leaves the native executable and inference unchanged;
#'   formatting costs remain. Per-variant draws cannot be recovered from this run.
#' @param dry_run Logical. Validate/fingerprint inputs and return a command plan
#'   without launching the native program. Default FALSE.
#'
#' @return On execution, standard PIP output and aligned posterior mean effects,
#'   with scale, checksums, upstream version, full modeled variant count, and
#'   native logs. Native credible sets are not synthesized. A dry run returns
#'   its command plan with \code{native_executed = FALSE}.
#' @export
run_sbayesrc_native <- function(gwas_file, ld_dir, annotation_file, output_dir,
                                scale, provenance, gctb_path = "gctb",
                                expected_version = "2.5.5", variant_ids = NULL,
                                n_iter = 3000L, burn_in = 1000L, seed = 22L,
                                threads = 1L, dry_run = FALSE, save_mcmc = FALSE,
                                thin = 10L, mcmc_scope = c("all", "core")) {
  scale <- match.arg(scale, c("genome_wide", "chromosome", "subset"))
  mcmc_scope <- match.arg(mcmc_scope)
  .sbayesrc_native_integer(n_iter, "n_iter", 1L)
  .sbayesrc_native_integer(burn_in, "burn_in", 0L)
  .sbayesrc_native_integer(seed, "seed", 1L)
  .sbayesrc_native_integer(threads, "threads", 1L)
  .sbayesrc_native_integer(thin, "thin", 1L)
  if (!is.logical(save_mcmc) || length(save_mcmc) != 1L || is.na(save_mcmc)) {
    stop("save_mcmc must be TRUE or FALSE.", call. = FALSE)
  }
  if (mcmc_scope == "core" && !save_mcmc) {
    stop("mcmc_scope = core requires save_mcmc = TRUE.", call. = FALSE)
  }
  if (burn_in >= n_iter) stop("burn_in must leave retained MCMC iterations.", call. = FALSE)
  if (!is.logical(dry_run) || length(dry_run) != 1L || is.na(dry_run)) {
    stop("dry_run must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(expected_version) || length(expected_version) != 1L ||
      is.na(expected_version) || !nzchar(expected_version)) {
    stop("expected_version must identify a pinned GCTB version.", call. = FALSE)
  }
  if (mcmc_scope == "core" && expected_version != "2.5.5") {
    stop("Core text retention has been audited only for GCTB 2.5.5.", call. = FALSE)
  }
  required_provenance <- c("gwas_source", "ld_source", "annotation_source",
                           "genome_build", "ancestry", "panel_id", "preprocessing")
  if (!is.list(provenance) || any(vapply(required_provenance, function(k) {
    x <- provenance[[k]]
    !is.character(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(x))
  }, logical(1)))) {
    stop("provenance must describe GWAS, LD, annotations, build, ancestry, panel and preprocessing.", call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir)) stop("output_dir must be one directory path.", call. = FALSE)
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("output_dir must be new or empty; existing fits are never reused.", call. = FALSE)
  }

  prepared <- .sbayesrc_native_validate_inputs(gwas_file, ld_dir, annotation_file, scale)
  reference <- prepared$reference
  if (is.null(variant_ids)) variant_ids <- reference$SNP
  if (!is.character(variant_ids) || !length(variant_ids) || anyNA(variant_ids) ||
      anyDuplicated(variant_ids) || any(!variant_ids %in% reference$SNP)) {
    stop("variant_ids must be a unique nonempty subset of the modeled LD panel.", call. = FALSE)
  }
  output_dir <- path.expand(output_dir)
  if (!grepl("^/|^[A-Za-z]:[/\\\\]", output_dir)) output_dir <- file.path(getwd(), output_dir)
  prefix <- file.path(output_dir, "sbayesrc_native")
  args <- c("--ldm-eigen", prepared$ld_dir, "--gwas-summary", prepared$gwas_file,
            "--sbayes", "RC", "--annot", prepared$annotation_file,
            "--chain-length", as.character(n_iter), "--burn-in", as.character(burn_in),
            "--seed", as.character(seed), "--thread", as.character(threads),
            "--thin", as.character(thin), "--out", prefix)
  if (save_mcmc) args <- c(args, "--write-mcmc-txt")
  manifest <- list(
    implementation = "official_gctb_sbayesrc", upstream_mode = "sbayes_RC",
    expected_version = expected_version,
    upstream_source = "https://github.com/jianzeng/GCTB/tree/cc7fa7d765c83a89c6375946cf77fe50ba1a317e",
    upstream_source_role = "reference source audit; release compiled commit not established",
    scale = scale, effect_scale = "original_per_GWAS_A1_allele",
    modeled_variant_count = nrow(reference),
    requested_variant_count = length(variant_ids), chromosomes = prepared$chromosomes,
    ld_block_count = prepared$ld_block_count, provenance = provenance,
    input_checksums = .sbayesrc_native_snapshot(prepared$files), checksum_algorithm = "md5",
    params = list(n_iter = n_iter, burn_in = burn_in, seed = seed, threads = threads,
                  save_mcmc = save_mcmc, thin = thin, mcmc_scope = mcmc_scope),
    text_retention = list(scope = mcmc_scope,
      discarded_arrays = if (mcmc_scope == "core") .sbayesrc_native_text_sinks() else character(),
      sink_target = if (mcmc_scope == "core") "/dev/null" else NULL,
      formatting_cpu_avoided = FALSE),
    native_executed = FALSE, args = args
  )
  if (dry_run) {
    return(list(command = gctb_path, args = args, output_prefix = prefix,
                manifest = manifest, native_executed = FALSE))
  }

  executable <- if (file.exists(gctb_path)) normalizePath(gctb_path, mustWork = TRUE) else
    unname(Sys.which(gctb_path))
  if (!nzchar(executable) || file.access(executable, 1) != 0) {
    stop("Official GCTB executable is unavailable; use dry_run for an unexecuted plan.", call. = FALSE)
  }
  # main.cpp prints its version before processing options; no invented
  # --version or --help flag is used. A no-argument probe normally exits 1.
  banner <- suppressWarnings(system2(executable, stdout = TRUE, stderr = TRUE))
  version_line <- grep("GCTB [0-9]+[.][0-9]+", banner, value = TRUE)
  version <- if (length(version_line)) {
    sub(".*GCTB ([0-9]+[.][0-9]+([.][0-9]+)?).*", "\\1", version_line[1L])
  } else NA_character_
  if (is.na(version) || !identical(version, expected_version)) {
    stop("GCTB banner does not match expected_version; refusing an unpinned run.", call. = FALSE)
  }
  manifest$actual_version <- version
  manifest$executable <- executable
  manifest$executable_checksum <- unname(tools::md5sum(executable))
  if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
    stop("output_dir changed while inputs were validated; refusing possible output reuse.", call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lock_dir <- file.path(output_dir, ".native_run_lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) stop("Another native run owns output_dir.", call. = FALSE)
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  if (mcmc_scope == "core") {
    manifest$text_retention$sinks <- .sbayesrc_native_prepare_text_sinks(prefix)
  }
  saveRDS(manifest, file.path(output_dir, "native_input_manifest.rds"))
  log_path <- file.path(output_dir, "native_process.log")
  started <- Sys.time()
  status <- system2(executable, args = vapply(args, shQuote, character(1)),
                    stdout = log_path, stderr = log_path)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (!identical(as.integer(status), 0L)) {
    stop("Native GCTB failed; inspect ", log_path, ". No result was imported.", call. = FALSE)
  }
  if (mcmc_scope == "core" &&
      !all(Sys.readlink(manifest$text_retention$sinks$path) == "/dev/null")) {
    stop("Native text retention sinks changed during execution; results are rejected.", call. = FALSE)
  }
  if (!identical(manifest$input_checksums, .sbayesrc_native_snapshot(prepared$files))) {
    stop("A frozen input changed during execution; results are rejected.", call. = FALSE)
  }
  if (!identical(manifest$executable_checksum, unname(tools::md5sum(executable)))) {
    stop("The pinned executable changed during execution; results are rejected.", call. = FALSE)
  }
  manifest$native_executed <- TRUE
  manifest$runtime_seconds <- elapsed
  manifest$completed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  native_result <- paste0(prefix, ".snpRes")
  imported <- .sbayesrc_native_import(native_result, reference, variant_ids)
  manifest$result_checksum <- unname(tools::md5sum(native_result))
  saveRDS(manifest, file.path(output_dir, "native_completed_manifest.rds"))
  list(
    pip = imported$pip, credible_sets = list(), method = "sbayesrc_native",
    input_type = "genome_wide_summary_files", params = manifest$params,
    runtime_seconds = elapsed,
    additional = list(
      posterior_mean_beta = imported$beta, posterior_sd_beta = imported$se,
      pip_event = "any_nonzero_mixture_component",
      variant_ids = variant_ids, effect_allele = imported$effect_allele,
      effects_oriented_to = "prepared_GWAS_A1", native_executed = TRUE,
      implementation = "official_gctb_sbayesrc", upstream_mode = "sbayes_RC",
      upstream_equivalent = TRUE,
      credible_set_semantics = "not_provided", convergence_claim = "not_established",
      manifest = manifest, process_log = log_path, native_result_file = native_result,
      native_parameter_file = if (file.exists(paste0(prefix, ".parRes"))) paste0(prefix, ".parRes") else NULL
    )
  )
}

# These filenames are the seven full-width arrays in the audited five-mixture
# RC implementation. CoreParameters and small annotation/block traces remain.
.sbayesrc_native_text_sinks <- function() {
  paste0(c("PIP", "PEP", paste0("DeltaPi", 1:5)), ".mcmcsamples.txt")
}

.sbayesrc_native_prepare_text_sinks <- function(prefix) {
  if (!identical(unname(Sys.info()["sysname"]), "Linux") ||
      !file.exists("/dev/null") || dir.exists("/dev/null")) {
    stop("Core text retention requires Linux /dev/null.", call. = FALSE)
  }
  directory <- paste0(prefix, ".mcmcsamples")
  if (!dir.create(directory, showWarnings = FALSE)) {
    stop("Native text trace directory must be new when creating retention sinks.", call. = FALSE)
  }
  paths <- file.path(directory, .sbayesrc_native_text_sinks())
  if (!all(file.symlink("/dev/null", paths)) ||
      !all(Sys.readlink(paths) == "/dev/null")) {
    stop("Could not create the exact audited native text retention sinks.", call. = FALSE)
  }
  data.frame(path = paths, target = "/dev/null", stringsAsFactors = FALSE)
}

.sbayesrc_native_integer <- function(x, name, lower) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) ||
      x < lower || x > .Machine$integer.max) {
    stop(name, " must be an integer >= ", lower, ".", call. = FALSE)
  }
}

.sbayesrc_native_read <- function(path, columns = NULL) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, select = columns, showProgress = FALSE,
                                         data.table = FALSE)))
  }
  if (file.info(path)$size > 128 * 1024^2) {
    stop("Install data.table to validate/import native panel files larger than 128 MiB.", call. = FALSE)
  }
  d <- utils::read.table(path, header = TRUE, stringsAsFactors = FALSE,
                         check.names = FALSE, comment.char = "", quote = "")
  if (is.null(columns)) d else d[, columns, drop = FALSE]
}

.sbayesrc_native_snapshot <- function(files) {
  hashes <- tools::md5sum(files)
  if (anyNA(hashes)) stop("Could not fingerprint every frozen input.", call. = FALSE)
  data.frame(path = files, size_bytes = unname(file.info(files)$size),
             checksum = unname(hashes), stringsAsFactors = FALSE)
}

.sbayesrc_native_validate_inputs <- function(gwas_file, ld_dir, annotation_file, scale) {
  for (path in c(gwas_file, annotation_file)) {
    if (length(path) != 1L || is.na(path) || !file.exists(path) || dir.exists(path)) {
      stop("GWAS and annotation inputs must be existing files.", call. = FALSE)
    }
  }
  if (length(ld_dir) != 1L || is.na(ld_dir) || !dir.exists(ld_dir)) {
    stop("ld_dir must contain a prepared native eigen-LD panel.", call. = FALSE)
  }
  gwas_file <- normalizePath(gwas_file, mustWork = TRUE)
  annotation_file <- normalizePath(annotation_file, mustWork = TRUE)
  ld_dir <- normalizePath(ld_dir, mustWork = TRUE)
  info_files <- file.path(ld_dir, c("snp.info", "ldm.info"))
  if (!all(file.exists(info_files))) stop("Native LD requires snp.info and ldm.info.", call. = FALSE)
  gwas <- .sbayesrc_native_read(gwas_file)
  required <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
  if (!identical(names(gwas), required) || !nrow(gwas)) {
    stop("GWAS columns must be exactly SNP A1 A2 freq b se p N, in COJO order.", call. = FALSE)
  }
  if (anyNA(gwas$SNP) || any(!nzchar(gwas$SNP)) || anyDuplicated(gwas$SNP)) {
    stop("GWAS SNP IDs must be unique and nonmissing.", call. = FALSE)
  }
  numeric_fields <- c("freq", "b", "se", "p", "N")
  if (any(vapply(gwas[numeric_fields], function(x) !is.numeric(x) || any(!is.finite(x)), logical(1))) ||
      any(gwas$freq <= 0 | gwas$freq >= 1) || any(gwas$se <= 0) ||
      any(gwas$p < 0 | gwas$p > 1) || any(gwas$N < 3 | gwas$N != floor(gwas$N))) {
    stop("GWAS frequency, effect, SE, p-value or sample size is invalid.", call. = FALSE)
  }
  panel <- .sbayesrc_native_read(info_files[1L])
  panel_columns <- c("Chrom", "ID", "Index", "GenPos", "PhysPos", "A1", "A2", "A1Freq", "N", "Block")
  if (!identical(names(panel), panel_columns) || !nrow(panel) ||
      anyNA(panel$ID) || anyDuplicated(panel$ID)) {
    stop("snp.info must use the official ten-column order with unique ID values.", call. = FALSE)
  }
  numeric_panel <- c("Chrom", "Index", "GenPos", "PhysPos", "A1Freq", "N", "Block")
  if (any(vapply(panel[numeric_panel], function(x) !is.numeric(x) || any(!is.finite(x)), logical(1))) ||
      any(panel$Chrom < 1 | panel$Chrom != floor(panel$Chrom)) ||
      !identical(as.numeric(panel$Index), as.numeric(seq_len(nrow(panel)) - 1L)) ||
      any(panel$A1Freq <= 0 | panel$A1Freq >= 1) || any(panel$N < 3 | panel$N != floor(panel$N))) {
    stop("snp.info has invalid chromosomes, zero-based indices, frequencies or sample sizes.", call. = FALSE)
  }
  ids <- as.character(panel$ID)
  if (!setequal(ids, as.character(gwas$SNP))) {
    stop("Prepared GWAS and LD SNP sets must match exactly; tidy/impute externally first.", call. = FALSE)
  }
  gwas <- gwas[match(ids, gwas$SNP), , drop = FALSE]
  valid_alleles <- function(a1, a2) !anyNA(a1) && !anyNA(a2) &&
    all(a1 %in% c("A", "C", "G", "T") & a2 %in% c("A", "C", "G", "T") & a1 != a2)
  if (!valid_alleles(gwas$A1, gwas$A2) || !valid_alleles(panel$A1, panel$A2)) {
    stop("Native inputs require explicit distinct A/C/G/T allele pairs.", call. = FALSE)
  }
  aligned <- panel$A1 == gwas$A1 & panel$A2 == gwas$A2
  reversed <- panel$A1 == gwas$A2 & panel$A2 == gwas$A1
  if (any(!(aligned | reversed))) stop("GWAS and LD alleles are incompatible.", call. = FALSE)
  chromosomes <- sort(unique(sub("^chr", "", tolower(as.character(panel$Chrom)))))
  if (anyNA(chromosomes) || any(!nzchar(chromosomes))) stop("LD chromosomes are missing.", call. = FALSE)
  if (scale == "genome_wide" && !all(as.character(1:22) %in% chromosomes)) {
    stop("A genome_wide declaration requires all 22 autosomes in the modeled panel.", call. = FALSE)
  }
  if (scale == "chromosome" && length(chromosomes) != 1L) {
    stop("A chromosome declaration requires exactly one chromosome.", call. = FALSE)
  }
  blocks <- .sbayesrc_native_read(info_files[2L])
  block_columns <- c("Block", "Chrom", "StartSnpIdx", "StartSnpID", "EndSnpIdx", "EndSnpID", "NumSnps")
  if (!identical(names(blocks), block_columns) || !nrow(blocks)) {
    stop("ldm.info must use the official seven-column GCTB format; convert legacy metadata externally.", call. = FALSE)
  }
  block_ids <- blocks$Block
  if (!is.numeric(block_ids) || any(!is.finite(block_ids)) || any(block_ids != floor(block_ids)) ||
      any(block_ids < 1) || anyDuplicated(block_ids)) stop("LD block identifiers are invalid.", call. = FALSE)
  if (!setequal(panel$Block, block_ids)) stop("snp.info and ldm.info block memberships differ.", call. = FALSE)
  for (i in seq_len(nrow(blocks))) {
    rows <- which(panel$Block == block_ids[i])
    if (anyNA(blocks[i, ]) || !length(rows) || !identical(rows, seq.int(min(rows), max(rows))) ||
        blocks$StartSnpIdx[i] != min(rows) - 1L || blocks$EndSnpIdx[i] != max(rows) - 1L ||
        blocks$NumSnps[i] != length(rows) || blocks$StartSnpID[i] != ids[min(rows)] ||
        blocks$EndSnpID[i] != ids[max(rows)] || any(panel$Chrom[rows] != blocks$Chrom[i])) {
      stop("LD block boundaries, chromosome or SNP counts do not match snp.info.", call. = FALSE)
    }
  }
  eigen_files <- file.path(ld_dir, paste0("block", block_ids, ".eigen.bin"))
  if (!all(file.exists(eigen_files)) || any(file.info(eigen_files)$size <= 0)) {
    stop("Every referenced native eigen-LD block must exist and be nonempty.", call. = FALSE)
  }
  annotation <- .sbayesrc_native_read(annotation_file)
  if (ncol(annotation) < 2L || !identical(names(annotation)[1:2], c("SNP", "Intercept")) ||
      anyDuplicated(names(annotation)) || anyNA(annotation$SNP) || anyDuplicated(annotation$SNP) ||
      !setequal(ids, as.character(annotation$SNP)) ||
      any(vapply(annotation[-1L], function(x) !is.numeric(x) || any(!is.finite(x)), logical(1))) ||
      any(annotation$Intercept != 1)) {
    stop("Annotations require the exact LD SNP set, a unit Intercept, and finite numeric categories.", call. = FALSE)
  }
  list(gwas_file = gwas_file, ld_dir = ld_dir, annotation_file = annotation_file,
       files = c(gwas_file, annotation_file, info_files, eigen_files),
       reference = data.frame(SNP = ids, A1 = gwas$A1, A2 = gwas$A2, stringsAsFactors = FALSE),
       chromosomes = chromosomes, ld_block_count = length(block_ids))
}

# Official GCTB output schema, as written by Data::outputSnpResults.
# Match by stable SNP ID, never by row position, and orient effects to GWAS A1.
.sbayesrc_native_import <- function(path, reference, variant_ids) {
  if (!is.character(variant_ids) || !length(variant_ids) || anyNA(variant_ids) ||
      anyDuplicated(variant_ids) || any(!variant_ids %in% reference$SNP)) {
    stop("variant_ids must be a unique subset of the native modeled panel.", call. = FALSE)
  }
  if (!file.exists(path)) stop("Native GCTB produced no .snpRes output.", call. = FALSE)
  d <- .sbayesrc_native_read(path)
  columns <- c("Name", "A1", "A2", "A1Effect", "SE", "PIP")
  if (!all(columns %in% names(d)) || anyNA(d$Name) || anyDuplicated(d$Name) ||
      !setequal(as.character(d$Name), reference$SNP)) {
    stop("Native result IDs/schema do not match the complete frozen modeled panel.", call. = FALSE)
  }
  d <- d[match(reference$SNP, d$Name), , drop = FALSE]
  same <- d$A1 == reference$A1 & d$A2 == reference$A2
  flip <- d$A1 == reference$A2 & d$A2 == reference$A1
  if (anyNA(same) || anyNA(flip) || any(!(same | flip)) ||
      any(vapply(d[c("A1Effect", "SE", "PIP")], function(x) !is.numeric(x) || any(!is.finite(x)), logical(1))) ||
      any(d$SE < 0) || any(d$PIP < 0 | d$PIP > 1)) {
    stop("Native result contains incompatible alleles or invalid PIPs/effects/SEs.", call. = FALSE)
  }
  idx <- match(variant_ids, reference$SNP)
  list(pip = as.numeric(d$PIP[idx]),
       beta = as.numeric(d$A1Effect[idx]) * ifelse(flip[idx], -1, 1),
       se = as.numeric(d$SE[idx]), effect_allele = reference$A1[idx])
}
