# =============================================================================
# sparsepro.R
#
# Wrapper for SparsePro (Zhang et al. 2023) fine-mapping.
#
# SparsePro is a Python CLI tool that performs variational fine-mapping using
# a sparse-projection prior on effect groups. Unlike Funmap (which exposes a
# Python module importable via reticulate), SparsePro is distributed as the
# script `sparsepro_zld.py` in the upstream GitHub repo, so this wrapper
# follows the BEATRICE / FINEMAP / PAINTOR pattern: write input files to a
# temp directory, invoke the script via system2(), parse the output files,
# and return results in the standardised format.
#
# The upstream repo is at https://github.com/zhwm/SparsePro. The user clones
# it once and passes the path via `sparsepro_dir`. See setup_sparsepro().
#
# This file provides:
#   - setup_sparsepro()       : verifies sparsepro_zld.py is reachable and that
#                                the Python env has numpy / scipy / pandas
#   - run_sparsepro()         : runs SparsePro on a single region (explicit
#                                inputs)
#   - run_sparsepro_region()  : adapter called by run_methods()
#
# Standard output format:
#   pip              Numeric vector (length p). Marginal PIPs from <prefix>.pip.
#   credible_sets    List of integer vectors (1-based). One per effect group
#                    in <prefix>.cs. Empty list if no credible sets reported.
#   method           Character. "sparsepro".
#   input_type       Character. Always "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. SparsePro-specific outputs (see run_sparsepro() docs).
#
# Reference:
#   Zhang W, Najafabadi H, Li Y (2023). SparsePro: an efficient fine-mapping
#   method integrating summary statistics and functional annotations.
#   PLoS Genetics, 19(12), e1011104.
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Set up SparsePro
#'
#' Verifies that the SparsePro script (\code{sparsepro_zld.py}) can be found
#' at \code{sparsepro_dir} and that the required Python packages
#' (\code{numpy}, \code{scipy}, \code{pandas}) are importable by the specified
#' Python executable.
#'
#' Install SparsePro by cloning the upstream repository and installing its
#' requirements into your Python environment:
#'
#' \preformatted{
#'   git clone https://github.com/zhwm/SparsePro
#'   cd SparsePro
#'   pip install -r requirements.txt
#' }
#'
#' SparsePro shares its Python dependencies (\code{numpy}, \code{scipy},
#' \code{pandas}) with BEATRICE and Funmap, so the existing benchmark conda
#' environment is sufficient — \code{pandas} may not be present in older
#' BEATRICE envs and can be added with \code{pip install pandas}.
#'
#' Then pass the repo directory and (if not on PATH) the Python executable:
#'
#' \preformatted{
#'   setup_sparsepro(
#'     sparsepro_dir = "~/SparsePro",
#'     python        = "~/anaconda3/envs/beatrice/bin/python"
#'   )
#' }
#'
#' @param sparsepro_dir Character. Path to the cloned SparsePro repository
#'   root (must contain \code{sparsepro_zld.py}).
#' @param python Character. Path to the Python executable to use.
#'   Default: \code{"python"} (searches PATH).
#'
#' @return Invisibly returns a named list with \code{sparsepro_script} (full
#'   path to \code{sparsepro_zld.py}) and \code{python} (resolved Python path).
#' @export
setup_sparsepro <- function(sparsepro_dir, python = "python") {

  sparsepro_dir <- path.expand(sparsepro_dir)

  script_path <- file.path(sparsepro_dir, "sparsepro_zld.py")
  if (!file.exists(script_path)) {
    stop(
      "sparsepro_zld.py not found in: ", sparsepro_dir, "\n\n",
      "Please clone the upstream repository:\n",
      "  git clone https://github.com/zhwm/SparsePro\n",
      "  cd SparsePro && pip install -r requirements.txt\n\n",
      "Then pass the path:\n",
      "  setup_sparsepro(sparsepro_dir = '/path/to/SparsePro')",
      call. = FALSE
    )
  }

  resolved_python <- if (file.exists(path.expand(python))) {
    .fm_python_path(python)
  } else {
    py <- Sys.which(python)
    if (nchar(py) == 0) {
      stop(
        "Python executable not found: '", python, "'\n\n",
        "Pass the full path to your Python:\n",
        "  setup_sparsepro(\n",
        "    sparsepro_dir = '/path/to/SparsePro',\n",
        "    python        = '~/anaconda3/envs/beatrice/bin/python'\n",
        "  )",
        call. = FALSE
      )
    }
    py
  }

  # Check required Python packages
  check_pkg <- function(pkg) {
    out <- system2(resolved_python,
                   args = c("-c", shQuote(paste0("import ", pkg))),
                   stdout = FALSE, stderr = FALSE)
    out == 0L
  }

  missing_pkgs <- Filter(Negate(check_pkg), c("numpy", "scipy", "pandas"))

  if (length(missing_pkgs) > 0) {
    stop(
      "The following Python packages are missing from '", resolved_python, "':\n",
      paste0("  ", missing_pkgs, collapse = "\n"), "\n\n",
      "Install SparsePro's requirements:\n",
      "  pip install -r ", file.path(sparsepro_dir, "requirements.txt"), "\n\n",
      "If you are using the BEATRICE conda env, you may need to add pandas:\n",
      "  conda run -n beatrice pip install pandas",
      call. = FALSE
    )
  }

  message("SparsePro ready.")
  message("  Script : ", script_path)
  message("  Python : ", resolved_python)

  invisible(list(sparsepro_script = script_path, python = resolved_python))
}


# =============================================================================
# Run SparsePro on a single region
# =============================================================================

#' Run SparsePro fine-mapping on a single region
#'
#' Writes the required input files to a temporary directory (z-score table,
#' LD matrix, and the --zld summary file SparsePro expects), invokes
#' \code{sparsepro_zld.py} via the specified Python executable, parses the
#' output \code{.pip} and \code{.cs} files, and returns results in the
#' standardised format. The temporary directory is deleted on exit.
#'
#' SparsePro fits one model per region (one entry in the --zld summary file).
#' The benchmark calls this wrapper one region at a time, so each call
#' produces a single-locus --zld file. This is less efficient than batching
#' many loci together, but keeps the wrapper interface symmetric with the
#' other per-region methods.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size (GWAS N).
#' @param variant_ids Character vector or NULL. Variant identifiers (length p).
#'   Required for SparsePro's output (the \code{.pip} and \code{.cs} files
#'   key on rsids). If NULL, synthetic IDs of the form \code{"snp1", "snp2",
#'   ...} are generated.
#'   Whitespace is normalized to underscores; the resulting IDs must be unique
#'   and must not contain comma, slash or quote delimiters.
#' @param sparsepro_dir Character. Path to the cloned SparsePro repository
#'   (must contain \code{sparsepro_zld.py}).
#' @param python Character. Path to the Python executable. Default:
#'   \code{"python"}.
#' @param K Integer. Maximum number of effect groups (causal signals)
#'   SparsePro will consider. Default: 5.
#' @param cthres Numeric. Coverage level for credible sets. Default: 0.95.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Posterior inclusion probabilities,
#'     re-ordered to match the input \code{variant_ids}.}
#'   \item{credible_sets}{List of integer vectors (1-based). One per effect
#'     group in the \code{.cs} file. Empty list if no credible sets were
#'     reported.}
#'   \item{method}{Character. Always \code{"sparsepro"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List of SparsePro-specific outputs:
#'     \describe{
#'       \item{cs_pip}{List of numeric vectors. Per-variant inclusion
#'         probabilities within each credible set (the \code{pip} column of
#'         the \code{.cs} file), in the same order as \code{credible_sets}.
#'         NULL if not parseable.}
#'       \item{cs_effect_size}{List of numeric vectors. Posterior effect
#'         sizes per variant per credible set (the \code{effect_size} column
#'         of the \code{.cs} file). NULL if not parseable.}
#'       \item{run_log}{Character. Captured stdout/stderr from the SparsePro
#'         invocation; useful for debugging unexpected failures.}
#'       \item{native_output_evidence}{On output-validation failure, bounded
#'         \code{.pip}/\code{.cs} text (or exact bytes when not UTF-8 text),
#'         with file sizes, decoding status and truncation flags. This is
#'         diagnostic evidence, never a substitute prediction.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if SparsePro failed.}
#' }
#'
#' @export
run_sparsepro <- function(z,
                          LD,
                          n,
                          variant_ids   = NULL,
                          sparsepro_dir,
                          python        = "python",
                          K             = 5,
                          cthres        = 0.95) {

  # --- Validate inputs --------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" =
      is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "n must be a positive integer" =
      is.numeric(n) && length(n) == 1 && n > 0,
    "K must be a positive integer" =
      is.numeric(K) && length(K) == 1 && K >= 1,
    "cthres must be in (0, 1)" =
      is.numeric(cthres) && length(cthres) == 1 &&
      cthres > 0 && cthres < 1
  )

  if (is.null(variant_ids)) variant_ids <- paste0("snp", seq_len(p))

  # VCF-derived variant_ids are "1 40023356 . A T" (embedded spaces).
  # SparsePro echoes the id back as the INDEX of its .pip output, and we
  # parse that file with read.table(sep = "") - i.e. split on ANY
  # whitespace - expecting exactly 3 fields (rsid, z, pip). A spaced id
  # yields 7 fields and the parse fails. Collapsing whitespace keeps each
  # id a single token on both the way out and the way back, and keeps
  # match(variant_ids, pip_df$rsid) consistent since both sides are
  # sanitised the same way.
  variant_ids <- gsub("\\s+", "_", variant_ids)

  sparsepro_dir <- path.expand(sparsepro_dir)
  python        <- if (file.exists(path.expand(python))) .fm_python_path(python) else python
  script_path   <- file.path(sparsepro_dir, "sparsepro_zld.py")

  params <- .sparsepro_params(n, K, cthres, sparsepro_dir, python)

  # Delimiters cannot be represented unambiguously in native set cells, and
  # whitespace normalization must remain a bijection over the input IDs.
  if (length(variant_ids) != p || anyNA(variant_ids) || any(!nzchar(variant_ids)) ||
      anyDuplicated(variant_ids) || any(grepl("[,/\"']", variant_ids)))
    return(.sparsepro_error_result(p, params, 0,
      "SparsePro requires unique nonempty normalized IDs without comma, slash or quote delimiters."))

  if (!file.exists(script_path)) {
    return(.sparsepro_error_result(
      p, params, 0,
      paste("sparsepro_zld.py not found at:", script_path,
            "\nRun setup_sparsepro() for instructions.")
    ))
  }

  # --- Set up temp working directory ------------------------------------------

  work_dir <- tempfile(pattern = "sparsepro_run_")
  out_dir  <- file.path(work_dir, "output")
  dir.create(work_dir, recursive = TRUE)
  dir.create(out_dir,  recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  region_label <- "region"
  zfile_name   <- paste0(region_label, ".z")
  ldfile_name  <- paste0(region_label, ".ld")

  # --- Write z-score file -----------------------------------------------------
  # Format: 2-column TSV, no header. Column 1 = variant ID, column 2 = z-score.

  z_path  <- file.path(work_dir, zfile_name)
  z_df    <- data.frame(rsid = variant_ids, z = as.numeric(z),
                        stringsAsFactors = FALSE)
  write.table(z_df, z_path, sep = "\t",
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # --- Write LD matrix file ---------------------------------------------------
  # SparsePro expects whitespace-separated Pearson correlations, no header.

  ld_path <- file.path(work_dir, ldfile_name)
  write.table(round(LD, 8), ld_path,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  # --- Write --zld summary file -----------------------------------------------
  # One row per locus: <zscore_file>\t<ld_file> (plus optional annotation
  # column, which we leave out). File paths are relative to --zdir, so we
  # use the bare filenames here.

  # SparsePro's sparsepro_zld.py reads this via
  #   pd.read_csv(args.zld, sep='\s+')
  # which treats the first row as a HEADER. Without a header line the
  # loop over ldlists iterates zero times, no .pip is written, and the
  # wrapper reports "produced no .pip output". Include the header explicitly.
  zld_path <- file.path(work_dir, "zld.txt")
  writeLines(c("z\tld",
               paste(zfile_name, ldfile_name, sep = "\t")),
             zld_path)

  # --- Build arguments --------------------------------------------------------

  args <- c(
    script_path,
    "--zld",    zld_path,
    "--zdir",   work_dir,
    "--N",      as.character(as.integer(n)),
    "--save",   out_dir,
    "--prefix", region_label,
    "--K",      as.character(as.integer(K)),
    "--cthres", as.character(cthres),
    "--verbose"
  )

  # --- Run SparsePro ----------------------------------------------------------

  start_time <- proc.time()

  run_output <- tryCatch({
    system2(python, args = args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    structure(conditionMessage(e), class = "sparsepro_error")
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  if (inherits(run_output, "sparsepro_error") ||
      (!is.null(attr(run_output, "status")) && attr(run_output, "status") != 0L)) {
    return(.sparsepro_error_result(p, params, elapsed, as.character(run_output),
                                    run_log = character(0)))
  }

  run_log <- if (is.character(run_output)) run_output else character(0)

  # --- Locate output files ----------------------------------------------------
  # SparsePro names .pip/.cs after the Z FILENAME from the --zld list, not
  # after --prefix (sparsepro_zld.py:254):
  #     z.to_csv(os.path.join(args.save, "{}.pip".format(zfile)), ...)
  # so with zfile = "region.z" it writes "region.z.pip" / "region.z.cs".
  # Looking for "<prefix>.pip" silently missed every successful run and
  # reported "produced no .pip output" even though SparsePro had converged
  # and detected effects. (--prefix only names the .h2/.wsep/.W* files.)

  pip_path <- file.path(out_dir, paste0(zfile_name, ".pip"))
  cs_path  <- file.path(out_dir, paste0(zfile_name, ".cs"))
  invalid_output <- function(error) .sparsepro_error_result(p, params, elapsed,
    error, run_log, native_output_evidence = .sparsepro_output_evidence(
      c(pip = pip_path, cs = cs_path)))

  if (!file.exists(pip_path)) {
    err_msg <- paste(
      c("SparsePro produced no .pip output.", run_log),
      collapse = "\n"
    )
    return(.sparsepro_error_result(p, params, elapsed, err_msg, run_log = run_log))
  }

  # R's text readers can silently truncate an embedded-NUL line. Validate
  # complete bytes first, so a valid prefix cannot hide damaged native output.
  for (path in c(pip_path, cs_path)) if (file.exists(path) && !.sparsepro_valid_output_text(path))
    return(invalid_output("SparsePro output contains NUL, invalid UTF-8 or unreadable text."))

  # --- Parse .pip file --------------------------------------------------------
  # Three columns: variant_id, z-score, pip. Whitespace-separated. The order
  # may differ from the input variant order; we re-index against variant_ids.

  pip_df <- tryCatch(
    utils::read.table(pip_path, header = FALSE, sep = "",
                      stringsAsFactors = FALSE,
                      colClasses = c("character", "numeric", "numeric"),
                      na.strings = character(0), quote = "", comment.char = "",
                      col.names = c("rsid", "z", "pip")),
    error = function(e) NULL
  )

  if (is.null(pip_df) || !"pip" %in% names(pip_df)) {
    return(invalid_output("Failed to parse SparsePro .pip output."))
  }

  # Complete predictions must align by unique IDs; never clip invalid output.
  ord <- match(variant_ids, pip_df$rsid)
  if (nrow(pip_df) != p || anyNA(ord) || anyNA(pip_df$rsid) || anyDuplicated(pip_df$rsid))
    return(invalid_output("SparsePro output IDs are incomplete or duplicated."))
  pip <- suppressWarnings(as.numeric(pip_df$pip[ord]))
  if (any(!is.finite(pip)) || any(pip < 0 | pip > 1))
    return(invalid_output("SparsePro output is not a finite probability vector."))

  # --- Parse .cs file ---------------------------------------------------------
  # SparsePro emits one row per effect group (credible set). Columns:
  #   cs            : slash-separated rsids in the pinned upstream version
  #   pip           : slash-separated within-group PIPs
  #   effect_size   : slash-separated effect sizes
  #
  # The exact intra-cell separator varies between SparsePro versions. We
  # also accept the comma/whitespace formats used by other versions. Splitting
  # all three columns the same way preserves each member's probability/effect.

  credible_sets  <- list()
  cs_pip         <- NULL
  cs_effect_size <- NULL

  if (file.exists(cs_path) && file.info(cs_path)$size > 0L) {
    cs_lines <- readLines(cs_path, warn = FALSE)
    cs_lines <- cs_lines[nchar(trimws(cs_lines)) > 0L]

    # A native ID such as "cs" is not a header. Require all three column names.
    header <- tolower(trimws(strsplit(cs_lines[1L], "\t", fixed = TRUE)[[1L]]))
    if (length(header) == 3L && header[1L] %in% c("cs", "set", "group") &&
        header[2L] %in% c("pip", "prob", "probability") &&
        header[3L] %in% c("beta", "effect", "effect_size")) {
      cs_lines <- cs_lines[-1L]
    }

    if (length(cs_lines) > 0L) {
      parsed <- lapply(cs_lines, function(line) {
        # Split top-level on tabs (between the three columns).
        # Append a sentinel before splitting: R otherwise drops a trailing
        # empty cell and could silently accept a missing effect column.
        cells <- head(strsplit(paste0(line, "\tEND"), "\t", fixed = TRUE)[[1L]], -1L)
        if (length(cells) != 3L || any(!nzchar(trimws(cells)))) return(NULL)

        # Robust intra-cell split.
        split_cell <- function(s) {
          s <- trimws(s)
          if (!nzchar(s)) return(character(0))
          comma <- grepl(",", s, fixed = TRUE); slash <- grepl("/", s, fixed = TRUE)
          if (comma && slash) return(character(0))
          if (comma || slash) {
            delimiter <- if (comma) "," else "/"
            tokens <- trimws(head(strsplit(paste0(s, delimiter, "END"),
                                           delimiter, fixed = TRUE)[[1L]], -1L))
          } else tokens <- strsplit(s, "\\s+")[[1L]]
          if (any(!nzchar(tokens))) return(character(0))
          tokens
        }

        rsids  <- split_cell(cells[1L])
        pips   <- if (length(cells) >= 2L) suppressWarnings(as.numeric(split_cell(cells[2L]))) else NA_real_
        effs   <- if (length(cells) >= 3L) suppressWarnings(as.numeric(split_cell(cells[3L]))) else NA_real_

        idx <- match(rsids, variant_ids)
        if (!length(idx) || anyNA(idx) || anyDuplicated(idx)) return(NULL)
        if (length(cells) >= 2L && (length(pips) != length(idx) ||
            any(!is.finite(pips)) || any(pips < 0 | pips > 1))) return(NULL)
        if (length(cells) >= 3L && (length(effs) != length(idx) || any(!is.finite(effs)))) return(NULL)
        position <- order(idx)
        list(idx = idx[position],
             pip = if (length(cells) >= 2L) pips[position] else NULL,
             effect = if (length(cells) >= 3L) effs[position] else NULL)
      })

      if (any(vapply(parsed, is.null, logical(1))))
        return(invalid_output(
          "SparsePro native sets contain invalid IDs or misaligned probabilities/effects."))

      if (length(parsed) > 0L) {
        credible_sets  <- lapply(parsed, `[[`, "idx")
        cs_pip         <- lapply(parsed, `[[`, "pip")
        cs_effect_size <- lapply(parsed, `[[`, "effect")
      }
    }
  }

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "sparsepro",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      cs_pip         = cs_pip,
      cs_effect_size = cs_effect_size,
      run_log        = run_log
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run SparsePro on a single region from simulation data structures
#'
#' Thin adapter that extracts the appropriate inputs from the simulation's
#' \code{region_geno} and \code{region_pheno} objects and calls
#' \code{\link{run_sparsepro}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n}, and optionally \code{variant_ids}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}.
#' @param ... Additional arguments passed to \code{\link{run_sparsepro}}
#'   (e.g. \code{sparsepro_dir}, \code{python}, \code{K}, \code{cthres}).
#'
#' @return The output of \code{\link{run_sparsepro}}.
#' @export
run_sparsepro_region <- function(region_geno, region_pheno, ...) {
  run_sparsepro(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    variant_ids = region_geno$variant_ids,
    ...
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

.sparsepro_params <- function(n, K, cthres, sparsepro_dir, python) {
  list(
    n             = n,
    K             = K,
    cthres        = cthres,
    sparsepro_dir = sparsepro_dir,
    python        = python
  )
}

.sparsepro_valid_output_text <- function(path) tryCatch({
  bytes <- suppressWarnings(readBin(path, "raw", n = file.info(path)$size))
  !any(bytes == as.raw(0L)) && validUTF8(rawToChar(bytes))
}, error = function(e) FALSE)

.sparsepro_output_evidence <- function(paths, max_bytes = 262144L) {
  stats::setNames(lapply(paths, function(path) tryCatch({
    if (!file.exists(path)) return(list(present = FALSE))
    size <- unname(file.info(path)$size)
    bytes <- suppressWarnings(readBin(path, "raw", n = min(size, max_bytes)))
    text <- tryCatch(rawToChar(bytes), error = function(e) NULL)
    if (!is.null(text) && !validUTF8(text)) text <- NULL
    if (!is.null(text)) Encoding(text) <- "UTF-8"
    list(present = TRUE, byte_size = size, captured_bytes = length(bytes),
         truncated = size > max_bytes, text = text,
         decoding_status = if (is.null(text)) "non_text_bytes" else "utf8_text",
         bytes = if (is.null(text)) bytes else NULL)
  }, error = function(e) list(present = file.exists(path),
                              capture_error = conditionMessage(e)))), names(paths))
}

.sparsepro_error_result <- function(p, params, elapsed, error_msg,
                                     run_log = character(0), native_output_evidence = NULL) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "sparsepro",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      cs_pip         = NULL,
      cs_effect_size = NULL,
      run_log        = run_log,
      native_output_evidence = native_output_evidence
    ),
    error           = error_msg
  )
}
