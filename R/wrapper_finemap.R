# =============================================================================
# finemap.R
#
# Wrapper for FINEMAP (Benner et al. 2016) fine-mapping.
#
# FINEMAP is an external C++ binary. This wrapper handles writing the required
# input files, calling the binary, and parsing the output back into the
# standardised result format.
#
# This file provides:
#   - setup_finemap()        : verifies the FINEMAP binary is accessible
#   - run_finemap()          : runs FINEMAP on a single region (explicit inputs)
#   - run_finemap_region()   : adapter called by run_methods(); extracts inputs
#                              from simulation data structures
#
# Standard output format:
#   pip              Numeric vector (length p). Marginal PIPs.
#   credible_sets    List of integer vectors (variant indices). Credible sets
#                    from the most probable number of causal variants k*.
#                    Empty list if k* = 0; invalid/missing native sets fail.
#   method           Character. "finemap".
#   input_type       Character. Always "summary".
#   params           List. Hyperparameters used.
#   runtime_seconds  Numeric. Wall-clock time.
#   additional       List. FINEMAP-specific outputs (see run_finemap() docs).
# =============================================================================


# =============================================================================
# Setup
# =============================================================================

#' Set up the FINEMAP binary
#'
#' Checks that the FINEMAP binary can be found at the supplied path or on the
#' system PATH. If not found and \code{download = TRUE} (the default), the
#' binary is downloaded automatically for your OS and cached locally — no
#' manual installation required.
#'
#' The downloaded binary is saved to your R user cache directory
#' (see \code{tools::R_user_dir("fmbenchmark", "cache")}). The resolved path
#' is returned invisibly, so you can capture it for use in \code{run_methods()}:
#'
#' \preformatted{
#'   fp  <- setup_finemap()
#'   out <- run_methods(sim, methods = "finemap",
#'                      method_args = list(finemap = list(finemap_path = fp)))
#' }
#'
#' @param finemap_path Character. Path to the FINEMAP binary, or \code{"finemap"}
#'   to search the system PATH. Default: \code{"finemap"}.
#' @param download Logical. If \code{TRUE} (default) and the binary is not
#'   found, attempt to download it automatically. Set to \code{FALSE} to
#'   disable downloading and receive manual installation instructions instead.
#'
#' @return Invisibly returns the resolved path to the binary.
#' @export
setup_finemap <- function(finemap_path = "finemap", download = TRUE) {

  # Full path supplied: check directly
  if (file.exists(finemap_path)) {
    message("FINEMAP binary found: ", normalizePath(finemap_path))
    return(invisible(normalizePath(finemap_path)))
  }

  # Search PATH
  resolved <- Sys.which(finemap_path)
  if (nchar(resolved) > 0) {
    message("FINEMAP binary found on PATH: ", resolved)
    return(invisible(resolved))
  }

  # Not found — try to download, or explain how to install manually
  if (download) {
    path <- .download_finemap()
    return(invisible(path))
  }

  sysname <- Sys.info()["sysname"]
  os_hint <- if (.Platform$OS.type == "windows") {
    "  FINEMAP does not provide a Windows binary. Consider using WSL (Linux)."
  } else if (sysname == "Darwin") {
    "  chmod +x finemap_v1.4.1_MacOSX && mv finemap_v1.4.1_MacOSX /usr/local/bin/finemap"
  } else {
    "  chmod +x finemap_v1.4.1_x86_64 && mv finemap_v1.4.1_x86_64 /usr/local/bin/finemap"
  }

  stop(
    "FINEMAP binary not found (searched for: '", finemap_path, "').\n\n",
    "To download automatically, run:\n",
    "  setup_finemap()   # download = TRUE is the default\n\n",
    "To install manually:\n",
    "  1. Go to: http://www.christianbenner.com\n",
    "  2. Download the binary for your OS.\n",
    "  3. ", os_hint, "\n",
    call. = FALSE
  )
}


#' Download the FINEMAP binary for the current OS
#'
#' Downloads a pre-compiled FINEMAP v1.4.1 binary from the echofinemap GitHub
#' releases page, saves it to the R user cache directory, makes it executable,
#' and returns the path.
#'
#' Called automatically by \code{\link{setup_finemap}} when the binary is not
#' found. Can also be called directly to force a fresh download.
#'
#' @return Character. Path to the downloaded binary.
#' @export
download_finemap <- function() {
  .download_finemap()
}


# Internal implementation
.download_finemap <- function() {

  sysname <- Sys.info()["sysname"]

  if (.Platform$OS.type == "windows") {
    stop(
      "Automatic FINEMAP download is not supported on Windows as no Windows\n",
      "binary is available. Consider using WSL (Linux subsystem).",
      call. = FALSE
    )
  }

  fname <- if (sysname == "Darwin") {
    "finemap_v1.4.1_MacOSX.tgz"
  } else {
    "finemap_v1.4.1_x86_64.tgz"
  }

  url <- paste0(
    "https://github.com/RajLabMSSM/echofinemap/releases/download/latest/",
    fname
  )

  # Cache directory
  cache_dir <- tools::R_user_dir("fmbenchmark", which = "cache")
  save_dir  <- file.path(cache_dir, "finemap")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  tgz_path <- file.path(save_dir, fname)

  message("Downloading FINEMAP binary from:\n  ", url)
  tryCatch(
    utils::download.file(url, destfile = tgz_path, mode = "wb", quiet = FALSE),
    error = function(e) {
      stop(
        "Download failed: ", conditionMessage(e), "\n",
        "Check your internet connection or download manually from:\n",
        "  http://www.christianbenner.com",
        call. = FALSE
      )
    }
  )

  # Extract
  utils::untar(tgz_path, exdir = save_dir)

  # Locate the extracted binary (strip .tgz suffix)
  binary_name <- sub("\\.tgz$", "", fname)
  binary_path <- file.path(save_dir, binary_name, binary_name)

  # Fallback: search for any file named like the binary in save_dir
  if (!file.exists(binary_path)) {
    candidates <- list.files(save_dir, pattern = binary_name,
                             full.names = TRUE, recursive = TRUE)
    candidates <- candidates[!grepl("\\.tgz$", candidates)]
    if (length(candidates) == 0) {
      stop(
        "Downloaded and extracted the archive but could not locate the binary.\n",
        "Archive saved at: ", tgz_path,
        call. = FALSE
      )
    }
    binary_path <- candidates[[1]]
  }

  # Make executable
  Sys.chmod(binary_path, "0755")

  if (sysname == "Darwin") {
    arch <- system("uname -m", intern = TRUE)
    if (identical(trimws(arch), "arm64")) {
      message(
        "Note: you are on Apple Silicon (arm64). The downloaded binary is\n",
        "x86_64 and will run via Rosetta 2. Two prerequisites:\n",
        "  1. Rosetta must be installed: softwareupdate --install-rosetta\n",
        "  2. libzstd must be at /usr/local/lib/libzstd.1.dylib (x86_64 build).\n",
        "     If missing, install x86_64 Homebrew and run:\n",
        "       arch -x86_64 /bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\"\n",
        "       arch -x86_64 /usr/local/bin/brew install zstd"
      )
    }
  }

  message("FINEMAP binary ready at:\n  ", binary_path)
  message(
    "\nTo use it, pass the path to run_methods():\n",
    "  run_methods(sim, methods = \"finemap\",\n",
    "              method_args = list(finemap = list(finemap_path = \"",
    binary_path, "\")))"
  )

  binary_path
}


# =============================================================================
# Run FINEMAP on a single region
# =============================================================================

# This is an executable-specific serialization contract, not a version-label
# or manual guarantee. The author download called 1.4.2 has a stale 1.4.1 log
# header; its executable console reports 1.4.2. Neither label selects bounds.
# Inspection of these exact executable bytes verified SNPs::writeLine(prob)
# -> floatToString(double, 6), a fresh default-format stringstream. The active
# floating-point rounding mode was not established, so retain a full adjacent
# decimal last place on each side, covering nearest and directed rounding.
.finemap_probability_writer_contract <- function() list(
  schema = 1L, method = "finemap", source = "native_snp_prob_literal",
  executable_sha256 = "41e0af31c18e2e92d2a5fd1fbf35c9832d62ae737746ecd0aa4629a5925c9bcc",
  writer = "SNPs::writeLine.prob -> floatToString(double,6)",
  significant_digits = 6L, numeric_format = "default_general",
  rounding_mode = "not_verified",
  envelope = "one_adjacent_decimal_last_place_nearest_or_directed",
  null_basis = "structural_zero_excluded_by_prior",
  null_event = "empty_configuration_in_nonnull_conditioned_support")

.finemap_binary_sha256 <- function(path) {
  resolved <- if (file.exists(path)) path else unname(Sys.which(path))
  if (!length(resolved) || !nzchar(resolved) || !file.exists(resolved)) return(NULL)
  .challenge_sha256(resolved)
}

.finemap_probability_token_bounds <- function(tokens, significant_digits = 6L) {
  if (!is.character(tokens) || !length(tokens) || anyNA(tokens) ||
      any(!grepl("^(0|[1-9][0-9]*)([.][0-9]+)?(e[+-]?[0-9]+)?$", tokens)))
    stop("Invalid native probability literal tokens")
  value <- suppressWarnings(as.numeric(tokens))
  if (any(!is.finite(value)) || any(value < 0 | value > 1))
    stop("Invalid native probability literal values")
  mantissa <- sub("e.*$", "", tokens)
  digits <- sub("^0+", "", gsub("[.]", "", mantissa))
  # General formatting suppresses trailing zeros but cannot create more than
  # the declared significant precision. A serialized zero is an exact double
  # zero; nonzero doubles have a nonzero general-format scientific token.
  if (any(nchar(digits) > significant_digits) ||
      any(value == 0 & tokens != "0")) stop("Native probability token exceeds its writer precision")
  lower_step <- upper_step <- rep(0, length(value))
  nonzero <- value > 0
  # Determine the decimal order and exact powers from the literal, avoiding
  # log10 boundary error and binary comparisons to an approximate power.
  for (i in which(nonzero)) {
    exponent <- if (grepl("e", tokens[i], fixed = TRUE))
      as.integer(sub("^.*e", "", tokens[i])) else 0L
    before <- strsplit(mantissa[i], ".", fixed = TRUE)[[1L]][1L]
    raw <- gsub("[.]", "", mantissa[i])
    first <- regexpr("[1-9]", raw)[1L]
    order <- nchar(before) - first + exponent
    is_power <- grepl("^10*$", substring(raw, first))
    upper_step[i] <- 10^(order - significant_digits + 1L)
    lower_step[i] <- 10^(order - significant_digits + if (is_power) 0L else 1L)
    # Decimal spacings may underflow for subnormal tokens. At those scales
    # one binary subnormal quantum is a conservative representable bound.
    quantum <- .Machine$double.xmin * .Machine$double.eps
    lower_step[i] <- max(lower_step[i], quantum)
    upper_step[i] <- max(upper_step[i], quantum)
  }
  list(value = value, lower = pmax(0, value - lower_step),
       upper = pmin(1, value + upper_step))
}

.finemap_precision_record_sha <- function(record) {
  path <- tempfile("finemap_precision_")
  on.exit(unlink(path), add = TRUE)
  saveRDS(record, path, version = 2L)
  .challenge_sha256(path)
}

.finemap_probability_precision <- function(tokens, variant_ids, native_variant_ids,
                                          executable_sha256, k_prior, k_posterior) {
  basis <- .finemap_probability_writer_contract()
  if (!identical(executable_sha256, basis$executable_sha256) ||
      is.null(k_prior) || !identical(unname(k_prior["0"]), 0) ||
      !identical(unname(k_posterior["0"]), 0)) return(NULL)
  bounds <- .finemap_probability_token_bounds(tokens, basis$significant_digits)
  record <- c(basis, list(variant_ids = variant_ids, native_variant_ids = native_variant_ids,
    tokens = tokens, probability_lower = bounds$lower, probability_upper = bounds$upper))
  record$record_sha256 <- .finemap_precision_record_sha(record)
  record
}

# Validate all provenance before applying any envelope. Metadata is an audited
# adapter declaration, sealed by campaign artifacts; it is not authentication
# of arbitrary public submissions. Bounds are never substituted for PIPs.
.finemap_validate_probability_precision <- function(pred) {
  record <- pred[["native_probability_precision"]]
  if (is.null(record)) return(NULL)
  basis <- .finemap_probability_writer_contract()
  fields <- c(names(basis), "variant_ids", "native_variant_ids", "tokens",
              "probability_lower", "probability_upper", "record_sha256")
  if (!is.list(record) || !identical(names(record), fields) ||
      !identical(pred$implementation, "finemap") ||
      !all(vapply(names(basis), function(k) identical(record[[k]], basis[[k]]), TRUE)))
    stop("Invalid native probability precision basis")
  declared_sha <- record$record_sha256
  plain <- record; plain$record_sha256 <- NULL
  if (!is.character(declared_sha) || length(declared_sha) != 1L ||
      is.na(declared_sha) || !identical(declared_sha, .finemap_precision_record_sha(plain)))
    stop("Native probability precision record changed")
  ids <- record$variant_ids
  if (!is.character(ids) || length(ids) != length(pred$pip) || anyNA(ids) ||
      anyDuplicated(ids) || !setequal(ids, pred$variant_id) ||
      !identical(record$native_variant_ids, paste0("variant_", seq_along(ids))))
    stop("Native probability precision variant IDs do not match")
  diagnostic <- pred[["native_diagnostics"]]
  if (!is.list(diagnostic) ||
      !identical(diagnostic[["variant_ids"]], ids) ||
      !identical(diagnostic[["native_variant_ids"]], record$native_variant_ids) ||
      !identical(diagnostic[["native_executable_sha256"]], basis$executable_sha256) ||
      !identical(diagnostic[["snp_probability_tokens"]], record$tokens) ||
      !identical(diagnostic[["native_probability_precision"]], record) ||
      !identical(diagnostic[["native_null_support"]], "excluded_by_prior") ||
      !identical(diagnostic[["posterior_conditioning"]], "at_least_one_modeled_cause") ||
      !identical(diagnostic[["null_probability_event"]], basis$null_event) ||
      !identical(diagnostic[["null_probability"]], 0) ||
      !identical(unname(diagnostic[["k_prior"]]["0"]), 0) ||
      !identical(unname(diagnostic[["k_posterior"]]["0"]), 0) ||
      !identical(pred[["native_null_support"]], "excluded_by_prior") ||
      !identical(pred[["posterior_conditioning"]], "at_least_one_modeled_cause") ||
      !identical(pred[["null_probability_event"]], basis$null_event) ||
      !identical(pred[["probability_target_aligned"]], FALSE) ||
      !identical(pred[["null_model_in_support"]], FALSE) ||
      !identical(pred[["null_probability"]], 0))
    stop("Native probability precision diagnostics or null basis do not match")
  bounds <- .finemap_probability_token_bounds(record$tokens, basis$significant_digits)
  ordering <- match(pred$variant_id, ids)
  if (!identical(unname(pred$pip), unname(bounds$value[ordering])) ||
      !identical(record$probability_lower, bounds$lower) ||
      !identical(record$probability_upper, bounds$upper))
    stop("Native probability tokens, PIPs, or bounds do not match")
  list(lower = bounds$lower[ordering], upper = bounds$upper[ordering])
}

.finemap_model_size_probabilities <- function(path, max_causal = NULL,
    header_label = "Post-Pr(# of causal SNPs is k)") {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path, warn = FALSE)
  header <- grep(header_label, lines, fixed = TRUE)
  if (!length(header) || tail(header, 1L) >= length(lines)) return(NULL)
  tail_lines <- lines[seq.int(tail(header, 1L) + 1L, length(lines))]
  pattern <- "^[[:space:]]*[(]?[0-9]+[[:space:]]*->[[:space:]]*[0-9.eE+-]+[)]?[[:space:]]*$"
  values <- character()
  for (line in tail_lines) {
    if (!nzchar(trimws(line)) && !length(values)) next
    if (!grepl(pattern, line)) break
    values <- c(values, line)
  }
  if (!length(values)) return(NULL)
  values <- gsub("[()]", "", values)
  pieces <- strsplit(trimws(values), "[[:space:]]*->[[:space:]]*")
  sizes <- vapply(pieces, function(x) x[1L], "")
  prob <- suppressWarnings(as.numeric(vapply(pieces, function(x) x[2L], "")))
  if (anyDuplicated(sizes) || any(!is.finite(prob)) || any(prob < 0 | prob > 1) ||
      !identical(as.integer(sizes), seq.int(0L, if (is.null(max_causal)) max(as.integer(sizes)) else max_causal)))
    return(NULL)
  # The tested 1.4.2 release reports log probabilities at three
  # significant digits, suppressing trailing zeros: .548+.404+.0476=.9996.
  # Validate the rounding interval, never replace the complete posterior by
  # sums of a truncated .config table or silently renormalize reported values.
  tokens <- vapply(pieces, function(x) x[2L], "")
  significant_digits <- vapply(tokens, function(token) {
    mantissa <- sub("[eE].*$", "", token)
    digits <- sub("^0+", "", gsub("[^0-9]", "", mantissa))
    max(3L, nchar(digits))
  }, 1L)
  half_unit <- ifelse(prob == 0, 0, .5 * 10^(floor(log10(pmax(prob, .Machine$double.xmin))) - significant_digits + 1L))
  # Under the existing half-unit log policy, an exact power of ten has ten
  # times finer preceding spacing. A symmetric lower envelope at 1 would
  # wrongly admit 1 + .0006 under that policy.
  power <- vapply(tokens, function(token) {
    mantissa <- sub("[eE].*$", "", token)
    digits <- sub("^0+", "", gsub("[^0-9]", "", mantissa))
    grepl("^10*$", digits)
  }, TRUE)
  lower <- pmax(0, prob - half_unit * ifelse(power, .1, 1)); upper <- pmin(1, prob + half_unit)
  if (sum(lower) > 1 + 1e-12 || sum(upper) < 1 - 1e-12) return(NULL)
  result <- stats::setNames(prob, sizes)
  attr(result, "reported_rounding") <- list(
    minimum_significant_digits = 3L, reported_tokens = stats::setNames(tokens, sizes),
    reported_sum = sum(prob),
    sum_interval = c(lower = sum(lower), upper = sum(upper)),
    probability_lower = stats::setNames(lower, sizes),
    probability_upper = stats::setNames(upper, sizes), renormalized = FALSE)
  result
}

.finemap_component_sets <- function(path, ids, k) {
  if (!file.exists(path)) return(NULL)
  table <- tryCatch(utils::read.table(path, header = TRUE, stringsAsFactors = FALSE,
    fill = TRUE, na.strings = c("", "NA"), comment.char = "#"), error = function(e) NULL)
  columns <- paste0("cred", seq_len(k))
  if (is.null(table) || !all(columns %in% names(table))) return(NULL)
  sets <- probabilities <- vector("list", k)
  for (j in seq_len(k)) {
    probability <- paste0("prob", j)
    if (!probability %in% names(table)) return(NULL)
    raw <- as.character(table[[columns[j]]])
    keep <- !is.na(raw) & nzchar(trimws(raw))
    raw <- raw[keep]
    idx <- match(raw, ids)
    prob <- suppressWarnings(as.numeric(table[[probability]][keep]))
    if (!length(idx) || anyNA(idx) || anyDuplicated(idx) ||
        any(!is.finite(prob)) || any(prob < 0 | prob > 1)) return(NULL)
    ordering <- order(idx)
    sets[[j]] <- as.integer(idx[ordering])
    probabilities[[j]] <- prob[ordering]
  }
  list(sets = sets, probabilities = probabilities)
}

#' Run FINEMAP fine-mapping on a single region
#'
#' Writes the required FINEMAP input files to a temporary directory, calls the
#' FINEMAP binary, parses the output, and returns results in the standardised
#' format. The temporary directory is deleted on exit.
#'
#' @param z Numeric vector. Marginal z-scores (length p).
#' @param LD Matrix. LD (correlation) matrix (p x p).
#' @param n Integer. Sample size.
#' @param beta Numeric vector or NULL. Marginal effect size estimates (length
#'   p). If NULL, derived as \code{z * se}. When both \code{beta} and
#'   \code{se} are NULL, \code{se} is set to 1 and \code{beta = z}.
#' @param se Numeric vector or NULL. Standard errors (length p). See above.
#' @param maf Numeric vector or NULL. Minor allele frequencies (length p).
#'   Used in the FINEMAP input file. If NULL, set to 0.3 for all variants.
#' @param variant_ids Character vector or NULL. Variant identifiers (length p).
#'   Used as rsids in FINEMAP input/output. If NULL, set to
#'   \code{"SNP_1", "SNP_2", ...}.
#' @param finemap_path Character. Path to the FINEMAP binary or name on PATH.
#'   Default: \code{"finemap"}.
#' @param n_causal Integer. Maximum number of causal variants considered.
#'   Default: 5.
#' @param n_iter Integer. Number of stochastic shotgun search iterations.
#'   Default: 100000.
#' @param prior_std Numeric. Prior standard deviation for effect sizes.
#'   Default: 0.05.
#' @param coverage Numeric. Coverage level for credible sets. Default: 0.95.
#'
#' @return A list with the standardised fine-mapping output:
#' \describe{
#'   \item{pip}{Numeric vector (length p). Marginal posterior inclusion
#'     probabilities.}
#'   \item{credible_sets}{List of integer vectors. Credible sets derived from
#'     the most probable number of causal variants k*. Each element is the
#'     set of variant indices for one causal signal. Empty list if k* = 0.}
#'   \item{method}{Character. Always \code{"finemap"}.}
#'   \item{input_type}{Character. Always \code{"summary"}.}
#'   \item{params}{List. Hyperparameters used.}
#'   \item{runtime_seconds}{Numeric. Wall-clock time in seconds.}
#'   \item{additional}{List of FINEMAP-specific outputs:
#'     \describe{
#'       \item{log10bf}{Numeric vector (length p). Log10 Bayes factor per
#'         variant.}
#'       \item{posterior_mean}{Numeric vector (length p). Posterior mean
#'         effect size.}
#'       \item{posterior_sd}{Numeric vector (length p). Posterior SD of
#'         effect size.}
#'       \item{k_posterior}{Named numeric vector. Posterior probability that
#'         exactly k variants are causal (names are k = 0, 1, 2, ...).}
#'       \item{best_k}{Integer. Most probable number of causal variants.}
#'       \item{configs}{Data frame. Top causal configurations with their
#'         posterior probabilities, as returned in FINEMAP's .config file.}
#'       \item{snp_probability_tokens}{Literal native probability strings,
#'         aligned to input variants. Numerical PIPs are not renormalized.}
#'       \item{native_probability_precision}{Serialization envelopes for the
#'         exact audited executable only, with six significant digits and an
#'         unverified rounding mode. Used solely to validate probability-event
#'         feasibility under the structural non-null prior; never as scores.}
#'       \item{null_probability_event}{Identifies the candidate model's empty
#'         configuration event. Its structural zero under the non-null prior
#'         is not a calibrated probability of the generating whole-locus null.}
#'     }
#'   }
#'   \item{error}{Character or NULL. Error message if FINEMAP failed.}
#' }
#'
#' @export
run_finemap <- function(z,
                        LD,
                        n,
                        beta         = NULL,
                        se           = NULL,
                        maf          = NULL,
                        variant_ids  = NULL,
                        finemap_path = "finemap",
                        n_causal     = 5,
                        n_iter       = 100000,
                        prior_std    = 0.05,
                        coverage     = 0.95) {

  # Expand ~ so system2 can exec the binary. Users often pass "~/tools/..."
  # and R's default tilde-expansion doesn't apply in system2's first arg.
  finemap_path <- path.expand(finemap_path)

  # --- Validate inputs --------------------------------------------------------

  p <- length(z)

  stopifnot(
    "LD must be a p x p matrix" = is.matrix(LD) && nrow(LD) == p && ncol(LD) == p,
    "n must be a positive integer" = is.numeric(n) && length(n) == 1 && n > 0
  )

  if (is.null(variant_ids)) variant_ids <- paste0("SNP_", seq_len(p))

  if (!is.character(variant_ids) || length(variant_ids) != p || anyNA(variant_ids) ||
      anyDuplicated(variant_ids) || any(!nzchar(variant_ids))) stop("Invalid variant IDs")
  original_variant_ids <- variant_ids
  # A bijective token map avoids whitespace parsing and normalization collisions.
  variant_ids <- paste0("variant_", seq_len(p))
  if (is.null(maf))         maf         <- rep(0.3, p)

  # Derive beta/se if not supplied
  if (is.null(se) && is.null(beta)) {
    se   <- rep(1.0, p)
    beta <- z
  } else if (is.null(beta)) {
    beta <- z * se
  } else if (is.null(se)) {
    # se = beta / z; guard against z = 0
    se <- ifelse(abs(z) > 1e-10, abs(beta / z), 1.0)
  }

  params <- list(
    n_causal     = n_causal,
    n_iter       = n_iter,
    prior_std    = prior_std,
    coverage     = coverage,
    finemap_path = finemap_path
  )

  # --- Set up temp working directory ------------------------------------------

  work_dir <- tempfile(pattern = "finemap_run_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  prefix <- file.path(work_dir, "fm")

  # --- Write .z file ----------------------------------------------------------

  z_file <- paste0(prefix, ".z")
  z_df <- data.frame(
    rsid       = variant_ids,
    chromosome = 1L,
    position   = seq_len(p),
    allele1    = "A",
    allele2    = "T",
    maf        = round(maf, 6),
    beta       = beta,
    se         = se,
    stringsAsFactors = FALSE
  )
  write.table(z_df, z_file,
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = " ")

  # --- Write .ld file ---------------------------------------------------------

  ld_file <- paste0(prefix, ".ld")
  write.table(round(LD, 6), ld_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  # --- Write .master file -----------------------------------------------------

  snp_file    <- paste0(prefix, ".snp")
  config_file <- paste0(prefix, ".config")
  cred_base   <- paste0(prefix, ".cred")   # FINEMAP appends the k (e.g. .cred1)
  # File suffix must equal the master-file header column name ("log"),
  # else FINEMAP v1.4.2 errors with 'Extension X of file Y ... cannot be
  # found in the header of the master file'. Previously ".log_sss" was
  # used and no longer accepted.
  log_file    <- paste0(prefix, ".log")

  master_file <- file.path(work_dir, "master")
  writeLines(
    c(
      "z;ld;snp;config;cred;log;n_samples",
      paste(z_file, ld_file, snp_file, config_file, cred_base, log_file, n,
            sep = ";")
    ),
    master_file
  )

  # --- Run FINEMAP ------------------------------------------------------------

  # FINEMAP v1.4.2 CLI: --n-iterations was renamed to --n-iter, and
  # --credible-config-value was replaced by --prob-cred-set. --log directs
  # FINEMAP to write to the log path specified in the master's 'log' column.
  args <- c(
    "--sss",
    "--in-files",       master_file,
    "--n-causal-snps",  as.character(n_causal),
    "--n-iter",         as.character(n_iter),
    "--prior-std",      as.character(prior_std),
    "--prob-cred-set",  as.character(coverage),
    "--log"
  )

  executable_sha256 <- .finemap_binary_sha256(finemap_path)
  # Custody checks belong to the task/PBS cost; the method's existing native
  # elapsed-time contract measures system2 execution without binary hashing.
  start_time <- proc.time()

  run_output <- tryCatch({
    system2(finemap_path, args = args, stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    structure(conditionMessage(e), class = "finemap_error")
  })

  elapsed <- as.numeric((proc.time() - start_time)["elapsed"])

  # Check for binary-level failure
  if (inherits(run_output, "finemap_error") ||
      (!is.null(attr(run_output, "status")) && attr(run_output, "status") != 0L)) {
    return(.finemap_error_result(p, params, elapsed, run_output))
  }

  # Check that the primary output file was produced
  if (!file.exists(snp_file)) {
    err_msg <- paste(
      c("FINEMAP produced no output.", run_output),
      collapse = "\n"
    )
    return(.finemap_error_result(p, params, elapsed, err_msg))
  }

  # --- Parse .snp file --------------------------------------------------------

  snp_df <- tryCatch(
    utils::read.table(snp_file, header = TRUE, stringsAsFactors = FALSE,
                      colClasses = c(prob = "character")),
    error = function(e) NULL
  )

  if (is.null(snp_df) || !all(c("rsid", "prob") %in% names(snp_df))) {
    return(.finemap_error_result(p, params, elapsed,
                                 "Failed to parse FINEMAP .snp output file."))
  }

  # Re-order rows to match the original variant order
  ord <- match(variant_ids, snp_df$rsid)
  if (nrow(snp_df) != p || anyNA(ord) || anyNA(snp_df$rsid) || anyDuplicated(snp_df$rsid))
    return(.finemap_error_result(p, params, elapsed, "FINEMAP output IDs do not match the complete input."))

  probability_tokens <- snp_df$prob[ord]
  pip            <- suppressWarnings(as.numeric(probability_tokens))
  if (any(!is.finite(pip)) || any(pip < 0 | pip > 1))
    return(.finemap_error_result(p, params, elapsed, "FINEMAP must return finite probabilities in [0,1]."))
  posterior_mean <- snp_df$mean[ord]
  posterior_sd   <- snp_df$sd[ord]
  log10bf        <- snp_df$log10bf[ord]

  # --- Parse .config file to get posterior over k -----------------------------

  config_df  <- tryCatch(
    utils::read.table(config_file, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  # .config contains only the top configurations; summing it need not recover
  # model-size probabilities. --log writes the complete native posterior.
  log_candidates <- c(paste0(log_file, "_sss"), log_file)
  actual_log <- log_candidates[file.exists(log_candidates)][1L]
  k_posterior <- if (!is.na(actual_log)) .finemap_model_size_probabilities(actual_log, n_causal) else NULL
  if (is.null(k_posterior)) return(.finemap_error_result(p, params, elapsed,
    "FINEMAP produced no valid complete model-size posterior in its log."))
  best_k <- as.integer(names(which.max(k_posterior)))
  rounding <- attr(k_posterior, "reported_rounding")
  k_posterior <- stats::setNames(as.numeric(k_posterior), names(k_posterior))
  prior <- .finemap_model_size_probabilities(actual_log, n_causal,
    header_label = "Prior-Pr(# of causal SNPs is k)")
  k_prior <- if (is.null(prior)) NULL else stats::setNames(as.numeric(prior), names(prior))
  if (!identical(executable_sha256, .finemap_binary_sha256(finemap_path)))
    return(.finemap_error_result(p, params, elapsed, "FINEMAP executable changed during execution."))
  precision <- tryCatch(.finemap_probability_precision(probability_tokens,
    original_variant_ids, variant_ids, executable_sha256, k_prior, k_posterior),
    error = identity)
  if (inherits(precision, "error")) return(.finemap_error_result(p, params, elapsed,
    paste("FINEMAP probability serialization contract failed:", conditionMessage(precision))))

  # --- Parse .cred<k*> file for credible sets ---------------------------------

  credible_sets <- list()
  cs_pip <- NULL

  if (best_k > 0L) {
    cred_k_file <- paste0(cred_base, best_k)

    native_sets <- .finemap_component_sets(cred_k_file, variant_ids, best_k)
    if (is.null(native_sets)) return(.finemap_error_result(p, params, elapsed,
      "FINEMAP selected-model credible sets are missing or contain invalid IDs/probabilities."))
    credible_sets <- native_sets$sets
    cs_pip <- native_sets$probabilities
  }

  # --- Return -----------------------------------------------------------------

  list(
    pip             = pip,
    credible_sets   = credible_sets,
    method          = "finemap",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      log10bf        = log10bf,
      posterior_mean = posterior_mean,
      posterior_sd   = posterior_sd,
      k_posterior    = k_posterior,
      k_posterior_rounding = rounding,
      k_prior = k_prior,
      native_executable_sha256 = executable_sha256,
      snp_probability_tokens = probability_tokens,
      native_probability_precision = precision,
      best_k         = best_k,
      configs        = config_df,
      cs_pip         = cs_pip,
      credible_set_semantics = "native_configuration_conditional_component_sets",
      k_posterior_source = "complete_native_log",
      best_k_selection = "maximum_reported_native_posterior",
      native_null_support = if (is.null(k_prior)) "not_established_from_log" else
        if (k_prior["0"] == 0) "excluded_by_prior" else "prior_includes_null",
      posterior_conditioning = if (is.null(k_prior)) "not_established_from_log" else
        "at_least_one_modeled_cause",
      null_probability = unname(k_posterior["0"]),
      null_probability_event = if (!is.null(k_prior) && k_prior["0"] == 0)
        "empty_configuration_in_nonnull_conditioned_support" else
          "separate_native_count_posterior_null_not_verified_jointly_with_snp_pips",
      null_probability_interpretation = if (!is.null(k_prior) && k_prior["0"] == 0)
        "structural_zero; operational_against_generating_null_not_event_calibration" else
          "cause_count_posterior_with_separate_null_prior; joint_pip_coherence_not_established",
      native_seed = NULL,
      native_rng_reproducibility = "not_established",
      variant_ids = original_variant_ids,
      native_variant_ids = variant_ids
    )
  )
}


# =============================================================================
# Region adapter (called by run_methods)
# =============================================================================

#' Run FINEMAP on a single region from simulation data structures
#'
#' Thin adapter that extracts the appropriate inputs from the simulation's
#' \code{region_geno} and \code{region_pheno} objects and calls
#' \code{\link{run_finemap}}. This is the function registered in the method
#' registry and called by \code{\link{run_methods}}.
#'
#' @param region_geno List. One element of \code{simulation$genotypes},
#'   containing \code{LD}, \code{n}, and optionally \code{maf} and
#'   \code{variant_ids}.
#' @param region_pheno List. One element of a scenario's \code{regions},
#'   containing \code{z}, \code{beta_hat}, and \code{se}.
#' @param ... Additional arguments passed to \code{\link{run_finemap}}
#'   (e.g. \code{n_causal}, \code{prior_std}, \code{finemap_path}).
#'
#' @return The output of \code{\link{run_finemap}}.
#' @export
run_finemap_region <- function(region_geno, region_pheno, ...) {
  run_finemap(
    z           = region_pheno$z,
    LD          = region_geno$LD,
    n           = region_geno$n,
    beta        = region_pheno$beta_hat,
    se          = region_pheno$se,
    maf         = region_geno$maf,
    variant_ids = region_geno$variant_ids,
    ...
  )
}


# =============================================================================
# Internal helpers
# =============================================================================

# Construct a standardised error result when FINEMAP fails, so the batch
# runner can continue rather than abort.
.finemap_error_result <- function(p, params, elapsed, error_msg) {
  list(
    pip             = rep(NA_real_, p),
    credible_sets   = list(),
    method          = "finemap",
    input_type      = "summary",
    params          = params,
    runtime_seconds = elapsed,
    additional      = list(
      log10bf        = rep(NA_real_, p),
      posterior_mean = rep(NA_real_, p),
      posterior_sd   = rep(NA_real_, p),
      k_posterior    = NULL,
      best_k         = NA_integer_,
      configs        = NULL
    ),
    error = error_msg
  )
}
