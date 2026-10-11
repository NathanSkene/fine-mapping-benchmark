# =============================================================================
# wrapper_fb_joint.R
#
# CROSS-REGION JOINT (Functional) BEATRICE  -  Iteration 003 model changes.
#
# Two methods that train ONE shared annotation prior head JOINTLY across all
# regions of a scenario (not the sequential warm-start chain, not per-locus):
#
#   fb_pooled   -> prior_head = 'linear'    (user idea #1, flagship):
#                  a shared logistic annotation->prior map, funmap/PolyFun-style.
#   fb_xregion  -> prior_head = 'lassonet'  (user idea #2):
#                  the LassoNet prior, shared and jointly trained across regions.
#
# Objective  L(phi,{psi_r}) = sum_r ELBO_r( psi_r ; p0_r = f_phi(v_r) ), with the
# shared head phi updated ONCE per step from the equal-weighted mean of the
# per-region gradients (EQUAL WEIGHT PER REGION; see docs/autoresearch/iteration-003.md).
#
# Mechanism. run_methods() calls a method's `run_<method>_scenario_setup(genotypes,
# regions, user_args)` ONCE per scenario (all regions visible) and merges the
# returned list into every region's args. We use that hook to run the joint Python
# trainer (BEATRICE_annot_sparse/beatrice_joint.py) over ALL regions in one process,
# parse each region's pip.csv/credible_set.txt, and return a cache keyed by a
# complete ordered-input digest. The per-region function retrieves its result.
#
# The annotation-free arm explicitly runs the locus-only BEATRICE ablation.
# Failed annotated joint fits remain failures; they never change estimator.
# =============================================================================


# --- complete ordered observed input identity, shared with PAINTOR -----------
.fb_fingerprint <- function(z, LD = NULL, annotations = NULL, n = NULL,
                             variant_ids = NULL) {
  .fm_pooled_input_key(z, LD, annotations, n, variant_ids)
}


# --- write one region's .z/.ld/.annot into work_dir (mirrors FB wrapper I/O) --
.fb_joint_write_region <- function(work_dir, idx, z, LD, annotations, variant_ids) {
  p <- length(z)
  identity <- .beatrice_variant_identity(variant_ids, p)
  variant_ids <- identity$native_ids

  z_path     <- file.path(work_dir, sprintf("region%d.z", idx))
  ld_path    <- file.path(work_dir, sprintf("region%d.ld", idx))
  annot_path <- file.path(work_dir, sprintf("region%d.annot", idx))
  target_dir <- file.path(work_dir, sprintf("out%d", idx))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

  writeLines(paste(variant_ids, z, sep = " "), z_path)
  utils::write.table(round(LD, 8), ld_path,
                     quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")
  annot_df <- cbind(variant_ids, as.data.frame(annotations))
  utils::write.table(annot_df, annot_path,
                     quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

  list(z = z_path, LD = ld_path, annot = annot_path, target = target_dir,
       variant_ids = variant_ids, original_ids = identity$original_ids,
       variant_alignment = identity$alignment, p = p)
}


# --- parse one region's joint output (same contract as the FB wrapper) --------
.fb_joint_requires_cs_diagnostics <- function(beatrice_dir) {
  marker <- file.path(beatrice_dir, "scripts", "credible_set_contract.json")
  if (!file.exists(marker)) return(FALSE)
  value <- jsonlite::fromJSON(marker, simplifyVector = FALSE)
  if (!identical(value$format, "beatrice-conditional-credible-set-diagnostics-v1") ||
      !is.numeric(value$index_base) || length(value$index_base) != 1L ||
      !is.finite(value$index_base) || value$index_base != 0 ||
      !identical(value$default_stopping_rule, "conditional_union"))
    stop("Invalid declared BEATRICE credible-set contract")
  TRUE
}

.fb_joint_read_cs_diagnostics <- function(path, variant_ids, native_sets, pip) {
  value <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  probability <- function(x) scalar(x) && x >= 0 && x <= 1
  indices <- function(x, empty = FALSE) {
    if (!is.list(x) || !is.null(names(x)) || (!empty && !length(x)) ||
        !all(vapply(x, scalar, TRUE))) stop("Invalid BEATRICE diagnostic indices")
    if (!length(x)) return(integer())
    result <- unlist(x, use.names = FALSE)
    if (any(result != floor(result) | result < 0 | result >= length(variant_ids)) ||
        anyDuplicated(result)) stop("Invalid BEATRICE diagnostic indices")
    as.integer(result)
  }
  if (!is.list(value) ||
      !identical(value$format, "beatrice-conditional-credible-set-diagnostics-v1") ||
      !scalar(value$index_base) || value$index_base != 0 ||
      !scalar(value$variant_count) || value$variant_count != length(variant_ids) ||
      !is.list(value$variant_names) || !is.null(names(value$variant_names)) ||
      !all(vapply(value$variant_names, function(x) is.character(x) && length(x) == 1L, TRUE)) ||
      !identical(unlist(value$variant_names, use.names = FALSE), variant_ids) ||
      !identical(value$posterior_memo_format, "beatrice-log-bf-final-prior-v1") ||
      !scalar(value$posterior_support_size) || value$posterior_support_size < 1 ||
      value$posterior_support_size != floor(value$posterior_support_size) ||
      !identical(value$null_probability_conditioning, "explored_support_under_frozen_final_prior") ||
      !probability(value$null_probability) || any(pip > 1 - value$null_probability + 1e-8) ||
      1 - value$null_probability > sum(pip) + 1e-8 ||
      !is.character(value$stopping_rule) || length(value$stopping_rule) != 1L ||
      !value$stopping_rule %in% c("conditional_union", "legacy_membership_sum") ||
      !is.list(value$sets) || !is.null(names(value$sets)) ||
      length(value$sets) != length(native_sets)) stop("Invalid BEATRICE credible-set diagnostics")
  for (i in seq_along(value$sets)) {
    s <- value$sets[[i]]
    set <- indices(s$set_indices_zero_based)
    given <- indices(s$condition_indices_zero_based, empty = TRUE)
    if (!identical(set, as.integer(native_sets[[i]])) || length(intersect(set, given)) ||
        !identical(s$event, "at_least_one_set_member_is_causal") ||
        !identical(s$posterior_conditioning, "given_variants_are_causal_within_explored_support") ||
        !probability(s$conditional_union_probability) ||
        !probability(s$target_probability) || s$target_probability <= 0 ||
        !is.logical(s$target_met) || length(s$target_met) != 1L || is.na(s$target_met) ||
        !identical(s$target_met, s$conditional_union_probability >= s$target_probability) ||
        !scalar(s$legacy_membership_sum) || s$legacy_membership_sum < 0 ||
        s$legacy_membership_sum > length(set) + 1e-8 ||
        s$conditional_union_probability + 1e-8 < s$legacy_membership_sum / length(set) ||
        s$conditional_union_probability > s$legacy_membership_sum + 1e-8 ||
        !probability(s$selection_probability_threshold) ||
        !probability(s$purity_threshold)) stop("Invalid BEATRICE credible-set event evidence")
    if (!is.null(s$stopping_rule) && !identical(s$stopping_rule, value$stopping_rule))
      stop("Inconsistent BEATRICE credible-set stopping rule")
    if (!length(given) &&
        (s$conditional_union_probability + 1e-8 < max(pip[set + 1L]) ||
         s$conditional_union_probability > 1 - value$null_probability + 1e-8 ||
         abs(s$legacy_membership_sum - sum(pip[set + 1L])) > 1e-8))
      stop("BEATRICE unconditional event disagrees with marginal probabilities")
    if (length(given) && any(pip[given + 1L] <= 0))
      stop("BEATRICE diagnostic conditions on an impossible variant event")
    value$sets[[i]]$set_indices <- sort(set + 1L)
    value$sets[[i]]$condition_indices <- given + 1L
  }
  value
}

.beatrice_posterior_metadata <- function(parsed) {
  d <- parsed$cs_diagnostics
  list(credible_set_semantics = if (is.null(d)) "native_semantics_not_verified" else
      if (identical(d$stopping_rule, "conditional_union")) "conditional_union_within_explored_support" else
        "legacy_conditional_membership_mass_candidates",
    credible_set_diagnostics = d,
    posterior_conditioning = if (is.null(d)) "not_declared" else "within_explored_configuration_support",
    native_null_support = if (is.null(d)) "not_declared" else "included_in_explored_support",
    null_model_in_support = if (is.null(d)) NA else TRUE,
    null_probability_not_provided = is.null(d),
    configuration_count = if (is.null(d)) NA_integer_ else d$posterior_support_size,
    null_probability_conditioning = if (is.null(d)) "not_declared" else d$null_probability_conditioning,
    null_probability = if (is.null(d)) NULL else d$null_probability)
}

.fb_joint_parse_output <- function(out_dir, variant_ids, p, require_diagnostics = FALSE) {
  pip_path  <- file.path(out_dir, "pip.csv")
  cred_path <- file.path(out_dir, "credible_set.txt")
  cond_path <- file.path(out_dir, "conditional_credible_variants_probability.txt")
  if (!file.exists(pip_path)) return(NULL)
  if (require_diagnostics && (!file.exists(cred_path) || !file.exists(cond_path))) return(NULL)

  pip_df <- tryCatch(utils::read.csv(pip_path, stringsAsFactors = FALSE),
                     error = function(e) NULL)
  if (is.null(pip_df) || !all(c("pip", "variant_names") %in% names(pip_df))) return(NULL)
  if (anyNA(pip_df$variant_names) || anyDuplicated(pip_df$variant_names)) return(NULL)
  ord <- match(variant_ids, pip_df$variant_names)
  if (anyNA(ord) || nrow(pip_df) != p) return(NULL)
  if ("variant_index" %in% names(pip_df) &&
      !identical(as.numeric(pip_df$variant_index[ord]), as.numeric(seq_len(p) - 1L))) return(NULL)
  pip <- suppressWarnings(as.numeric(pip_df$pip[ord]))
  if (length(pip) != p || any(!is.finite(pip)) || any(pip < 0 | pip > 1)) return(NULL)

  credible_sets <- list()
  native_sets <- list()
  set_orders <- list()
  if (file.exists(cred_path)) {
    cred_lines <- readLines(cred_path, warn = FALSE)
    cred_lines <- cred_lines[nchar(trimws(cred_lines)) > 0]
    if (length(cred_lines) > 0) {
      for (line in cred_lines) {
        tokens <- strsplit(trimws(line), "\\s+")[[1L]]
        if (any(!grepl("^[0-9]+$", tokens))) return(NULL)
        native_idx <- suppressWarnings(as.numeric(tokens))
        if (any(!is.finite(native_idx)) || any(native_idx < 0 | native_idx >= p) ||
            anyDuplicated(native_idx)) return(NULL)
        idx <- as.integer(native_idx) + 1L  # native zero-based indices
        native_sets[[length(native_sets) + 1L]] <- as.integer(native_idx)
        set_orders[[length(set_orders) + 1L]] <- order(idx)
        credible_sets[[length(credible_sets) + 1L]] <- sort(idx)
      }
    }
  }

  cs_pip <- NULL
  if (file.exists(cond_path)) {
    cond_lines <- readLines(cond_path, warn = FALSE)
    cond_lines <- cond_lines[nchar(trimws(cond_lines)) > 0]
    if (length(cond_lines) != length(credible_sets)) return(NULL)
    cs_pip <- vector("list", length(cond_lines))
    for (i in seq_along(cond_lines)) {
      values <- suppressWarnings(as.numeric(strsplit(trimws(cond_lines[i]), "\\s+")[[1L]]))
      if (length(values) != length(credible_sets[[i]]) || any(!is.finite(values)) ||
          any(values < 0 | values > 1)) return(NULL)
      cs_pip[[i]] <- values[set_orders[[i]]]
    }
  }

  diagnostic_path <- file.path(out_dir, "credible_set_diagnostics.json")
  diagnostics <- NULL
  if (file.exists(diagnostic_path)) {
    diagnostics <- tryCatch(.fb_joint_read_cs_diagnostics(diagnostic_path,
      variant_ids, native_sets, pip), error = function(e) NULL)
    if (is.null(diagnostics)) return(NULL)
    # The human-readable conditional probabilities are rounded to three
    # decimals. Check consistency using their documented half-unit allowance.
    if (!is.null(cs_pip)) for (i in seq_along(cs_pip)) {
      s <- diagnostics$sets[[i]]; q <- cs_pip[[i]]
      if (abs(s$legacy_membership_sum - sum(q)) > .0005 * length(q) + 1e-8 ||
          s$conditional_union_probability + .0005 + 1e-8 < max(q)) return(NULL)
    }
  } else if (require_diagnostics) return(NULL)
  list(pip = pip, credible_sets = credible_sets, cs_pip = cs_pip,
       cs_diagnostics = diagnostics)
}


# --- the shared scenario-setup core ------------------------------------------
# Returns list(.fb_joint_cache = <fingerprint -> parsed result>,
#              .fb_joint_annotated = TRUE/FALSE, .fb_joint_head = head).
# Empty annotated-arm cache means the joint run failed, not another estimator.
.fb_joint_scenario_setup <- function(genotypes, regions, user_args, prior_head) {
  n_regions <- length(regions)
  annotation_inputs <- lapply(seq_len(n_regions), function(i) {
    A <- genotypes[[i]]$annotations_matrix
    if (is.null(A)) A <- regions[[i]]$annotations_matrix
    A
  })
  if (all(vapply(annotation_inputs, is.null, logical(1)))) {
    return(list(.fb_joint_annotated = FALSE, .fb_joint_head = prior_head))
  }

  # collect regions that have an annotation matrix (joint prior needs annotations)
  specs <- vector("list", n_regions)
  keys  <- character(n_regions)
  m_ann <- NULL
  for (i in seq_len(n_regions)) {
    A <- annotation_inputs[[i]]
    z <- regions[[i]]$z
    if (is.null(A) || !is.matrix(A) || is.null(z) || nrow(A) != length(z)) {
      stop("Joint Functional BEATRICE requires a valid annotation matrix for every region.", call. = FALSE)
    }
    if (is.null(m_ann)) m_ann <- ncol(A) else if (ncol(A) != m_ann) {
      stop("Joint Functional BEATRICE annotation columns must agree across regions.", call. = FALSE)
    }
    specs[[i]] <- list(z = z, LD = genotypes[[i]]$LD, A = A,
                       variant_ids = genotypes[[i]]$variant_ids, n = genotypes[[i]]$n)
    keys[i] <- .fb_fingerprint(z, specs[[i]]$LD, A, specs[[i]]$n, specs[[i]]$variant_ids)
  }

  beatrice_dir <- user_args$beatrice_dir
  python       <- user_args$python
  script       <- file.path(beatrice_dir, "beatrice_joint.py")
  if (is.null(beatrice_dir) || is.null(python) || !file.exists(script)) {
    stop("Joint Functional BEATRICE needs Python and beatrice_joint.py; no fallback is performed.", call. = FALSE)
  }
  py <- if (file.exists(path.expand(python))) .fm_python_path(python) else python

  work_dir <- tempfile(pattern = "fbjoint_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  # write regions + manifest
  written <- vector("list", n_regions)
  man <- c("z\tLD\tannot\ttarget\tN")
  for (i in seq_len(n_regions)) {
    w <- .fb_joint_write_region(work_dir, i, specs[[i]]$z, specs[[i]]$LD,
                                specs[[i]]$A, specs[[i]]$variant_ids)
    written[[i]] <- w
    man <- c(man, sprintf("%s\t%s\t%s\t%s\t%d",
                          w$z, w$LD, w$annot, w$target, as.integer(specs[[i]]$n)))
  }
  manifest_path <- file.path(work_dir, "manifest.tsv")
  writeLines(man, manifest_path)

  gv <- function(k, d) { v <- user_args[[k]]; if (is.null(v)) d else v }
  args <- c(
    script,
    "--manifest",             manifest_path,
    "--prior_head",           prior_head,
    "--max_iter",             as.character(as.integer(gv("max_iter", 1500))),
    "--n_caus",               as.character(as.integer(gv("n_caus", 5))),
    "--sigma_sq",             as.character(gv("sigma_sq", 0.05)),
    "--gamma_coverage",       as.character(gv("gamma_coverage", 0.95)),
    "--sparse_concrete",      as.character(as.integer(gv("sparse_concrete", 50))),
    "--prior_regularisation", as.character(gv("prior_regularisation", 1.0)),
    "--lambda_l1",            as.character(gv("lambda_l1", 0.01)),
    "--hierarchy_M",          as.character(gv("hierarchy_M", 10.0))
  )

  started <- proc.time()
  run_output <- tryCatch(
    system2(py, args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) structure(conditionMessage(e), class = "fb_joint_error"))
  elapsed <- as.numeric((proc.time() - started)["elapsed"])

  if (inherits(run_output, "fb_joint_error") || .beatrice_process_failed(run_output)) {
    stop(paste(c(sprintf("Joint Functional BEATRICE (%s) execution failed:", prior_head),
                 tail(as.character(run_output), 20L)), collapse = "\n"), call. = FALSE)
  }

  # parse per-region outputs into the fingerprint-keyed cache
  cache <- list()
  for (i in seq_len(n_regions)) {
    parsed <- .fb_joint_parse_output(written[[i]]$target,
                                     written[[i]]$variant_ids, written[[i]]$p,
                                     .fb_joint_requires_cs_diagnostics(beatrice_dir))
    if (!is.null(parsed)) {
      parsed$variant_ids <- written[[i]]$original_ids
      parsed$variant_alignment <- written[[i]]$variant_alignment
      parsed$runtime_seconds <- elapsed / n_regions
      parsed$pooled_runtime_seconds <- elapsed
      parsed$n_loci_pooled <- n_regions
      cache[[keys[i]]] <- parsed
    }
  }
  if (length(cache) != n_regions) {
    stop("Joint Functional BEATRICE did not return valid output for every region.", call. = FALSE)
  }
  list(.fb_joint_cache = cache, .fb_joint_annotated = TRUE, .fb_joint_head = prior_head)
}


# --- public scenario-setup hooks (looked up by name in run_methods) ----------
#' Prepare a jointly learned Functional BEATRICE annotation prior
#' @param genotypes List of per-region LD, variant and annotation inputs.
#' @param regions List of corresponding summary-statistic inputs; ordinary runs omit truth.
#' @param user_args Prespecified arguments for the Functional BEATRICE Python runner.
#' @return Shared fit cache and declared annotation-prior head for region dispatch.
#' @details The pooled variant uses a shared linear head; xregion uses a shared
#'   LassoNet head. Annotation learning on supplied loci is transductive and
#'   must be declared when comparing methods. See the scientific contract.
#' @export
run_fb_pooled_scenario_setup <- function(genotypes, regions, user_args) {
  .fb_joint_scenario_setup(genotypes, regions, user_args, prior_head = "linear")
}
#' @rdname run_fb_pooled_scenario_setup
#' @export
run_fb_xregion_scenario_setup <- function(genotypes, regions, user_args) {
  .fb_joint_scenario_setup(genotypes, regions, user_args, prior_head = "lassonet")
}


# --- the shared per-region core ----------------------------------------------
.fb_joint_region <- function(region_geno, region_pheno, method_name,
                             .fb_joint_cache = NULL, .fb_joint_annotated = TRUE,
                             .fb_joint_head = NULL, ...) {
  p <- length(region_pheno$z)
  A <- region_geno$annotations_matrix %||% region_pheno$annotations_matrix
  key <- .fb_fingerprint(region_pheno$z, region_geno$LD, A,
                         region_geno$n, region_geno$variant_ids)
  entry <- if (!is.null(.fb_joint_cache)) .fb_joint_cache[[key]] else NULL

  if (!is.null(entry)) {
    return(list(
      pip             = entry$pip,
      credible_sets   = entry$credible_sets,
      method          = method_name,
      input_type      = "summary",
      params          = list(prior_head = .fb_joint_head, joint = TRUE),
      runtime_seconds = entry$runtime_seconds,
      additional      = c(list(cs_pip = entry$cs_pip, joint_fallback = FALSE,
                             variant_ids = entry$variant_ids,
                             variant_alignment = entry$variant_alignment,
                             pooled_runtime_seconds = entry$pooled_runtime_seconds,
                             n_loci_pooled = entry$n_loci_pooled,
                             pooled_across_loci = TRUE), .beatrice_posterior_metadata(entry))
    ))
  }

  if (!identical(.fb_joint_annotated, FALSE)) {
    return(list(pip = rep(NA_real_, p), credible_sets = list(), method = method_name,
                input_type = "summary", params = list(prior_head = .fb_joint_head, joint = TRUE),
                runtime_seconds = NA_real_, additional = list(joint_fallback = FALSE),
                error = "Annotated joint fit has no valid cached prediction; no fallback is performed."))
  }

  # The prespecified annotation-free ablation is ordinary locus-only BEATRICE.
  fb <- tryCatch(
    run_functional_beatrice_region(region_geno = region_geno,
                                   region_pheno = region_pheno, ...),
    error = function(e) list(pip = rep(NA_real_, p), credible_sets = list(),
                             method = method_name, input_type = "summary",
                             params = list(), runtime_seconds = NA_real_,
                             additional = list(), error = conditionMessage(e)))
  fb$method <- method_name
  if (is.null(fb$additional)) fb$additional <- list()
  fb$additional$joint_fallback <- FALSE
  fb$additional$annotation_free_equivalent <- "beatrice_locus_only"
  fb$additional$joint <- FALSE
  fb
}

#' Return a region prediction from a joint Functional BEATRICE fit
#' @param region_geno Region LD, variant and annotation inputs.
#' @param region_pheno Region summary statistics.
#' @param ... Shared scenario cache or arguments for the per-region fallback
#'   when annotations are absent.
#' @return Standard method prediction with PIPs, native sets and fit metadata.
#' @export
run_fb_pooled_region <- function(region_geno, region_pheno, ...) {
  .fb_joint_region(region_geno, region_pheno, method_name = "fb_pooled", ...)
}
#' @rdname run_fb_pooled_region
#' @export
run_fb_xregion_region <- function(region_geno, region_pheno, ...) {
  .fb_joint_region(region_geno, region_pheno, method_name = "fb_xregion", ...)
}
