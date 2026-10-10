# Frozen campaign orchestration helpers. The public package API remains the
# simulation/method/evaluation functions; scripts/campaign/ uses these helpers.

.campaign_atomic_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, version = 2)
  if (!file.rename(tmp, path)) stop("Atomic checkpoint rename failed: ", path)
  invisible(path)
}

.campaign_object_sha <- function(object) {
  path <- tempfile()
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 2)
  .challenge_sha256(path)
}

.campaign_snapshot <- function(repo) {
  paths <- c(list.files(file.path(repo, "R"), "\\.R$", full.names = TRUE),
             list.files(file.path(repo, "inst/extdata"), full.names = TRUE),
             list.files(file.path(repo, "BEATRICE_annot_sparse"), "\\.py$",
                        full.names = TRUE, recursive = TRUE),
             list.files(file.path(repo, "scripts/campaign"), "\\.(R|py|sh)$", full.names = TRUE,
                        recursive = TRUE))
  paths <- sort(paths[!dir.exists(paths)], method = "radix")
  data.frame(path = substring(paths, nchar(repo) + 2L),
             sha256 = unname(vapply(paths, .challenge_sha256, "")), row.names = NULL,
             stringsAsFactors = FALSE)
}

.campaign_input_manifest <- function(path, directories) {
  if (nzchar(path)) {
    data <- utils::read.csv(path, stringsAsFactors = FALSE)
    if (!all(c("path", "sha256") %in% names(data)) || anyDuplicated(data$path) ||
        any(!grepl("^[a-f0-9]{64}$", data$sha256))) stop("Invalid input manifest")
    data$path <- normalizePath(data$path, mustWork = TRUE)
    data <- data[order(data$path, method = "radix"), c("path", "sha256")]
    rownames(data) <- NULL
    return(data)
  }
  paths <- unique(unlist(lapply(directories[nzchar(directories)], function(d)
    list.files(d, recursive = TRUE, full.names = TRUE))))
  paths <- sort(paths[!dir.exists(paths)], method = "radix")
  data.frame(path = normalizePath(paths, mustWork = TRUE),
              sha256 = unname(vapply(paths, .challenge_sha256, "")), row.names = NULL,
              stringsAsFactors = FALSE)
}

.campaign_check_inputs <- function(paths, manifest) {
  paths <- normalizePath(paths, mustWork = TRUE)
  expected <- manifest$sha256[match(paths, manifest$path)]
  if (anyNA(expected)) stop("A used input is absent from the frozen data manifest")
  for (i in seq_along(paths))
    if (!identical(.challenge_sha256(paths[i]), expected[i])) stop("Frozen input changed: ", paths[i])
  invisible(TRUE)
}

.campaign_runtime_fingerprint <- function(method_args) {
  versions <- utils::installed.packages()[, c("Package", "Version"), drop = FALSE]
  versions <- versions[order(versions[, 1L], versions[, 2L], method = "radix"), , drop = FALSE]
  rownames(versions) <- NULL
  paths <- unique(unlist(lapply(method_args, function(args)
     unlist(args[intersect(names(args), c("finemap_path", "paintor_path", "caviar_path",
                         "beatrice_dir", "sparsepro_dir", "finemap_inf_dir", "python"))]))))
  artifacts <- list()
  for (path in paths) {
    resolved <- if (file.exists(path)) normalizePath(path) else unname(Sys.which(path))
    if (!nzchar(resolved)) { artifacts[[path]] <- "missing"; next }
    if (dir.exists(resolved)) {
      files <- sort(list.files(resolved, "\\.py$", full.names = TRUE, recursive = TRUE), method = "radix")
      artifacts[[path]] <- data.frame(path = substring(files, nchar(resolved) + 2L),
                                    sha256 = unname(vapply(files, .challenge_sha256, "")), row.names = NULL)
    } else artifacts[[path]] <- .challenge_sha256(resolved)
  }
  py <- unique(vapply(method_args, function(a) a$python %||% "", ""))
  py <- py[nzchar(py)]
  py_versions <- lapply(py, function(exe) tryCatch(system2(exe, c("-c", shQuote(
      "import importlib.metadata as m; print('\\n'.join(sorted(d.metadata['Name']+'=='+d.version for d in m.distributions())))")),
      stdout = TRUE, stderr = FALSE), error = function(e) "python_metadata_unavailable"))
  list(R = R.version.string, package_versions = versions, artifacts = artifacts,
       python_versions = stats::setNames(py_versions, py),
       threading = Sys.getenv(c("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")))
}

.campaign_method_fingerprint <- function(env, method) {
  # Reuse is allowed only for exactly equal reachable authored functions and
  # globals. Dynamic evaluation and private-truth readers conservatively opt out.
  name <- env$.FM_REGISTRY[[method]]
  hook <- paste0("run_", method, "_scenario_setup")
  pending <- c(name, if (exists(hook, envir = env, inherits = FALSE)) hook)
  functions <- globals <- list(); seen <- character()
  while (length(pending)) {
    nm <- pending[1L]; pending <- pending[-1L]
    if (nm %in% seen) next
    seen <- c(seen, nm)
    fun <- get(nm, envir = env, inherits = FALSE)
    body <- paste(deparse(fun, width.cutoff = 500L), collapse = "\n")
    found <- codetools::findGlobals(fun, merge = FALSE)
    if (any(found$functions %in% c("get", "eval", "evalq", "parse", "source", "sys.source", "load")) ||
        (method != "polyfun_oracle" && grepl("truth", body, fixed = TRUE)))
      return(list(eligible = FALSE, reason = "dynamic lookup or private-truth reader"))
    functions[[nm]] <- body
    for (variable in found$variables) if (exists(variable, envir = env, inherits = FALSE)) {
      value <- get(variable, envir = env, inherits = FALSE)
      if (is.environment(value) || is.function(value) || typeof(value) == "externalptr")
        return(list(eligible = FALSE, reason = "mutable authored global"))
      globals[[variable]] <- value
    }
    pending <- c(pending, found$functions[vapply(found$functions,
       function(f) exists(f, envir = env, mode = "function", inherits = FALSE), TRUE)])
  }
  functions <- functions[order(names(functions), method = "radix")]
  if (length(globals)) globals <- globals[order(names(globals), method = "radix")]
  list(eligible = TRUE, fingerprint = .campaign_object_sha(list(functions = functions, globals = globals)))
}

.campaign_assert_snapshot <- function(repo, expected) {
  actual <- .campaign_snapshot(repo)
  if (!identical(actual, expected))
    stop("Campaign source changed; create a new campaign rather than reuse checkpoints")
  invisible(TRUE)
}

.campaign_default_args <- function(repo, capacity = c("original", "common")) {
  capacity <- match.arg(capacity)
  tool <- Sys.getenv("FMB_TOOL_ROOT", path.expand("~/tools"))
  py <- Sys.getenv("FMB_PYTHON", "python")
  fb <- Sys.getenv("FMB_BEATRICE_DIR", file.path(repo, "BEATRICE_annot_sparse"))
  neural <- list(beatrice_dir = fb, python = py, max_iter = 2000L, n_caus = 10,
                 sigma_sq = .0899, gamma_coverage = .95, sparse_concrete = 51L)
  functional <- c(neural, list(prior_regularisation = .1255, lambda_l1 = .1223,
                             hierarchy_M = 10.15))
  list(susie = list(L = 10L, coverage = .95),
       susie_inf = list(L = 10L), abf = list(prior_variance = .04),
       carma = list(num.causal = 10L), marginal_z = list(coverage = .95),
       polyfun_oracle = list(L = 10L), polyfun_est = list(L = 10L),
       polyfun_ldsc = list(L = 10L),
       sbayesrc = list(n_iter = 300L, burn_in = 150L, gamma_update_every = 10L),
       finemap = list(finemap_path = Sys.getenv("FMB_FINEMAP_PATH",
                         file.path(tool, "finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64")),
                      n_causal = if (capacity == "original") 5L else 10L,
                      n_iter = 100000L, prior_std = .05, coverage = .95),
       paintor = list(paintor_path = Sys.getenv("FMB_PAINTOR_PATH",
                         file.path(tool, "PAINTOR_V3.0/PAINTOR")),
                      max_causal = 2L, mcmc = FALSE, coverage = .95),
       caviar = list(caviar_path = Sys.getenv("FMB_CAVIAR_PATH",
                         file.path(tool, "caviar/CAVIAR-C++/CAVIAR")),
                     max_causal = 2L, rho = .95),
       finimom = list(maxsize = 10L, niter = 12500L, burnin = 2500L, cs_level = .95),
       finemap_inf = list(finemap_inf_dir = Sys.getenv("FMB_FINEMAP_INF_DIR",
                         file.path(tool, "fine-mapping-inf")), python = py,
                         num_sparse_effects = 10L, coverage = .95, purity = .5),
       sparsepro = list(sparsepro_dir = Sys.getenv("FMB_SPARSEPRO_DIR",
                         file.path(tool, "SparsePro")), python = py, cthres = .95),
       funmap = list(python = py, max_iter = 100L, tol = 5e-5),
       beatrice = neural, functional_beatrice = functional,
       fb_pooled = functional, fb_xregion = functional)
}

.campaign_load_source <- function(repo) {
  env <- new.env(parent = globalenv())
  for (path in list.files(file.path(repo, "R"), "\\.R$", full.names = TRUE))
    sys.source(path, envir = env)
  env
}

.campaign_legacy_index <- function(root, methods) {
  sims <- list.files(root, "^sim\\.rds$", recursive = TRUE, full.names = TRUE)
  if (!length(sims)) stop("No frozen legacy sim.rds found under --legacy-root")
  rows <- list(); j <- 0L
  for (path in sort(sims)) {
    sim <- readRDS(path)
    sim_sha <- .challenge_sha256(path)
    genotype_sha <- .campaign_object_sha(sim$genotypes)
    if (is.null(sim$genotypes) || !length(sim$scenarios)) stop("Malformed legacy simulation: ", path)
    for (sc in seq_along(sim$scenarios)) {
      sd <- file.path(dirname(path), sprintf("scenario_%03d", sc))
      candidates <- file.path(sd, c("results.rds", "results_supp.rds"))
      candidates <- candidates[file.exists(candidates)]
      j <- j + 1L
      rows[[j]] <- data.frame(dataset_id = sprintf("legacy_%05d", j),
          family = "legacy", row = NA_integer_, scenario = sc, sim_file = normalizePath(path),
          sim_sha256 = sim_sha,
          genotypes_sha256 = genotype_sha,
          result_files = paste(normalizePath(candidates), collapse = "|"),
          result_sha256 = paste(vapply(candidates, .challenge_sha256, ""), collapse = "|"),
          original_label = basename(dirname(path)), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

.campaign_read_frozen <- function(path, sha256) {
  if (!file.exists(path) || !identical(.challenge_sha256(path), sha256))
    stop("Frozen artifact missing or changed: ", path)
  readRDS(path)
}

.campaign_merge_archived <- function(paths, hashes) {
  files <- strsplit(paths, "|", fixed = TRUE)[[1L]]
  sha <- strsplit(hashes, "|", fixed = TRUE)[[1L]]
  if (!length(files) || !nzchar(files[1L])) return(NULL)
  if (length(files) != length(sha)) stop("Archived result manifest is malformed")
  merged <- list()
  for (i in seq_along(files)) {
    z <- .campaign_read_frozen(files[i], sha[i])
    for (m in z$methods_run) merged[[m]] <- z[[m]]
  }
  merged$methods_run <- setdiff(names(merged), "methods_run")
  merged
}

.campaign_one_scenario <- function(sim, scenario) {
  sim$scenarios <- sim$scenarios[scenario]
  sim$scenarios[[1L]]$scenario_id <- 1L
  sim
}

.campaign_corrected_phenotype <- function(sim, scenario, seed) {
  mini <- .campaign_one_scenario(sim, scenario)
  sc <- mini$scenarios[[1L]]; pa <- sim$params
  generated <- simulate_phenotypes(sim$genotypes, S = sc$S, phi = sc$phi,
      model = sc$model %||% pa$model, p_causal = sc$p_causal %||% .5,
      inf_model = pa$inf_model %||% "beatrice",
      effect_distribution = pa$effect_distribution %||% "normal",
      effect_variance = pa$effect_variance %||% .36,
      annotations = pa$annotation_type %||% "none", n_annotations = pa$n_annotations %||% 3L,
      annotation_proportions = pa$annotation_proportions,
      enrichment = pa$enrichment_actual %||% pa$enrichment,
      annotation_correlation = pa$annotation_correlation %||% 0,
      seed = seed, verbose = FALSE)
  fields <- c("y", "z", "beta_hat", "se", "statistic_type", "degrees_of_freedom",
              "annotations_matrix", "truth")
  mini$scenarios[[1L]]$regions <- lapply(generated, function(x) x[fields])
  mini$params$campaign_phenotype_seed <- seed
  mini$params$campaign_input_contrast <- "corrected phenotype; frozen original genotypes"
  mini
}

.campaign_panel_simulation <- function(panel, indices) {
  xs <- panel$inputs[indices]
  regions <- lapply(xs, function(x) list(z = x$z, beta_hat = x$beta_hat, se = x$se,
                                        annotations_matrix = x$annotations))
  geno <- lapply(xs, function(x) list(LD = x$LD, n = x$n, p = length(x$z),
      maf = rep(.3, length(x$z)), variant_ids = x$variant_id,
      annotations_matrix = x$annotations))
  list(genotypes = geno, scenarios = list(list(scenario_id = 1L, regions = regions)),
       params = list(annotation_type = "synthetic_continuous", input_type = "normal_score"))
}

.campaign_fit_to_prediction <- function(fit, input) {
  event <- .campaign_event_metadata(fit$method %||% "unknown", fit$additional %||% list())
  if (!is.null(fit$error)) return(list(task_id = input$task_id, status = "failed",
      error = fit$error, runtime_seconds = fit$runtime_seconds,
      pip_event = event$pip_event, probability_target_aligned = event$probability_target_aligned,
      credible_set_semantics = fit$additional$credible_set_semantics %||% "unspecified_method_native_sets"))
  list(task_id = input$task_id, variant_id = input$variant_id, pip = fit$pip,
       status = "ok", credible_sets = lapply(fit$credible_sets, function(s) input$variant_id[s]),
       converged = fit$additional$converged %||% NA,
       runtime_seconds = fit$runtime_seconds,
       pip_event = event$pip_event, probability_target_aligned = event$probability_target_aligned,
       credible_set_semantics = fit$additional$credible_set_semantics %||% "unspecified_method_native_sets",
       implementation = fit$additional$implementation %||% fit$method)
}

.campaign_event_metadata <- function(method, additional = list()) {
  event <- additional$pip_event %||% if (method == "sbayesrc") "any_nonzero_mixture_component" else
    if (method == "marginal_z") "normalized_ranking_score" else "method_model_sparse_inclusion"
  aligned <- if (method == "marginal_z") NA else if (method == "sbayesrc") FALSE else TRUE
  list(score_target = "sparse_component_membership", pip_event = event,
       probability_target_aligned = aligned,
       probability_score_interpretation = if (method == "sbayesrc")
         "operational_sparse_membership_diagnostic; event mismatch" else
         if (method == "marginal_z") "not_a_probability_score" else "operational_sparse_membership_score")
}

.campaign_exact_log_loss <- function(y, score) {
  # Index the realized class to avoid 0 * Inf at valid endpoint probabilities.
  y <- as.logical(y)
  -mean(c(log(score[y]), log1p(-score[!y])))
}

.campaign_ap_failure_zero <- function(ap, success, n_observed_causal) {
  ifelse(n_observed_causal > 0, ifelse(success, ap, 0), NA_real_)
}

.campaign_set_metrics <- function(sets, causal, semantics, success = TRUE) {
  native_global <- identical(semantics, "native_all_causal_containment_set")
  reported <- success && length(sets) > 0L
  n <- length(causal)
  covered <- if (n && reported) any(vapply(sets, function(s) all(causal %in% s), TRUE)) else FALSE
  list(credible_set_semantics = semantics,
       cs_distinct_recall = if (native_global || !n) NA_real_ else
         .challenge_distinct_hits(sets, causal) / n,
       native_all_causal_coverage_reported = if (native_global && n && reported)
         as.numeric(covered) else NA_real_,
       native_all_causal_recovery_all_attempts = if (native_global && n)
         as.numeric(covered) else NA_real_,
       cs_causal_variant_recall_all_attempts = if (n)
         sum(causal %in% unlist(sets, use.names = FALSE)) / n else NA_real_,
       cs_support_size = length(unique(unlist(sets, use.names = FALSE))))
}

.campaign_challenge_task_metrics <- function(table, predictions, panel, method) {
  extras <- lapply(seq_len(nrow(table)), function(i) {
    id <- table$task_id[i]; input <- panel$inputs[[id]]; truth <- panel$truth[[id]]
    pred <- predictions[[match(id, vapply(predictions, function(x) x$task_id, ""))]]
    ok <- isTRUE(table$success[i])
    if (ok) pred <- .challenge_validate_prediction(pred, input)
    semantics <- pred$credible_set_semantics %||% if (method == "caviar")
      "native_all_causal_containment_set" else "unspecified_method_native_sets"
    values <- .campaign_set_metrics(if (ok) pred$credible_sets %||% list() else list(),
                                    truth$causal_variant_id, semantics, ok)
    values <- c(values, .campaign_event_metadata(method, pred))
    values$n_observed_causal <- sum(input$variant_id %in% truth$causal_variant_id)
    values$log_loss <- if (ok && method != "marginal_z")
      .campaign_exact_log_loss(input$variant_id %in% truth$causal_variant_id, pred$pip) else NA_real_
    total <- (truth$beta_sparse + truth$beta_background)[input$variant_id]
    values$n_any_nonzero_total_effect <- sum(total != 0)
    values$brier_any_nonzero_total_effect <- if (ok && method == "sbayesrc")
      mean((pred$pip - (total != 0))^2) else NA_real_
    values$log_loss_any_nonzero_total_effect <- if (ok && method == "sbayesrc")
      .campaign_exact_log_loss(total != 0, pred$pip) else NA_real_
    as.data.frame(values, stringsAsFactors = FALSE)
  })
  extras <- do.call(rbind, extras)
  for (field in names(extras)) table[[field]] <- extras[[field]]
  table
}

.campaign_challenge_predict <- function(panel, indices, method, args) {
  xs <- panel$inputs[indices]
  if (method %in% c("exact", "forced_one"))
    return(run_challenge(xs, stats::setNames(list(if (method == "exact") challenge_exact
                                                else challenge_forced_one), method))[[1L]])
  if (method == "polyfun_oracle") stop("Enrichment oracle is not a public-input challenge predictor")
  unsupported <- vapply(xs, function(x) !is.null(x$cohorts) || !is.null(x$n_per_variant), TRUE)
  pred <- lapply(xs, function(x) list(task_id = x$task_id, status = "failed",
              error = "unsupported_input: regional single-cohort scalar-N method"))
  good <- which(!unsupported)
  if (!length(good)) return(pred)
  sim <- .campaign_panel_simulation(panel, indices[good])
  # Exact mixture/empirical surrogate knows these are normal score statistics.
  if (method %in% c("sbayesrc", "susie", "susie_inf", "polyfun_est", "polyfun_ldsc")) {
    args$statistic_type <- "normal_score"
    args$residual_variance <- 1
    if (method != "sbayesrc") args$estimate_residual_variance <- FALSE
  }
  if (method == "finemap_inf") args$statistic_type <- "normal_score"
  out <- run_methods(sim, method, stats::setNames(list(args), method), verbose = FALSE)
  out <- .campaign_reject_substitution(out, sim, method)
  for (k in seq_along(good)) pred[[good[k]]] <- .campaign_fit_to_prediction(
      out[[method]]$results[[k]], xs[[good[k]]])
  pred
}

.campaign_task_table <- function(result, sim, method, metadata) {
  fits <- .annotate_fits_with_truth(.complete_fit_grid(result[[method]]$results, sim, method), sim)
  do.call(rbind, lapply(fits, function(f) {
    ok <- is.null(f$error); y <- seq_len(f$n_variants) %in% f$causal_indices
    score <- if (ok) f$pip else rep(NA_real_, f$n_variants)
    called <- if (ok) score >= .95 else rep(FALSE, f$n_variants)
    cs <- if (ok) f$credible_sets else list()
    semantics <- f$additional$credible_set_semantics %||% if (method == "caviar")
      "native_all_causal_containment_set" else "unspecified_method_native_sets"
    set_metrics <- .campaign_set_metrics(cs, which(y), semantics, ok)
    total <- sim$scenarios[[1L]]$regions[[f$region_id]]$truth$beta_total
    has_total <- length(total) == f$n_variants && all(is.finite(total))
    tp <- sum(called & y); fp <- sum(called & !y)
    data.frame(dataset_id = metadata$dataset_id, family = metadata$family,
        stage = metadata$stage, method = method, region_id = f$region_id,
        cluster_id = paste(metadata$input_panel_id, f$region_id, sep = ":"),
        panel_cluster_id = metadata$input_panel_id,
        source_locus_cluster = sim$genotypes[[f$region_id]]$vcf_source %||%
            paste(metadata$input_panel_id, f$region_id, sep = ":"),
        S = f$S %||% sum(y), phi = f$phi %||% NA_real_,
        p_causal = f$p_causal %||% NA_real_, iter = f$iter %||% NA_integer_,
        annotation_type = sim$params$annotation_type %||% NA_character_,
        enrichment_fold = if (length(sim$params$enrichment_actual %||% sim$params$enrichment))
           max(sim$params$enrichment_actual %||% sim$params$enrichment) else NA_real_,
        n_gwas = sim$genotypes[[f$region_id]]$n,
        n_ref = sim$params$n_ref %||% NA_integer_,
        nominal_p = f$region_size, actual_p = f$n_variants, success = ok,
        n_causal = sum(y), n_observed_causal = sum(y),
        tp = tp, fp = fp, recall = if (sum(y)) tp / sum(y) else NA_real_,
        any_call = any(called), ap = if (ok) .challenge_ap(y, score) else NA_real_,
        brier = if (ok && method != "marginal_z") mean((score - y)^2) else NA_real_,
        log_loss = if (ok && method != "marginal_z")
          .campaign_exact_log_loss(y, score) else NA_real_,
        n_any_nonzero_total_effect = if (has_total) sum(total != 0) else NA_integer_,
        brier_any_nonzero_total_effect = if (ok && has_total && method == "sbayesrc")
          mean((score - (total != 0))^2) else NA_real_,
        log_loss_any_nonzero_total_effect = if (ok && has_total && method == "sbayesrc")
          .campaign_exact_log_loss(total != 0, score) else NA_real_,
        sum_pip = if (ok) sum(score) else NA_real_,
        cs_count = length(cs),
        set_metrics,
        .campaign_event_metadata(method, f$additional %||% list()),
        runtime_seconds = f$runtime_seconds %||% NA_real_,
        converged = f$additional$converged %||% NA,
        joint_fallback = isTRUE(f$additional$joint_fallback),
        error = f$error %||% "", stringsAsFactors = FALSE)
  }))
}

.campaign_assert_plan <- function(plan, repo) {
  declared <- plan$plan_sha256
  plain <- plan; plain$plan_sha256 <- NULL
  if (!identical(declared, .campaign_object_sha(plain))) stop("Campaign plan was changed")
  .campaign_assert_snapshot(repo, plan$corrected_source)
  if (nzchar(plan$original_source))
    .campaign_assert_snapshot(plan$original_source, plan$original_source_hashes)
  if (!identical(plan$runtime_fingerprint, .campaign_runtime_fingerprint(plan$method_args)))
    stop("Frozen runtime dependencies changed; create a new campaign plan")
  invisible(TRUE)
}

.campaign_fit_cache_identity <- function(sim, method, args, seed, source_env, runtime,
                                         genotypes_sha256 = NULL) {
  authored <- .campaign_method_fingerprint(source_env, method)
  if (!isTRUE(authored$eligible)) return(authored)
  observed <- lapply(sim$scenarios[[1L]]$regions, function(r)
     if (method == "polyfun_oracle") r else r[setdiff(names(r), "truth")])
  input <- list(genotypes_sha256 = genotypes_sha256 %||% .campaign_object_sha(sim$genotypes),
                 regions = observed)
  list(eligible = TRUE, key = .campaign_object_sha(list(method = method, authored = authored$fingerprint,
          input = input, args = args, seed = seed, runtime = runtime)),
       authored_sha256 = authored$fingerprint, observed_input_sha256 = .campaign_object_sha(input))
}

.campaign_artifact_manifest <- function(paths) {
  paths <- sort(unique(normalizePath(paths, mustWork = TRUE)), method = "radix")
  data.frame(path = unname(paths), sha256 = unname(vapply(paths, .challenge_sha256, "")),
             row.names = NULL, stringsAsFactors = FALSE)
}

.campaign_validate_artifacts <- function(artifacts) {
  if (!is.data.frame(artifacts) || !all(c("path", "sha256") %in% names(artifacts)) ||
      !is.character(artifacts$path) || !is.character(artifacts$sha256) ||
      anyNA(artifacts[c("path", "sha256")]) || anyDuplicated(artifacts$path) ||
      any(!nzchar(artifacts$path)) || any(!grepl("^[a-f0-9]{64}$", artifacts$sha256)))
    stop("Malformed sealed artifact manifest")
  for (i in seq_len(nrow(artifacts))) {
    path <- artifacts$path[i]
    if (!file.exists(path)) stop("Sealed artifact is missing: ", path)
    if (!identical(.challenge_sha256(path), artifacts$sha256[i]))
      stop("Sealed artifact bytes changed: ", path)
  }
  invisible(TRUE)
}

.campaign_validate_checkpoint <- function(path, plan, task, validate_artifacts = TRUE) {
  if (!file.exists(path)) return(FALSE)
  z <- readRDS(path)
  if (!identical(z$plan_sha256, plan$plan_sha256) ||
      !identical(z$task, task) || !isTRUE(z$terminal))
    stop("Checkpoint is incompatible or incomplete: ", path)
  if (isTRUE(validate_artifacts) && !is.null(z$artifacts)) .campaign_validate_artifacts(z$artifacts)
  TRUE
}

.campaign_check_receipt <- function(directory, checkpoint, plan_sha256, id,
                                     receipt_name = "terminal.json") {
  path <- file.path(directory, receipt_name)
  if (!file.exists(path)) return(FALSE)
  value <- tryCatch(jsonlite::fromJSON(path), error = function(e)
    stop("Malformed terminal receipt: ", path))
  if (!identical(value$schema, 1L) || !isTRUE(value$terminal) ||
      !identical(value$id, id) || !identical(value$plan_sha256, plan_sha256) ||
      !is.character(value$checkpoint_sha256) || length(value$checkpoint_sha256) != 1L ||
      !grepl("^[a-f0-9]{64}$", value$checkpoint_sha256))
    stop("Foreign or malformed terminal receipt: ", path)
  if (!file.exists(checkpoint)) stop("Receipted checkpoint is missing: ", checkpoint)
  if (!identical(.challenge_sha256(checkpoint), value$checkpoint_sha256))
    stop("Receipted checkpoint bytes changed: ", checkpoint)
  TRUE
}

.campaign_claim <- function(directory) {
  dir.create(dirname(directory), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(directory, showWarnings = FALSE)) stop("Task has an existing lock: ", directory)
  writeLines(c(paste("pid", Sys.getpid()), paste("host", Sys.info()[["nodename"]]),
               paste("created", format(Sys.time(), tz = "UTC", usetz = TRUE))),
             file.path(directory, "owner.txt"))
  invisible(directory)
}

.campaign_start_attempt <- function(directory, plan, task, dataset) {
  token <- paste0(dataset$dataset_id, "-", Sys.getpid(), "-", format(Sys.time(), "%Y%m%d%H%M%OS6"))
  .campaign_atomic_rds(list(plan_sha256 = plan$plan_sha256, task = task, dataset = dataset,
        started = Sys.time(), pid = Sys.getpid(), host = Sys.info()[["nodename"]]),
        file.path(directory, "attempts", paste0(token, ".rds")))
}

.campaign_receipt <- function(directory, checkpoint, plan_sha256, id,
                               receipt_name = "terminal.json") {
  if (!grepl("^[A-Za-z0-9_.-]+$", id) || !grepl("^[a-f0-9]{64}$", plan_sha256))
    stop("Invalid terminal receipt identity")
  if (!grepl("^[A-Za-z0-9_.-]+$", receipt_name)) stop("Invalid receipt filename")
  path <- file.path(directory, receipt_name)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  writeLines(sprintf(paste0('{"schema":1,"terminal":true,"id":"%s",',
       '"plan_sha256":"%s","checkpoint_sha256":"%s"}'),
       id, plan_sha256, .challenge_sha256(checkpoint)), tmp)
  if (!file.rename(tmp, path)) stop("Terminal receipt rename failed")
  invisible(path)
}

.campaign_analysis_manifest <- function(directory, plan, status, summary = list()) {
  files <- list.files(directory, full.names = TRUE, recursive = TRUE)
  files <- files[!basename(files) %in% c("analysis-manifest.rds", "analysis-manifest.json")]
  manifest <- list(schema = 1L, plan_sha256 = plan$plan_sha256, status = status,
       collected_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
       source = plan$corrected_source, runtime_fingerprint = plan$runtime_fingerprint,
       summary = summary, artifacts = .campaign_artifact_manifest(files))
  path <- file.path(directory, "analysis-manifest.rds")
  .campaign_atomic_rds(manifest, path)
  json_path <- file.path(directory, "analysis-manifest.json")
  json_tmp <- paste0(json_path, ".tmp.", Sys.getpid())
  jsonlite::write_json(list(schema = 1L, plan_sha256 = plan$plan_sha256, status = status,
       collected_utc = manifest$collected_utc, summary = summary,
       manifest_sha256 = .challenge_sha256(path)), json_tmp,
       auto_unbox = TRUE, pretty = TRUE)
  if (!file.rename(json_tmp, json_path)) stop("Analysis manifest rename failed")
  invisible(manifest)
}

.campaign_reject_substitution <- function(result, sim, method) {
  if (!method %in% c("fb_pooled", "fb_xregion")) return(result)
  for (i in seq_along(result[[method]]$results)) {
    fit <- result[[method]]$results[[i]]
    geno <- sim$genotypes[[fit$region_id]]
    annotated <- !is.null(geno$annotations_matrix) ||
       !is.null(sim$scenarios[[fit$scenario_id]]$regions[[fit$region_id]]$annotations_matrix)
    if (annotated && isTRUE(fit$additional$joint_fallback)) {
      fit$error <- "Joint annotated method substituted a per-region fit; rejected by campaign contract"
      result[[method]]$results[[i]] <- fit
    }
  }
  result[[method]]$n_failed <- sum(vapply(result[[method]]$results,
                                          function(f) !is.null(f$error), TRUE))
  result
}


.campaign_dependency_inventory <- function(method_args) {
  packages <- c("susieR", "CARMA", "finimom", "nnet", "reticulate", "sim1000G")
  r <- data.frame(dependency = packages, kind = "R_package", available = vapply(packages,
       requireNamespace, TRUE, quietly = TRUE), version = vapply(packages, function(p)
       if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "missing", ""))
  fields <- list(finemap = "finemap_path", paintor = "paintor_path", caviar = "caviar_path",
      sparsepro = "sparsepro_dir", finemap_inf = "finemap_inf_dir", beatrice = "beatrice_dir")
  z <- do.call(rbind, lapply(names(fields), function(m) {
    path <- method_args[[m]][[fields[[m]]]]
    resolved <- if (file.exists(path)) normalizePath(path) else unname(Sys.which(path))
    data.frame(dependency = m, kind = "external_path", available = nzchar(resolved),
               version = if (nzchar(resolved)) resolved else path)
  }))
  rbind(r, z)
}
