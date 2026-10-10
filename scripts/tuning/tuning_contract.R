# Dependency-free helpers shared by the tuning CLI and scientific tests.
.tuning_integer <- function(x, name, minimum = 0L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) || x < minimum)
    stop(name, " must be one integer >= ", minimum, call. = FALSE)
  as.integer(x)
}

subset_tuning_scenario <- function(sim, scenario, n_regions = 0L,
                                   rotate_regions = FALSE) {
  scenario <- .tuning_integer(scenario, "scenario", 1L)
  n_regions <- .tuning_integer(n_regions, "regions")
  if (isTRUE(rotate_regions))
    stop("One-region-per-scenario rotation confounds architecture with locus size; ",
         "use the same --regions subset in every scenario.", call. = FALSE)
  n_sc <- length(sim$scenarios)
  n_reg <- length(sim$genotypes)
  if (!n_sc || !n_reg || scenario > n_sc)
    stop("Scenario is outside the cached simulation grid.", call. = FALSE)
  if (n_regions > n_reg)
    stop("Requested region count exceeds the cached simulation.", call. = FALSE)
  if (length(sim$scenarios[[scenario]]$regions) != n_reg)
    stop("Cached scenario and genotype region counts disagree.", call. = FALSE)
  keep <- seq_len(if (n_regions == 0L) n_reg else n_regions)
  mini <- sim
  mini$genotypes <- sim$genotypes[keep]
  mini$scenarios <- sim$scenarios[scenario]
  mini$scenarios[[1]]$scenario_id <- 1L
  mini$scenarios[[1]]$regions <- sim$scenarios[[scenario]]$regions[keep]
  # Region IDs are list indices throughout run_methods/evaluate_methods.
  for (k in seq_along(keep)) {
    mini$genotypes[[k]]$region_id <- k
    mini$scenarios[[1]]$regions[[k]]$region_id <- k
  }
  for (nm in c("p", "p_actual", "vcf_files_used")) {
    values <- sim$params[[nm]]
    if (!is.null(values)) {
      if (length(values) == 1L && nm != "vcf_files_used") values <- rep(values, n_reg)
      if (length(values) != n_reg)
        stop("Per-region metadata ", nm, " does not match the cached loci.", call. = FALSE)
      mini$params[[nm]] <- values[keep]
    }
  }
  mini$params$n_regions <- length(keep)
  mini$params$n_scenarios <- 1L
  mini$params$source_scenario_index <- scenario
  mini$params$source_region_indices <- keep
  mini
}

# No pooled-AP fallback: old evaluators used a different estimand and must be
# upgraded rather than silently changing this study's objective.
tuning_attempted_ap <- function(global) {
  if (is.null(global$ap_macro_attempted))
    stop("Evaluator does not provide ap_macro_attempted; upgrade fmbenchmark.", call. = FALSE)
  value <- global$ap_macro_attempted
  if (length(value) != 1L || (!is.na(value) && (!is.finite(value) || value < 0 || value > 1)))
    stop("Invalid all-attempt locus-macro average precision.", call. = FALSE)
  value
}
