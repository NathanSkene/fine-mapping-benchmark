#!/usr/bin/env Rscript
# Pure contracts and a tiny sealed, synthetic metadata fixture. No native
# execution, full-campaign fits, genotype matrix or scientific truth is read.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
source(file.path(repo, "scripts/reporting/resource_cost_report.R"))
if (dir.exists(file.path(repo, ".local-lib"))) .libPaths(c(file.path(repo, ".local-lib"), .libPaths()))
checks <- 0L
ok <- function(x, message) {
  checks <<- checks + 1L
  if (!isTRUE(x)) stop("FAILED: ", message)
}
reject <- function(expr, pattern, message) {
  error <- tryCatch({ force(expr); NULL }, error = conditionMessage)
  if (!is.character(error) || !grepl(pattern, error))
    cat("Expected rejection pattern: ", pattern, "; observed: ", error %||rc% "no rejection", "\n", sep = "")
  ok(is.character(error) && grepl(pattern, error), message)
}
near <- function(x, y) isTRUE(all.equal(x, y, tolerance = 1e-12, check.attributes = FALSE))

requested <- rep(c(500L, 750L, 1000L, 1500L, 2000L), each = 2L)
realized <- c(498L, 498L, 749L, 750L, 997L, 998L, 1496L, 1496L, 1998L, 1999L)
audit_fixture <- function(req = requested, actual = realized, n = 5000L, nref = NA_integer_, known = TRUE,
                          source_sha = if (known) .rc_known_filter_sha else strrep("a", 64L)) {
  basis <- if (known) .rc_filter_basis else "unverified_filter_origin_no_dimension_loss_permitted"
  list(schema = 1L, source_file = "R/simulate_genotypes.R", source_file_sha256 = source_sha,
    known_post_simulation_filter = known, design_semantics = "nominal_requested_variant_count",
    variant_universe_semantics = "retained_polymorphic_gwas_and_reference_if_present",
    post_simulation_filter_basis = basis,
    loci = data.frame(region_id = seq_along(req), requested_p = req, realized_p = actual,
      dropped_p = req - actual, retained_fraction = actual / req, dropped_fraction = (req - actual) / req,
      n_gwas = rep(n, length(req)), n_reference = rep(nref, length(req)),
      post_simulation_filter_basis = rep(basis, length(req)), stringsAsFactors = FALSE),
    warnings = data.frame(index = integer(), message = character(), classification = character()))
}
audit <- audit_fixture()
census <- .rc_dimension_census(audit, requested, 5000L, full_probe = TRUE)
ok(max(census$realized_p) == 1999L, "legitimate largest realized locus below2000 is accepted")
ok(identical(census$realized_p, realized), "the full realized ten-locus vector is retained")
ok(all(census$requested_p == requested), "nominal counts remain distinct")
ok(all(census$design_label == "nominal full-size filtered ten-locus design"), "filtered full-size scope is explicit")
ok(all(census$N_basis == census$N_basis[1L]) && grepl("actual simulated", census$N_basis[1L]), "N is actual simulated rows")
sizes <- .rc_size_summary(requested, realized)
ok(sizes$realized_p_vector == paste(realized, collapse = "|"), "size vector is exported")
ok(near(sizes$realized_to_requested_sum_p, sum(realized) / sum(requested)), "p sensitivity ratio")
ok(near(sizes$realized_to_requested_sum_p2, sum(realized^2) / sum(requested^2)), "p squared sensitivity ratio")
ok(near(sizes$realized_to_requested_sum_p3, sum(realized^3) / sum(requested^3)), "p cubed sensitivity ratio")
ok(grepl("not runtime", sizes$scaling_interpretation), "dimension sensitivity is not a cost bound")
for (field in c("requested_p", "realized_p", "dropped_p", "n_gwas", "n_reference", "retained_fraction", "dropped_fraction")) {
  bad <- audit; bad$loci[[field]][1L] <- 42
  reject(.rc_dimension_census(bad, requested, 5000L, full_probe = TRUE), "Contradictory|Invalid", paste("contradictory audit", field))
}
bad <- audit; bad$source_file_sha256 <- strrep("b", 64L); bad$known_post_simulation_filter <- FALSE
bad$post_simulation_filter_basis <- "unverified_filter_origin_no_dimension_loss_permitted"
bad$loci$post_simulation_filter_basis <- rep(bad$post_simulation_filter_basis, 10L)
reject(.rc_dimension_census(bad, requested, 5000L, full_probe = TRUE), "Unknown origin", "unknown filter cannot explain reduction")
for (warning in c("Region 1: requested 2000 SNPs but only 1500 variants available after MAF filtering.",
                  "Capping p at 500 for regions where p > 500")) {
  bad <- audit; bad$warnings <- data.frame(index = 1L, message = warning,
                          classification = "other_original_warning_preserved")
  reject(.rc_dimension_census(bad, requested, 5000L, full_probe = TRUE), "undersupply", "known filter does not excuse undersupply")
}
bad <- audit; bad$loci$region_id[2L] <- 1L
reject(.rc_dimension_census(bad, requested, 5000L), "Contradictory", "duplicate locus cannot pass")
bad <- audit; bad$post_simulation_filter_basis <- "assumed"
reject(.rc_dimension_census(bad, requested, 5000L), "filter basis", "filter basis is literal")
reject(.rc_dimension_census(audit, requested, 1000L, full_probe = TRUE), "full probe", "N1000 cannot be relabelled as N5000 probe")

design_fixture <- function(n = 5000L, nref = NA_integer_, p = requested, S = 2L, phi = .1,
                           annotation = "binary", dataset = "test_dataset") {
  data.frame(family = "legacy_reconstruction", nominal_p_vector = paste(p, collapse = "|"),
    n = n, n_ref = nref, model = "sparse", annotation = annotation, enrichment = "1|2|10.8",
    annotation_correlation = NA_real_, p_causal = NA_real_, S = S, phi = phi,
    replicates = 1L, n_loci = length(p), dataset_id = dataset, stringsAsFactors = FALSE)
}
legacy_table <- function(design = design_fixture(), census = census, method = "abf", stage = "corrected_refit") {
  data.frame(dataset_id = design$dataset_id, family = design$family, stage = stage, method = method,
    region_id = census$region_id, nominal_p = census$requested_p, actual_p = census$realized_p,
    n_gwas = census$n_gwas, n_ref = census$n_reference, S = design$S, phi = design$phi,
    p_causal = design$p_causal, annotation_type = design$annotation, success = TRUE,
    runtime_seconds = seq_len(nrow(census)) / 10, error = "", stringsAsFactors = FALSE)
}
d <- design_fixture(); t <- legacy_table(d, census)
loci <- .rc_locus_census(t, d, census, "abf", "corrected_refit")
support <- .rc_execution_support(loci)
ok(support$new_inference_anchor_eligible, "successful filtered largest loci can provide a measured anchor")
ok(support$n_applicable_largest_loci == 2L && support$n_successful_largest_loci == 2L,
   "largest two nominal loci are both checked")
ok(near(loci$fit_runtime_seconds, t$runtime_seconds), "recorded fit times are unchanged")
for (field in c("actual_p", "nominal_p", "n_gwas", "n_ref", "S", "phi", "p_causal")) {
  bad <- t; bad[[field]][1L] <- 42
  reject(.rc_locus_census(bad, d, census, "abf", "corrected_refit"), "disagrees", paste("table contradicts", field))
}
bad <- t; bad$dataset_id[1L] <- "foreign"
reject(.rc_locus_census(bad, d, census, "abf", "corrected_refit"), "identity", "foreign table identity")
bad <- t; bad$success[10L] <- FALSE; bad$error[10L] <- "unexpected_failure"; bad$runtime_seconds[10L] <- .001
failed_loci <- .rc_locus_census(bad, d, census, "abf", "corrected_refit")
failed_support <- .rc_execution_support(failed_loci)
ok(!failed_support$new_inference_anchor_eligible, "a cheap largest-locus failure cannot clear support")
ok(failed_support$n_finite_failed_fit_runtime == 1L, "cheap failed timing is retained")
ok(failed_support$inference_support_reason == "largest_applicable_locus_failed", "largest-locus reason is explicit")
bad <- t; bad$success[1L] <- FALSE; bad$error[1L] <- "unexpected_failure"
ok(!.rc_execution_support(.rc_locus_census(bad, d, census, "abf", "corrected_refit"))$new_inference_anchor_eligible,
   "unexpected smaller-locus failure also blocks the complete-input anchor")
ok(!.rc_execution_support(loci, TRUE)$new_inference_anchor_eligible, "reused fits are not new measured inference")
none <- design_fixture(annotation = "none")
reject_table <- legacy_table(none, census, "funmap")
reject_table$success <- FALSE
reject_table$error <- "No annotations available for this region. Funmap requires functional annotations. Re-run simulation with annotations = 'binary' or annotations = 'continuous'."
reject_table$runtime_seconds <- .001
rejections <- .rc_locus_census(reject_table, none, census, "funmap", "corrected_refit")
reject_support <- .rc_execution_support(rejections)
ok(reject_support$n_expected_unsupported == 10L, "expected unsupported outcomes remain in the census")
ok(reject_support$n_finite_failed_fit_runtime == 10L, "all cheap reject times remain recorded")
ok(!reject_support$new_inference_anchor_eligible, "all-rejected tasks cannot provide inference feasibility")
ok(reject_support$inference_support_reason == "all_attempts_are_expected_unsupported_rejections", "unsupported scope is explicit")
bad <- reject_table; bad$error <- paste0(bad$error, " Additional unexpected failure.")
ok(!any(.rc_expected_rejection(bad, "funmap", "legacy_reconstruction", "none")), "error substrings cannot excuse an unexpected failure")

mechanisms <- c("matched_prior", "null", "single_sparse", "multi_sparse", "many_sparse", "opposing_effects",
 "sparse_inf", "ld_reference_small", "ld_reference_shift", "annotation_informative", "annotation_irrelevant",
 "annotation_reversal", "missing_causal", "perfect_ld", "allele_flip", "allele_recoding", "selected_null",
 "selected_signal", "meta_missingness", "meta_known_overlap", "multi_ancestry_shared", "multi_ancestry_heterogeneous")
challenge_design <- data.frame(family = "challenge", nominal_p_vector = "1000", n = 10000L, n_ref = NA_real_,
 model = "conditional_Gaussian", annotation = "mixed_mechanisms", enrichment = NA_character_,
 annotation_correlation = NA_real_, p_causal = NA_real_, S = NA_real_, phi = NA_real_,
 replicates = 10L, n_loci = 220L, dataset_id = "stress_panel", stringsAsFactors = FALSE)
challenge_table <- data.frame(dataset_id = "stress_panel", family = "challenge", stage = "challenge", method = "abf",
  task_id = sprintf("stress_%03d", 1:220), scenario = rep(mechanisms, each = 10L), p = 1000L, n = 10000L,
  actual_p = 1000L - as.integer(rep(mechanisms, each = 10L) == "missing_causal"), success = TRUE,
  runtime_seconds = .1, error = "", stringsAsFactors = FALSE)
unsupported <- challenge_table$scenario %in% mechanisms[19:22]
challenge_table$success[unsupported] <- FALSE
challenge_table$error[unsupported] <- "unsupported_input: regional single-cohort scalar-N method"
stress <- .rc_locus_census(challenge_table, challenge_design, NULL, "abf", "challenge")
stress_support <- .rc_execution_support(stress)
ok(nrow(stress) == 220L, "all220 Gaussian attempts retained")
ok(sum(stress$dropped_p) == 10L, "ten intentional missing-candidate reductions retained")
ok(sum(stress$classification == "expected_unsupported_rejection") == 40L, "Gaussian scalar-N unsupported classes explicit")
ok(stress_support$new_inference_anchor_eligible, "supported applicable Gaussian batch with declared unsupported mechanisms")
ok(all(grepl("design parameter", stress$N_basis)), "Gaussian N is not called individual rows")
small_reference <- challenge_table; small_reference$method <- "exact"; small_reference$p <- 200L
small_reference$actual_p <- 200L - as.integer(small_reference$scenario == "missing_causal")
small_reference$success <- small_reference$scenario != "meta_missingness"
small_reference$error <- ifelse(small_reference$success, "", "Variant-wise n requires sample-overlap metadata for this adapter")
small_design <- challenge_design; small_design$nominal_p_vector <- "200"
for (reference in c("exact", "forced_one")) {
  small_reference$method <- reference
  ref_loci <- .rc_locus_census(small_reference, small_design, NULL, reference, "challenge")
  ok(sum(ref_loci$classification == "expected_unsupported_rejection") == 10L,
     paste(reference, "missing overlap metadata is a declared absence"))
  bad <- small_reference; bad$scenario[bad$scenario == "meta_missingness"] <- "meta_known_overlap"
  ok(!any(.rc_expected_rejection(bad, reference, "challenge", "mixed_mechanisms")),
     paste(reference, "overlap absence is not inferred for other mechanisms"))
}
bad <- challenge_table; bad$method <- "exact"
reject(.rc_locus_census(bad, challenge_design, NULL, "exact", "challenge"), "p<=200", "small-p reference cannot claim measured large-p stress")
bad <- challenge_table; bad$actual_p[1L] <- 999L
reject(.rc_locus_census(bad, challenge_design, NULL, "abf", "challenge"), "Unexplained", "unplanned candidate reduction rejected")
bad <- challenge_table; bad$actual_p[121L] <- 1000L
reject(.rc_locus_census(bad, challenge_design, NULL, "abf", "challenge"), "Unexplained", "missing-cause reduction cannot be hidden")
bad <- challenge_table[-1L, ]
reject(.rc_locus_census(bad, challenge_design, NULL, "abf", "challenge"), "Missing", "incomplete stress batch rejected")
bad <- challenge_table; bad$scenario[1L] <- "unknown"
reject(.rc_locus_census(bad, challenge_design, NULL, "abf", "challenge"), "Incomplete", "unknown mechanism rejected")
bad <- challenge_table; bad$task_id[2L] <- bad$task_id[1L]
reject(.rc_locus_census(bad, challenge_design, NULL, "abf", "challenge"), "identity", "duplicate stress identity rejected")
bad <- challenge_table; bad$n <- 500L
reject(.rc_locus_census(bad, challenge_design, NULL, "abf", "challenge"), "dimensions", "small integration cannot masquerade as stress")

resource <- list(schema = 1L, n_new = 1L, n_datasets = 1L, script_elapsed_seconds = 9,
 prelude_seconds = 2, input_io_seconds = 1, compaction_seconds = 1, data_processing_seconds = 5, overhead_seconds = 4)
timing <- .rc_worker_timing(resource, 1L, 1L, 1)
ok(timing$complete_new_worker_timing, "complete new worker timing qualifies")
ok(identical(timing$script_elapsed_seconds, 9), "actual script timing unchanged")
ok(!.rc_worker_timing(NULL, 1L, 1L, 1)$complete_new_worker_timing, "missing timing unresolved")
partial <- resource; partial$n_new <- 0L
ok(!.rc_worker_timing(partial, 1L, 0L, 1)$complete_new_worker_timing, "resumed worker is not a fresh complete timing anchor")
bad <- resource; bad$overhead_seconds <- 10
reject(.rc_worker_timing(bad, 1L, 1L, 1), "Contradictory", "timing components must reconcile")
bad <- resource; bad$input_io_seconds <- -1
reject(.rc_worker_timing(bad, 1L, 1L, 1), "Invalid", "negative timing rejected")
reject(.rc_worker_timing(resource, 1L, 1L, 6), "Contradictory", "dataset elapsed cannot exceed whole data-processing measurement")

tail_p <- rep(c(100L, 200L, 400L, 500L, 1000L), each = 2L)
tail_d <- do.call(rbind, lapply(c(1L, 2L, 3L, 5L, 10L), function(S) design_fixture(
  n = 1000L, nref = 200L, p = tail_p, S = S, phi = .4, dataset = paste0("tail_S", S))))
batch <- .rc_batch_design(tail_d, 1:5, "corrected_refit", "abf")
carma <- .rc_batch_design(tail_d, 5L, "corrected_refit", "carma")
ok(batch$n == 1000L && batch$n_ref == 200L, "N1000 finite-ref tail accepted as its own stratum")
ok(batch$S_pack == "1|2|3|5|10" && batch$dataset_records_per_task == 5L, "packed tail S values retained")
ok(carma$S_pack == "10" && carma$dataset_records_per_task == 1L, "CARMA S10-only distinct worker")
ref500 <- tail_d; ref500$n_ref <- 500L
in_sample <- tail_d; in_sample$n_ref <- NA_integer_
ok(.rc_batch_design(ref500, 1:5, "corrected_refit", "abf")$workload_key != batch$workload_key, "ref500 andref200 are distinct")
ok(.rc_batch_design(in_sample, 1:5, "corrected_refit", "abf")$workload_key != batch$workload_key, "in-sample andfinite-ref are distinct")
ok(.rc_batch_design(tail_d, 1:5, "original_refit", "abf")$workload_key != batch$workload_key, "stages cannot pool cost anchors")
bad <- tail_d; bad$n_ref[1L] <- 500L
reject(.rc_batch_design(bad, 1:5, "corrected_refit", "abf"), "Mixed-design", "a mixed-reference batch cannot hide heterogeneity")
workload <- batch; workload$planned_task_units <- 10L; workload$planned_dataset_records <- 50L; workload$planned_locus_attempts <- 500L
measured <- cbind(data.frame(task_complete = TRUE, new_inference_anchor_eligible = TRUE,
  complete_new_worker_timing = TRUE, script_elapsed_seconds = 36, primary_retention_bytes = 1000,
  inference_support_reason = "observed_applicable_execution_supported"), batch)
envelope <- .rc_envelope(workload, measured)
ok(near(envelope$conditional_R_worker_hours_median, .1), "exact-batch conditional R envelope has hand-calculated cost")
ok(!envelope$feasibility_cleared, "even measured anchors are not feasibility guarantees")
cheap <- measured; cheap$new_inference_anchor_eligible <- FALSE; cheap$script_elapsed_seconds <- .001
cheap$inference_support_reason <- "all_attempts_are_expected_unsupported_rejections"
ok(is.na(.rc_envelope(workload, cheap)$conditional_R_worker_hours_median), "cheap rejection cannot lower an inference envelope")
empty <- measured[0L, ]
ok(.rc_envelope(workload, empty)$unresolved == "unmeasured_exact_worker_design", "missing stress/other design is explicitly unresolved")
ok(is.na(.rc_envelope(workload, empty)$conditional_R_worker_hours_median), "missing evidence does not become zero")

helper <- .rc_load_helpers(repo)
before_source <- helper$.campaign_snapshot(repo)
tmp <- tempfile("resource-report-fixture-"); dir.create(tmp)
small_metadata <- file.path(tmp, "small-metadata.json"); writeLines("{}", small_metadata)
ok(identical(.rc_bounded_file(small_metadata, 3L, "Tail selection"), normalizePath(small_metadata)),
   "metadata read bound accepts a bounded existing file")
reject(.rc_bounded_file(small_metadata, 2L, "Tail selection"), "exceeds explicit read bound",
   "tail metadata size refused before parsing")
reject(.rc_bounded_file(tmp, 10L, "Tail selection"), "existing file", "metadata read route cannot be a directory")
make_plan <- function(root, full = FALSE) {
  dir.create(root)
  root <- normalizePath(root)
  grid <- data.frame(n_regions = 10L, n = 5000L, p_values = paste(requested, collapse = "|"), n_iter = 1L,
    S_values = "2", phi_values = "0.1", model = "sparse", annotation_type = "binary",
    enrichment_values = "1|2|10.8", annotation_correlation = NA_real_, p_causal = NA_real_, n_ref = NA_integer_)
  if (full) {
    tail_grid <- grid; tail_grid$n <- 1000L; tail_grid$p_values <- paste(tail_p, collapse = "|")
    tail_grid$S_values <- "1|2|3|5|10"; tail_grid$phi_values <- "0.4"; tail_grid$n_ref <- 200L
    grid <- rbind(grid, tail_grid)
  }
  generation <- data.frame(generation_id = 1L, family = "legacy_reconstruction", row = 1L,
    p = NA_integer_, n = 5000L, replicates = NA_integer_, seed = 1001L)
  if (full) generation <- rbind(generation, data.frame(generation_id = 2L,
    family = "legacy_reconstruction", row = 2L, p = NA_integer_, n = 1000L, replicates = NA_integer_, seed = 1002L))
  gi <- nrow(generation) + 1L
  generation <- rbind(generation, data.frame(generation_id = gi, family = "challenge", row = 1L,
    p = 1000L, n = 10000L, replicates = 10L, seed = 2001L))
  if (full) generation <- rbind(generation, data.frame(generation_id = gi + 1L, family = "challenge", row = 1L,
    p = 200L, n = 1000L, replicates = 10L, seed = 2002L))
  dataset <- function(id, family, row, scenario, gen) data.frame(dataset_id = id, family = family, row = row,
    scenario = scenario, sim_file = file.path(root, "datasets", sprintf("generation_%05d", gen),
      if (family == "challenge") "panel.rds" else "simulation.rds"), sim_sha256 = "generated",
    result_files = "", result_sha256 = "", genotypes_sha256 = "generated", original_label = "fixture")
  datasets <- dataset("fixture_legacy", "legacy_reconstruction", 1L, 1L, 1L)
  if (full) for (Sindex in 1:5) datasets <- rbind(datasets, dataset(paste0("fixture_tail", Sindex),
    "legacy_reconstruction", 2L, Sindex, 2L))
  datasets <- rbind(datasets, dataset("stress_panel", "challenge", gi, NA_integer_, gi))
  if (full) datasets <- rbind(datasets, dataset("unmeasured_exact", "challenge", gi + 1L, NA_integer_, gi + 1L))
  task_rows <- list(list(ids = "fixture_legacy", stage = "corrected_refit", method = "abf"))
  if (full) task_rows <- c(task_rows, list(list(ids = paste0("fixture_tail", 1:5, collapse = "|"), stage = "corrected_refit", method = "abf"),
    list(ids = "fixture_tail5", stage = "corrected_refit", method = "carma")))
  task_rows <- c(task_rows, list(list(ids = "stress_panel", stage = "challenge", method = "abf")))
  if (full) task_rows <- c(task_rows, list(list(ids = "unmeasured_exact", stage = "challenge", method = "exact")))
  tasks <- do.call(rbind, lapply(seq_along(task_rows), function(i) data.frame(task_id = sprintf("task_%07d", i),
    index = i, dataset_ids = task_rows[[i]]$ids, stage = task_rows[[i]]$stage, method = task_rows[[i]]$method, seed = 3000L + i)))
  plan <- list(profile = if (full) "full" else "resource_probe", original_source = repo,
    corrected_source = helper$.campaign_snapshot(repo), original_source_hashes = helper$.campaign_snapshot(repo),
    runtime_fingerprint = helper$.campaign_runtime_fingerprint(list()), method_args = list(),
    tasks = tasks, datasets = datasets, generation = generation, grid = grid,
    challenge = list(variants = if (full) c(200L, 1000L) else 1000L, n = if (full) c(1000L, 10000L) else 10000L,
      mixed_mechanism_batch_replicates = 10L))
  plan$plan_sha256 <- helper$.campaign_object_sha(plan)
  saveRDS(plan, file.path(root, "plan.rds"), version = 2)
  write.csv(grid, file.path(root, "legacy-grid.csv"), row.names = FALSE)
  jsonlite::write_json(list(schema = 1L, fixture = TRUE), file.path(root, "dispatch-manifest.json"), auto_unbox = TRUE)
  plan
}
probe_root <- file.path(tmp, "probe"); full_root <- file.path(tmp, "full")
probe <- make_plan(probe_root); full <- make_plan(full_root, TRUE)
helper$.campaign_assert_plan(probe, repo); helper$.campaign_assert_plan(full, repo)
full_designs <- .rc_scenario_designs(full)
cache_workload <- .rc_workload_census(full, full_designs)
ok(sum(cache_workload$planned_task_units) == nrow(full$tasks), "cached census retains every planned worker")
ok(sum(cache_workload$planned_dataset_records) == 9L, "cached packing count is exact")
large_plan <- full
large_plan$tasks <- full$tasks[rep(seq_len(nrow(full$tasks)), length.out = 649220L), , drop = FALSE]
large_plan$tasks$task_id <- sprintf("task_%07d", seq_len(nrow(large_plan$tasks)))
large_plan$tasks$index <- seq_len(nrow(large_plan$tasks))
began <- proc.time()[["elapsed"]]
large_workload <- .rc_workload_census(large_plan, full_designs)
large_seconds <- proc.time()[["elapsed"]] - began
ok(sum(large_workload$planned_task_units) == 649220L, "649220-worker metadata census preserves all workers")
ok(nrow(large_workload) == 5L, "repeated metadata is cached by exact batch design")
ok(large_seconds < 30, "cached repeated 649220-worker metadata fixture completes within30seconds")
cat(sprintf("Cached repeated metadata fixture: 649220 tasks, %d design groups, %.3f seconds.\n", nrow(large_workload), large_seconds))
rm(large_plan, large_workload)
seal_generation <- function(plan, root, gi) {
  g <- plan$generation[gi, , drop = FALSE]; directory <- file.path(root, "datasets", sprintf("generation_%05d", gi))
  dir.create(directory, recursive = TRUE)
  artifact <- plan$datasets$sim_file[match(gi, ifelse(plan$datasets$family == "challenge", plan$datasets$row,
    vapply(plan$datasets$sim_file, function(p) as.integer(sub("generation_", "", basename(dirname(p)))), 1L)))]
  saveRDS(list(fixture = "metadata only; no scientific truth or genotype rows"), artifact, version = 2)
  dimensions <- NULL
  if (g$family != "challenge") {
    gr <- plan$grid[g$row, , drop = FALSE]; p <- .rc_values(gr$p_values, "fixture p", TRUE)
    dimensions <- audit_fixture(p, p, gr$n, gr$n_ref, known = FALSE,
      source_sha = helper$.challenge_sha256(file.path(repo, "R/simulate_genotypes.R")))
  }
  z <- list(terminal = TRUE, plan_sha256 = plan$plan_sha256, task = g,
    started = as.POSIXct("2026-10-11 12:00:00", tz = "UTC"), completed = as.POSIXct("2026-10-11 12:00:01", tz = "UTC"),
    generation_dimension_audit = dimensions, artifacts = NULL)
  artifacts <- artifact
  if (g$family == "challenge") for (scope in c("public", "assessor")) {
    release <- file.path(directory, scope); dir.create(release)
    files <- if (scope == "public") c("inputs.rds", "metadata.rds") else c("truth.rds", "generator.rds")
    for (file in files) writeLines(paste("opaque metadata fixture", scope, file), file.path(release, file))
    manifest <- data.frame(file = files, sha256 = vapply(file.path(release, files), .rc_file_sha, ""))
    write.csv(manifest, file.path(release, "manifest.csv"), row.names = FALSE)
    artifacts <- c(artifacts, file.path(release, "manifest.csv"))
  }
  z$artifacts <- helper$.campaign_artifact_manifest(artifacts)
  saveRDS(z, file.path(directory, "generated.rds"), version = 2)
  helper$.campaign_receipt(directory, file.path(directory, "generated.rds"), plan$plan_sha256, paste0("generation_", gi))
}
seal_task <- function(plan, root, ti) {
  task <- plan$tasks[ti, , drop = FALSE]; ids <- strsplit(task$dataset_ids, "|", fixed = TRUE)[[1L]]
  directory <- file.path(root, "tasks", task$task_id); dir.create(directory, recursive = TRUE)
  designs <- .rc_scenario_designs(plan)
  for (id in ids) {
    dataset <- plan$datasets[match(id, plan$datasets$dataset_id), , drop = FALSE]
    design <- designs[match(id, designs$dataset_id), , drop = FALSE]
    if (design$family == "challenge") {
      table <- challenge_table; table$method <- task$method
    } else {
      req <- .rc_values(design$nominal_p_vector, "fixture p", TRUE)
      a <- .rc_dimension_census(audit_fixture(req, req, design$n, design$n_ref, FALSE), req, design$n, design$n_ref)
      table <- legacy_table(design, a, task$method, task$stage)
    }
    z <- list(terminal = TRUE, plan_sha256 = plan$plan_sha256, task = list(task = task, dataset = dataset),
      method_task = task, metadata = as.list(dataset), reuse = list(reused = FALSE), table = table, elapsed_seconds = 1)
    helper$.campaign_start_attempt(directory, plan, task, dataset)
    file <- file.path(directory, paste0(id, ".rds")); saveRDS(z, file, version = 2)
    helper$.campaign_receipt(directory, file, plan$plan_sha256, paste(task$task_id, id, sep = "."), paste0(id, ".terminal.json"))
  }
  helper$.campaign_compact_task(directory, plan, task, list(n_new = length(ids)))
  receipt <- file.path(directory, "terminal.json"); v <- jsonlite::fromJSON(receipt)
  timing <- resource; timing$n_datasets <- timing$n_new <- length(ids)
  timing$data_processing_seconds <- 5 * length(ids); timing$script_elapsed_seconds <- timing$data_processing_seconds + timing$overhead_seconds
  v$resource_measurement <- timing
  jsonlite::write_json(v, receipt, auto_unbox = TRUE, pretty = TRUE, digits = NA)
}
for (gi in 1:2) seal_generation(probe, probe_root, gi)
for (ti in 1:2) seal_task(probe, probe_root, ti)
seal_generation(full, full_root, 2L)
for (ti in 2:3) seal_task(full, full_root, ti)
selection <- file.path(tmp, "tail-selection.json")
jsonlite::write_json(list(schema = 1L, source_commit = "0e4a7ead774b2e8d3f8ab2f48bc61eadb31c6e91",
  plan_sha256 = full$plan_sha256, plan_rds_sha256 = helper$.challenge_sha256(file.path(full_root, "plan.rds")),
  projection_manifest_sha256 = helper$.challenge_sha256(file.path(full_root, "dispatch-manifest.json")),
  legacy_grid_csv_sha256 = helper$.challenge_sha256(file.path(full_root, "legacy-grid.csv")),
  original_scenario_source_sha256 = helper$.challenge_sha256(file.path(repo, "R/run_simulation.R")),
  method_task_indices = c(2L, 3L), generation_indices = 2L), selection, auto_unbox = TRUE, pretty = TRUE)
args <- c("--source", repo, "--probe-plan", file.path(probe_root, "plan.rds"), "--full-plan", file.path(full_root, "plan.rds"),
          "--tail-selections", selection, "--out", file.path(tmp, "resource-cost-fixture-v1"))
input_paths <- list.files(tmp, recursive = TRUE, full.names = TRUE)
input_paths <- input_paths[!dir.exists(input_paths)]
before_inputs <- helper$.campaign_artifact_manifest(input_paths)
report <- .rc_main(args, file.path(repo, "scripts/reporting/resource_cost_report.R"))
ok(report$complete && report$source_plans_and_receipts_unchanged, "sealed synthetic metadata fixture completes")
ok(report$measured_task_rows == 4L, "two probe andtwo tail tasks retained")
ok(report$planned_method_task_units == 5L, "every full planned task counted")
ok(report$planned_dataset_records == 9L, "packed worker dataset census exact")
ok(report$planned_locus_attempts == 510L, "all full planned attempts include unmeasured exact scope")
ok(!report$full_R_cost_established && !report$inference_feasibility_cleared, "report cannot claim full cost or feasibility")
ok(report$PBS_actual_attempts_observed == 0L, "absence of PBS evidence explicit")
ok(identical(before_source, helper$.campaign_snapshot(repo)), "source snapshot bytes unchanged")
ok(identical(before_inputs, helper$.campaign_artifact_manifest(input_paths)), "all fixture plans/receipts/checkpoints unchanged")
cost <- read.csv(file.path(tmp, "resource-cost-fixture-v1", "conditional-cost-envelope.csv"), stringsAsFactors = FALSE)
ok(nrow(cost) == 5L, "exact worker designs are distinct")
ok(cost$unresolved[cost$method == "exact"] == "unmeasured_exact_worker_design", "unmeasured exact reference unresolved")
ok(is.na(cost$conditional_R_worker_hours_median[cost$method == "exact"]), "unmeasured exact cost stays NA")
ok(cost$S_pack[cost$method == "carma"] == "10", "CARMA output does not inherit pooled S pack")
ok(any(cost$n == 1000L & cost$n_ref == 200L, na.rm = TRUE), "N1000ref200 strata survive CLI")
attempts <- read.csv(file.path(tmp, "resource-cost-fixture-v1", "started-attempts.csv"), stringsAsFactors = FALSE)
ok(nrow(attempts) == 8L, "all started attempts retained including packed datasets")
ok(!any(attempts$duration_available), "start ledger does not invent duration")
reject(.rc_main(args), "new versioned", "immutable output directory cannot be overwritten")

# Missing stress evidence remains valid unresolved reporting, not a fake zero
# or an installation-only success masquerading as full stress feasibility.
missing_probe_root <- file.path(tmp, "missing-probe"); missing_probe <- make_plan(missing_probe_root)
seal_generation(missing_probe, missing_probe_root, 1L); seal_task(missing_probe, missing_probe_root, 1L)
missing_args <- args; missing_args[match("--probe-plan", args) + 1L] <- file.path(missing_probe_root, "plan.rds")
missing_args[match("--out", args) + 1L] <- file.path(tmp, "resource-cost-missing-stress-v1")
missing_report <- .rc_main(missing_args, file.path(repo, "scripts/reporting/resource_cost_report.R"))
missing_cost <- read.csv(file.path(tmp, "resource-cost-missing-stress-v1", "conditional-cost-envelope.csv"))
ok(is.na(missing_cost$conditional_R_worker_hours_median[missing_cost$family == "challenge" & missing_cost$method == "abf"]),
   "zero completed stress samples cannot produce a Gaussian cost")
ok(!missing_report$inference_feasibility_cleared, "missing stress cannot clear feasibility")
# A dependency fingerprint alone does not bind scalar search/training/prior
# settings. Use a valid resealed plan with unchanged runtime to exercise that.
altered_full <- full; altered_full$method_args <- list(abf = list(prior_variance = .1))
ok(identical(helper$.campaign_runtime_fingerprint(altered_full$method_args), full$runtime_fingerprint),
   "non-path setting mutation really preserves the dependency fingerprint")
altered_full$plan_sha256 <- NULL; altered_full$plan_sha256 <- helper$.campaign_object_sha(altered_full)
altered_root <- file.path(tmp, "altered-full"); dir.create(altered_root)
altered_path <- file.path(altered_root, "plan.rds"); saveRDS(altered_full, altered_path, version = 2)
helper$.campaign_assert_plan(altered_full, repo)
settings_args <- args; settings_args[match("--full-plan", args) + 1L] <- altered_path
settings_args[match("--tail-selections", args) + 1L] <- ""
settings_args[match("--out", args) + 1L] <- file.path(tmp, "resource-cost-settings-mismatch-v1")
reject(.rc_main(settings_args), "method settings", "valid foreign method settings cannot inherit probe costs")
# An independently valid original-source lineage is also incompatible. Copy
# only snapshot-listed bytes to a disposable directory, change one comment,
# and reseal its plan so this tests anchor compatibility rather than a broken
# source snapshot. Nothing in the actual inference source is edited or run.
old_copy <- file.path(tmp, "different-original-source"); dir.create(old_copy)
for (relative in full$original_source_hashes$path) {
  destination <- file.path(old_copy, relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(repo, relative), destination)) stop("Could not prepare original-source fixture")
}
old_copy <- normalizePath(old_copy)
cat("\n# disposable original-lineage fixture\n", file = file.path(old_copy, "R/simulate_genotypes.R"), append = TRUE)
lineage_full <- full; lineage_full$original_source <- old_copy
lineage_full$original_source_hashes <- helper$.campaign_snapshot(old_copy)
lineage_full$plan_sha256 <- NULL; lineage_full$plan_sha256 <- helper$.campaign_object_sha(lineage_full)
lineage_root <- file.path(tmp, "different-original-full"); dir.create(lineage_root)
lineage_path <- file.path(lineage_root, "plan.rds"); saveRDS(lineage_full, lineage_path, version = 2)
helper$.campaign_assert_plan(lineage_full, repo)
ok(identical(lineage_full$runtime_fingerprint, full$runtime_fingerprint) &&
   !identical(lineage_full$original_source_hashes, full$original_source_hashes),
   "valid original-source fixture changes only byte lineage relevant to compatibility")
lineage_args <- settings_args; lineage_args[match("--full-plan", args) + 1L] <- lineage_path
lineage_args[match("--out", args) + 1L] <- file.path(tmp, "resource-cost-original-lineage-mismatch-v1")
reject(.rc_main(lineage_args), "original-generator byte lineage", "valid foreign original-source lineage cannot inherit probe costs")
empty_root <- file.path(tmp, "empty-probe"); empty_probe <- make_plan(empty_root)
empty_args <- args; empty_args[match("--probe-plan", args) + 1L] <- file.path(empty_root, "plan.rds")
empty_args[match("--tail-selections", args) + 1L] <- ""
empty_args[match("--out", args) + 1L] <- file.path(tmp, "resource-cost-all-unmeasured-v1")
empty_report <- .rc_main(empty_args, file.path(repo, "scripts/reporting/resource_cost_report.R"))
ok(is.na(empty_report$conditional_qualified_subset_R_hours_median), "zero qualified groups have unknown total cost, not zero")
ok(empty_report$qualified_exact_design_groups == 0L, "no qualified groups are manufactured")

# Opaque child payloads are only hashed/stat'ed; a plain-text assessor payload
# would fail if the utility attempted to deserialize scientific truth.
stress_dir <- file.path(probe_root, "datasets/generation_00002")
gen <- readRDS(file.path(stress_dir, "generated.rds"))
retention <- .rc_generation_retention(gen, stress_dir, "challenge", function(path) .rc_file_sha(path), helper)
expected_files <- file.path(stress_dir, c("panel.rds", "public/manifest.csv", "assessor/manifest.csv",
  "public/inputs.rds", "public/metadata.rds", "assessor/truth.rds", "assessor/generator.rds"))
ok(retention$retained_input_files == 7L, "all seven Gaussian input/release files counted")
ok(retention$retained_input_bytes == sum(file.info(expected_files)$size), "retained release payload bytes are exact")
child <- file.path(stress_dir, "assessor/truth.rds")
saved_child <- readBin(child, "raw", n = file.info(child)$size)
writeLines("changed opaque child", child)
reject(.rc_generation_retention(gen, stress_dir, "challenge", function(path) .rc_file_sha(path), helper),
  "artifact bytes changed", "changed assessor child refused without deserializing it")
writeBin(saved_child, child)
release_manifest <- file.path(stress_dir, "public/manifest.csv")
saved_manifest <- readBin(release_manifest, "raw", n = file.info(release_manifest)$size)
m <- read.csv(release_manifest); m$file[1L] <- "../escape.rds"
write.csv(m, release_manifest, row.names = FALSE)
reject(.rc_generation_retention(gen, stress_dir, "challenge", function(path) .rc_file_sha(path), helper),
  "artifact bytes changed", "changed release manifest cannot redirect payload reads")
writeBin(saved_manifest, release_manifest)

# A method checkpoint with no compatible generator receipt is a custody error,
# distinct from an as-yet unmeasured method in a valid plan.
gen_receipt <- file.path(stress_dir, "terminal.json")
saved_receipt <- readBin(gen_receipt, "raw", n = file.info(gen_receipt)$size)
unlink(gen_receipt)
no_gen_args <- args; no_gen_args[match("--out", args) + 1L] <- file.path(tmp, "resource-cost-missing-generator-v1")
reject(.rc_main(no_gen_args), "Gaussian fit lacks compatible sealed generation", "sealed stress task cannot qualify without input custody")
ok(!file.exists(file.path(tmp, "resource-cost-missing-generator-v1/report-manifest.json")), "missing generator cannot produce completed report")
writeBin(saved_receipt, gen_receipt)

# Strict raw PBS binding; no scheduler command is executed in this test.
script <- file.path(tmp, "job-fixture.sh")
script_lines <- c("#!/bin/bash", "#PBS -N fm_fixture", vapply(2:3, function(i) paste("Rscript --vanilla",
  file.path(repo, "scripts/campaign/run_task.R"), "--plan", file.path(full_root, "plan.rds"), "--index", i), ""))
writeLines(script_lines, script)
qstat <- file.path(tmp, "qstat-fixture.json")
jobs <- stats::setNames(list(list(Job_Name = "fm_fixture", Submit_arguments = script, job_state = "F", Exit_status = 1L,
  resources_used = list(walltime = "00:01:01", cput = "00:00:50", mem = "2048kb"))), "123.pbs-7")
jsonlite::write_json(list(Jobs = jobs), qstat, auto_unbox = TRUE, pretty = TRUE)
pbs_path <- file.path(tmp, "pbs-bindings.json")
pbs_mapping <- list(schema = 1L, kind = "fmbenchmark-resource-pbs-bindings-v1", attempts = list(list(
 plan_sha256 = full$plan_sha256, kind = "methods", indices = c(2L, 3L), job_id = "123.pbs-7",
 script_path = script, script_sha256 = .rc_file_sha(script), qstat_path = qstat, qstat_sha256 = .rc_file_sha(qstat))))
jsonlite::write_json(pbs_mapping, pbs_path, auto_unbox = TRUE, pretty = TRUE)
plans <- stats::setNames(list(list(plan = full, path = file.path(full_root, "plan.rds"))), full$plan_sha256)
track <- function(path) invisible(.rc_file_sha(path))
observed <- .rc_pbs_observations(pbs_path, plans, repo, track)
ok(nrow(observed) == 1L && observed$exit_status == 1L, "actual failed PBS attempt retained")
ok(observed$wall_seconds == 61 && observed$cpu_seconds == 50, "PBS times read in actual units")
ok(observed$peak_memory_bytes == 2048 * 1024, "PBS kb converted explicitly")
ok(grepl("whole PBS job", observed$resource_scope), "batch resources not charged per method")
for (prefix in c("# ", "echo ", "true; ")) {
  changed <- script_lines; changed[3:4] <- paste0(prefix, changed[3:4]); writeLines(changed, script)
  bad <- pbs_mapping; bad$attempts[[1L]]$script_sha256 <- .rc_file_sha(script)
  jsonlite::write_json(bad, pbs_path, auto_unbox = TRUE)
  reject(.rc_pbs_observations(pbs_path, plans, repo, track), "exact source", paste("PBS comment/echo/compound prefix rejected", prefix))
}
writeLines(script_lines, script)
bad <- pbs_mapping; bad$attempts[[1L]]$indices <- 2L
jsonlite::write_json(bad, pbs_path, auto_unbox = TRUE)
reject(.rc_pbs_observations(pbs_path, plans, repo, track), "exact source", "batch-to-unit mismatch rejected")
bad <- pbs_mapping; bad$attempts[[1L]]$script_sha256 <- strrep("c", 64L)
jsonlite::write_json(bad, pbs_path, auto_unbox = TRUE)
reject(.rc_pbs_observations(pbs_path, plans, repo, track), "Changed PBS", "changed submitted script rejected")
bad <- pbs_mapping; bad$attempts <- c(bad$attempts, bad$attempts)
jsonlite::write_json(bad, pbs_path, auto_unbox = TRUE)
reject(.rc_pbs_observations(pbs_path, plans, repo, track), "Duplicate", "PBS attempt cannot be double counted")
bad_jobs <- jobs; bad_jobs[[1L]]$job_state <- "R"
jsonlite::write_json(list(Jobs = bad_jobs), qstat, auto_unbox = TRUE)
bad <- pbs_mapping; bad$attempts[[1L]]$qstat_sha256 <- .rc_file_sha(qstat)
jsonlite::write_json(bad, pbs_path, auto_unbox = TRUE)
reject(.rc_pbs_observations(pbs_path, plans, repo, track), "terminal", "running PBS job is not terminal measurement")

# Changed receipt/source custody is refused without touching real source.
receipt <- file.path(probe_root, "tasks/task_0000001/terminal.json")
v <- jsonlite::fromJSON(receipt); v$checkpoint_sha256 <- strrep("d", 64L)
jsonlite::write_json(v, receipt, auto_unbox = TRUE)
bad_args <- args; bad_args[match("--out", args) + 1L] <- file.path(tmp, "resource-cost-changed-receipt-v1")
reject(.rc_main(bad_args), "checkpoint bytes changed", "changed receipt cannot be blessed by reporter")
ok(!file.exists(file.path(tmp, "resource-cost-changed-receipt-v1/report-manifest.json")), "failed custody never emits a completion manifest")
source_copy <- file.path(tmp, "tampered-source"); dir.create(file.path(source_copy, "R"), recursive = TRUE)
dir.create(file.path(source_copy, "scripts/campaign"), recursive = TRUE)
invisible(file.copy(file.path(repo, "scripts/campaign/common.R"), file.path(source_copy, "scripts/campaign/common.R")))
invisible(file.copy(file.path(repo, "R/campaign.R"), file.path(source_copy, "R/campaign.R")))
expected <- data.frame(path = "R/campaign.R", sha256 = .rc_file_sha(file.path(repo, "R/campaign.R")))
cat("\n# fixture tamper\n", file = file.path(source_copy, "R/campaign.R"), append = TRUE)
reject(.rc_load_helpers(source_copy, expected), "before helper loading", "helper source bytes checked before any source execution")
invisible(file.copy(file.path(repo, "R/campaign.R"), file.path(source_copy, "R/campaign.R"), overwrite = TRUE))
marker <- file.path(tmp, "extra-helper-executed")
writeLines(paste0("writeLines('executed', ", deparse(marker), ")"), file.path(source_copy, "R/zzz_unlisted.R"))
reject(.rc_load_helpers(source_copy, expected), "Unlisted.*before helper loading", "unlisted helper is refused before source execution")
ok(!file.exists(marker), "unlisted source file never executed")
ok(identical(before_source, helper$.campaign_snapshot(repo)), "real scientific source still unchanged after custody failures")
unlink(tmp, recursive = TRUE)
cat(sprintf("PASS: %d resource reporting assertions; no inference or scheduler calls.\n", checks))
