#!/usr/bin/env Rscript
# A design-only prefix for a full immutable campaign. No predictions are read.
.campaign_priority_design <- function(plan) {
  if (!identical(plan$profile, "full")) stop("Priority audit requires a full scientific plan")
  ds <- plan$datasets; tasks <- plan$tasks; gen <- plan$generation; grid <- plan$grid
  if (anyDuplicated(tasks$index) || anyDuplicated(gen$generation_id) ||
      anyDuplicated(ds$dataset_id)) stop("Duplicate frozen design identities")
  if (any(ds$family != "challenge" & ds$family != "legacy_reconstruction"))
    stop("This priority design requires reconstructed grid/scenario identities, not inferred archives")
  latest <- which(!startsWith(grid$label, "iter003_control_"))
  if (length(latest) != 45L || any(grid$n_regions != 10L))
    stop("Expected the preserved 45-row original design and ten whole loci per row")
  near <- function(x, values) vapply(x, function(v)
    is.finite(v) && any(abs(v - values) < 1e-10), TRUE)
  controls <- which(startsWith(grid$label, "iter003_control_") &
      (grid$model == "sparse" | (grid$model == "sparse_inf" & near(grid$p_causal, c(.5, .9)))) &
      (grid$annotation_type == "none" |
        (grid$annotation_type %in% c("binary", "continuous") & near(grid$enrichment_fold, 5.4))))
  if (any(startsWith(grid$label, "iter003_control_")) && length(controls) != 27L)
    stop("Expected 27 fixed control audit rows (three backgrounds, three annotations, three LD regimes)")
  audit_rows <- c(latest, controls)
  design <- vector("list", nrow(grid))
  for (r in seq_len(nrow(grid))) {
    g <- grid[r, ]
    split_num <- function(s) as.numeric(strsplit(s, "|", fixed = TRUE)[[1L]])
    if (!identical(split_num(g$S_values), c(1, 2, 3, 5, 10)) || g$n_iter < 2L)
      stop("Unrecognized S/iteration design in row ", r)
    # Upstream run_simulation.R:205-220 uses expand.grid(S, phi, [scalar
    # p_causal], iter). S varies fastest; a scalar p_causal does not add cells.
    sc <- expand.grid(S = split_num(g$S_values), phi = split_num(g$phi_values),
                      iter = seq_len(g$n_iter), KEEP.OUT.ATTRS = FALSE)
    found <- ds[ds$family == "legacy_reconstruction" & ds$row == r, ]
    if (nrow(found) != nrow(sc) || anyDuplicated(found$scenario) ||
        !setequal(found$scenario, seq_len(nrow(sc)))) stop("Incomplete scenario mapping for row ", r)
    sc$dataset_id <- found$dataset_id[match(seq_len(nrow(sc)), found$scenario)]
    sc$row <- r; sc$scenario <- seq_len(nrow(sc))
    sc$audit_requested <- r %in% audit_rows & sc$iter %in% c(1L, 2L) &
      near(sc$phi, if (r %in% latest) c(.05, .2, .6) else c(.0075, .1, .4))
    if (r %in% audit_rows && sum(sc$audit_requested) != 30L)
      stop("Expected all five S levels, three fixed phi levels and two iterations in audit row ", r)
    design[[r]] <- sc
  }
  design <- do.call(rbind, design)
  requested <- design$dataset_id[design$audit_requested]
  # Parse each distinct task chunk once. Hundreds of thousands of method tasks
  # share these strings; repeated splitting would make planning needlessly slow.
  chunks <- unique(tasks$dataset_ids)
  chunk_ids <- strsplit(chunks, "|", fixed = TRUE)
  lengths <- lengths(chunk_ids)
  flat <- unlist(chunk_ids, use.names = FALSE)
  flat_chunks <- rep(seq_along(chunks), lengths)
  pos <- match(flat, ds$dataset_id)
  if (anyNA(pos)) stop("Task references an unknown frozen dataset")
  selected_chunks <- unique(flat_chunks[ds$dataset_id[pos] %in% requested])
  task_chunks <- match(tasks$dataset_ids, chunks)
  audit_indices <- which(task_chunks %in% selected_chunks &
      tasks$stage %in% c("original_refit", "corrected_refit", "corrected_resimulation"))
  challenge_indices <- which(tasks$stage == "challenge")
  if (!length(challenge_indices) || !length(audit_indices)) stop("Empty challenge or paired audit prefix")
  methods <- c(plan$methods, "exact", "forced_one")
  first_dataset <- match(vapply(chunk_ids, `[`, "", 1L), ds$dataset_id)
  order_tasks <- function(indices, stages) indices[order(first_dataset[task_chunks[indices]],
      match(tasks$stage[indices], stages), match(tasks$method[indices], methods),
      tasks$index[indices], method = "radix")]
  challenge_indices <- order_tasks(challenge_indices, "challenge")
  audit_indices <- order_tasks(audit_indices,
        c("original_refit", "corrected_refit", "corrected_resimulation"))
  prefix <- c(challenge_indices, audit_indices)
  actual_ids <- unique(flat[flat_chunks %in% unique(task_chunks[audit_indices])])
  actual <- design[match(actual_ids, design$dataset_id), ]
  if (anyNA(actual$dataset_id)) stop("Audit chunk crossed outside the legacy family")
  actual <- actual[order(actual$row, actual$scenario, method = "radix"), ]
  actual$selected_by_chunk_boundary <- !actual$audit_requested
  audit_generation <- gen$generation_id[gen$family == "legacy_reconstruction" & gen$row %in% actual$row]
  generation_prefix <- c(gen$generation_id[gen$family == "challenge"], audit_generation)
  if (anyDuplicated(generation_prefix) || anyDuplicated(tasks$index[prefix]))
    stop("Duplicate priority identity")
  whole_loci <- ifelse(ds$family[pos] == "challenge",
       gen$replicates[match(ds$row[pos], gen$generation_id)] * 22L,
       grid$n_regions[ds$row[pos]])
  chunk_loci <- as.numeric(rowsum(whole_loci, flat_chunks, reorder = FALSE))
  counts <- data.frame(group = c("Gaussian challenge", "paired legacy audit", "full remainder"),
    method_tasks = c(length(challenge_indices), length(audit_indices), nrow(tasks) - length(prefix)),
    locus_method_attempts = c(sum(chunk_loci[task_chunks[challenge_indices]]),
       sum(chunk_loci[task_chunks[audit_indices]]), sum(chunk_loci[task_chunks[-prefix]])))
  list(method_task_indices = as.integer(tasks$index[prefix]),
       generation_indices = as.integer(generation_prefix),
       audit_datasets = actual, counts = counts,
       summary = list(planned_method_tasks = nrow(tasks), planned_generation_units = nrow(gen),
         priority_method_tasks = length(prefix), priority_generation_units = length(generation_prefix),
         audit_grid_rows = length(unique(actual$row)),
         audit_requested_scenarios = length(requested), audit_actual_scenarios = nrow(actual),
         audit_boundary_extra_scenarios = sum(actual$selected_by_chunk_boundary),
         audit_loci_per_scenario = 10L,
         schedule_uses_results = FALSE, full_remainder_retained = TRUE))
}

.campaign_priority_main <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1L]
  script <- normalizePath(sub("^--file=", "", file_arg))
  source(file.path(dirname(script), "common.R"))
  opt <- .campaign_options(list(plan = "", `worker-repo` = "", out = ""))
  if (!nzchar(opt$plan) || !nzchar(opt$out)) stop("Supply --plan and a new external --out JSON path")
  path <- normalizePath(opt$plan); root <- dirname(path); plan <- readRDS(path)
  worker_repo <- normalizePath(if (nzchar(opt$`worker-repo`)) opt$`worker-repo` else plan$repo)
  # The builder may be newer than the immutable scientific workers. Validate
  # using THEIR source helpers so adding this external script cannot rebind it.
  worker <- .campaign_load_source(worker_repo)
  worker$.campaign_assert_plan(plan, worker_repo)
  manifest_path <- file.path(root, "dispatch-manifest.json")
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  files <- c("tasks.csv", "datasets.csv", "generation.csv")
  if (!identical(manifest$schema, 1L) ||
      !identical(manifest$plan_sha256, plan$plan_sha256) ||
      !identical(readLines(file.path(root, "plan.sha256"), warn = FALSE), plan$plan_sha256) ||
      !identical(manifest$plan_rds_sha256, .challenge_sha256(path)) ||
      !setequal(names(manifest$files), files)) stop("Incomplete or foreign dispatch commitment")
  for (name in files) if (!identical(manifest$files[[name]], .challenge_sha256(file.path(root, name))))
    stop("Changed frozen projection: ", name)
  original <- worker$.campaign_load_source(plan$original_source)
  body_text <- gsub("[[:space:]]", "", paste(deparse(body(original$run_simulation)), collapse = ""))
  if (!grepl("expand.grid(S=S,phi=phi,iter=seq_len(n_iter)", body_text, fixed = TRUE) ||
      !grepl("expand.grid(S=S,phi=phi,p_causal=p_causal,iter=seq_len(n_iter)", body_text, fixed = TRUE))
    stop("Unrecognized upstream scenario expansion; do not infer scenario numbering")
  result <- .campaign_priority_design(plan)
  out <- path.expand(opt$out)
  stem <- sub("[.]json$", "", out)
  reports <- paste0(stem, c("-counts.csv", "-audit-datasets.csv"))
  if (any(file.exists(c(out, reports)))) stop("Priority artifacts are immutable; choose a new --out")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  write.csv(result$counts, reports[1L], row.names = FALSE)
  write.csv(result$audit_datasets, reports[2L], row.names = FALSE)
  value <- list(schema = 1L, plan_sha256 = plan$plan_sha256,
    plan_rds_sha256 = .challenge_sha256(path),
    dispatch_manifest_sha256 = .challenge_sha256(manifest_path),
    method_task_indices = as.list(result$method_task_indices),
    generation_indices = as.list(result$generation_indices),
    design = "full Gaussian first; fixed paired legacy audit; untouched full remainder",
    audit = list(original_rows = "all 45", control_rows = "sparse + sparse_inf .5/.9; none/binary5.4/continuous5.4; all three LD references",
      S = "all five", original_phi = list(.05, .2, .6), control_phi = list(.0075, .1, .4),
      iterations = list(1L, 2L), locus_scope = "all ten loci, including entire overlapping task chunks"),
    summary = result$summary,
    builder_sha256 = .challenge_sha256(script), worker_repo = worker_repo,
    reports = stats::setNames(as.list(vapply(reports, .challenge_sha256, "")), basename(reports)))
  tmp <- paste0(out, ".tmp.", Sys.getpid())
  on.exit(unlink(tmp), add = TRUE)
  jsonlite::write_json(value, tmp, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  if (!file.rename(tmp, out)) stop("Cannot atomically seal priority manifest")
  print(result$counts)
  cat("External priority manifest:", normalizePath(out), "\n")
  cat("Every remaining full-plan task is still scheduled. No result was inspected.\n")
}

if (sys.nframe() == 0L) .campaign_priority_main()
