#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
opt <- .campaign_options(list(out = "results/full-campaign", `legacy-root` = "",
  `original-source` = "", `method-args` = "", profile = "full",
  `inputs-manifest` = "",
  methods = paste(names(.FM_REGISTRY), collapse = ","),
  stages = "archived_rescore,original_refit,corrected_refit,corrected_resimulation",
  `vcf-dir` = file.path(.campaign_repo, "data/vcf"),
  `genetic-map-dir` = file.path(.campaign_repo, "data/genetic_maps"),
  `challenge-p` = "40,200,1000", `challenge-n` = "1000,10000",
  `challenge-replicates` = "100", `challenge-batch-replicates` = "10",
  `legacy-chunk` = "5", `iter003-controls` = "true", seed = "20261011"))
if (!opt$profile %in% c("full", "integration", "canary", "resource_probe")) stop("Unknown profile")
if (dir.exists(opt$out)) stop("A campaign plan is immutable; choose a new --out")
methods <- trimws(strsplit(opt$methods, ",", fixed = TRUE)[[1L]])
if (anyDuplicated(methods) || any(!methods %in% names(.FM_REGISTRY))) stop("Unknown/duplicate method")
stages <- strsplit(opt$stages, ",", fixed = TRUE)[[1L]]
valid_stages <- c("archived_rescore", "original_refit", "corrected_refit", "corrected_resimulation")
if (any(!stages %in% valid_stages)) stop("Unknown legacy stages")
original <- if (nzchar(opt$`original-source`)) normalizePath(opt$`original-source`, mustWork = TRUE) else ""
if ("original_refit" %in% stages && !nzchar(original)) stop("Original refits require frozen --original-source")
method_args <- if (nzchar(opt$`method-args`)) readRDS(opt$`method-args`) else .campaign_default_args(.campaign_repo)
if (!is.list(method_args) || any(!methods %in% names(method_args))) stop("Missing method arguments")
seed <- as.integer(opt$seed)
if (is.na(seed) || seed < 0) stop("Invalid seed")
out <- opt$out
dir.create(out, recursive = TRUE)
out <- normalizePath(out)

# Preserve the complete Iteration-004 design. Integration is solely an
# installation/contract run and is never presented as the scientific campaign.
grid_path <- file.path(out, "legacy-grid.csv")
status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(file.path(.campaign_repo, "scripts/hpc/generate_params_grid.R")),
      shQuote(grid_path)), stdout = file.path(out, "grid-generation.log"),
      stderr = file.path(out, "grid-generation-error.log"))
if (status) stop("Original grid generation failed")
grid <- read.csv(grid_path, stringsAsFactors = FALSE)
if (identical(opt$`iter003-controls`, "true") && opt$profile == "full") {
  controls <- grid[rep(seq_len(nrow(grid)), each = 3L), ]
  controls$n <- 1000L
  controls$p_values <- paste(rep(c(100L, 200L, 400L, 500L, 1000L), each = 2L), collapse = "|")
  controls$phi_values <- "0.0075|0.05|0.1|0.2|0.4"
  # The sparse_inf values are changed by architectural index rather than their
  # floating labels, preserving all four original-003 background fractions.
  controls$p_causal <- ifelse(controls$model == "sparse", NA_real_,
      c(`0.2` = .5, `0.4` = .7, `0.6` = .9, `0.8` = 1)[as.character(controls$p_causal)])
  controls$n_ref <- rep(c(NA_integer_, 500L, 200L), nrow(grid))
  controls$label <- paste0("iter003_control_", controls$label, "_",
                        ifelse(is.na(controls$n_ref), "insample", paste0("ref", controls$n_ref)))
  controls$job_id <- seq.int(nrow(grid) + 1L, nrow(grid) + nrow(controls))
  grid <- rbind(grid, controls)
  write.csv(grid, grid_path, row.names = FALSE)
}
if (opt$profile == "integration") {
  grid <- grid[c(1L, 2L), ]
  grid$n <- 200L; grid$n_regions <- 2L; grid$p_values <- "24|40"
  grid$n_iter <- 1L; grid$S_values <- "1|2"; grid$phi_values <- ".1"
  write.csv(grid, grid_path, row.names = FALSE)
}
if (opt$profile %in% c("canary", "resource_probe")) {
  grid <- grid[c(1L, 2L, 10L, 11L), ]
  grid$n_iter <- 1L; grid$S_values <- "2"; grid$phi_values <- ".1"
  write.csv(grid, grid_path, row.names = FALSE)
}
generation <- list(); datasets <- list(); gi <- di <- 0L
add_dataset <- function(row) { di <<- di + 1L; datasets[[di]] <<- row }
if (nzchar(opt$`legacy-root`)) {
  archived <- .campaign_legacy_index(normalizePath(opt$`legacy-root`, mustWork = TRUE), methods)
  for (i in seq_len(nrow(archived))) add_dataset(archived[i, ])
} else {
  # A replay reconstructs the original generator under a frozen old snapshot;
  # it is not called an archived prediction or exact historical replication.
  if (!nzchar(original)) stop("Without original artifacts, original generation requires --original-source")
  for (i in seq_len(nrow(grid))) {
    gi <- gi + 1L
    generation[[gi]] <- data.frame(generation_id = gi, family = "legacy_reconstruction",
        row = i, p = NA_integer_, n = grid$n[i], replicates = NA_integer_,
        seed = 1000L + grid$job_id[i], stringsAsFactors = FALSE)
    nsc <- length(strsplit(grid$S_values[i], "|", fixed = TRUE)[[1L]]) *
           length(strsplit(grid$phi_values[i], "|", fixed = TRUE)[[1L]]) * grid$n_iter[i]
    for (sc in seq_len(nsc)) add_dataset(data.frame(
        dataset_id = sprintf("reconstructed_row%03d_sc%03d", i, sc),
        family = "legacy_reconstruction", row = i, scenario = sc,
        sim_file = file.path(out, "datasets", sprintf("generation_%05d", gi), "simulation.rds"),
        sim_sha256 = "generated", result_files = "", result_sha256 = "",
        genotypes_sha256 = "generated",
        original_label = grid$label[i], stringsAsFactors = FALSE))
  }
  stages <- setdiff(stages, "archived_rescore")
}
cp <- as.integer(strsplit(opt$`challenge-p`, ",", fixed = TRUE)[[1L]])
cn <- as.integer(strsplit(opt$`challenge-n`, ",", fixed = TRUE)[[1L]])
replicates <- as.integer(opt$`challenge-replicates`)
batch <- as.integer(opt$`challenge-batch-replicates`)
if (opt$profile == "integration") { cp <- 24L; cn <- 500L; replicates <- 2L; batch <- 1L }
if (opt$profile == "canary") { cp <- 200L; cn <- 10000L; replicates <- 1L; batch <- 1L }
if (opt$profile == "resource_probe") { cp <- 1000L; cn <- 10000L; replicates <- 10L; batch <- 10L }
if (anyNA(c(cp, cn, replicates, batch)) || any(cp < 12L) || any(cn < 3L) ||
    replicates < 1L || batch < 1L || replicates %% batch) stop("Invalid/divisible challenge dimensions")
for (p in cp) for (n in cn) for (b in seq_len(replicates / batch)) {
  gi <- gi + 1L
  generation[[gi]] <- data.frame(generation_id = gi, family = "challenge", row = b,
       p = p, n = n, replicates = batch, seed = seed + 100000L + gi,
       stringsAsFactors = FALSE)
  add_dataset(data.frame(dataset_id = sprintf("challenge_p%d_n%d_batch%03d", p, n, b),
       family = "challenge", row = gi, scenario = NA_integer_,
       sim_file = file.path(out, "datasets", sprintf("generation_%05d", gi), "panel.rds"),
       sim_sha256 = "generated", result_files = "", result_sha256 = "",
       genotypes_sha256 = "generated",
       original_label = "Gaussian score stress experiment", stringsAsFactors = FALSE))
}
datasets <- do.call(rbind, datasets)
generation <- if (length(generation)) do.call(rbind, generation) else data.frame()
tasks <- list(); ti <- 0L
chunk <- as.integer(opt$`legacy-chunk`)
if (is.na(chunk) || chunk < 1L) stop("Invalid legacy chunk")
add_task <- function(ids, stage, method, packed = FALSE) {
  values <- if (packed) ids else paste(ids, collapse = "|")
  indices <- ti + seq_along(values)
  ti <<- ti + length(values)
  tasks[[length(tasks) + 1L]] <<- data.frame(task_id = sprintf("task_%07d", indices), index = indices,
      dataset_ids = values, stage = stage, method = method,
      seed = seed + 1000000L + indices, stringsAsFactors = FALSE)
}
legacy <- datasets[datasets$family != "challenge", ]
for (stage in stages) for (method in methods) {
  for (file in unique(legacy$sim_file)) {
    ids <- legacy$dataset_id[legacy$sim_file == file]
    # CARMA measured hours per whole scenario at large p; do not batch those
    # into tasks that cannot finish. Other methods amortize matrix loading.
    actual_chunk <- if (method == "carma") 1L else chunk
    packed <- vapply(split(seq_along(ids), ceiling(seq_along(ids) / actual_chunk)),
                      function(indices) paste(ids[indices], collapse = "|"), "")
    add_task(packed, stage, method, packed = TRUE)
  }
}
for (i in which(datasets$family == "challenge")) {
  for (method in setdiff(methods, "polyfun_oracle")) add_task(datasets$dataset_id[i], "challenge", method)
  p <- generation$p[generation$generation_id == datasets$row[i]]
  if (p <= 200L) for (method in c("exact", "forced_one")) add_task(datasets$dataset_id[i], "challenge", method)
}
tasks <- do.call(rbind, tasks)
manifest <- list(specification = "fmbenchmark-full-campaign-v1", profile = opt$profile,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), seed = seed,
  repo = .campaign_repo, original_source = original,
  corrected_source = .campaign_snapshot(.campaign_repo),
  original_source_hashes = if (nzchar(original)) .campaign_snapshot(original) else NULL,
  methods = methods, method_args = method_args, stages = stages,
  datasets = datasets, generation = generation, tasks = tasks, grid = grid,
  data_paths = list(vcf_dir = opt$`vcf-dir`, genetic_map_dir = opt$`genetic-map-dir`),
  input_file_manifest = .campaign_input_manifest(opt$`inputs-manifest`,
       c(opt$`vcf-dir`, opt$`genetic-map-dir`)),
  runtime_fingerprint = .campaign_runtime_fingerprint(method_args),
  interpretation = if (nzchar(opt$`legacy-root`)) "frozen archived input attribution" else
      "original-generator reconstruction; unavailable historical predictions cannot be rescored",
  challenge = list(variants = cp, n = cn, replicates_per_cell = replicates,
      mixed_mechanism_batch_replicates = batch,
      prediction_scope = "regional scalar-N methods reject multicohort and variable-N inputs",
      exact_reference_scope = "p<=200; two-effect working model; not universal oracle"))
manifest$plan_sha256 <- .campaign_object_sha(manifest)
.campaign_atomic_rds(manifest, file.path(out, "plan.rds"))
writeLines(manifest$plan_sha256, file.path(out, "plan.sha256"))
write.csv(tasks, file.path(out, "tasks.csv"), row.names = FALSE)
write.csv(datasets, file.path(out, "datasets.csv"), row.names = FALSE)
write.csv(generation, file.path(out, "generation.csv"), row.names = FALSE)
dispatch_files <- file.path(out, c("tasks.csv", "datasets.csv", "generation.csv"))
jsonlite::write_json(list(schema = 1L, plan_sha256 = manifest$plan_sha256,
    plan_rds_sha256 = .challenge_sha256(file.path(out, "plan.rds")),
    files = stats::setNames(as.list(vapply(dispatch_files, .challenge_sha256, "")),
                           basename(dispatch_files))),
    file.path(out, "dispatch-manifest.json"), auto_unbox = TRUE, pretty = TRUE)
write.csv(.campaign_dependency_inventory(method_args), file.path(out, "dependencies.csv"), row.names = FALSE)
dput(manifest[c("specification", "profile", "seed", "methods", "stages", "interpretation", "challenge", "plan_sha256")],
     file.path(out, "design.txt"))
cat(sprintf("Frozen %s plan: %d generation jobs, %d datasets, %d method tasks.\n%s\n",
    opt$profile, nrow(generation), nrow(datasets), nrow(tasks), out))
