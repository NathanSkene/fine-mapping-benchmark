#!/usr/bin/env Rscript
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
opt <- .native_options(list(out = "results/native-campaign", profile = "full", seed = "20261012",
  `vcf-dir` = "data/gwfm_vcf", `genetic-map-dir` = "data/genetic_maps",
  `regions-file` = file.path(.native_repo, "inst/extdata/gwfm_regions.csv"),
  `trace-contract` = file.path(.native_repo, "scripts/native/campaign/trace-contract.json"),
  `gctb-path` = "gctb", methods = paste(names(.FM_REGISTRY), collapse = ",")))
if (dir.exists(opt$out)) stop("Use a new output directory for an immutable native plan")
if (!requireNamespace("posterior", quietly = TRUE)) stop("Install posterior before freezing a native campaign")
seed <- as.integer(opt$seed); if (is.na(seed) || seed < 1L) stop("Invalid seed")
design <- .native_design(opt$profile, seed)
region_file <- normalizePath(opt$`regions-file`, mustWork = TRUE)
regions <- read.csv(region_file, stringsAsFactors = FALSE)
if (!all(c("region_id", "chrom", "start", "end") %in% names(regions)) ||
    anyDuplicated(regions$region_id) || anyNA(regions[c("region_id", "chrom", "start", "end")]) ||
    any(regions$start < 1 | regions$end <= regions$start)) stop("Invalid frozen region coordinates")
if (design$config$regions == 128L && (nrow(regions) != 128L ||
    !setequal(as.integer(sub("^chr", "", regions$chrom)), 1:22)))
  stop("Full-size native inputs require 128 unique source regions across all 22 autosomes")
if (nrow(regions) < design$config$regions) stop("Not enough regions for the frozen native profile")
if (opt$profile != "full") regions <- regions[seq_len(design$config$regions), ]
vcf <- normalizePath(opt$`vcf-dir`, mustWork = TRUE)
maps <- normalizePath(opt$`genetic-map-dir`, mustWork = TRUE)
files <- c(region_file, file.path(vcf, paste0(regions$region_id, ".vcf.gz")),
           list.files(maps, recursive = TRUE, full.names = TRUE))
files <- sort(unique(files[!dir.exists(files)]), method = "radix")
if (any(!file.exists(files))) stop("Prepare every required source VCF before freezing the native plan")
methods <- strsplit(opt$methods, ",", fixed = TRUE)[[1L]]
if (anyDuplicated(methods) || any(!methods %in% names(.FM_REGISTRY))) stop("Invalid regional methods")
regional <- expand.grid(input_id = design$inputs$input_id, method = methods, stringsAsFactors = FALSE)
regional$regional_id <- sprintf("regional_%06d", seq_len(nrow(regional)))
binary <- if (file.exists(opt$`gctb-path`)) normalizePath(opt$`gctb-path`) else unname(Sys.which(opt$`gctb-path`))
if (!nzchar(binary)) stop("Pin an existing official GCTB executable")
trace_contract <- jsonlite::fromJSON(opt$`trace-contract`)
.native_validate_trace_contract(trace_contract, .challenge_sha256(binary), opt$profile)
method_args <- .campaign_default_args(.native_repo)[methods]
plan <- c(list(specification = "native-sbayesrc-comparison-v1", profile = opt$profile, seed = seed,
  source = .native_source_snapshot(.native_repo), regions = regions, regions_file = region_file,
  vcf_dir = vcf, map_dir = maps,
  data_snapshot = data.frame(path = unname(files), sha256 = unname(vapply(files, .challenge_sha256, "")),
                             row.names = NULL, stringsAsFactors = FALSE),
  executable = binary, executable_sha256 = .challenge_sha256(binary), expected_version = "2.5.5",
  trace_contract = trace_contract,
  regional = regional, method_args = method_args,
  runtime_fingerprint = .campaign_runtime_fingerprint(method_args),
  posterior_version = as.character(utils::packageVersion("posterior")),
  thresholds = list(rhat = 1.01, ess = 400),
  extension = list(eligible = "all four chains of any input failing core diagnostics; no accuracy-dependent selection",
     n_iter = 100000L, burn_in = 20000L, thin = 10L,
     disposition = "fresh independent extension chains, separately frozen extension plan; keep initial outcomes"),
  interpretation = "empirical-EUR-template representative subset; shared phenotype; block LD approximation",
  created = Sys.time()), design)
plan$plan_sha256 <- .campaign_object_sha(plan)
dir.create(opt$out, recursive = TRUE)
.campaign_atomic_rds(plan, file.path(opt$out, "plan.rds"))
writeLines(plan$plan_sha256, file.path(opt$out, "plan.sha256"))
for (name in c("draws", "inputs", "chains", "regional"))
  write.csv(plan[[name]], file.path(opt$out, paste0(name, ".csv")), row.names = FALSE)
dispatch_files <- file.path(opt$out, c("draws.csv", "inputs.csv", "chains.csv", "regional.csv"))
jsonlite::write_json(list(schema = 1L, plan_sha256 = plan$plan_sha256,
    plan_rds_sha256 = .challenge_sha256(file.path(opt$out, "plan.rds")),
    files = stats::setNames(as.list(unname(vapply(dispatch_files, .challenge_sha256, ""))),
                           basename(dispatch_files))),
    file.path(opt$out, "dispatch-manifest.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("Frozen %s native plan: %d genotype draws, %d inputs, %d native chains, %d regional fits\n",
  opt$profile, nrow(plan$draws), nrow(plan$inputs), nrow(plan$chains), nrow(plan$regional)))
