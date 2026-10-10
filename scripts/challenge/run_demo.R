#!/usr/bin/env Rscript
# Run with Rscript --vanilla so a source-only challenge does not restore renv.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "../.."))
setwd(repo)
if (dir.exists(".local-lib")) .libPaths(c(normalizePath(".local-lib"), .libPaths()))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
args <- commandArgs(TRUE)
opts <- list(replicates = "10", variants = "24", seed = "20261010", out = "results/challenge-demo")
if (length(args) %% 2L) stop("Arguments must be --name value pairs")
if (length(args)) for (i in seq(1, length(args), 2)) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opts)) stop("Unknown option ", key)
  opts[[key]] <- args[i + 1]
}
if (dir.exists(opts$out)) stop("Output already exists; choose a new --out for a frozen run")
dir.create(opts$out, recursive = TRUE)
panel <- generate_challenge(n_replicates = as.integer(opts$replicates),
                            p = as.integer(opts$variants), seed = as.integer(opts$seed))
write_challenge_release(panel, file.path(opts$out, "public"), file.path(opts$out, "assessor"))
public <- read_challenge_release(file.path(opts$out, "public"))
methods <- list(exact = challenge_exact, forced_one = challenge_forced_one,
                annotation_positive = function(x) challenge_exact(x, annotation_coefficient = 1),
                abf = challenge_region_adapter(run_abf_region))
if (requireNamespace("susieR", quietly = TRUE)) {
  methods$susie <- challenge_region_adapter(run_susie_region, list(L = 10, coverage = .95))
  if (utils::packageVersion("susieR") >= "0.15.0")
    methods$susie_inf <- challenge_region_adapter(run_susie_inf_region, list(L = 10, coverage = .95))
}
message("Evaluating ", length(panel$inputs), " development tasks with ", length(methods), " predictors")
predictions <- run_challenge(public$inputs, methods)
evaluation <- evaluate_challenge(panel, predictions)
saveRDS(predictions, file.path(opts$out, "predictions.rds"), version = 2)
saveRDS(evaluation, file.path(opts$out, "evaluation.rds"), version = 2)
utils::write.csv(evaluation$summary, file.path(opts$out, "summary.csv"), row.names = FALSE)
utils::write.csv(evaluation$tasks, file.path(opts$out, "tasks.csv"), row.names = FALSE)
comparisons <- do.call(rbind, lapply(unique(evaluation$tasks$scenario), function(s) {
  z <- challenge_compare(evaluation, "exact", "forced_one", "brier", s, n_boot = 500)
  as.data.frame(c(list(scenario = s, metric = "brier", contrast = "exact - forced_one"), z))
}))
utils::write.csv(comparisons, file.path(opts$out, "paired-comparisons.csv"), row.names = FALSE)
provenance <- list(specification = panel$config$specification, dataset = "development-v1",
  run = basename(opts$out), seed = opts$seed, replicates = opts$replicates, variants = opts$variants,
  R = R.version.string, session = utils::sessionInfo(),
  source_sha256 = vapply(list.files("R", full.names = TRUE), .challenge_sha256, ""),
  commit = if (file.exists(".git")) system2("git", c("rev-parse", "HEAD"), stdout = TRUE) else
    Sys.getenv("FMB_SOURCE_COMMIT", unset = NA_character_),
  dirty = if (file.exists(".git")) length(system2("git", c("status", "--porcelain"), stdout = TRUE)) > 0 else NA,
  public_manifest_sha256 = .challenge_sha256(file.path(opts$out, "public/manifest.csv")))
saveRDS(provenance, file.path(opts$out, "provenance.rds"), version = 2)
dput(provenance[c("specification", "dataset", "run", "seed", "replicates", "variants", "R", "commit", "dirty")],
     file = file.path(opts$out, "provenance.txt"))

# Exploratory plots; intervals are in the paired-comparisons table and remain
# mechanism-specific. No whole-suite rank is manufactured.
grDevices::pdf(file.path(opts$out, "report.pdf"), width = 12, height = 8)
for (metric in c("success_rate", "brier_successful", "recall_all_attempts", "null_false_call_rate_successful")) {
  d <- evaluation$summary
  wide <- matrix(NA_real_, nrow = length(methods), ncol = length(unique(d$scenario)),
                 dimnames = list(names(methods), sort(unique(d$scenario))))
  for (i in seq_len(nrow(d))) wide[d$method[i], d$scenario[i]] <- d[[metric]][i]
  graphics::par(mar = c(12, 10, 4, 2))
  graphics::image(seq_len(ncol(wide)), seq_len(nrow(wide)), t(wide),
    axes = FALSE, xlab = "", ylab = "", main = metric,
    col = grDevices::hcl.colors(40, "YlOrRd", rev = TRUE))
  graphics::axis(1, at = seq_len(ncol(wide)), labels = colnames(wide), las = 2, cex.axis = .8)
  graphics::axis(2, at = seq_len(nrow(wide)), labels = rownames(wide), las = 2, cex.axis = .8)
  graphics::mtext("Development diagnostic; read availability and task definition before comparing", side = 3, cex = .7)
}
grDevices::dev.off()
message("Completed: ", normalizePath(opts$out), "; see summary.csv, report.pdf and provenance.rds")
