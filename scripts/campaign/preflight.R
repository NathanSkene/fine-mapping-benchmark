#!/usr/bin/env Rscript
# A successful dependency import is insufficient. Every requested regional
# implementation must complete an actual observed-input native/method fit.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])), "common.R"))
main <- function() {
opt <- .campaign_options(list(plan = "results/full-campaign/plan.rds", index = "1"))
plan <- readRDS(opt$plan); .campaign_assert_plan(plan, .campaign_repo)
index <- as.integer(opt$index)
if (is.na(index) || index < 1L || index > length(plan$methods)) stop("Invalid method index")
method <- plan$methods[index]
directory <- file.path(dirname(normalizePath(opt$plan)), "preflight", method)
checkpoint <- file.path(directory, "completed.rds")
if (.campaign_validate_checkpoint(checkpoint, plan, method)) {
  .campaign_check_receipt(directory, checkpoint, plan$plan_sha256, method)
  z <- readRDS(checkpoint)
  if (!z$success) stop("An existing failed preflight is evidence; use a new plan after repairing dependencies")
  .campaign_receipt(directory, checkpoint, plan$plan_sha256, method)
  cat("Successful method preflight already frozen\n"); quit(save = "no", status = 0L)
}
.campaign_claim(paste0(directory, ".lock"))
on.exit(unlink(paste0(directory, ".lock"), recursive = TRUE), add = TRUE)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
set.seed(873)
genotypes <- lapply(1:2, function(i) {
  X <- scale(matrix(rnorm(200 * 24), 200L, 24L))
  list(X = X, X_raw = X, n = 200L, p = 24L, LD = cor(X), maf = rep(.3, 24L),
       variant_ids = paste0("preflight:", i, ":", 1:24), region_id = i,
       annotations_matrix = cbind(category = rep(0:1, 12L), decoy = rnorm(24L)))
})
z <- simulate_phenotypes(genotypes, S = 2L, phi = .2, annotations = "binary",
                         n_annotations = 2L, enrichment = c(2, 1), seed = 913L, verbose = FALSE)
sim <- list(genotypes = genotypes, scenarios = list(list(scenario_id = 1L, S = 2L, phi = .2,
        iter = 1L, model = "sparse", regions = lapply(z, function(x) x[c("y", "z", "beta_hat", "se",
                           "annotations_matrix", "truth", "statistic_type")]))),
        params = list(p = rep(24L, 2L), annotation_type = "binary"))
started <- Sys.time(); error <- NULL
fits <- tryCatch(run_methods(sim, method, stats::setNames(list(plan$method_args[[method]]), method),
                             verbose = FALSE), error = function(e) { error <<- conditionMessage(e); NULL })
if (!is.null(fits)) {
  fits <- .campaign_reject_substitution(fits, sim, method)
  errors <- vapply(fits[[method]]$results, function(f) f$error %||% "", "")
  success <- length(errors) == 2L && all(!nzchar(errors))
  if (!success) error <- paste(unique(errors[nzchar(errors)]), collapse = "\n")
} else success <- FALSE
z <- list(terminal = TRUE, plan_sha256 = plan$plan_sha256, task = method, success = success,
    error = error, fits = fits, started = started, completed = Sys.time(),
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    session = utils::sessionInfo(), dependency_inventory = .campaign_dependency_inventory(plan$method_args),
    role = "installation/integration gate; not an accuracy or convergence experiment")
.campaign_atomic_rds(z, checkpoint)
if (!success) stop("Method preflight failed for ", method, ": ", error)
.campaign_receipt(directory, checkpoint, plan$plan_sha256, method)
cat(sprintf("Actual %s execution completed both annotated preflight loci in %.1fs\n", method, z$elapsed_seconds))
}
main()
