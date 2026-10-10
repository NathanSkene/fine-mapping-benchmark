test_that("campaign checkpoints reject changed source, data and task identity", {
  repo <- tempfile(); dir.create(repo); dir.create(file.path(repo, "R"))
  writeLines("x <- 1", file.path(repo, "R/x.R"))
  snapshot <- .campaign_snapshot(repo)
  expect_silent(.campaign_assert_snapshot(repo, snapshot))
  writeLines("x <- 2", file.path(repo, "R/x.R"))
  expect_error(.campaign_assert_snapshot(repo, snapshot), "source changed")
  path <- file.path(repo, "frozen.rds")
  .campaign_atomic_rds(list(value = 1), path)
  sha <- .challenge_sha256(path)
  expect_identical(.campaign_read_frozen(path, sha), list(value = 1))
  .campaign_atomic_rds(list(value = 2), path)
  expect_error(.campaign_read_frozen(path, sha), "missing or changed")
  plan <- list(plan_sha256 = "plan-a")
  task <- data.frame(index = 1L)
  .campaign_atomic_rds(list(terminal = TRUE, plan_sha256 = "plan-a", task = task), path)
  expect_true(.campaign_validate_checkpoint(path, plan, task))
  expect_error(.campaign_validate_checkpoint(path, plan, data.frame(index = 2L)), "incompatible")
  expect_error(.campaign_validate_checkpoint(path, list(plan_sha256 = "plan-b"), task), "incompatible")
})

test_that("campaign configuration reproduces the original neural and FINEMAP settings", {
  args <- .campaign_default_args(tempdir())
  expect_equal(args$beatrice$n_caus, 10)
  expect_equal(args$beatrice$sigma_sq, .0899)
  expect_equal(args$beatrice$max_iter, 2000L)
  expect_equal(args$beatrice$sparse_concrete, 51L)
  expect_equal(args$functional_beatrice$prior_regularisation, .1255)
  expect_equal(args$functional_beatrice$lambda_l1, .1223)
  expect_equal(args$functional_beatrice$hierarchy_M, 10.15)
  expect_equal(args$finemap$n_causal, 5L)
  expect_equal(.campaign_default_args(tempdir(), "common")$finemap$n_causal, 10L)
})

test_that("a failed joint method cannot be counted as its local fallback", {
  fit <- list(scenario_id = 1L, region_id = 1L, pip = c(.2, .8),
              additional = list(joint_fallback = TRUE))
  result <- list(fb_pooled = list(results = list(fit)))
  sim <- list(genotypes = list(list(annotations_matrix = matrix(1, 2L, 1L))),
              scenarios = list(list(regions = list(list()))))
  rejected <- .campaign_reject_substitution(result, sim, "fb_pooled")
  expect_match(rejected$fb_pooled$results[[1L]]$error, "substituted")
  expect_equal(rejected$fb_pooled$n_failed, 1L)
  sim$genotypes[[1L]]$annotations_matrix <- NULL
  expect_null(.campaign_reject_substitution(result, sim, "fb_pooled")$fb_pooled$results[[1L]]$error)
})

test_that("challenge campaign preserves out-of-scope attempts and observed-only inputs", {
  panel <- generate_challenge(scenarios = c("single_sparse", "meta_known_overlap"),
                               p = 12L, n_replicates = 1L, seed = 721L)
  pred <- .campaign_challenge_predict(panel, seq_along(panel$inputs), "abf", list())
  unsupported <- vapply(panel$inputs, function(x) !is.null(x$n_per_variant), TRUE)
  expect_equal(sum(unsupported), 1L)
  expect_equal(pred[[which(unsupported)]]$status, "failed")
  expect_match(pred[[which(unsupported)]]$error, "unsupported_input")
  expect_equal(pred[[which(!unsupported)]]$status, "ok")
  sim <- .campaign_panel_simulation(panel, which(!unsupported))
  expect_null(sim$scenarios[[1L]]$regions[[1L]]$truth)
  expect_null(sim$genotypes[[1L]]$X)
})

test_that("regenerated phenotypes retain fixed genotypes and vary only the intended contrast", {
  set.seed(59)
  X <- scale(matrix(rnorm(90), 30L, 3L))
  geno <- list(list(X = X, n = 30L, p = 3L, LD = cor(X), maf = rep(.3, 3L)))
  truth <- simulate_phenotypes(geno, S = 1L, phi = .2, seed = 23L, verbose = FALSE)
  sim <- list(genotypes = geno, scenarios = list(list(S = 1L, phi = .2, model = "sparse",
           scenario_id = 1L, regions = truth, iter = 1L)),
           params = list(annotation_type = "none", model = "sparse"))
  regenerated <- .campaign_corrected_phenotype(sim, 1L, 182L)
  expect_identical(regenerated$genotypes, sim$genotypes)
  expect_identical(regenerated$scenarios[[1L]]$S, sim$scenarios[[1L]]$S)
  expect_false(identical(regenerated$scenarios[[1L]]$regions[[1L]]$y, truth[[1L]]$y))
  expect_identical(regenerated, .campaign_corrected_phenotype(sim, 1L, 182L))
})

test_that("fit reuse requires exact reachable code, inputs, settings and randomness", {
  a <- new.env(parent = baseenv()); b <- new.env(parent = baseenv())
  for (env in list(a, b)) {
    env$.FM_REGISTRY <- list(toy = "run_toy_region")
    evalq({helper <- function(z) z * 2
           run_toy_region <- function(region_geno, region_pheno) helper(region_pheno$z)}, env)
  }
  sim <- list(genotypes = list(list(LD = diag(2L), n = 20L)),
               scenarios = list(list(regions = list(list(z = c(1, 2), truth = list(secret = 1))))))
  x <- .campaign_fit_cache_identity(sim, "toy", list(), 7L, a, list(R = "same"))
  y <- .campaign_fit_cache_identity(sim, "toy", list(), 7L, b, list(R = "same"))
  expect_true(x$eligible)
  expect_identical(x$key, y$key)
  evalq(helper <- function(z) z * 3, b)
  expect_false(identical(x$key, .campaign_fit_cache_identity(sim, "toy", list(), 7L, b, list(R = "same"))$key))
  expect_false(identical(x$key, .campaign_fit_cache_identity(sim, "toy", list(), 8L, a, list(R = "same"))$key))
  sim$scenarios[[1L]]$regions[[1L]]$z[1L] <- 5
  expect_false(identical(x$key, .campaign_fit_cache_identity(sim, "toy", list(), 7L, a, list(R = "same"))$key))
  evalq(run_toy_region <- function(region_geno, region_pheno) region_pheno$truth$secret, b)
  expect_false(.campaign_method_fingerprint(b, "toy")$eligible)
  evalq(run_toy_region <- function(region_geno, region_pheno) get("helper")(region_pheno$z), b)
  expect_false(.campaign_method_fingerprint(b, "toy")$eligible)
})

test_that("prepared input manifests reject absent and altered source files", {
  dir <- tempfile(); dir.create(dir)
  path <- file.path(dir, "source.vcf"); writeLines("verified alleles", path)
  manifest <- .campaign_input_manifest("", dir)
  expect_silent(.campaign_check_inputs(path, manifest))
  writeLines("changed alleles", path)
  expect_error(.campaign_check_inputs(path, manifest), "input changed")
  other <- file.path(dir, "other.vcf"); writeLines("unplanned", other)
  expect_error(.campaign_check_inputs(other, manifest), "absent")
})

test_that("exact log loss preserves wrong endpoint probabilities", {
  expect_equal(.campaign_exact_log_loss(c(TRUE, FALSE), c(1, 0)), 0)
  expect_equal(.campaign_exact_log_loss(c(TRUE, FALSE), c(0, 1)), Inf)
  expect_equal(.campaign_exact_log_loss(c(TRUE, FALSE), c(.5, .5)), log(2))
})

test_that("native all-causal sets keep their containment event distinct from per-signal recovery", {
  global <- .campaign_set_metrics(list(1:8), c(2L, 5L), "native_all_causal_containment_set")
  expect_true(is.na(global$cs_distinct_recall))
  expect_equal(global$native_all_causal_coverage_reported, 1)
  expect_equal(global$cs_causal_variant_recall_all_attempts, 1)
  expect_equal(global$cs_support_size, 8L)
  partial <- .campaign_set_metrics(list(1:3), c(2L, 5L), "native_all_causal_containment_set")
  expect_equal(partial$native_all_causal_coverage_reported, 0)
  expect_equal(partial$cs_causal_variant_recall_all_attempts, .5)
  local <- .campaign_set_metrics(list(1:3, 4:6), c(2L, 5L), "per_signal")
  expect_equal(local$cs_distinct_recall, 1)
  expect_true(is.na(local$native_all_causal_coverage_reported))
  failed <- .campaign_set_metrics(list(), c(2L, 5L), "native_all_causal_containment_set", FALSE)
  expect_true(is.na(failed$native_all_causal_coverage_reported))
  expect_equal(failed$native_all_causal_recovery_all_attempts, 0)
})

test_that("the surrogate mixture event is preserved as a distinct probability target", {
  event <- .campaign_event_metadata("sbayesrc", list(pip_event = "any_nonzero_mixture_component"))
  expect_equal(event$pip_event, "any_nonzero_mixture_component")
  expect_false(event$probability_target_aligned)
  expect_match(event$probability_score_interpretation, "diagnostic")
  panel <- generate_challenge(scenarios = "sparse_inf", p = 12L, n_replicates = 1L, seed = 922L)
  input <- panel$inputs[[1L]]
  pred <- list(task_id = input$task_id, variant_id = input$variant_id, pip = rep(1, 12L),
               status = "ok", credible_sets = list(), pip_event = "any_nonzero_mixture_component")
  ev <- evaluate_challenge(panel, list(sbayesrc = stats::setNames(list(pred), input$task_id)))
  scores <- .campaign_challenge_task_metrics(ev$tasks, list(pred), panel, "sbayesrc")
  expect_false(scores$probability_target_aligned)
  expect_equal(scores$brier_any_nonzero_total_effect, 0)
  expect_equal(scores$log_loss_any_nonzero_total_effect, 0)
  expect_equal(scores$log_loss, Inf)
})

test_that("null failures cannot change the AP eligibility denominator", {
  ap <- c(NA_real_, .8, NA_real_, NA_real_)
  n_positive <- c(0L, 1L, 1L, 0L)
  scores <- .campaign_ap_failure_zero(ap, c(TRUE, TRUE, FALSE, FALSE), n_positive)
  expect_equal(mean(scores, na.rm = TRUE), .4)
  changed_null_status <- .campaign_ap_failure_zero(ap, c(FALSE, TRUE, FALSE, TRUE), n_positive)
  expect_identical(scores, changed_null_status)
})

test_that("terminal receipts distinguish missing commits from altered valid checkpoints", {
  directory <- tempfile(); dir.create(directory)
  checkpoint <- file.path(directory, "prediction.rds")
  plan <- list(plan_sha256 = strrep("a", 64L)); task <- list(dataset = "example")
  .campaign_atomic_rds(list(terminal = TRUE, plan_sha256 = plan$plan_sha256,
                           task = task, pip = c(.1, .9)), checkpoint)
  expect_false(.campaign_check_receipt(directory, checkpoint, plan$plan_sha256, "example"))
  .campaign_receipt(directory, checkpoint, plan$plan_sha256, "example")
  expect_true(.campaign_check_receipt(directory, checkpoint, plan$plan_sha256, "example"))
  expect_error(.campaign_check_receipt(directory, checkpoint, plan$plan_sha256, "foreign"), "Foreign")
  .campaign_atomic_rds(list(terminal = TRUE, plan_sha256 = plan$plan_sha256,
                           task = task, pip = c(.8, .2)), checkpoint)
  expect_true(.campaign_validate_checkpoint(checkpoint, plan, task))
  expect_error(.campaign_check_receipt(directory, checkpoint, plan$plan_sha256, "example"), "bytes changed")
  unlink(checkpoint)
  expect_error(.campaign_check_receipt(directory, checkpoint, plan$plan_sha256, "example"), "is missing")
  writeLines("invalid JSON", file.path(directory, "terminal.json"))
  expect_error(.campaign_check_receipt(directory, checkpoint, plan$plan_sha256, "example"), "Malformed")
})

test_that("completed tasks seal inner artifacts and compact readers explicitly select their evidence", {
  directory <- tempfile(); dir.create(directory)
  inner <- file.path(directory, "inner.rds"); checkpoint <- file.path(directory, "completed.rds")
  .campaign_atomic_rds(list(pip = c(.2, .8)), inner)
  artifacts <- .campaign_artifact_manifest(inner)
  expect_silent(.campaign_validate_artifacts(artifacts))
  plan <- list(plan_sha256 = "plan-a"); task <- list(index = 1L)
  .campaign_atomic_rds(list(terminal = TRUE, plan_sha256 = plan$plan_sha256,
                           task = task, artifacts = artifacts), checkpoint)
  expect_true(.campaign_validate_checkpoint(checkpoint, plan, task))
  .campaign_atomic_rds(list(pip = c(.7, .3)), inner)
  expect_error(.campaign_validate_checkpoint(checkpoint, plan, task), "bytes changed")
  expect_true(.campaign_validate_checkpoint(checkpoint, plan, task, validate_artifacts = FALSE))
  unlink(inner)
  expect_error(.campaign_validate_checkpoint(checkpoint, plan, task), "is missing")
  expect_error(.campaign_validate_artifacts(data.frame(path = inner, sha256 = "bad")), "Malformed")
})

test_that("analysis manifests seal derived outputs and expose incomplete collection", {
  directory <- tempfile(); dir.create(directory)
  writeLines("fresh derived evidence", file.path(directory, "metrics.csv"))
  plan <- list(plan_sha256 = strrep("a", 64L), corrected_source = list(source = "frozen"),
               runtime_fingerprint = list(version = "same"))
  manifest <- .campaign_analysis_manifest(directory, plan, "incomplete", list(locus_rows = 3L))
  expect_equal(manifest$status, "incomplete")
  expect_equal(manifest$summary$locus_rows, 3L)
  expect_identical(manifest$source, plan$corrected_source)
  expect_identical(manifest$runtime_fingerprint, plan$runtime_fingerprint)
  expect_silent(.campaign_validate_artifacts(manifest$artifacts))
  receipt <- jsonlite::fromJSON(file.path(directory, "analysis-manifest.json"))
  expect_identical(receipt$manifest_sha256, .challenge_sha256(file.path(directory, "analysis-manifest.rds")))
  writeLines("altered derived evidence", file.path(directory, "metrics.csv"))
  expect_error(.campaign_validate_artifacts(manifest$artifacts), "bytes changed")
})
