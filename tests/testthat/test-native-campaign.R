helpers <- c("scripts/native/campaign/helpers.R", "../../scripts/native/campaign/helpers.R")
helpers <- helpers[file.exists(helpers)]
if (length(helpers)) source(helpers[1L], local = TRUE)

.native_fixture <- function(p = 12L, model = "sparse_inf") {
  set.seed(1101)
  gs <- lapply(1:2, function(i) {
    raw <- matrix(rbinom(200 * p, 2, .3), 200L, p)
    list(X_raw = raw, X = scale(raw), n = 200L, p = p, LD = cor(raw), counted_allele = "ALT",
      annotations_matrix = cbind(rep(0:1, length.out = p), rep(c(0, 0, 1), length.out = p),
                                rep(c(1, 0, 0, 1), length.out = p)),
      vcf_source = paste0("fixture", i, ".vcf.gz"), region_label = paste0("block", i),
      variant_metadata = data.frame(chrom = i, position = i * 100000 + seq_len(p),
        source_id = paste0("source", i, ":", seq_len(p)), ref = "A", alt = "G"))
  })
  sparse <- lapply(1:2, function(i) c(.3, rep(0, p - 1L)))
  diffuse <- lapply(1:2, function(i) if (model == "sparse_inf") c(0, rep(.002, p - 1L)) else rep(0, p))
  y <- Reduce("+", lapply(1:2, function(i) as.numeric(gs[[i]]$X %*% (sparse[[i]] + diffuse[[i]])))) + rnorm(200)
  regions <- lapply(1:2, function(i) {
    s <- compute_summary_statistics(gs[[i]]$X, y)
    c(s, list(y = y, annotations_matrix = gs[[i]]$annotations_matrix,
      truth = list(causal_indices = 1L, beta_true = sparse[[i]], beta_sparse = sparse[[i]],
        beta_infinitesimal = diffuse[[i]], beta_total = sparse[[i]] + diffuse[[i]])))
  })
  list(genotypes = gs, scenarios = list(list(scenario_id = 1L, pi = .03, h2 = .2,
       iter = 1L, model = model, regions = regions)), params = list(annotation_type = "binary"))
}

test_that("native full design preserves 540 inputs and four independent chains", {
  skip_if_not(length(helpers) > 0L)
  x <- .native_design("full", 222L)
  expect_equal(nrow(x$draws), 10L)
  expect_equal(nrow(x$inputs), 540L)
  expect_equal(nrow(x$chains), 2160L)
  expect_true(all(table(x$chains$input_id) == 4L))
  expect_equal(anyDuplicated(x$chains$seed), 0L)
  expect_equal(x$config$replicates, 3L)
  expect_equal(x$config$pi, c(.001, .01))
  expect_equal(x$config$n_iter, 25000L)
  expect_equal(.native_profile("integration")$n, 200L)
  expect_equal(.native_profile("canary")$regions, 128L)
})

test_that("native empirical eigen layout reconstructs PSD LD and singular blocks", {
  skip_if_not(length(helpers) > 0L)
  path <- tempfile()
  R <- outer(1:12, 1:12, function(i, j) .8^abs(i - j))
  q <- .native_write_eigen(R, path)
  expect_equal(.native_read_eigen(path), R, tolerance = 2e-6)
  expect_equal(q$retained_rank, 12L)
  R[, 2] <- R[, 1]; R[2, ] <- R[1, ]; R[2, 2] <- 1
  q <- .native_write_eigen(R, path)
  expect_equal(q$retained_rank, 11L)
  expect_equal(.native_read_eigen(path), R, tolerance = 2e-6)
})

test_that("native preparation reverses dosage scale and keeps annotation exposures paired", {
  skip_if_not(length(helpers) > 0L)
  sim <- .native_fixture()
  design <- .native_design("integration", 322L)
  design$config$regions <- 2L
  design$config$p <- 12L
  directory <- tempfile(); dir.create(directory)
  plan <- c(list(plan_sha256 = "fixture-plan"), design)
  prepared <- .native_prepare_draw(sim, design$draws[1L, ], plan, directory)
  expect_equal(nrow(prepared$inputs), 3L)
  map <- read.csv(file.path(directory, "variant-map.csv"), stringsAsFactors = FALSE)
  expect_equal(map$A1, rep("G", 24L))
  expect_equal(map$A2, rep("A", 24L))
  expect_equal(anyDuplicated(map$SNP), 0L)
  gwas <- read.table(prepared$inputs$gwas_file[1L], header = TRUE)
  bhat <- unlist(lapply(sim$scenarios[[1L]]$regions, `[[`, "beta_hat"), use.names = FALSE)
  se <- unlist(lapply(sim$scenarios[[1L]]$regions, `[[`, "se"), use.names = FALSE)
  expect_equal(gwas$b * map$genotype_sd, bhat, tolerance = 1e-12)
  expect_equal(gwas$se * map$genotype_sd, se, tolerance = 1e-12)
  expect_equal(gwas$b / gwas$se, bhat / se, tolerance = 1e-12)
  controls <- readRDS(file.path(directory, "annotation-controls.rds"))
  for (i in 1:2) {
    expect_equal(colSums(controls$permuted[[i]]), colSums(controls$informative[[i]]))
    expect_equal(sort(apply(controls$permuted[[i]], 1L, paste, collapse = ":")),
                 sort(apply(controls$informative[[i]], 1L, paste, collapse = ":")))
  }
  exposed <- .native_expose_annotations(sim, 1L, "none", controls)
  expect_null(exposed$genotypes[[1L]]$annotations_matrix)
  expect_identical(exposed$scenarios[[1L]]$regions[[1L]]$truth, sim$scenarios[[1L]]$regions[[1L]]$truth)
  expect_identical(exposed$scenarios[[1L]]$regions[[1L]]$y, sim$scenarios[[1L]]$regions[[1L]]$y)
  expect_silent(.sbayesrc_native_validate_inputs(prepared$inputs$gwas_file[1L],
    prepared$inputs$ld_dir[1L], prepared$inputs$annotation_file[1L], "subset"))
  sim$genotypes[[1L]]$counted_allele <- "unknown"
  expect_error(.native_variant_table(sim$genotypes), "ALT allele")
})

test_that("native event scoring distinguishes sparse membership from nonzero effects", {
  skip_if_not(length(helpers) > 0L)
  sim <- .native_fixture(); input <- .native_design("integration", 1L)$inputs[1L, ]
  rows <- .native_metric_rows(sim, input, "native", pip = rep(.8, 24L), success = TRUE,
                            attempted = TRUE, complete = TRUE, pip_event = "any_nonzero_effect")
  expect_equal(rows$n_positive[rows$event == "sparse_component_membership"], c(1L, 1L))
  expect_equal(rows$n_positive[rows$event == "any_nonzero_effect"], c(12L, 12L))
  expect_false(any(rows$probability_target_aligned[rows$event == "sparse_component_membership"]))
  expect_true(all(rows$probability_target_aligned[rows$event == "any_nonzero_effect"]))
  missing <- .native_metric_rows(sim, input, "native")
  attempted <- .native_metric_rows(sim, input, "native", attempted = TRUE, complete = TRUE)
  metrics <- .native_aggregate_draw(rbind(missing, attempted))
  expect_true(all(metrics$n_planned == 4L))
  expect_true(all(metrics$n_attempted == 2L))
  expect_true(all(metrics$n_complete == 2L))
  expect_true(all(metrics$recall_all_completed == 0))
  expect_true(all(is.na(metrics$brier_successful)))
  expect_equal(.native_event_metrics(c(0, 1), c(FALSE, TRUE))$log_loss, 0)
  expect_true(is.infinite(.native_event_metrics(c(1, 0), c(FALSE, TRUE))$log_loss))
  expect_true(is.finite(.native_event_metrics(c(1, 0), c(FALSE, TRUE))$log_loss_clipped))
  good <- rows[1L, ]; good$complete <- TRUE
  null_failure <- good; null_failure$success <- FALSE
  null_failure$n_positive <- 0L; null_failure$ap <- NA_real_; null_failure$recall <- NA_real_
  ap <- .native_aggregate_draw(rbind(good, null_failure))
  expect_equal(ap$ap_failure_zero_completed, good$ap)
  expect_equal(ap$n_ap_eligible_completed, 1L)
  expect_equal(ap$n_planned, 2L)
})

test_that("CoreParameters contains post-burn draws and rejects truncated traces", {
  skip_if_not(length(helpers) > 0L)
  path <- tempfile()
  write.table(data.frame(Vg = rnorm(40), Ve = rexp(40)), path, row.names = FALSE, quote = FALSE)
  expect_equal(nrow(.native_read_core(path, 40L)), 40L)
  expect_error(.native_read_core(path, 20L), "draw count")
  skip_if_not_installed("posterior")
  set.seed(52)
  paths <- vapply(1:4, function(i) {
    p <- tempfile(); write.table(data.frame(Vg = rnorm(2000), Ve = rexp(2000)), p,
                                  row.names = FALSE, quote = FALSE); p
  }, "")
  d <- .native_diagnostics(paths, 2000L)
  expect_true(all(d$assessed_pass))
  bad <- read.table(paths[4L], header = TRUE); bad$Vg <- bad$Vg + 4
  write.table(bad, paths[4L], row.names = FALSE, quote = FALSE)
  d <- .native_diagnostics(paths, 2000L)
  expect_false(d$assessed_pass[d$parameter == "Vg"])
  expect_error(.native_diagnostics(paths[1:3], 2000L), "Four")
  bad$Vg <- 0; write.table(bad, paths[4L], row.names = FALSE, quote = FALSE)
  d <- .native_diagnostics(paths, 2000L)
  expect_identical(d$diagnostic_status[d$parameter == "Vg"], "constant_unassessable")
  expect_equal(d$minimum_within_chain_variance[d$parameter == "Vg"], 0)
  expect_false(d$assessed_pass[d$parameter == "Vg"])
})

test_that("native uncertainty uses paired complete genotype draws", {
  skip_if_not(length(helpers) > 0L)
  a <- data.frame(draw_id = paste0("d", 1:5), ap_successful = 1:5 / 10,
                  n_complete = 3L, n_planned = 3L)
  b <- a; b$ap_successful <- b$ap_successful - .03
  result <- .native_paired_draws(a, b, "ap_successful", 100L)
  expect_equal(result$difference, .03)
  expect_equal(result$n_paired_genotype_draws, 5L)
  a$n_complete[1L] <- 2L
  expect_equal(.native_paired_draws(a, b, "ap_successful", 100L)$n_paired_genotype_draws, 4L)
})

test_that("native CLI retains partial coverage and does not claim short-chain convergence", {
  skip_if_not(length(helpers) > 0L)
  skip_if_not_installed("posterior")
  original_repo <- normalizePath(file.path(dirname(helpers[1L]), "../../.."))
  root <- tempfile(); dir.create(root)
  # A source-frozen fixture checkout prevents concurrent production edits from
  # changing this test's source commitment between subprocess invocations.
  repo <- file.path(root, "source"); dir.create(repo)
  file.copy(file.path(original_repo, "R"), repo, recursive = TRUE)
  dir.create(file.path(repo, "scripts/native"), recursive = TRUE)
  file.copy(file.path(original_repo, "scripts/native/campaign"),
            file.path(repo, "scripts/native"), recursive = TRUE)
  repo <- normalizePath(repo)
  if (dir.exists(file.path(original_repo, ".local-lib")))
    file.symlink(file.path(original_repo, ".local-lib"), file.path(repo, ".local-lib"))
  binary <- file.path(root, "fixture-gctb.py")
  # This executable tests adapter/orchestration contracts only. Its outputs
  # are deliberately synthetic and are never scientific method evidence.
  writeLines(c("#!/usr/bin/env python3", "import sys,csv,os,random",
    "if len(sys.argv)==1:", " print('GCTB 2.5.5'); sys.exit(1)",
    "args=[x for x in sys.argv[1:] if x!='--write-mcmc-txt']",
    "a=dict(zip(args[::2],args[1::2]))", "r=random.Random(int(a['--seed']))",
    "prefix=a['--out']", "os.makedirs(os.path.dirname(prefix),exist_ok=True)",
    "rows=list(csv.DictReader(open(a['--gwas-summary']),delimiter=' '))",
    "with open(prefix+'.snpRes','w') as f:",
    " f.write('Name A1 A2 A1Effect SE PIP\\n')",
    " for x in rows: f.write(f\"{x['SNP']} {x['A1']} {x['A2']} 0 0.01 {0.2+r.random()*0.01}\\n\")",
    "os.makedirs(prefix+'.mcmcsamples',exist_ok=True)",
    "count=int(a['--chain-length'])//int(a['--thin'])-int(a['--burn-in'])//int(a['--thin'])",
    "with open(prefix+'.mcmcsamples/CoreParameters.mcmcsamples.txt','w') as f:",
    " f.write('Vg Ve\\n')",
    " for i in range(count): f.write(f'{r.gauss(0,1)} {r.expovariate(1)}\\n')"), binary)
  Sys.chmod(binary, "0755")
  # Only this source-frozen test checkout overrides biological simulation.
  writeLines(c(paste0(".native_cli_fixture <- ", paste(deparse(.native_fixture), collapse = "\n")),
    "simulate_gwfm_data <- function(...) .native_cli_fixture(24L, 'sparse')"),
    file.path(repo, "R/zzz_native_cli_fixture.R"))
  vcf <- file.path(root, "vcf"); maps <- file.path(root, "maps")
  dir.create(vcf); dir.create(maps)
  regions <- data.frame(region_id = paste0("fixture", 1:2), chrom = 1:2, start = 10000L, end = 20000L)
  region_file <- file.path(root, "regions.csv"); write.csv(regions, region_file, row.names = FALSE)
  for (id in regions$region_id) writeLines("frozen synthetic VCF fixture", file.path(vcf, paste0(id, ".vcf.gz")))
  writeLines("frozen synthetic map fixture", file.path(maps, "map.txt"))
  output_root <- file.path(root, "campaign"); plan_path <- file.path(output_root, "plan.rds")
  run <- function(script, args, include_plan = TRUE) {
    output <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla",
      shQuote(file.path(repo, "scripts/native/campaign", script)),
      if (include_plan) c("--plan", shQuote(plan_path)), args),
      stdout = TRUE, stderr = TRUE)
    expect_equal(attr(output, "status") %||% 0L, 0L, info = paste(output, collapse = "\n"))
  }
  run("create_plan.R", c("--out", shQuote(output_root), "--profile", "integration", "--seed", "143",
    "--vcf-dir", shQuote(vcf), "--genetic-map-dir", shQuote(maps), "--regions-file", shQuote(region_file),
    "--gctb-path", shQuote(binary), "--methods", "abf"), include_plan = FALSE)
  commitment <- jsonlite::fromJSON(file.path(output_root, "dispatch-manifest.json"))
  expect_equal(commitment$schema, 1L)
  expect_identical(commitment$plan_rds_sha256, .challenge_sha256(plan_path))
  expect_setequal(names(commitment$files), c("draws.csv", "inputs.csv", "chains.csv", "regional.csv"))
  for (name in names(commitment$files))
    expect_identical(commitment$files[[name]], .challenge_sha256(file.path(output_root, name)))
  run("prepare.R", c("--index", "1"))
  root <- output_root
  for (index in 1:4) run("run_chain.R", c("--index", index))
  run("run_regional.R", c("--input-id", "input_00001"))
  run("collect.R", c("--bootstrap-draws", "100"))
  coverage <- read.csv(file.path(root, "analysis/planned-chain-coverage.csv"))
  expect_equal(sum(coverage$checkpointed), 4L)
  expect_equal(nrow(coverage), 12L)
  diagnostic <- read.csv(file.path(root, "analysis/input-diagnostics.csv"))
  expect_false(any(diagnostic$core_diagnostics_pass))
  expect_true(all(diagnostic$extension_eligible))
  outcomes <- read.csv(file.path(root, "analysis/task-metrics.csv"))
  expect_true(all(outcomes$success[outcomes$input_id == "input_00001"]))
  expect_true(all(!outcomes$complete[outcomes$input_id != "input_00001"]))
  consensus <- readRDS(file.path(root, "analysis/consensus/input_00001.rds"))
  expect_equal(ncol(consensus$pip_chains), 4L)
  expect_false(consensus$core_diagnostics_pass)
  frozen_plan <- readRDS(plan_path)
  analysis_manifest <- readRDS(file.path(root, "analysis/analysis-manifest.rds"))
  expect_identical(analysis_manifest$collection_status, "incomplete")
  expect_false(analysis_manifest$campaign_complete)
  expect_identical(analysis_manifest$plan_sha256, frozen_plan$plan_sha256)
  expect_identical(analysis_manifest$source, frozen_plan$source)
  expect_true(.campaign_check_receipt(file.path(root, "analysis"),
    file.path(root, "analysis/analysis-manifest.rds"), frozen_plan$plan_sha256, "analysis"))
  expect_silent(.campaign_validate_artifacts(analysis_manifest$artifacts))
  expect_false(dir.exists(file.path(root, "analysis.lock")))
  # Resume must leave existing attempt bytes unchanged.
  file <- file.path(root, "chains", "chain_000001", "completed.rds")
  before <- .challenge_sha256(file)
  run("run_chain.R", c("--index", "1"))
  expect_identical(.challenge_sha256(file), before)
  # An interruption after checkpoint rename but before receipt rename must
  # resume without recomputation and regenerate the missing terminal receipt.
  resume <- list(
    list(script = "prepare.R", args = c("--index", "1"), path = file.path(root, "prepared/draw_001")),
    list(script = "run_chain.R", args = c("--index", "1"), path = file.path(root, "chains/chain_000001")),
    list(script = "run_regional.R", args = c("--input-id", "input_00001"), path = file.path(root, "regional/regional_000001")))
  for (unit in resume) {
    completed <- file.path(unit$path, "completed.rds")
    digest <- .challenge_sha256(completed)
    unlink(file.path(unit$path, "terminal.json"))
    run(unit$script, unit$args)
    expect_identical(.challenge_sha256(completed), digest)
    receipt <- jsonlite::fromJSON(file.path(unit$path, "terminal.json"))
    expect_true(receipt$terminal)
    expect_identical(receipt$checkpoint_sha256, digest)
  }
  # Collection sees a missing receipt as pending, even when its checkpoint
  # exists. It archives all prior derived outputs before rebuilding them.
  for (unit in resume) {
    receipt_path <- file.path(unit$path, "terminal.json")
    receipt_bytes <- readBin(receipt_path, "raw", file.info(receipt_path)$size)
    unlink(receipt_path)
    run("collect.R", c("--bootstrap-draws", "100"))
    fresh <- readRDS(file.path(root, "analysis/analysis-manifest.rds"))
    expect_identical(fresh$collection_status, "incomplete")
    expect_false(file.exists(file.path(root, "analysis/consensus/input_00001.rds")) &&
                 unit$script != "run_regional.R")
    if (unit$script == "prepare.R") {
      expect_equal(fresh$coverage$n_prepared, 0L)
      expect_false(file.exists(file.path(root, "analysis/task-metrics.csv")))
      expect_false(file.exists(file.path(root, "analysis/draw-metrics.csv")))
    } else if (unit$script == "run_chain.R") {
      current <- read.csv(file.path(root, "analysis/planned-chain-coverage.csv"))
      expect_false(current$checkpointed[current$chain_id == "chain_000001"])
    } else {
      current <- read.csv(file.path(root, "analysis/planned-regional-coverage.csv"))
      expect_false(current$checkpointed[current$regional_id == "regional_000001"])
    }
    writeBin(receipt_bytes, receipt_path)
  }
  expect_equal(length(list.dirs(file.path(root, "analysis-history"), recursive = FALSE)), 3L)
  expect_false(dir.exists(file.path(root, "analysis.lock")))
  fail_collection <- function(pattern) {
    rejected <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla",
      shQuote(file.path(repo, "scripts/native/campaign/collect.R")), "--plan", shQuote(plan_path),
      "--bootstrap-draws", "100"), stdout = TRUE, stderr = TRUE))
    expect_equal(attr(rejected, "status"), 1L)
    expect_true(any(grepl(pattern, rejected, fixed = TRUE)), info = paste(rejected, collapse = "\n"))
    manifest <- readRDS(file.path(root, "analysis/analysis-manifest.rds"))
    expect_identical(manifest$collection_status, "failed")
    expect_true(grepl(pattern, manifest$error, fixed = TRUE))
    expect_false(dir.exists(file.path(root, "analysis.lock")))
  }
  # A changed PIP checkpoint and a changed external parameter trace are both
  # rejected, rather than combined into a new apparently valid analysis.
  completed <- file.path(root, "chains/chain_000001/completed.rds")
  original_bytes <- readBin(completed, "raw", file.info(completed)$size)
  changed <- readRDS(completed); changed$fit$pip[1L] <- .9
  saveRDS(changed, completed, version = 2)
  fail_collection("Receipted checkpoint bytes changed")
  writeBin(original_bytes, completed)
  trace <- readRDS(completed)$core_trace
  original_trace <- readBin(trace, "raw", file.info(trace)$size)
  cat("0 0\n", file = trace, append = TRUE)
  fail_collection("Sealed artifact bytes changed")
  writeBin(original_trace, trace)
  run("collect.R", c("--bootstrap-draws", "100"))
  expect_equal(sum(read.csv(file.path(root, "analysis/planned-chain-coverage.csv"))$checkpointed), 4L)
  expect_true(file.exists(file.path(root, "analysis/consensus/input_00001.rds")))
  # Preparation validates the original source data; worker/collector compact
  # checks validate the prepared bytes that they actually use.
  writeLines("tampered fixture", file.path(vcf, "fixture1.vcf.gz"))
  rejected <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), c("--vanilla",
    shQuote(file.path(repo, "scripts/native/campaign/prepare.R")), "--plan", shQuote(plan_path), "--index", "1"),
    stdout = TRUE, stderr = TRUE))
  expect_equal(attr(rejected, "status"), 1L)
  expect_true(any(grepl("VCF/map/region bytes changed", rejected, fixed = TRUE)))
})
