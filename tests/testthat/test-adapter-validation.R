test_that("normal-score SuSiE-inf cannot silently disable its infinitesimal model", {
  skip_if_not_installed("susieR")
  expect_error(.fm_susie_rss(c(3, 1), diag(2), 100,
    statistic_type = "normal_score", unmappable_effects = "inf"), "unsupported_input")
  fit <- suppressMessages(run_susie_inf(c(3, 1), diag(2), 100,
    statistic_type = "normal_score", L = 1))
  expect_match(fit$error, "unsupported_input")
  expect_true(all(is.na(fit$pip)))
  expect_identical(fit$method, "susie_inf")
})

test_that("individual inputs cannot claim a Gaussian summary likelihood", {
  skip_if_not_installed("susieR")
  set.seed(82)
  X <- matrix(rnorm(90), 30, 3)
  y <- X[, 1] + rnorm(30)
  for (fn in list(run_susie, run_susie_inf)) {
    expect_error(fn(X = X, y = y, statistic_type = "normal_score"), "individual data")
    fit <- suppressMessages(fn(X = X, y = y, L = 1))
    expect_null(fit$error)
    expect_identical(fit$input_type, "individual")
    expect_true(fit$params$effective_estimate_residual_variance)
  }
})

test_that("OLS-t reconstruction matches independently centered sufficient statistics", {
  set.seed(917)
  n <- 37
  X <- matrix(rnorm(n * 4), n, 4)
  X <- sweep(X, 2, c(2, -3, 8, 1), "+")
  y <- 20 + 5 * X[, 1] - 2 * X[, 3] + rnorm(n, sd = 3)
  t <- vapply(seq_len(ncol(X)), function(j)
    coef(summary(lm(y ~ X[, j])))[2L, "t value"], numeric(1))
  Xc <- sweep(X, 2, colMeans(X), "-")
  Xpop <- sweep(Xc, 2, sqrt(colMeans(Xc^2)), "/")
  yc <- y - mean(y)
  for (meansq in c(1, mean(yc^2), 4)) {
    chosen_y <- yc * sqrt(meansq / mean(yc^2))
    exact_score <- as.numeric(crossprod(Xpop, chosen_y)) / sqrt(n)
    expect_equal(.fm_ols_t_to_score(t, n, meansq), exact_score, tolerance = 1e-12)
    expect_equal(crossprod(Xpop), n * cor(X), tolerance = 1e-12)
    expect_equal(sum(chosen_y^2), n * meansq, tolerance = 1e-12)
  }
  expect_equal(.fm_ols_t_to_score(c(0, 1e308, -1e308), n), c(0, sqrt(n), -sqrt(n)))
  expect_error(.fm_ols_t_to_score(t, 2), "n > 2")
  expect_error(.fm_ols_t_to_score(t, 37.5), "n > 2")
  expect_error(.fm_ols_t_to_score(t, n, 0), "positive meansq")
})

test_that("FINEMAP-inf types and writes the declared reconstructed likelihood", {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  writeLines("# fixture", file.path(path, "run_fine_mapping.py"))
  fn <- run_finemap_inf
  e <- new.env(parent = environment(fn))
  e$system2 <- function(command, args, ...) {
    z <- read.delim(args[match("--sumstats", args) + 1L])$Z
    expect_equal(z, .fm_ols_t_to_score(c(3, 1), 100, 4), tolerance = 1e-12)
    expect_equal(as.numeric(args[match("--meansq", args) + 1L]), 4)
    prefix <- args[match("--output-prefix", args) + 1L]
    write.table(data.frame(SNP = c("a", "b"), prob = c(.8, .2)),
      paste0(prefix, ".finemapinf.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    "success"
  }
  environment(fn) <- e
  fit <- fn(c(3, 1), diag(2), 100, variant_ids = c("a", "b"), finemap_inf_dir = path, meansq = 4)
  expect_null(fit$error)
  expect_identical(fit$params$statistic_type, "marginal_ols_t")
  expect_identical(fit$params$genotype_effect_units, "population_sd_standardized")
  expect_identical(fit$params$phenotype_effect_units, "centered_trait_with_supplied_mean_square")
  expect_identical(fit$params$likelihood, "reconstructed_centered_ols_sufficient_statistics")
  unsupported <- fn(c(3, 1), diag(2), 100, finemap_inf_dir = path, statistic_type = "normal_score")
  expect_match(unsupported$error, "unsupported_input")
  expect_true(all(is.na(unsupported$pip)))
  from_region <- run_finemap_inf_region(list(LD = diag(2), n = 100),
    list(z = c(3, 1), statistic_type = "normal_score"), finemap_inf_dir = path)
  expect_match(from_region$error, "unsupported_input")
})

test_that("joint FB sets reject corruption and retain member probability alignment", {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  write.csv(data.frame(variant_names = c("a", "b", "c"), pip = c(.8, .2, .1)),
    file.path(path, "pip.csv"), row.names = FALSE)
  parse <- function(set, probabilities = NULL) {
    writeLines(set, file.path(path, "credible_set.txt"))
    conditional <- file.path(path, "conditional_credible_variants_probability.txt")
    if (is.null(probabilities)) unlink(conditional) else writeLines(probabilities, conditional)
    .fb_joint_parse_output(path, c("a", "b", "c"), 3)
  }
  valid <- parse("1 0", ".2 .8")
  expect_identical(valid$credible_sets, list(c(1L, 2L)))
  expect_equal(valid$cs_pip, list(c(.8, .2)))
  for (bad in c("0 broken", "0 0", "0 -1", "0 3", "0 1.5", "2147483647"))
    expect_null(parse(bad))
  for (bad in c(".8", ".8 NaN", ".8 1.2")) expect_null(parse("0 1", bad))
  write.csv(data.frame(pip = c(.8, .2, .1)), file.path(path, "pip.csv"), row.names = FALSE)
  expect_null(parse("0 1"))
})

test_that("joint FB returns measured and amortized pooled inference time", {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  writeLines("# fixture", file.path(path, "beatrice_joint.py"))
  fn <- .fb_joint_scenario_setup
  e <- new.env(parent = environment(fn))
  e$system2 <- function(command, args, ...) {
    manifest <- read.delim(args[match("--manifest", args) + 1L])
    for (i in seq_len(nrow(manifest))) {
      ids <- read.table(manifest$z[i], stringsAsFactors = FALSE)[[1L]]
      write.csv(data.frame(variant_names = ids, pip = rep(.2 * i, length(ids))),
        file.path(manifest$target[i], "pip.csv"), row.names = FALSE)
    }
    "success"
  }
  tick <- 0
  e$proc.time <- function() { tick <<- tick + 4; c(elapsed = tick) }
  environment(fn) <- e
  geno <- list(LD = diag(2), annotations_matrix = matrix(1, 2, 1),
               n = 100, variant_ids = c("a", "b"))
  other <- geno; other$LD[1, 2] <- other$LD[2, 1] <- .3
  genotypes <- list(geno, other)
  regions <- list(list(z = c(1, 2)), list(z = c(1, 2)))
  setup <- fn(genotypes, regions, list(beatrice_dir = path, python = "python"), "linear")
  fits <- lapply(seq_along(regions), function(i) do.call(.fb_joint_region,
    c(list(region_geno = genotypes[[i]], region_pheno = regions[[i]], method_name = "fb_pooled"), setup)))
  expect_equal(fits[[1]]$pip, c(.2, .2))
  expect_equal(fits[[2]]$pip, c(.4, .4))
  expect_equal(vapply(fits, `[[`, numeric(1), "runtime_seconds"), c(2, 2))
  expect_equal(sum(vapply(fits, `[[`, numeric(1), "runtime_seconds")), 4)
  for (fit in fits) {
    expect_equal(fit$additional$pooled_runtime_seconds, 4)
    expect_equal(fit$additional$n_loci_pooled, 2)
    expect_true(fit$additional$pooled_across_loci)
  }
})

test_that("pooled caches distinguish complete ordered scientific inputs", {
  z <- c(1, 2, 3, 4)
  LD <- diag(4)
  A <- matrix(c(1, 0, 1, 0), 4, 1)
  ids <- letters[1:4]
  for (fingerprint in list(.fb_fingerprint, .paintor_fingerprint)) {
    key <- fingerprint(z, LD, A, 100L, ids)
    expect_match(key, "^[a-f0-9]{64}$")
    expect_identical(key, fingerprint(as.integer(z), LD, A, 100, ids))
    expect_false(identical(key, fingerprint(c(1, 3, 2, 4), LD, A, 100, ids)))
    other_ld <- LD; other_ld[1, 2] <- other_ld[2, 1] <- .3
    expect_false(identical(key, fingerprint(z, other_ld, A, 100, ids)))
    expect_false(identical(key, fingerprint(z, LD, 1 - A, 100, ids)))
    expect_false(identical(key, fingerprint(z, LD, A, 200, ids)))
    expect_false(identical(key, fingerprint(z, LD, A, 100, rev(ids))))
    expect_false(identical(key, fingerprint(z, LD, matrix(A, 2, 2), 100, ids)))
  }
  geno <- list(LD = LD, n = 100, p = 4, annotations_matrix = A, variant_ids = ids)
  pheno <- list(z = z)
  key <- .paintor_fingerprint(z, LD, A, geno$n, ids)
  cached <- list(pip = c(.1, .2, .3, .4), credible_sets = list(), runtime_seconds = 2,
                 additional = list(pooled_runtime_seconds = 4, n_loci_pooled = 2))
  cache <- setNames(list(cached), key)
  expect_identical(run_paintor_region(geno, pheno, .paintor_cache = cache), cached)
  changed <- geno; changed$variant_ids <- rev(ids)
  expect_match(run_paintor_region(changed, pheno, .paintor_cache = cache)$error, "no single-locus fallback")
})

sparsepro_set_fixture <- function(set, variant_ids = c("a", "b"), pip_bytes = NULL) {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  writeLines("# fixture", file.path(path, "sparsepro_zld.py"))
  fn <- run_sparsepro
  e <- new.env(parent = environment(fn))
  e$system2 <- function(command, args, ...) {
    out <- args[match("--save", args) + 1L]
    write.table(data.frame(id = gsub("\\s+", "_", variant_ids), z = c(3, 1), pip = c(.8, .2)),
      file.path(out, "region.z.pip"), row.names = FALSE, col.names = FALSE, quote = FALSE)
    if (!is.null(pip_bytes)) writeBin(pip_bytes, file.path(out, "region.z.pip"))
    if (is.raw(set)) writeBin(set, file.path(out, "region.z.cs"))
    else writeLines(set, file.path(out, "region.z.cs"))
    "success"
  }
  environment(fn) <- e
  fn(c(3, 1), diag(2), 100, variant_ids = variant_ids, sparsepro_dir = path)
}

test_that("SparsePro rejects damaged sets instead of discarding bad members", {
  valid <- sparsepro_set_fixture("b,a\t0.2,0.8\t-0.3,0.7")
  expect_null(valid$error)
  expect_identical(valid$credible_sets, list(c(1L, 2L)))
  expect_equal(valid$additional$cs_pip, list(c(.8, .2)))
  expect_equal(valid$additional$cs_effect_size, list(c(.7, -.3)))
  for (bad in c("a,missing\t0.8,0.2\t0.7,0.3", "a,a\t0.8,0.2\t0.7,0.3",
                "a,b\t0.8\t0.7,0.3", "a,b\t0.8,NaN\t0.7,0.3",
                "a,b\t0.8,1.2\t0.7,0.3", "a,b\t0.8,0.2\t0.7,Inf")) {
    fit <- sparsepro_set_fixture(bad)
    expect_match(fit$error, "native sets contain invalid")
    expect_true(all(is.na(fit$pip)))
  }
})

test_that("SparsePro rejects corrupted bytes rather than accepting a valid prefix", {
  cs_prefix <- charToRaw("a/b\t0.8/0.2\t0.7/0.3")
  corrupt_cs <- c(cs_prefix, as.raw(0L), charToRaw("JUNK"))
  non_utf8_cs <- c(cs_prefix, as.raw(255L))
  pip_bytes <- c(charToRaw("a 3 0.8\nb 1 0.2"), as.raw(0L), charToRaw("JUNK"))
  for (cs in list(corrupt_cs, non_utf8_cs)) {
    fit <- sparsepro_set_fixture(cs)
    expect_match(fit$error, "NUL, invalid UTF-8")
    expect_true(all(is.na(fit$pip)))
    expect_identical(fit$additional$native_output_evidence$cs$bytes, cs)
  }
  fit <- sparsepro_set_fixture("a/b\t0.8/0.2\t0.7/0.3", pip_bytes = pip_bytes)
  expect_match(fit$error, "NUL, invalid UTF-8")
  expect_true(all(is.na(fit$pip)))
  expect_identical(fit$additional$native_output_evidence$pip$bytes, pip_bytes)
})

test_that("SparsePro parses the upstream slash-separated multi-member sets", {
  fit <- sparsepro_set_fixture("cs\tpip\tbeta\nb/a\t0.2/0.8\t-0.3/0.7")
  expect_null(fit$error)
  expect_identical(fit$credible_sets, list(c(1L, 2L)))
  expect_equal(fit$additional$cs_pip, list(c(.8, .2)))
  expect_equal(fit$additional$cs_effect_size, list(c(.7, -.3)))
  expect_equal(fit$pip, c(.8, .2))
  for (bad in c("a/missing\t0.8/0.2\t0.7/0.3", "a/a\t0.8/0.2\t0.7/0.3",
                "a/b\t0.8\t0.7/0.3", "a/b\t0.8/NaN\t0.7/0.3",
                "a/b\t0.8/1.2\t0.7/0.3", "a/b\t0.8/0.2\t0.7/Inf",
                "a/b/\t0.8/0.2/\t0.7/0.3/", "a/b\t0.8/0.2\t",
                "a/b\t0.8/0.2\t0.7/0.3\tJUNK", "a//b\t0.8/0.2\t0.7/0.3",
                "a,b/\t0.8,0.2\t0.7,0.3", "a,b,\t0.8,0.2,\t0.7,0.3,")) {
    rejected <- sparsepro_set_fixture(bad)
    expect_match(rejected$error, "native sets contain invalid")
    expect_true(all(is.na(rejected$pip)))
    expect_identical(rejected$additional$native_output_evidence$cs$text,
                     paste0(bad, "\n"))
    expect_false(rejected$additional$native_output_evidence$cs$truncated)
  }
})

test_that("SparsePro retains singleton IDs matching a header prefix and rejects ambiguous IDs", {
  singleton <- sparsepro_set_fixture("cs\t0.8\t0.7", c("cs", "b"))
  expect_null(singleton$error)
  expect_identical(singleton$credible_sets, list(1L))
  spaced <- sparsepro_set_fixture("b / a\t0.2 / 0.8\t-0.3 / 0.7")
  expect_null(spaced$error)
  expect_identical(spaced$credible_sets, list(c(1L, 2L)))
  for (ids in list(c("a b", "a_b"), c("a/b", "b"), c("a,b", "b"), c('a"b', "b"))) {
    rejected <- sparsepro_set_fixture("a\t0.8\t0.7", ids)
    expect_match(rejected$error, "unique nonempty normalized IDs")
    expect_true(all(is.na(rejected$pip)))
  }
})

test_that("SparsePro failure evidence is bounded and survives challenge conversion", {
  path <- tempfile(); on.exit(unlink(path))
  writeBin(charToRaw(paste(rep("x", 100L), collapse = "")), path)
  evidence <- .sparsepro_output_evidence(c(cs = path), max_bytes = 16L)
  expect_identical(evidence$cs$captured_bytes, 16L)
  expect_identical(evidence$cs$text, paste(rep("x", 16L), collapse = ""))
  expect_true(evidence$cs$truncated)
  binary <- as.raw(c(97L, 0L, 98L))
  writeBin(binary, path)
  binary_evidence <- .sparsepro_output_evidence(c(cs = path))
  expect_null(binary_evidence$cs$text)
  expect_identical(binary_evidence$cs$bytes, binary)
  expect_identical(binary_evidence$cs$decoding_status, "non_text_bytes")
  expect_false(binary_evidence$cs$truncated)
  unreadable <- .sparsepro_output_evidence(c(cs = tempdir()))
  expect_true(unreadable$cs$present)
  expect_type(unreadable$cs$capture_error, "character")
  failed <- sparsepro_set_fixture("a/missing\t0.8/0.2\t0.7/0.3")
  prediction <- .campaign_fit_to_prediction(failed, list(task_id = "fixture"))
  expect_identical(prediction$status, "failed")
  expect_identical(prediction$native_diagnostics, failed$additional)
  expect_null(prediction[["pip"]])
})
