# Adapter contract fixtures, not simulated upstream MCMC outputs or a native fit.
native_sbayesrc_fixture <- function(chromosomes = c(1, 1, 2, 2)) {
  root <- tempfile("native-sbayesrc-")
  dir.create(root)
  ld <- file.path(root, "ld panel")
  dir.create(ld)
  p <- length(chromosomes)
  ids <- paste0("rs", seq_len(p))
  gwas <- data.frame(SNP = ids, A1 = "A", A2 = "G", freq = .3,
                     b = seq_len(p) / 100, se = .05, p = .5, N = 500)
  panel <- data.frame(Chrom = chromosomes, ID = ids, Index = seq_len(p) - 1L,
                      GenPos = 0, PhysPos = seq_len(p) * 100,
                      A1 = "A", A2 = "G", A1Freq = .3, N = 500,
                      Block = match(chromosomes, unique(chromosomes)))
  blocks <- do.call(rbind, lapply(unique(panel$Block), function(b) {
    rows <- which(panel$Block == b)
    data.frame(Block = b, Chrom = chromosomes[rows[1]], StartSnpIdx = min(rows) - 1L,
                StartSnpID = ids[min(rows)], EndSnpIdx = max(rows) - 1L,
                EndSnpID = ids[max(rows)], NumSnps = length(rows))
  }))
  annotation <- data.frame(SNP = rev(ids), Intercept = 1, category = seq_len(p) %% 2)
  writer <- function(d, path) write.table(d, path, quote = FALSE, row.names = FALSE)
  gf <- file.path(root, "prepared gwas.ma")
  af <- file.path(root, "prepared annotation.txt")
  writer(gwas, gf); writer(panel, file.path(ld, "snp.info"))
  writer(blocks, file.path(ld, "ldm.info")); writer(annotation, af)
  # A nonempty placeholder tests presence checks only. It is not valid binary
  # eigen-LD and these fixtures are used only for dry-run/import validation.
  for (b in blocks$Block) writeBin(as.raw(1:8), file.path(ld, paste0("block", b, ".eigen.bin")))
  provenance <- as.list(setNames(rep("fixture; no empirical inference", 7),
                                 c("gwas_source", "ld_source", "annotation_source",
                                   "genome_build", "ancestry", "panel_id", "preprocessing")))
  result <- data.frame(Name = rev(ids), A1 = "A", A2 = "G",
                        A1Effect = rev(seq_len(p) / 10), SE = .02,
                        PIP = rev(seq_len(p) / (p + 1)))
  rf <- file.path(root, "fixture.snpRes")
  writer(result, rf)
  list(root = root, gwas_file = gf, ld_dir = ld, annotation_file = af,
        output_dir = file.path(root, "new outputs"), provenance = provenance,
        ids = ids, gwas = gwas, annotation = annotation, result = result,
        result_file = rf, writer = writer)
}

native_fixture_args <- function(f) {
  f[c("gwas_file", "ld_dir", "annotation_file", "output_dir", "provenance")]
}

test_that("native plan retains the whole frozen panel when importing a subset", {
  f <- native_sbayesrc_fixture()
  plan <- do.call(run_sbayesrc_native, c(native_fixture_args(f), list(
    scale = "subset", variant_ids = c("rs4", "rs1"), dry_run = TRUE)))
  expect_false(plan$native_executed)
  expect_false(dir.exists(f$output_dir))
  expect_equal(plan$manifest$modeled_variant_count, 4L)
  expect_equal(plan$manifest$requested_variant_count, 2L)
  expect_equal(plan$manifest$effect_scale, "original_per_GWAS_A1_allele")
  expect_equal(sum(plan$args == "--gwas-summary"), 1L)
  expect_equal(sum(plan$args == "--ldm-eigen"), 1L)
  expect_equal(plan$args[match("--sbayes", plan$args) + 1L], "RC")
  expect_equal(plan$manifest$checksum_algorithm, "md5")
  expect_equal(nrow(plan$manifest$input_checksums), 6L)
  expect_equal(plan$manifest$provenance, f$provenance)
  expect_false(any(c("--help", "--version", "--extract") %in% plan$args))
})

test_that("native scope, identity and preprocessing declarations are enforced", {
  f <- native_sbayesrc_fixture()
  args <- native_fixture_args(f)
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "genome_wide", dry_run = TRUE))),
               "all 22 autosomes")
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "chromosome", dry_run = TRUE))),
               "exactly one chromosome")
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", variant_ids = c("rs1", "rs1"),
                                                        dry_run = TRUE))), "unique nonempty subset")
  args$provenance$genome_build <- NULL
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", dry_run = TRUE))), "provenance")
  all_autosomes <- native_sbayesrc_fixture(1:22)
  plan <- do.call(run_sbayesrc_native, c(native_fixture_args(all_autosomes),
                                        list(scale = "genome_wide", dry_run = TRUE)))
  expect_equal(plan$manifest$scale, "genome_wide")
  expect_equal(plan$manifest$modeled_variant_count, 22L)
  # An all-autosome declaration is not a claim of high-density coverage.
})

test_that("native adapter rejects stale outputs and invalid MCMC controls", {
  f <- native_sbayesrc_fixture()
  args <- native_fixture_args(f)
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", n_iter = 50,
                                                        burn_in = 50, dry_run = TRUE))), "retained")
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", seed = .5,
                                                        dry_run = TRUE))), "seed")
  trace <- do.call(run_sbayesrc_native, c(args, list(scale = "subset", save_mcmc = TRUE,
                                                   thin = 5L, dry_run = TRUE)))
  expect_true("--write-mcmc-txt" %in% trace$args)
  expect_equal(trace$args[match("--thin", trace$args) + 1L], "5")
  expect_true(trace$manifest$params$save_mcmc)
  no_text <- do.call(run_sbayesrc_native, c(args, list(scale = "subset", thin = 5L,
                                                     dry_run = TRUE)))
  expect_equal(no_text$args[match("--thin", no_text$args) + 1L], "5")
  expect_false("--write-mcmc-txt" %in% no_text$args)
  core <- do.call(run_sbayesrc_native, c(args, list(scale = "subset", save_mcmc = TRUE,
                                                  mcmc_scope = "core", dry_run = TRUE)))
  expect_equal(core$manifest$text_retention$discarded_arrays,
               paste0(c("PIP", "PEP", paste0("DeltaPi", 1:5)), ".mcmcsamples.txt"))
  expect_equal(core$manifest$text_retention$sink_target, "/dev/null")
  expect_false(core$manifest$text_retention$formatting_cpu_avoided)
  expect_false(dir.exists(f$output_dir))
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", mcmc_scope = "core",
                                                        dry_run = TRUE))), "requires save_mcmc")
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", save_mcmc = TRUE,
    mcmc_scope = "core", expected_version = "9.9.9", dry_run = TRUE))), "audited only")
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", thin = 0,
                                                        dry_run = TRUE))), "thin")
  dir.create(f$output_dir)
  writeLines("stale", file.path(f$output_dir, ".old-fit"))
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", dry_run = TRUE))), "new or empty")
})

test_that("native core retention uses only the audited sinks in a fresh directory", {
  prefix <- file.path(tempdir(), paste0("native-retention-", sample.int(1e8, 1)))
  if (!identical(unname(Sys.info()["sysname"]), "Linux")) {
    expect_error(.sbayesrc_native_prepare_text_sinks(prefix), "requires Linux")
    expect_false(dir.exists(paste0(prefix, ".mcmcsamples")))
    return(invisible(NULL))
  }
  sinks <- .sbayesrc_native_prepare_text_sinks(prefix)
  on.exit(unlink(paste0(prefix, ".mcmcsamples"), recursive = TRUE), add = TRUE)
  expect_equal(basename(sinks$path), .sbayesrc_native_text_sinks())
  expect_true(all(Sys.readlink(sinks$path) == "/dev/null"))
  expect_false(file.exists(file.path(paste0(prefix, ".mcmcsamples"), "CoreParameters.mcmcsamples.txt")))
  expect_error(.sbayesrc_native_prepare_text_sinks(prefix), "must be new")
})

test_that("native GWAS and annotations must cover the exact native panel", {
  f <- native_sbayesrc_fixture()
  args <- native_fixture_args(f)
  f$writer(f$gwas[-1, ], f$gwas_file)
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", dry_run = TRUE))), "match exactly")
  f$writer(f$gwas, f$gwas_file)
  bad <- f$gwas; bad$A2[1] <- "T"
  f$writer(bad, f$gwas_file)
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", dry_run = TRUE))), "incompatible")
  f$writer(f$gwas, f$gwas_file)
  bad <- f$annotation; bad$Intercept[1] <- 0
  f$writer(bad, f$annotation_file)
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", dry_run = TRUE))), "unit Intercept")
  f$writer(f$annotation, f$annotation_file)
  unlink(file.path(f$ld_dir, "block1.eigen.bin"))
  expect_error(do.call(run_sbayesrc_native, c(args, list(scale = "subset", dry_run = TRUE))), "nonempty")
})

test_that("native result import aligns IDs and orients effects to prepared GWAS A1", {
  f <- native_sbayesrc_fixture()
  reference <- .sbayesrc_native_validate_inputs(f$gwas_file, f$ld_dir, f$annotation_file, "subset")$reference
  bad_order <- f$result
  bad_order$A1[1] <- "G"; bad_order$A2[1] <- "A"
  f$writer(bad_order, f$result_file)
  got <- .sbayesrc_native_import(f$result_file, reference, c("rs4", "rs1"))
  expect_equal(got$pip, c(.8, .2))
  expect_equal(got$beta, c(-.4, .1))
  expect_equal(got$se, c(.02, .02))
  expect_equal(got$effect_allele, c("A", "A"))
})

test_that("native full-panel output is validated before selecting common SNPs", {
  f <- native_sbayesrc_fixture()
  reference <- .sbayesrc_native_validate_inputs(f$gwas_file, f$ld_dir, f$annotation_file, "subset")$reference
  bad <- f$result; bad$PIP[bad$Name == "rs2"] <- 1.01
  f$writer(bad, f$result_file)
  expect_error(.sbayesrc_native_import(f$result_file, reference, "rs1"), "invalid PIPs")
  bad <- f$result; bad$Name[1] <- bad$Name[2]
  f$writer(bad, f$result_file)
  expect_error(.sbayesrc_native_import(f$result_file, reference, "rs1"), "complete frozen")
  f$writer(f$result[-1, ], f$result_file)
  expect_error(.sbayesrc_native_import(f$result_file, reference, "rs1"), "complete frozen")
})

test_that("fingerprints include changed native LD bytes", {
  f <- native_sbayesrc_fixture()
  args <- c(native_fixture_args(f), list(scale = "subset", dry_run = TRUE))
  first <- do.call(run_sbayesrc_native, args)
  writeBin(as.raw(2:9), file.path(f$ld_dir, "block1.eigen.bin"))
  second <- do.call(run_sbayesrc_native, args)
  expect_false(identical(first$manifest$input_checksums, second$manifest$input_checksums))
})

test_that("native positional LD metadata and block boundaries are validated", {
  f <- native_sbayesrc_fixture()
  args <- c(native_fixture_args(f), list(scale = "subset", dry_run = TRUE))
  panel_file <- file.path(f$ld_dir, "snp.info")
  panel <- read.table(panel_file, header = TRUE)
  f$writer(panel[, c(2, 1, 3:10)], panel_file)
  expect_error(do.call(run_sbayesrc_native, args), "ten-column order")
  f$writer(panel, panel_file)
  blocks_file <- file.path(f$ld_dir, "ldm.info")
  blocks <- read.table(blocks_file, header = TRUE)
  blocks$EndSnpIdx[1] <- blocks$EndSnpIdx[1] + 1L
  f$writer(blocks, blocks_file)
  expect_error(do.call(run_sbayesrc_native, args), "boundaries")
})

test_that("an unavailable real executable is an explicit failure", {
  f <- native_sbayesrc_fixture()
  expect_error(do.call(run_sbayesrc_native, c(native_fixture_args(f), list(
    scale = "subset", gctb_path = file.path(f$root, "missing-gctb")))), "executable is unavailable")
  expect_false(dir.exists(f$output_dir))
})
