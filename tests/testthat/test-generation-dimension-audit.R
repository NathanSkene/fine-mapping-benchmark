generation_audit_fixture <- function(reference = FALSE, requested = 4L) {
  raw <- matrix(c(0, 1, 2, 0, 1, 2, 0, 0, 1, 1, 2, 2, 2, 1, 0, 2, 1, 0), 6L, 3L)
  X <- scale(raw)
  ids <- paste0("variant", 1:3)
  geno <- list(X = X, X_raw = raw, n = 6L, p = 3L, variant_ids = ids,
               region_id = 1L, LD = cor(X), LD_true = cor(X))
  if (reference) {
    geno$X_ref_raw <- raw[1:4, , drop = FALSE]
    geno$X_ref <- scale(geno$X_ref_raw)
    geno$n_ref <- 4L
    geno$LD <- cor(geno$X_ref)
  }
  list(genotypes = list(geno), params = list(p = requested, p_actual = 3L, n = 6L,
       n_ref = if (reference) 4L else NULL),
       scenarios = list(list(regions = list(list(y = seq_len(6L), z = c(1, 2, 3),
         truth = list(beta_true = c(0, .1, 0)))))))
}

generation_audit_known_sha <- "651dbb1640b63187028c12406a82801a676f10be1f9a52b635aedad5a49be139"
audit_generation_fixture <- function(sim = generation_audit_fixture(),
    requested = sim$params$p, warnings = character(), sha = generation_audit_known_sha,
    n = 6L, n_ref = sim$params$n_ref) {
  .campaign_generation_dimension_audit(sim, requested, n, n_ref, sha, warnings)
}

test_that("known upstream post-simulation filtering has a nominal versus retained audit", {
  sim <- generation_audit_fixture(TRUE)
  set.seed(992L)
  rng <- .Random.seed
  payload <- serialize(sim, NULL, version = 2L)
  audit <- audit_generation_fixture(sim)
  expect_identical(serialize(sim, NULL, version = 2L), payload)
  expect_identical(.Random.seed, rng)
  expect_true(audit$known_post_simulation_filter)
  expect_identical(audit$source_file_sha256, generation_audit_known_sha)
  expect_identical(audit$loci$requested_p, 4L)
  expect_identical(audit$loci$realized_p, 3L)
  expect_identical(audit$loci$dropped_p, 1L)
  expect_identical(audit$loci$retained_fraction, .75)
  expect_identical(audit$loci$dropped_fraction, .25)
  expect_identical(audit$loci$n_gwas, 6L)
  expect_identical(audit$loci$n_reference, 4L)
  expect_match(audit$variant_universe_semantics, "retained_polymorphic")
  expect_match(audit$post_simulation_filter_basis, "known_f053")
  enriched <- sim; enriched$params$generation_dimension_audit <- audit
  expect_identical(enriched$genotypes, sim$genotypes)
  expect_identical(enriched$scenarios, sim$scenarios)
})

test_that("warning capture preserves every message and native RNG consumption", {
  make <- function() {
    warning("original first warning", call. = FALSE)
    sim <- generation_audit_fixture()
    sim$scenarios[[1L]]$regions[[1L]]$y <- rnorm(6L)
    warning("original second warning", call. = FALSE)
    warning("original first warning", call. = FALSE)
    sim
  }
  set.seed(821L)
  expected <- suppressWarnings(make())
  expected_rng <- .Random.seed
  set.seed(821L)
  delivered <- character()
  captured <- withCallingHandlers(.campaign_capture_generation(make), warning = function(w) {
    delivered <<- c(delivered, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_identical(captured$value, expected)
  expect_identical(.Random.seed, expected_rng)
  expect_identical(delivered, c("original first warning", "original second warning", "original first warning"))
  expect_identical(captured$warnings, delivered)
  audit <- audit_generation_fixture(captured$value, warnings = captured$warnings)
  expect_identical(audit$warnings$message, delivered)
  expect_identical(audit$warnings$classification, rep("other_original_warning_preserved", 3L))
})

test_that("heterogeneous loci retain their own requested and realized dimensions", {
  sim <- generation_audit_fixture()
  second <- sim$genotypes[[1L]]; second$region_id <- 2L
  ids <- second$variant_ids
  for (field in c("X", "X_raw")) colnames(second[[field]]) <- ids
  for (field in c("LD", "LD_true")) dimnames(second[[field]]) <- list(ids, ids)
  sim$genotypes[[2L]] <- second
  sim$params$p <- c(4L, 3L); sim$params$p_actual <- c(3L, 3L)
  audit <- audit_generation_fixture(sim)
  expect_identical(audit$loci$region_id, 1:2)
  expect_identical(audit$loci$requested_p, c(4L, 3L))
  expect_identical(audit$loci$realized_p, c(3L, 3L))
  expect_identical(audit$loci$dropped_p, c(1L, 0L))
  expect_identical(audit$loci$retained_fraction, c(.75, 1))
  expect_identical(audit$loci$dropped_fraction, c(.25, 0))
})

test_that("source undersupply and unverified dimension losses remain fatal", {
  shortage <- "Region 1: requested 4 SNPs but only 3 variants available after MAF filtering. Using all 3."
  expect_error(audit_generation_fixture(warnings = shortage), "VCF undersupply remains fatal")
  expect_error(audit_generation_fixture(generation_audit_fixture(requested = 3L), warnings = shortage),
               "VCF undersupply remains fatal")
  cap <- paste0("The bundled sim1000G example VCF contains ~567 variants. ",
    "Capping p at 500 for regions where p > 500. Use custom VCF files for larger regions.")
  expect_error(audit_generation_fixture(warnings = cap), "VCF undersupply remains fatal")
  expect_error(audit_generation_fixture(generation_audit_fixture(requested = 3L), warnings = cap),
               "VCF undersupply remains fatal")
  expect_error(audit_generation_fixture(sha = strrep("a", 64L)), "Unknown original-generator origin")
  exact <- audit_generation_fixture(generation_audit_fixture(requested = 3L), sha = strrep("a", 64L))
  expect_false(exact$known_post_simulation_filter)
  expect_identical(exact$loci$dropped_p, 0L)
  expect_true(is.na(exact$loci$n_reference))
  expect_error(audit_generation_fixture(sha = "unverified"), "source-file SHA256")
  expect_error(audit_generation_fixture(warnings = NA_character_), "warnings")
})

test_that("requested, realized and parameter dimensions must agree exactly", {
  sim <- generation_audit_fixture()
  expect_error(audit_generation_fixture(requested = 4.5), "requested variant")
  expect_error(audit_generation_fixture(requested = 0L), "requested variant")
  extra <- sim; extra$params$p <- 2L
  expect_error(audit_generation_fixture(extra), "extra variants")
  bad <- sim; bad$genotypes[[1L]]$p <- 3.5
  expect_error(audit_generation_fixture(bad), "genotype p")
  bad <- sim; bad$genotypes[[1L]]$p <- 4L
  expect_error(audit_generation_fixture(bad), "Contradictory generation genotype")
  bad <- sim; bad$params$p_actual <- 2L
  expect_error(audit_generation_fixture(bad), "params.p_actual")
  bad <- sim; bad$params$p_actual <- 3.5
  expect_error(audit_generation_fixture(bad), "params.p_actual")
  bad <- sim; bad$params$p <- 5L
  expect_error(audit_generation_fixture(bad, requested = 4L), "design parameters")
  bad <- sim; bad$params$n <- 5L
  expect_error(audit_generation_fixture(bad), "design parameters")
  expect_error(audit_generation_fixture(n = 6.5), "GWAS sample")
  bad <- sim; bad$genotypes <- c(bad$genotypes, bad$genotypes)
  expect_error(audit_generation_fixture(bad), "requested variant")
  bad <- sim; bad$genotypes[[1L]]$region_id <- 2L
  expect_error(audit_generation_fixture(bad), "region ID")
})

test_that("finite variable genotype columns and GWAS sample sizes are required", {
  sim <- generation_audit_fixture()
  for (field in c("X", "X_raw")) {
    bad <- sim; bad$genotypes[[1L]][[field]][1L, 1L] <- NA_real_
    expect_error(audit_generation_fixture(bad), "dimensions or values")
    bad <- sim; bad$genotypes[[1L]][[field]][, 1L] <- 1
    expect_error(audit_generation_fixture(bad), "nonvariable")
    bad <- sim; bad$genotypes[[1L]][[field]] <- bad$genotypes[[1L]][[field]][1:5, , drop = FALSE]
    expect_error(audit_generation_fixture(bad), "dimensions or values")
  }
  bad <- sim; bad$genotypes[[1L]]$X_raw <- NULL
  expect_error(audit_generation_fixture(bad), "X_raw")
  bad <- sim; bad$genotypes[[1L]]$n <- 5L
  expect_error(audit_generation_fixture(bad), "genotype dimensions")
  bad <- sim; bad$genotypes[[1L]]$X <- matrix(numeric(), 6L, 0L)
  expect_error(audit_generation_fixture(bad), "realized variant")
})

test_that("reference panels, variant identities and LD matrices cannot contradict the audit", {
  sim <- generation_audit_fixture(TRUE)
  for (field in c("X_ref", "X_ref_raw")) {
    bad <- sim; bad$genotypes[[1L]][[field]][, 1L] <- 0
    expect_error(audit_generation_fixture(bad), "nonvariable")
    bad <- sim; bad$genotypes[[1L]][[field]] <- bad$genotypes[[1L]][[field]][1:3, , drop = FALSE]
    expect_error(audit_generation_fixture(bad), "dimensions or values")
  }
  bad <- sim; bad$genotypes[[1L]]$n_ref <- 3L
  expect_error(audit_generation_fixture(bad), "reference dimensions")
  bad <- sim; bad$params$n_ref <- 3L
  expect_error(audit_generation_fixture(bad, n_ref = 4L), "reference sample size")
  expect_error(audit_generation_fixture(sim, n_ref = NULL), "Unexpected generation reference")
  bad <- sim; bad$genotypes[[1L]]$variant_ids <- c("a", "a", "b")
  expect_error(audit_generation_fixture(bad), "variant IDs")
  bad <- sim; bad$genotypes[[1L]]$variant_ids <- c("a", "b")
  expect_error(audit_generation_fixture(bad), "variant IDs")
  bad <- sim; colnames(bad$genotypes[[1L]]$X) <- c("wrong", "order", "ids")
  expect_error(audit_generation_fixture(bad), "variant IDs disagree")
  for (field in c("LD", "LD_true")) {
    bad <- sim; bad$genotypes[[1L]][[field]] <- diag(2L)
    expect_error(audit_generation_fixture(bad), "dimensions or values")
    bad <- sim; bad$genotypes[[1L]][[field]][1L, 2L] <- Inf
    expect_error(audit_generation_fixture(bad), "dimensions or values")
    bad <- sim; rownames(bad$genotypes[[1L]][[field]]) <- c("wrong", "order", "ids")
    expect_error(audit_generation_fixture(bad), "variant IDs disagree")
  }
})
