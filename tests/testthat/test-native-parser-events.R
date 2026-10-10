finemap_parser_fixture <- function(pip = c(.8, .2), ids = c("variant_1", "variant_2"),
                                   sets = c("variant_1", "variant_2"), status = NULL) {
  fn <- run_finemap
  e <- new.env(parent = environment(fn))
  e$system2 <- function(command, args, ...) {
    master <- readLines(args[match("--in-files", args) + 1L])
    paths <- stats::setNames(strsplit(master[2L], ";", fixed = TRUE)[[1L]],
                            strsplit(master[1L], ";", fixed = TRUE)[[1L]])
    expect_identical(read.table(paths[["z"]], header = TRUE)$rsid,
                     c("variant_1", "variant_2"))
    write.table(data.frame(rsid = ids, prob = pip, mean = 0, sd = .1, log10bf = 1),
      paths[["snp"]], row.names = FALSE, quote = FALSE)
    # Intentionally truncated: summing this table would give k=1 probability .01.
    write.table(data.frame(config = "variant_1", prob = .01, k = 1),
      paths[["config"]], row.names = FALSE, quote = FALSE)
    writeLines(c("- Prior-Pr(# of causal SNPs is k) :", " (0 -> 0)", " 1 -> .5", " 2 -> .5",
      "- Post-Pr(# of causal SNPs is k) :", " (0 -> 0)", " 1 -> .8", " 2 -> .2"),
      paste0(paths[["log"]], "_sss"))
    writeLines(c("# Post-Pr(# of causal SNPs is 1) = .8", "#log10bf 2 NA",
      "index cred1 prob1", paste(1:2, sets, c(.8, .2))), paste0(paths[["cred"]], "1"))
    out <- "native fixture"; if (!is.null(status)) attr(out, "status") <- status
    out
  }
  environment(fn) <- e
  fn(c(3, 1), diag(2), 100, variant_ids = c("a b", "a_b"))
}

test_that("FINEMAP uses the complete posterior and reads native set columns", {
  fit <- finemap_parser_fixture()
  expect_null(fit$error)
  expect_equal(fit$pip, c(.8, .2))
  expect_equal(fit$additional$k_posterior, c(`0` = 0, `1` = .8, `2` = .2))
  expect_identical(fit$additional$best_k, 1L)
  expect_equal(fit$credible_sets, list(c(1L, 2L)))
  expect_equal(fit$additional$cs_pip, list(c(.8, .2)))
  expect_identical(fit$additional$variant_ids, c("a b", "a_b"))
  expect_identical(fit$additional$k_posterior_source, "complete_native_log")
})

test_that("FINEMAP rejects damaged native evidence even when output exists", {
  expect_match(finemap_parser_fixture(ids = c("variant_1", "variant_1"))$error, "IDs")
  expect_match(finemap_parser_fixture(pip = c(1.2, .2))$error, "finite probabilities")
  expect_match(finemap_parser_fixture(sets = c("variant_1", "foreign"))$error, "invalid IDs")
  expect_false(is.null(finemap_parser_fixture(status = 1L)$error))
})

r_native_parser_fixture <- function(method, pip = c(.9, .2, .8, .1),
                                   sets = list(c(1L, 2L), 3L)) {
  fn <- get(paste0("run_", method))
  e <- new.env(parent = environment(fn))
  e$requireNamespace <- function(...) TRUE
  e$getExportedValue <- function(...) function(...) {
    if (method == "carma") list(list(PIPs = pip, `Credible set` = list(NULL, sets)))
    else list(pip = pip, sets = sets)
  }
  environment(fn) <- e
  fn(c(4, 1, 3, 0), diag(4), 100)
}

test_that("CARMA preserves native groups and labels their probability construction", {
  fit <- r_native_parser_fixture("carma")
  expect_null(fit$error)
  expect_equal(fit$credible_sets, list(c(1L, 2L), 3L))
  expect_identical(fit$additional$credible_set_semantics, "ld_clustered_pip_mass_atleastone")
  expect_false(fit$additional$nominal_joint_coverage_guarantee)
})

test_that("R native adapters cannot silently clip PIPs or repair set membership", {
  for (method in c("carma", "finimom")) {
    expect_null(r_native_parser_fixture(method)$error)
    expect_match(r_native_parser_fixture(method, pip = c(1.2, .2, .8, .1))$error, "finite probability|finite probabilities")
    expect_match(r_native_parser_fixture(method, pip = c(.9, NaN, .8, .1))$error, "finite probability|finite probabilities")
    for (sets in list(list(c(1L, 1L)), list(c(1L, 5L)), list(c(1, 1.5))))
      expect_match(r_native_parser_fixture(method, sets = sets)$error, "invalid native set")
  }
})
