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
  fn(c(3, 1), diag(2), 100, variant_ids = c("a b", "a_b"), n_causal = 2L)
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
  expect_identical(fit$additional$posterior_conditioning, "at_least_one_modeled_cause")
  expect_identical(fit$additional$native_null_support, "excluded_by_prior")
  expect_identical(fit$additional$null_probability, 0)
})

test_that("FINEMAP validates actual rounded complete posteriors without renormalization", {
  path <- tempfile(); on.exit(unlink(path), add = TRUE)
  write_log <- function(prob, sizes = 0:5) writeLines(c("- Post-Pr(# of causal SNPs is k) :",
    paste(sizes, "->", prob), "- Run time : 0 seconds"), path)
  for (prob in list(c("0", ".548", ".404", ".0476", "0", "0"),
                   c("0", ".536", ".32", ".145", "0", "0"))) {
    write_log(prob)
    got <- .finemap_model_size_probabilities(path, 5L)
    expect_equal(as.numeric(got), as.numeric(prob))
    expect_false(attr(got, "reported_rounding")$renormalized)
    expect_lte(attr(got, "reported_rounding")$sum_interval["lower"], 1)
    expect_gte(attr(got, "reported_rounding")$sum_interval["upper"], 1)
  }
  write_log(c("0", ".6", ".3", ".2", "0", "0"))
  expect_null(.finemap_model_size_probabilities(path, 5L))
  write_log(c("0", ".8", ".2"), 0:2)
  expect_null(.finemap_model_size_probabilities(path, 5L))
  write_log(c("0", ".8", ".2", "0", "0", "0"), c(0L, 1L, 1L, 3L, 4L, 5L))
  expect_null(.finemap_model_size_probabilities(path, 5L))
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

test_that("audited PAINTOR enumeration advances past null and counts nonempty support", {
  # Independently follow the upstream next_combo state machine for small p.
  # Do not use the wrapper's binomial-count formula to generate configurations.
  enumerate <- function(p) {
    state <- c(0L, 0L)
    support <- list()
    repeat {
      advanced <- FALSE
      for (i in 2:1) {
        state[i] <- state[i] + 1L
        if (state[i] <= p + 1L - (3L - i)) {
          if (i == 1L) state[2L] <- state[1L] + 1L
          advanced <- TRUE
          break
        }
      }
      if (!advanced) return(support)
      support[[length(support) + 1L]] <- state[state > 0L]
    }
  }
  for (p in c(2L, 3L, 5L, 12L)) {
    support <- enumerate(p)
    expected <- c(as.list(seq_len(p)), utils::combn(p, 2L, simplify = FALSE))
    keys <- function(x) sort(vapply(x, paste, "", collapse = ","))
    expect_identical(keys(support), keys(expected))
    expect_false(any(lengths(support) == 0L))
    expect_equal(.paintor_probability_contract(p, 2, FALSE)$configuration_count,
                 length(support))
  }
  expect_identical(.paintor_probability_contract(3, 2, FALSE)$posterior_conditioning,
                   "at_least_one_modeled_causal_variant")
  expect_identical(.paintor_probability_contract(3, 2, TRUE)$posterior_conditioning,
                   "null_including_model_space")
  expect_true(is.na(.paintor_probability_contract(3, 2, TRUE)$configuration_count))
})

paintor_conditioning_fixture <- function(pooled = FALSE, mcmc = FALSE) {
  fn <- if (pooled) run_paintor_scenario_setup else run_paintor
  e <- new.env(parent = environment(fn))
  e$system2 <- function(command, args, ...) {
    arg <- function(name) args[match(name, args) + 1L]
    for (locus in readLines(arg("-input"))) {
      ids <- read.table(file.path(arg("-in"), locus), header = TRUE)$rsid
      write.table(data.frame(rsid = ids, Posterior_Prob = c(.8, .2)),
        file.path(arg("-out"), paste0(locus, ".results")),
        row.names = FALSE, quote = FALSE)
    }
    "native parser fixture"
  }
  environment(fn) <- e
  if (!pooled) return(list(fn(c(3, 1), diag(2), mcmc = mcmc)))
  genotypes <- list(list(LD = diag(2), n = 100, variant_ids = c("a", "b")),
                    list(LD = diag(2), n = 100, variant_ids = c("c", "d")))
  regions <- list(list(z = c(3, 1)), list(z = c(2, 0)))
  unname(fn(genotypes, regions, list(mcmc = mcmc))$.paintor_cache)
}

test_that("PAINTOR labels conditioning for successful single and pooled fits without null fabrication", {
  for (pooled in c(FALSE, TRUE)) for (mcmc in c(FALSE, TRUE)) {
    fits <- paintor_conditioning_fixture(pooled, mcmc)
    expect_length(fits, if (pooled) 2L else 1L)
    for (fit in fits) {
      expect_null(fit$error)
      expect_equal(fit$pip, c(.8, .2))
      expect_identical(fit$additional$posterior_conditioning,
        if (mcmc) "null_including_model_space" else "at_least_one_modeled_causal_variant")
      expect_identical(fit$additional$null_model_in_support, mcmc)
      expect_true(fit$additional$null_probability_not_provided)
      expect_null(fit$additional$null_prob)
      expect_null(fit$additional$null_probability)
      expect_equal(fit$additional$configuration_count, if (mcmc) NA_real_ else 3)
    }
  }
})
