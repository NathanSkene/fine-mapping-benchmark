#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "../.."))
setwd(repo)
if (dir.exists(".local-lib")) .libPaths(c(normalizePath(".local-lib"), .libPaths()))
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
panel <- generate_challenge(n_replicates = 1, p = 12)
pred <- run_challenge(panel$inputs)
evaluation <- evaluate_challenge(panel, pred)
stopifnot(nrow(evaluation$tasks) == nrow(challenge_scenarios()) * 2,
          all(evaluation$tasks$recall[evaluation$tasks$scenario == "missing_causal"] == 0),
          all(!evaluation$tasks$success[evaluation$tasks$scenario == "meta_missingness"]),
          all(evaluation$tasks$success[evaluation$tasks$scenario == "meta_known_overlap"]))
if (requireNamespace("testthat", quietly = TRUE)) {
  for (f in c("test-challenge.R", "test-evaluate-validity.R", "test-simulation-validity.R",
              "test-sbayesrc-native.R", "test-tuning-contract.R")) {
    testthat::test_file(file.path("tests/testthat", f), stop_on_failure = TRUE, reporter = "summary")
  }
} else stop("Install testthat to run the complete scientific invariant suite")
message("Scientific invariant suite passed")
