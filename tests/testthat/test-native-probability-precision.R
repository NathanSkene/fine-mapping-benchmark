precision_validate <- function(pred, input) .challenge_validate_prediction(pred, input,
  allow_native_probability_precision = TRUE)

precision_prediction_fixture <- function(tokens, ids = paste0("v", seq_along(tokens))) {
  basis <- .finemap_probability_writer_contract()
  prior <- posterior <- c(`0` = 0, `1` = 1)
  native_ids <- paste0("variant_", seq_along(tokens))
  precision <- .finemap_probability_precision(tokens, ids, native_ids,
    basis$executable_sha256, prior, posterior)
  additional <- list(variant_ids = ids, native_variant_ids = native_ids,
    native_executable_sha256 = basis$executable_sha256,
    snp_probability_tokens = tokens, native_probability_precision = precision,
    native_null_support = "excluded_by_prior",
    posterior_conditioning = "at_least_one_modeled_cause",
    null_probability_event = basis$null_event, k_prior = prior,
    k_posterior = posterior, null_probability = 0)
  fit <- list(method = "finemap", pip = as.numeric(tokens), credible_sets = list(),
              runtime_seconds = 0, additional = additional)
  input <- list(task_id = "precision_fixture", variant_id = ids, z = rep(0, length(ids)))
  list(fit = fit, input = input, prediction = .campaign_fit_to_prediction(fit, input))
}

test_that("literal decimal envelopes respect adjacent powers and directed rounding", {
  bounds <- .finemap_probability_token_bounds(c("1", "0.1", "0.01", "1e-05", "0.100001", "0"))
  expect_equal(bounds$lower, c(1 - 1e-6, .1 - 1e-7, .01 - 1e-8, 1e-5 - 1e-11, .1, 0), tolerance = 1e-15)
  expect_equal(bounds$upper, c(1, .1 + 1e-6, .01 + 1e-7, 1e-5 + 1e-10, .100002, 0), tolerance = 1e-15)
  # The power's lower spacing is one tenth of its upper spacing, not a
  # symmetric interval. These independently chosen negative examples expose
  # a too-wide lower envelope that could accept an impossible event.
  expect_gt(bounds$lower[2], .1 - 2e-7)
  expect_lt(bounds$upper[2], .1 + 2e-6)
  expect_equal(.finemap_probability_token_bounds("0.999999")$upper, 1)
  subnormal <- .finemap_probability_token_bounds("4.94066e-324")
  expect_lte(subnormal$lower, 0)
  expect_gte(subnormal$upper, .Machine$double.xmin * .Machine$double.eps)
  for (token in c("NA", "NaN", "Inf", "-0.1", "1.01", ".1", "0e-9", "0.1234567"))
    expect_error(.finemap_probability_token_bounds(token), "native probability|writer precision")
})

test_that("actual archived FINEMAP deficits validate without changing any point values", {
  # Exact literals from the two failing author-binary outputs, in input order.
  # They are retained regression evidence, not benchmark performance results.
  actual <- list(
    c("0.0279124", "0.0298004", "0.041207", "0.0310035", "0.0303409", "0.029256",
      "0.026872", "0.0268027", "0.0389764", "0.0994689", "0.0279357", "0.0552624",
      "0.034529", "0.0364924", "0.0271742", "0.0350059", "0.029788", "0.0459373",
      "0.0289729", "0.0423501", "0.0627792", "0.117508", "0.0371783", "0.0374462"),
    c("0.03881", "0.0329409", "0.0360661", "0.0862981", "0.0464306", "0.0329839",
      "0.0427008", "0.0374994", "0.0436348", "0.0320576", "0.0317941", "0.0321575",
      "0.0319227", "0.0420286", "0.0378533", "0.101944", "0.0406588", "0.0316427",
      "0.0326706", "0.0482516", "0.0436726", "0.0322112", "0.0317296", "0.03204"))
  for (i in seq_along(actual)) {
    x <- precision_prediction_fixture(actual[[i]])
    raw <- x$prediction$pip
    expect_lt(abs((1 - sum(raw)) - c(2e-7, 5e-7)[i]), 1e-14)
    checked <- precision_validate(x$prediction, x$input)
    expect_identical(checked$pip, raw)
    expect_identical(checked$null_probability, 0)
    expect_identical(checked$native_diagnostics$snp_probability_tokens, actual[[i]])
    expect_false(checked$probability_target_aligned)
    expect_match(checked$probability_score_interpretation, "not unconditional calibration")
    expect_identical(checked$null_probability_event,
      "empty_configuration_in_nonnull_conditioned_support")
    before <- mean((raw - c(TRUE, rep(FALSE, length(raw) - 1)))^2)
    after <- mean((checked$pip - c(TRUE, rep(FALSE, length(raw) - 1)))^2)
    expect_identical(after, before)
    # The same archived point vector without verified serialization provenance
    # must retain the old strict union check, rather than a blanket tolerance.
    strict <- x$prediction; strict$native_probability_precision <- NULL
    expect_error(precision_validate(strict, x$input), "union bound")
    shuffled <- x$prediction; ordering <- rev(seq_along(raw))
    shuffled$variant_id <- shuffled$variant_id[ordering]; shuffled$pip <- raw[ordering]
    expect_identical(precision_validate(shuffled, x$input)$pip, raw)
  }
})

test_that("precision provenance cannot admit unrelated or contradictory predictions", {
  x <- precision_prediction_fixture(c("0.499999", "0.5"))
  expect_identical(precision_validate(x$prediction, x$input)$pip, x$fit$pip)
  for (field in c("tokens", "variant_ids", "native_variant_ids", "probability_upper", "record_sha256")) {
    bad <- x$prediction
    if (field == "probability_upper") bad$native_probability_precision[[field]][1] <- 1
    else bad$native_probability_precision[[field]][1] <- "changed"
    expect_error(precision_validate(bad, x$input), "precision record|precision basis")
  }
  for (field in c("executable_sha256", "writer", "significant_digits", "rounding_mode", "null_basis")) {
    bad <- x$prediction; bad$native_probability_precision[[field]] <- "changed"
    expect_error(precision_validate(bad, x$input), "precision basis")
  }
  bad <- x$prediction; bad$pip[1] <- bad$pip[1] + 1e-7
  expect_error(precision_validate(bad, x$input), "tokens, PIPs")
  bad <- x$prediction; bad$implementation <- "another_predictor"
  expect_error(precision_validate(bad, x$input), "precision basis")
  bad <- x$prediction; bad$native_diagnostics$native_executable_sha256 <- paste(rep("0", 64), collapse = "")
  expect_error(precision_validate(bad, x$input), "diagnostics or null basis")
  bad <- x$prediction; bad$null_probability <- 1e-6
  expect_error(precision_validate(bad, x$input), "diagnostics or null basis")
  bad <- x$prediction; bad$native_diagnostics$k_prior["0"] <- .01
  expect_error(precision_validate(bad, x$input), "diagnostics or null basis")
  bad <- x$prediction; bad$native_diagnostics <- NULL
  expect_error(precision_validate(bad, x$input), "diagnostics or null basis")
  for (field in c("posterior_conditioning", "native_null_support", "null_probability_event",
                   "probability_target_aligned", "null_model_in_support")) {
    bad <- x$prediction; bad[[field]] <- "changed"
    expect_error(precision_validate(bad, x$input), "diagnostics or null basis")
  }
  impossible <- precision_prediction_fixture(c("0.49", "0.5"))
  expect_error(precision_validate(impossible$prediction, impossible$input), "union bound")
  strict <- x$prediction; strict$native_probability_precision <- NULL
  strict$implementation <- "another_predictor"
  expect_error(precision_validate(strict, x$input), "union bound")
  expect_null(.finemap_probability_precision(c("0.499999", "0.5"), x$input$variant_id,
    c("variant_1", "variant_2"), paste(rep("0", 64), collapse = ""),
    c(`0` = 0, `1` = 1), c(`0` = 0, `1` = 1)))
  expect_null(.finemap_probability_precision(c("0.499999", "0.5"), x$input$variant_id,
    c("variant_1", "variant_2"), .finemap_probability_writer_contract()$executable_sha256,
    c(`0` = .1, `1` = .9), c(`0` = .1, `1` = .9)))
})

test_that("public submissions remain strict until the assessor explicitly checks custody", {
  panel <- generate_challenge("null", n_replicates = 1L, p = 12L)
  input <- panel$inputs[[1L]]
  x <- precision_prediction_fixture(c("0.499999", "0.5", rep("0", 10)), input$variant_id)
  x$prediction$task_id <- input$task_id
  # A copied or client-created self-consistent record is not authentication.
  expect_error(.challenge_validate_prediction(x$prediction, input), "union bound")
  ran <- run_challenge(panel$inputs, list(client = function(input) x$prediction))
  expect_identical(ran$client[[1L]]$status, "failed")
  expect_match(ran$client[[1L]]$error, "union bound")
  predictions <- list(finemap = setNames(list(x$prediction), input$task_id))
  strict <- evaluate_challenge(panel, predictions)
  expect_false(strict$tasks$success)
  expect_match(strict$tasks$error, "union bound")
  audited <- evaluate_challenge(panel, predictions, allow_native_probability_precision = TRUE)
  expect_true(audited$tasks$success)
  expect_identical(audited$tasks$brier, mean(x$prediction$pip^2))
  # The candidate non-null support's structural P(empty)=0 is scored against
  # a generating-null locus operationally; it is not same-event calibration.
  expect_identical(audited$tasks$null_brier, 1)
  expect_error(evaluate_challenge(panel, predictions,
    allow_native_probability_precision = NA), "assessor option")
})

test_that("FINEMAP parser preserves literal precision and ties envelopes to executable bytes", {
  parser <- function(known = TRUE, changes = FALSE) {
    fn <- run_finemap; e <- new.env(parent = environment(fn)); calls <- 0L
    e$.finemap_binary_sha256 <- function(path) {
      calls <<- calls + 1L
      if (known && (!changes || calls == 1L)) .finemap_probability_writer_contract()$executable_sha256
      else paste(rep("0", 64), collapse = "")
    }
    e$system2 <- function(command, args, ...) {
      master <- readLines(args[match("--in-files", args) + 1L])
      paths <- setNames(strsplit(master[2L], ";", fixed = TRUE)[[1L]],
                        strsplit(master[1L], ";", fixed = TRUE)[[1L]])
      # Deliberately reverse output order, with literal trailing zero omission.
      writeLines(c("rsid prob mean sd log10bf", "variant_2 0.5 0 .1 1",
                   "variant_1 0.499999 0 .1 1"), paths[["snp"]])
      writeLines(c("config prob", "variant_1 1"), paths[["config"]])
      writeLines(c("- Prior-Pr(# of causal SNPs is k) :", "0 -> 0", "1 -> 1",
        "- Post-Pr(# of causal SNPs is k) :", "0 -> 0", "1 -> 1"), paste0(paths[["log"]], "_sss"))
      writeLines(c("index cred1 prob1", "1 variant_1 0.5", "2 variant_2 0.5"), paste0(paths[["cred"]], "1"))
      "native parser fixture"
    }
    environment(fn) <- e
    fn(c(0, 1), diag(2), 100, variant_ids = c("a b", "a_b"), n_causal = 1L)
  }
  fit <- parser()
  expect_null(fit$error)
  expect_identical(fit$pip, c(.499999, .5))
  expect_identical(fit$additional$snp_probability_tokens, c("0.499999", "0.5"))
  input <- list(task_id = "parser", variant_id = c("a b", "a_b"), z = c(0, 1))
  pred <- .campaign_fit_to_prediction(fit, input)
  expect_identical(precision_validate(pred, input)$pip, fit$pip)
  expect_identical(pred$native_diagnostics$snp_probability_tokens, fit$additional$snp_probability_tokens)
  expect_false("configs" %in% names(pred$native_diagnostics))
  unknown <- parser(FALSE)
  expect_null(unknown$error)
  expect_null(unknown$additional$native_probability_precision)
  expect_identical(unknown$pip, fit$pip)
  expect_error(precision_validate(.campaign_fit_to_prediction(unknown, input), input), "union bound")
  expect_match(parser(changes = TRUE)$error, "executable changed")
})
