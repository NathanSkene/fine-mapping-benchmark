beatrice_event_fixture <- function(mutator = identity, rule = "conditional_union",
                                   diagnostic = TRUE, sets = list(c(0L, 1L))) {
  path <- tempfile(); dir.create(path)
  ids <- paste0("variant_", 1:3)
  write.csv(data.frame(variant_names = ids, variant_index = 0:2,
    pip = c(.7, .6, .1)), file.path(path, "pip.csv"), row.names = FALSE)
  writeLines(vapply(sets, paste, "", collapse = " "), file.path(path, "credible_set.txt"))
  writeLines(vapply(sets, function(set) paste(c(.7, .6, .1)[set + 1L], collapse = " "), ""),
    file.path(path, "conditional_credible_variants_probability.txt"))
  # This posterior is .4{0,1} + .3{0} + .2{1} + .1{2}.
  # CS {0,1} has membership sum1.3 but union probability .9.
  d <- list(format = "beatrice-conditional-credible-set-diagnostics-v1",
    index_base = 0L, variant_count = 3L, variant_names = as.list(ids),
    posterior_memo_format = "beatrice-log-bf-final-prior-v1",
    posterior_support_size = 5L, null_probability = 0,
    null_probability_conditioning = "explored_support_under_frozen_final_prior",
    stopping_rule = rule, sets = lapply(sets, function(set) list(
      set_indices_zero_based = as.list(set), condition_indices_zero_based = list(),
      event = "at_least_one_set_member_is_causal",
      posterior_conditioning = "given_variants_are_causal_within_explored_support",
      conditional_union_probability = .9, target_probability = .95,
      target_met = FALSE, legacy_membership_sum = 1.3,
      selection_probability_threshold = .001, purity_threshold = .5)))
  if (diagnostic) jsonlite::write_json(mutator(d),
    file.path(path, "credible_set_diagnostics.json"), auto_unbox = TRUE)
  list(path = path, ids = ids)
}

test_that("BEATRICE imports event coverage without replacing it by summed marginals", {
  f <- beatrice_event_fixture(); on.exit(unlink(f$path, recursive = TRUE))
  fit <- .fb_joint_parse_output(f$path, f$ids, 3L, require_diagnostics = TRUE)
  expect_equal(fit$pip, c(.7, .6, .1))
  expect_identical(fit$credible_sets, list(c(1L, 2L)))
  d <- fit$cs_diagnostics$sets[[1L]]
  expect_equal(d$conditional_union_probability, .9)
  expect_equal(d$legacy_membership_sum, 1.3)
  expect_false(d$target_met)
  expect_identical(d$set_indices, c(1L, 2L))
  expect_identical(d$condition_indices, integer())
  metadata <- .beatrice_posterior_metadata(fit)
  expect_identical(metadata$credible_set_semantics, "conditional_union_within_explored_support")
  expect_equal(metadata$null_probability, 0)
  expect_true(metadata$null_model_in_support)
  expect_false(metadata$null_probability_not_provided)
  expect_equal(metadata$configuration_count, 5)
  expect_identical(metadata$posterior_conditioning, "within_explored_configuration_support")
})

test_that("BEATRICE separates undeclared external outputs and the legacy mass diagnostic", {
  legacy <- beatrice_event_fixture(rule = "legacy_membership_sum")
  old <- beatrice_event_fixture(diagnostic = FALSE)
  on.exit(unlink(c(legacy$path, old$path), recursive = TRUE))
  fit <- .fb_joint_parse_output(legacy$path, legacy$ids, 3L, TRUE)
  expect_identical(.beatrice_posterior_metadata(fit)$credible_set_semantics,
    "legacy_conditional_membership_mass_candidates")
  external <- .fb_joint_parse_output(old$path, old$ids, 3L)
  expect_identical(.beatrice_posterior_metadata(external)$credible_set_semantics,
    "native_semantics_not_verified")
  expect_null(.beatrice_posterior_metadata(external)$null_probability)
  expect_null(.fb_joint_parse_output(old$path, old$ids, 3L, TRUE))
})

test_that("BEATRICE requires aligned, finite, coherent event diagnostics", {
  damage <- list(
    function(d) { d$variant_names <- rev(d$variant_names); d },
    function(d) { d$variant_count <- 2L; d },
    function(d) { d$index_base <- 1L; d },
    function(d) { d$null_probability <- .4; d },
    function(d) { d$posterior_support_size <- 0L; d },
    function(d) { d$posterior_memo_format <- "unknown"; d },
    function(d) { d$sets[[1L]]$set_indices_zero_based <- list(1L, 0L); d },
    function(d) { d$sets[[1L]]$set_indices_zero_based <- list(0L, 3L); d },
    function(d) { d$sets[[1L]]$condition_indices_zero_based <- list(0L); d },
    function(d) { d$sets[[1L]]$conditional_union_probability <- 1.1; d },
    function(d) { d$sets[[1L]]$conditional_union_probability <- NULL; d },
    function(d) { d$sets[[1L]]$target_met <- TRUE; d },
    function(d) { d$sets[[1L]]$target_met <- 0L; d },
    function(d) { d$sets[[1L]]$legacy_membership_sum <- .8; d },
    function(d) { d$sets[[1L]]$conditional_union_probability <- .5; d },
    function(d) { d$sets[[1L]]$legacy_membership_sum <- 1.4; d },
    function(d) { d$sets[[1L]]$stopping_rule <- "legacy_membership_sum"; d },
    function(d) { d$sets[[1L]]$event <- "all_causes_are_in_set"; d },
    function(d) { d$sets[[1L]]$purity_threshold <- -.1; d },
    function(d) { d$stopping_rule <- "invented"; d })
  for (mutate in damage) {
    f <- beatrice_event_fixture(mutator = mutate)
    expect_null(.fb_joint_parse_output(f$path, f$ids, 3L, TRUE))
    unlink(f$path, recursive = TRUE)
  }
})

test_that("BEATRICE declared fork contract makes diagnostics mandatory", {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  expect_false(.fb_joint_requires_cs_diagnostics(path))
  dir.create(file.path(path, "scripts"))
  marker <- file.path(path, "scripts", "credible_set_contract.json")
  jsonlite::write_json(list(format = "beatrice-conditional-credible-set-diagnostics-v1",
    index_base = 0L, default_stopping_rule = "conditional_union"), marker, auto_unbox = TRUE)
  expect_true(.fb_joint_requires_cs_diagnostics(path))
  jsonlite::write_json(list(format = "beatrice-conditional-credible-set-diagnostics-v1",
    index_base = FALSE, default_stopping_rule = "conditional_union"), marker, auto_unbox = TRUE)
  expect_error(.fb_joint_requires_cs_diagnostics(path), "Invalid declared")
  jsonlite::write_json(list(format = "unknown"), marker, auto_unbox = TRUE)
  expect_error(.fb_joint_requires_cs_diagnostics(path), "Invalid declared")
})
