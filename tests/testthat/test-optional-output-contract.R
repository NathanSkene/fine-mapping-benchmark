test_that("failed annotated joint fits cannot become locus-only estimates", {
  geno <- list(LD = diag(3), variant_ids = c("a", "b", "c"), n = 100,
               annotations_matrix = matrix(1, 3, 1))
  pheno <- list(z = c(4, 0, -1))
  x <- .fb_joint_region(geno, pheno, "fb_pooled", .fb_joint_cache = list())
  expect_match(x$error, "no fallback")
  expect_true(all(is.na(x$pip)))
  expect_false(x$additional$joint_fallback)
  expect_error(.fb_joint_scenario_setup(list(geno), list(pheno), list(), "linear"),
               "needs Python")
  no_annotation <- geno
  no_annotation$annotations_matrix <- NULL
  expect_error(.fb_joint_scenario_setup(list(geno, no_annotation), list(pheno, pheno), list(), "linear"),
               "every region")
  setup <- .fb_joint_scenario_setup(list(no_annotation), list(pheno), list(), "linear")
  expect_false(setup$.fb_joint_annotated)
})

test_that("joint output parser rejects ID corruption and invalid probabilities", {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  file <- file.path(path, "pip.csv")
  write.csv(data.frame(variant_names = c("b", "a"), pip = c(.2, .8)), file, row.names = FALSE)
  expect_equal(.fb_joint_parse_output(path, c("a", "b"), 2)$pip, c(.8, .2))
  write.csv(data.frame(variant_names = c("a", "a"), pip = c(.2, .8)), file, row.names = FALSE)
  expect_null(.fb_joint_parse_output(path, c("a", "b"), 2))
  write.csv(data.frame(variant_names = c("a", "b"), pip = c(1.2, -.1)), file, row.names = FALSE)
  expect_null(.fb_joint_parse_output(path, c("a", "b"), 2))
  write.csv(data.frame(variant_names = "a", pip = .8), file, row.names = FALSE)
  expect_null(.fb_joint_parse_output(path, c("a", "b"), 2))
})

fminf_contract_fixture <- function(table, status = NULL) {
  path <- tempfile(); dir.create(path)
  writeLines("# fixture", file.path(path, "run_fine_mapping.py"))
  fn <- run_finemap_inf
  e <- new.env(parent = environment(fn))
  e$system2 <- function(command, args, ...) {
    prefix <- args[match("--output-prefix", args) + 1L]
    write.table(table, paste0(prefix, ".finemapinf.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    x <- "fixture output"
    if (!is.null(status)) attr(x, "status") <- status
    x
  }
  environment(fn) <- e
  x <- fn(z = c(3, 1), LD = diag(2), n = 100, variant_ids = c("a", "b"), finemap_inf_dir = path)
  unlink(path, recursive = TRUE)
  x
}

test_that("FINEMAP-inf retains native events and keeps PIP mass descriptive", {
  x <- fminf_contract_fixture(data.frame(SNP = c("b", "a"), prob = c(.3, .8)))
  expect_equal(x$pip, c(.8, .3))
  expect_length(x$credible_sets, 0)
  expect_equal(x$additional$pip_mass_set, c(1L, 2L))
  expect_equal(x$additional$credible_set_semantics, "not_provided")
  y <- fminf_contract_fixture(data.frame(SNP = c("b", "a"), prob = c(.3, .8), cs = c(-1, 1)))
  expect_equal(y$credible_sets, list(1L))
  expect_null(y$additional$pip_mass_set)
  expect_equal(y$additional$credible_set_semantics, "native_component_sets")
})

test_that("FINEMAP-inf rejects failed executions and damaged complete outputs", {
  expect_match(fminf_contract_fixture(data.frame(prob = c(.8, .3)))$error, "SNP IDs are missing")
  expect_match(fminf_contract_fixture(data.frame(SNP = c("a", "a"), prob = c(.8, .3)))$error,
               "SNP IDs")
  expect_match(fminf_contract_fixture(data.frame(SNP = c("a", "b"), prob = c(1.2, .3)))$error,
               "finite probability")
  expect_match(fminf_contract_fixture(data.frame(SNP = c("a", "b"), prob = c(.8, .3)), 1L)$error,
               "exited unsuccessfully")
})

beatrice_contract_fixture <- function(method, table = NULL, status = NULL,
                                     set = "1 0", conditional = ".2 .8",
                                     feature = NULL, ids = c("a b", "a_b")) {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  writeLines("# fixture", file.path(path, "beatrice_annot.py"))
  fn <- get(paste0("run_", method))
  e <- new.env(parent = environment(fn))
  captured <- NULL
  e$system2 <- function(command, args, ...) {
    input <- read.table(args[match("--z", args) + 1L], stringsAsFactors = FALSE)
    captured <<- list(z = input)
    if ("--annot" %in% args)
      captured$annotations <<- read.table(args[match("--annot", args) + 1L], stringsAsFactors = FALSE)
    out <- args[match("--target", args) + 1L]
    native_table <- table
    if (is.null(native_table)) native_table <- data.frame(
      variant_names = rev(input[[1L]]), variant_index = c(1L, 0L), pip = c(.2, .8))
    write.csv(native_table, file.path(out, "pip.csv"), row.names = FALSE)
    if (!is.null(set)) writeLines(set, file.path(out, "credible_set.txt"))
    if (!is.null(conditional)) writeLines(conditional,
      file.path(out, "conditional_credible_variants_probability.txt"))
    if (!is.null(feature)) write.csv(feature, file.path(out, "feature_importance.csv"), row.names = FALSE)
    result <- "fixture output"
    if (!is.null(status)) attr(result, "status") <- status
    result
  }
  environment(fn) <- e
  args <- list(z = c(3, 1), LD = diag(2), n = 100, variant_ids = ids,
               beatrice_dir = path, max_iter = 500)
  if (method == "functional_beatrice") args$annotations <- matrix(c(1, 0), 2, 1)
  result <- do.call(fn, args)
  result$fixture_input <- captured
  result
}

test_that("BEATRICE wrappers preserve IDs and align native probabilities and set members", {
  for (method in c("beatrice", "functional_beatrice")) {
    fit <- beatrice_contract_fixture(method)
    expect_null(fit$error)
    expect_equal(fit$pip, c(.8, .2))
    expect_identical(fit$credible_sets, list(c(1L, 2L)))
    expect_equal(fit$additional$cs_pip, list(c(.8, .2)))
    expect_identical(fit$additional$variant_ids, c("a b", "a_b"))
    expect_identical(fit$additional$variant_alignment$variant_id, c("a b", "a_b"))
    expect_identical(fit$fixture_input$z[[1L]], c("variant_1", "variant_2"))
    expect_equal(fit$fixture_input$z[[2L]], c(3, 1))
    expect_equal(ncol(fit$fixture_input$z), 2)
    if (method == "functional_beatrice") {
      expect_identical(fit$fixture_input$annotations[[1L]], c("variant_1", "variant_2"))
      expect_equal(fit$fixture_input$annotations[[2L]], c(1, 0))
    }
    no_sets <- beatrice_contract_fixture(method, set = NULL, conditional = NULL, status = 0L)
    expect_null(no_sets$error)
    expect_length(no_sets$credible_sets, 0)
    expect_null(no_sets$additional$cs_pip)
    expect_match(beatrice_contract_fixture(method, status = 1L)$error, "exited unsuccessfully")
  }
})

test_that("BEATRICE wrappers reject missing, duplicate, unknown and partial IDs and invalid PIPs", {
  valid <- data.frame(variant_names = c("variant_1", "variant_2"), pip = c(.8, .2))
  damaged <- list(valid["pip"], valid[1L, ], rbind(valid, valid[1L, ]))
  for (ids in list(c("variant_1", "variant_1"), c("variant_1", "missing"))) {
    x <- valid; x$variant_names <- ids; damaged[[length(damaged) + 1L]] <- x
  }
  for (pip in c(-.1, 1.1, NA_real_, NaN, Inf)) {
    x <- valid; x$pip[1L] <- pip; damaged[[length(damaged) + 1L]] <- x
  }
  x <- valid; x$variant_index <- c(1L, 0L); damaged[[length(damaged) + 1L]] <- x
  for (method in c("beatrice", "functional_beatrice")) {
    for (table in damaged) {
      fit <- beatrice_contract_fixture(method, table = table)
      expect_match(fit$error, "invalid variant IDs, finite probabilities, or credible-set evidence")
      expect_true(all(is.na(fit$pip)))
      expect_length(fit$credible_sets, 0L)
    }
    for (ids in list(c("a", "a"), c("a", NA_character_), "a", c("a", " ")))
      expect_error(beatrice_contract_fixture(method, ids = ids), "unique nonempty")
  }
})

test_that("BEATRICE wrappers reject damaged credible sets and conditional probabilities", {
  for (method in c("beatrice", "functional_beatrice")) {
    for (set in c("0 broken", "0 0", "-1", "2", "0.5", "99999999999"))
      expect_false(is.null(beatrice_contract_fixture(method, set = set)$error))
    for (conditional in c(".8", ".8 NaN", ".8 1.2", ".8 -0.1"))
      expect_false(is.null(beatrice_contract_fixture(method, conditional = conditional)$error))
    expect_false(is.null(beatrice_contract_fixture(method, set = NULL)$error))
  }
  valid <- data.frame(annotation_index = 0, importance = .3)
  expect_null(beatrice_contract_fixture("functional_beatrice", feature = valid)$error)
  for (feature in list(valid["importance"], data.frame(annotation_index = 1, importance = .3),
                       data.frame(annotation_index = 0, importance = NaN),
                       data.frame(annotation_index = 0, importance = -.1)))
    expect_match(beatrice_contract_fixture("functional_beatrice", feature = feature)$error,
                 "feature importance")
})

test_that("BEATRICE setup and runner select the same candidate and joint serialization is collision free", {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  writeLines("# fixture", file.path(path, "beatrice_annot.py"))
  fn <- setup_beatrice
  e <- new.env(parent = environment(fn)); e$system2 <- function(...) 0L
  environment(fn) <- e
  python <- file.path(path, "python")
  writeLines("# fixture", python)
  expect_identical(suppressMessages(fn(path, python))$beatrice_script, .beatrice_script_path(path))
  writeLines("# fixture", file.path(path, "beatrice.py"))
  expect_identical(suppressMessages(fn(path, python))$beatrice_script, file.path(path, "beatrice_annot.py"))
  unlink(file.path(path, "beatrice_annot.py"))
  expect_identical(suppressMessages(fn(path, python))$beatrice_script, file.path(path, "beatrice.py"))
  ids <- c("a\tb", "a_b")
  written <- .fb_joint_write_region(path, 1, c(3, 1), diag(2), matrix(c(1, 0), 2, 1), ids)
  expect_identical(read.table(written$z)[[1L]], c("variant_1", "variant_2"))
  expect_identical(read.table(written$annot)[[1L]], c("variant_1", "variant_2"))
  expect_identical(written$original_ids, ids)
  expect_identical(written$variant_alignment$variant_id, ids)
})

test_that("Python resolution preserves distinct symlinked virtual environments", {
  skip_if(.Platform$OS.type == "windows", "POSIX symlinked virtual environment regression")
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "Python3 not available for virtual environment regression")
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  prefixes <- character()
  for (name in c("first", "second")) {
    venv <- file.path(path, name)
    status <- system2(python, c("-m", "venv", "--without-pip", "--symlinks", shQuote(venv)))
    expect_equal(status, 0L)
    exe <- file.path(venv, "bin", "python")
    resolved <- .fm_python_path(exe)
    expect_identical(resolved, file.path(normalizePath(dirname(exe)), "python"))
    expect_false(identical(resolved, normalizePath(exe)))
    prefix <- system2(resolved, c("-c", shQuote("import sys; print(sys.prefix)")), stdout = TRUE)
    expect_identical(prefix, normalizePath(venv))
    prefixes <- c(prefixes, prefix)
  }
  expect_false(identical(prefixes[[1L]], prefixes[[2L]]))
})

test_that("PAINTOR output cannot silently reorder, clip or substitute fits", {
  valid <- data.frame(rsid = c("variant_2", "variant_1"), Posterior_Prob = c(.3, .8))
  expect_equal(.paintor_parse_probabilities(valid, c("variant_1", "variant_2"))$pip, c(.8, .3))
  damaged <- valid; damaged$rsid <- rep("variant_1", 2)
  expect_null(.paintor_parse_probabilities(damaged, c("variant_1", "variant_2")))
  damaged <- valid; damaged$Posterior_Prob[1] <- 1.2
  expect_null(.paintor_parse_probabilities(damaged, c("variant_1", "variant_2")))
  expect_null(.paintor_parse_probabilities(valid[1, ], c("variant_1", "variant_2")))
  fit <- run_paintor_region(list(LD = diag(2)), list(z = c(2, 1)), .paintor_cache = list())
  expect_match(fit$error, "no single-locus fallback")
  expect_true(all(is.na(fit$pip)))
})

external_contract_fixture <- function(method, table, status = NULL, set = character()) {
  path <- tempfile(); dir.create(path)
  on.exit(unlink(path, recursive = TRUE))
  writeLines("# fixture", file.path(path, "sparsepro_zld.py"))
  fn <- get(paste0("run_", method))
  e <- new.env(parent = environment(fn))
  e$system2 <- function(command, args, ...) {
    if (method == "caviar") {
      prefix <- args[match("-o", args) + 1L]
      write.table(table, paste0(prefix, "_post"), row.names = FALSE, quote = FALSE)
      writeLines(set, paste0(prefix, "_set"))
    } else {
      out <- args[match("--save", args) + 1L]
      write.table(table, file.path(out, "region.z.pip"), row.names = FALSE,
        col.names = FALSE, quote = FALSE)
    }
    x <- "fixture output"; if (!is.null(status)) attr(x, "status") <- status
    x
  }
  environment(fn) <- e
  args <- list(z = c(3, 1), LD = diag(2), variant_ids = c("a", "b"))
  if (method == "sparsepro") args <- c(args, list(n = 100, sparsepro_dir = path))
  do.call(fn, args)
}

test_that("CAVIAR and SparsePro require valid complete native outputs", {
  cv <- data.frame(SNP = c("b", "a"), Causal_Post._Prob = c(.3, .8))
  sp <- data.frame(SNP = c("b", "a"), z = c(1, 3), pip = c(.3, .8))
  expect_equal(external_contract_fixture("caviar", cv, set = "a")$pip, c(.8, .3))
  expect_equal(external_contract_fixture("sparsepro", sp)$pip, c(.8, .3))
  for (method in c("caviar", "sparsepro")) {
    table <- if (method == "caviar") cv else sp
    damaged <- table; damaged[[1]] <- c("a", "a")
    expect_match(external_contract_fixture(method, damaged)$error, "IDs")
    damaged <- table; damaged[[ncol(damaged)]][1] <- 1.2
    expect_match(external_contract_fixture(method, damaged)$error, "finite probability")
    expect_false(is.null(external_contract_fixture(method, table, status = 1L)$error))
  }
  expect_match(external_contract_fixture("caviar", cv, set = "wrong")$error, "invalid IDs")
})
