#!/usr/bin/env Rscript
# Plotting/parser checks only; no inference or bootstrap.
source("paper/read_selected_csv.R")
checks <- 0L
check <- function(ok) { stopifnot(ok); checks <<- checks + 1L }
path <- tempfile(fileext = ".csv")
x <- data.frame(id = 1:12, value = c(NA, Inf, -Inf, seq_len(9)),
                text = c("comma, text", "quote\"text", "embedded\nnewline", rep("value", 9)))
write.csv(x, path, row.names = FALSE)
observed <- integer()
selected <- .paper_read_selected_csv(path, function(chunk) {
  observed <<- c(observed, nrow(chunk)); chunk[chunk$id %% 2L == 0L, , drop = FALSE]
}, chunk_size = 3L)
expected <- read.csv(path, stringsAsFactors = FALSE)
expected <- expected[expected$id %% 2L == 0L, ]; rownames(expected) <- NULL
check(identical(selected, expected))
check(length(observed) == 4L && max(observed) == 3L)
all_rows <- .paper_read_selected_csv(path, chunk_size = 2L)
check(identical(all_rows, read.csv(path, stringsAsFactors = FALSE)))
check(is.infinite(all_rows$value[2L]) && is.infinite(all_rows$value[3L]))
check(identical(all_rows$text[3L], "embedded\nnewline"))
check(!nrow(.paper_read_selected_csv(path, function(chunk) chunk[FALSE, ], chunk_size = 3L)))
check(inherits(try(.paper_read_selected_csv(path, max_rows = 5L), silent = TRUE), "try-error"))
check(inherits(try(.paper_read_selected_csv(path, function(chunk) list()), silent = TRUE), "try-error"))
empty <- tempfile(fileext = ".csv"); write.csv(x[FALSE, ], empty, row.names = FALSE)
check(identical(names(.paper_read_selected_csv(empty)), names(x)))
large <- tempfile(fileext = ".csv")
write.csv(data.frame(id = 1:150000, method = rep(c("abf", "susie"), 75000)), large, row.names = FALSE)
max_chunk <- 0L; calls <- 0L
tail_rows <- .paper_read_selected_csv(large, function(chunk) {
  calls <<- calls + 1L; max_chunk <<- max(max_chunk, nrow(chunk))
  chunk[chunk$id > 149990L, , drop = FALSE]
}, chunk_size = 10000L)
check(nrow(tail_rows) == 10L && identical(tail_rows$id, 149991:150000))
check(calls == 15L && max_chunk == 10000L)
unlink(c(path, empty, large))
cat("Passed", checks, "paper CSV reader checks; no scientific outcomes recomputed.\n")
