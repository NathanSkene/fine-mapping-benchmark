# Retain only declared plotting rows; read.csv supplies the CSV quoting/parser
# contract, including embedded newlines. No scientific summaries are recomputed.
.paper_read_selected_csv <- function(path, select = identity,
                                     chunk_size = 10000L, max_rows = 100000L) {
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L ||
      length(max_rows) != 1L || is.na(max_rows) || max_rows < 1L)
    stop("Invalid figure CSV memory bound")
  connection <- file(path, "rt")
  on.exit(close(connection), add = TRUE)
  # scan's nmax=0 means unlimited, so nrows=0 must never be used for a header.
  chunk <- utils::read.csv(connection, colClasses = "character", nrows = chunk_size,
                          stringsAsFactors = FALSE, check.names = FALSE)
  header <- chunk[FALSE, , drop = FALSE]
  pieces <- list(); retained <- 0L; template <- NULL
  repeat {
    if (!nrow(chunk)) break
    # Character parsing preserves CSV values until this bounded chunk is typed,
    # using the same conversion as read.csv. Undefined and infinite scores stay.
    chunk[] <- lapply(chunk, utils::type.convert, as.is = TRUE)
    selected <- select(chunk)
    if (!is.data.frame(selected)) stop("Figure CSV selector must return a data frame")
    if (is.null(template)) template <- selected[FALSE, , drop = FALSE]
    retained <- retained + nrow(selected)
    if (retained > max_rows)
      stop("Figure CSV selection exceeds ", format(max_rows, big.mark = ","),
           " rows; select fewer settings/methods")
    if (nrow(selected)) pieces[[length(pieces) + 1L]] <- selected
    if (nrow(chunk) < chunk_size) break
    chunk <- utils::read.csv(connection, header = FALSE, col.names = names(header),
      colClasses = "character", nrows = chunk_size, stringsAsFactors = FALSE,
      check.names = FALSE)
  }
  if (!length(pieces)) return(if (is.null(template)) select(header) else template)
  result <- do.call(rbind, pieces)
  rownames(result) <- NULL
  result
}
