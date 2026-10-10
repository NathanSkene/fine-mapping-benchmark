# Cache identity is a digest of the complete ordered observed input. Normalize
# storage types and discard incidental dimnames; explicit variant/annotation
# labels and matrix dimensions remain part of the scientific input identity.
.fm_pooled_input_key <- function(z, LD = NULL, annotations = NULL, n = NULL,
                                  variant_ids = NULL) {
  matrix_input <- function(x, column_labels = FALSE) {
    if (is.null(x)) return(NULL)
    list(dim = as.integer(dim(x)), values = as.numeric(x),
         column_names = if (!column_labels || is.null(colnames(x))) NULL else as.character(colnames(x)))
  }
  input <- list(schema = 1L, z = as.numeric(z), LD = matrix_input(LD),
                annotations = matrix_input(annotations, column_labels = TRUE),
                n = if (is.null(n)) NULL else as.numeric(n),
                variant_ids = if (is.null(variant_ids)) NULL else as.character(variant_ids))
  path <- tempfile("pooled_input_")
  on.exit(unlink(path), add = TRUE)
  writeBin(serialize(input, NULL, version = 2), path)
  .challenge_sha256(path)
}
