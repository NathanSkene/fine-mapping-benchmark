file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.native_repo <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1L])), "../../.."))
if (nzchar(Sys.getenv("FMB_R_LIB"))) .libPaths(c(Sys.getenv("FMB_R_LIB"), .libPaths()))
if (dir.exists(file.path(.native_repo, ".local-lib"))) .libPaths(c(file.path(.native_repo, ".local-lib"), .libPaths()))
for (.path in list.files(file.path(.native_repo, "R"), "\\.R$", full.names = TRUE)) source(.path)
source(file.path(.native_repo, "scripts/native/campaign/helpers.R"))
.native_options <- function(defaults) {
  args <- commandArgs(TRUE)
  if (length(args) %% 2L) stop("Use --option value pairs")
  if (length(args)) for (i in seq(1L, length(args), 2L)) {
    key <- sub("^--", "", args[i]); if (!key %in% names(defaults)) stop("Unknown --", key)
    defaults[[key]] <- args[i + 1L]
  }
  defaults
}
