file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) .campaign_repo <- normalizePath(file.path(
  dirname(sub("^--file=", "", file_arg[1L])), "../..")) else .campaign_repo <- normalizePath(getwd())
if (nzchar(Sys.getenv("FMB_R_LIB"))) .libPaths(c(Sys.getenv("FMB_R_LIB"), .libPaths()))
if (dir.exists(file.path(.campaign_repo, ".local-lib")))
  .libPaths(c(file.path(.campaign_repo, ".local-lib"), .libPaths()))
for (.path in list.files(file.path(.campaign_repo, "R"), "\\.R$", full.names = TRUE)) source(.path)
.campaign_options <- function(defaults) {
  args <- commandArgs(TRUE)
  if (length(args) %% 2L) stop("Use --option value pairs")
  if (length(args)) for (i in seq(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i])
    if (!key %in% names(defaults)) stop("Unknown option --", key)
    defaults[[key]] <- args[i + 1L]
  }
  defaults
}
.campaign_main_source <- .campaign_repo

