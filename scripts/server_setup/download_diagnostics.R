#!/usr/bin/env Rscript
# Network-only preparation for diagnostics added after the main cache freeze.
root <- Sys.getenv("FMB_CAMPAIGN_ROOT")
stopifnot(nzchar(root))
cache <- file.path(root, "downloads/R")
repos <- "https://cloud.r-project.org"
ap <- available.packages(repos = repos)
packages <- unique(c("posterior", unlist(tools::package_dependencies(
  "posterior", ap, which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))
packages <- setdiff(packages, rownames(installed.packages(priority = c("base", "recommended"))))
files <- file.path(cache, paste0(packages, "_", ap[packages, "Version"], ".tar.gz"))
needed <- packages[!file.exists(files)]
if (length(needed)) download.packages(needed, destdir = cache, available = ap,
                                     repos = repos, type = "source")
write.csv(data.frame(package = packages, version = ap[packages, "Version"]),
  file.path(cache, "diagnostics-versions.csv"), row.names = FALSE)
