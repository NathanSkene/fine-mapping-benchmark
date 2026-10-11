#!/usr/bin/env Rscript
# Topological offline source installation on a scheduled compute node.
root <- Sys.getenv("FMB_CAMPAIGN_ROOT")
lib <- Sys.getenv("R_LIBS_USER")
stopifnot(nzchar(root), nzchar(lib))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
archives <- list.files(file.path(root, "downloads/R"), "\\.tar\\.gz$", full.names = TRUE)
staging <- file.path(root, "build/R-sources")
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
sources <- list()
for (archive in archives) {
  files <- utils::untar(archive, list = TRUE)
  desc <- files[grepl("^[^/]+/DESCRIPTION$", files)][1]
  utils::untar(archive, files = desc, exdir = staging)
  d <- read.dcf(file.path(staging, desc))
  sources[[d[1, "Package"]]] <- list(path = archive, description = d)
}
for (pkg in c("CARMA", "finimom", "susieR")) {
  path <- file.path(root, "tools", pkg)
  sources[[pkg]] <- list(path = path, description = read.dcf(file.path(path, "DESCRIPTION")))
}
dependencies <- function(d) {
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(d))
  setdiff(unique(unlist(lapply(fields, function(k)
    trimws(sub("\\s*\\(.*", "", strsplit(d[1, k], ",")[[1]]))))), "R")
}
remaining <- names(sources)
probe <- tempfile(fileext = ".R")
fresh_load <- function(pkg) {
  writeLines(sprintf("quit(status=if(requireNamespace(%s,quietly=TRUE)) 0L else 1L)",
                     dQuote(pkg, q = FALSE)), probe)
  system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(probe)),
          stdout = FALSE, stderr = FALSE) == 0L
}
while (length(remaining)) {
  ready <- remaining[vapply(remaining, function(pkg) {
    !length(intersect(dependencies(sources[[pkg]]$description), remaining))
  }, logical(1))]
  if (!length(ready)) stop("Unresolved dependency cycle: ", paste(remaining, collapse = ", "))
  for (pkg in ready) {
    item <- sources[[pkg]]
    installed <- installed.packages()
    # Git checkouts/compatibility patches can change without a package-version
    # bump. Always cleanly rebuild the three pinned source checkouts.
    if (pkg %in% c("CARMA", "finimom", "susieR") || !pkg %in% rownames(installed) ||
        installed[pkg, "Version"] != item$description[1, "Version"] || !fresh_load(pkg)) {
      install.packages(item$path, repos = NULL, type = "source", lib = lib,
                       INSTALL_opts = c("--preclean", "--clean"),
                       Ncpus = as.integer(Sys.getenv("PBS_NCPUS", "4")))
    }
    # A dependency may have been upgraded during this invocation. Probe in a
    # fresh process, avoiding an older namespace cached by the installer.
    if (!fresh_load(pkg)) stop("Installation failed: ", pkg)
  }
  remaining <- setdiff(remaining, ready)
}
write.csv(as.data.frame(installed.packages()[, c("Package", "Version", "LibPath")]),
          file.path(root, "R-versions.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(root, "R-session.txt"))
