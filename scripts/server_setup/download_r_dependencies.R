#!/usr/bin/env Rscript
# Downloads only: compilation belongs on a scheduled compute node.
root <- Sys.getenv("FMB_CAMPAIGN_ROOT")
stopifnot(nzchar(root))
cache <- file.path(root, "downloads/R")
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
repos <- "https://cloud.r-project.org"
ap <- available.packages(repos = repos)
roots <- c("susieR", "nnet", "reticulate", "data.table", "yaml", "testthat",
           "ggplot2", "gridExtra", "Rcpp", "RcppArmadillo", "RcppGSL", "roptim",
           "glmnet", "mvtnorm", "stringr", "readr", "knitr", "rmarkdown", "remotes", "posterior")
local_sources <- c(file.path(root, "tools/CARMA/DESCRIPTION"),
                   file.path(root, "tools/finimom/DESCRIPTION"),
                   file.path(root, "tools/susieR/DESCRIPTION"))
for (path in local_sources) {
  d <- read.dcf(path)
  for (field in intersect(c("Depends", "Imports", "LinkingTo"), colnames(d))) {
    pkgs <- trimws(sub("\\s*\\(.*", "", strsplit(d[1, field], ",")[[1]]))
    roots <- c(roots, setdiff(pkgs, "R"))
  }
}
roots <- unique(roots)
deps <- unique(c(roots, unlist(tools::package_dependencies(roots, ap, which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))
base <- rownames(installed.packages(priority = c("base", "recommended")))
deps <- setdiff(deps, base)
missing <- setdiff(deps, rownames(ap))
if (length(missing)) stop("CRAN cannot resolve: ", paste(missing, collapse = ", "))
dl <- download.packages(deps, destdir = cache, available = ap, repos = repos, type = "source")
write.csv(data.frame(package = deps, version = ap[deps, "Version"]), file.path(cache, "cran-versions.csv"), row.names = FALSE)
for (item in list(c("hapsim", "0.31"), c("sim1000G", "1.40"))) {
  url <- sprintf("https://cran.r-project.org/src/contrib/Archive/%s/%s_%s.tar.gz", item[1], item[1], item[2])
  download.file(url, file.path(cache, basename(url)), mode = "wb")
}
writeLines(capture.output(sessionInfo()), file.path(cache, "download-session.txt"))
