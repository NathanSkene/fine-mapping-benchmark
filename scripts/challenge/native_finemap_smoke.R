#!/usr/bin/env Rscript
# Tiny external-binary test, not a performance benchmark or a null-locus test.
# FINEMAP's n_causal=1 analysis here conditions on a non-null locus.
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
repo_default <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "../.."))
args <- commandArgs(TRUE)
opts <- list(binary = "finemap", out = "results/native-finemap-smoke", repo = repo_default)
if (length(args) %% 2L) stop("Arguments must be --name value pairs")
if (length(args)) for (i in seq(1, length(args), 2L)) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opts)) stop("Unknown option ", key)
  opts[[key]] <- args[i + 1L]
}
source(file.path(opts$repo, "R/wrapper_finemap.R"))
binary <- normalizePath(opts$binary, mustWork = TRUE)
if (dir.exists(opts$out)) stop("Output already exists; choose a fresh --out")
dir.create(opts$out, recursive = TRUE)
out <- normalizePath(opts$out)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# A declared deterministic summary-statistic fixture: independent variants,
# n=2000, standardised effect estimates/se consistent with the supplied z.
# IDs deliberately disagree with lexical and posterior order.
p <- 12L
n <- 2000L
ids <- paste0("fixture_", c(9, 2, 12, 4, 1, 11, 3, 8, 5, 10, 6, 7))
z <- c(.1, -1, 6, .3, -.7, 1.2, -.2, .9, -.5, .6, 0, -.4)
LD <- diag(p)
se <- rep(1 / sqrt(n), p)
beta <- z * se
maf <- rep(.25, p)
fixture <- list(variant_ids = ids, z = z, LD = LD, n = n, beta = beta,
                se = se, maf = maf, n_causal = 1L, prior_std = .05,
                declaration = "Synthetic independent-LD non-null singleton smoke fixture")
saveRDS(fixture, file.path(out, "fixture.rds"), version = 2)

# Direct native output is an independent alignment reference. The wrapper's
# temporary files are not read, and the second wrapper run permutes every input.
prefix <- file.path(out, "native")
write.table(data.frame(rsid = ids, chromosome = 1L, position = seq_len(p),
                       allele1 = "A", allele2 = "T", maf = maf, beta = beta, se = se),
            paste0(prefix, ".z"), quote = FALSE, row.names = FALSE, sep = " ")
write.table(LD, paste0(prefix, ".ld"), quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = " ")
master <- file.path(out, "master")
writeLines(c("z;ld;snp;config;cred;log;n_samples",
             paste(paste0(prefix, c(".z", ".ld", ".snp", ".config", ".cred", ".log")),
                   collapse = ";") |> paste(n, sep = ";")), master)
stdout <- file.path(out, "native.stdout.log")
start <- proc.time()[["elapsed"]]
exit_status <- system2(binary, c("--sss", "--in-files", shQuote(master),
                                "--n-causal-snps", "1", "--n-iter", "2000",
                                "--prior-std", ".05", "--prob-cred-set", ".95", "--log"),
                       stdout = stdout, stderr = stdout)
native_runtime <- proc.time()[["elapsed"]] - start
stopifnot("Native FINEMAP returned nonzero status" = identical(exit_status, 0L))
native <- read.table(paste0(prefix, ".snp"), header = TRUE, stringsAsFactors = FALSE)
stopifnot(setequal(native$rsid, ids), !anyDuplicated(native$rsid))
reference <- native$prob[match(ids, native$rsid)]
stopifnot(length(reference) == p, all(is.finite(reference)),
          all(reference >= 0 & reference <= 1), abs(sum(reference) - 1) < 1e-4)

run <- function(order) run_finemap(z = z[order], LD = LD[order, order], n = n,
                                  beta = beta[order], se = se[order], maf = maf[order],
                                  variant_ids = ids[order], finemap_path = binary,
                                  n_causal = 1L, n_iter = 2000L, prior_std = .05)
a <- run(seq_len(p))
perm <- c(12, 3, 1, 7, 5, 2, 10, 6, 4, 11, 8, 9)
b <- run(perm)
for (result in list(a, b)) {
  stopifnot(is.null(result$error), identical(result$method, "finemap"),
            length(result$pip) == p, all(is.finite(result$pip)),
            all(result$pip >= 0 & result$pip <= 1),
            is.finite(result$runtime_seconds), result$runtime_seconds >= 0,
            length(result$credible_sets) > 0L,
            all(unlist(result$credible_sets) %in% seq_len(p)))
}
aligned_b <- b$pip[match(ids, ids[perm])]
max_native_difference <- max(abs(a$pip - reference))
max_permutation_difference <- max(abs(a$pip - aligned_b))
stopifnot("Wrapper PIPs differ from directly read native rsid alignment" = max_native_difference < 1e-6,
          "Wrapper PIPs change when inputs are permuted together" = max_permutation_difference < 1e-6,
          "Strongest independent z does not map to the correct variant" = which.max(a$pip) == which.max(abs(z)))
saveRDS(list(original = a, permuted = b), file.path(out, "wrapper-results.rds"), version = 2)
write.csv(data.frame(variant_id = ids, z = z, native_pip = reference,
                     wrapper_pip = a$pip, permuted_aligned_pip = aligned_b),
          file.path(out, "aligned-pips.csv"), row.names = FALSE)
checksum <- function(path) {
  value <- system2("sha256sum", shQuote(path), stdout = TRUE)
  if (!length(value)) stop("Could not checksum ", path)
  strsplit(value[[1]], "[[:space:]]+")[[1]][1]
}
provenance <- list(status = "PASS", binary = binary, binary_sha256 = checksum(binary),
                   binary_banner = readLines(stdout, warn = FALSE)[seq_len(min(15L, length(readLines(stdout))))],
                   wrapper_sha256 = checksum(file.path(opts$repo, "R/wrapper_finemap.R")),
                   fixture_sha256 = checksum(file.path(out, "fixture.rds")),
                   exit_status = exit_status, native_runtime_seconds = native_runtime,
                   wrapper_runtime_seconds = a$runtime_seconds,
                   permutation_runtime_seconds = b$runtime_seconds,
                   max_native_difference = max_native_difference,
                   max_permutation_difference = max_permutation_difference,
                   R = R.version.string, host = Sys.info()[["nodename"]],
                   PBS_JOBID = Sys.getenv("PBS_JOBID"), utc = format(Sys.time(), tz = "UTC"))
saveRDS(provenance, file.path(out, "provenance.rds"), version = 2)
dput(provenance, file.path(out, "provenance.txt"))
message("PASS: native FINEMAP, wrapper PIP alignment and input permutation; ", out)
