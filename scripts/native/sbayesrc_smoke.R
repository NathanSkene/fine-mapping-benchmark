#!/usr/bin/env Rscript
# Real upstream executable smoke, chromosome-scale synthetic summary inputs.
# This is not a calibration experiment or a genome-wide scale demonstration.
# The eigen binary layout is verified against GCTB eigen.cpp at commit
# cc7fa7d765c83a89c6375946cf77fe50ba1a317e, readEigenMatrixBinaryFileAndMakeWandQ.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript sbayesrc_smoke.R ADAPTER_FILE GCTB_BINARY NEW_WORK_DIR")
source(normalizePath(args[1], mustWork = TRUE))
binary <- normalizePath(args[2], mustWork = TRUE)
work <- args[3]
if (dir.exists(work)) stop("Use a fresh smoke work directory")
dir.create(work, recursive = TRUE)
work <- normalizePath(work)
ld <- file.path(work, "ld")
dir.create(ld)
set.seed(4101)
p <- 256L
block_size <- 64L
n <- 10000L
ids <- paste0("synthetic", seq_len(p))
freq <- rep(.3, p)
category <- rep(c(0, 1), p / 2)
truth <- numeric(p)
causal <- sample(which(category == 1), 16L)
truth[causal] <- sample(c(-1, 1), 16L, replace = TRUE) * sqrt(.25 / (16 * 2 * .3 * .7))
se <- sqrt(1 / (n * 2 * freq * (1 - freq)))
effect <- truth + rnorm(p, sd = se)
gwas <- data.frame(SNP = ids, A1 = "A", A2 = "G", freq = freq,
                    b = effect, se = se, p = 2 * pnorm(-abs(effect / se)), N = n)
panel <- data.frame(Chrom = 1L, ID = ids, Index = seq_len(p) - 1L,
                     GenPos = 0, PhysPos = seq_len(p) * 1000L,
                     A1 = "A", A2 = "G", A1Freq = freq, N = n,
                     Block = rep(1:4, each = block_size))
blocks <- data.frame(Block = 1:4, Chrom = 1L,
                      StartSnpIdx = seq(0L, p - block_size, by = block_size),
                      StartSnpID = ids[seq(1L, p, by = block_size)],
                      EndSnpIdx = seq(block_size - 1L, p - 1L, by = block_size),
                      EndSnpID = ids[seq(block_size, p, by = block_size)], NumSnps = block_size)
annotation <- data.frame(SNP = ids, Intercept = 1L, synthetic_category = category)
writer <- function(x, path) write.table(x, path, quote = FALSE, row.names = FALSE)
writer(gwas, file.path(work, "gwas.ma"))
writer(annotation, file.path(work, "annotations.txt"))
writer(panel, file.path(ld, "snp.info"))
writer(blocks, file.path(ld, "ldm.info"))
# Binary header: int32 m, int32 k, float32 sum(lambda), float32 cutoff;
# then float32 eigenvalues and column-major m x k float32 eigenvectors.
# lambda=1 and U=I reconstruct the specified identity correlation exactly.
for (b in blocks$Block) {
  con <- file(file.path(ld, paste0("block", b, ".eigen.bin")), "wb")
  writeBin(c(block_size, block_size), con, size = 4L, endian = "little")
  writeBin(c(as.double(block_size), 1), con, size = 4L, endian = "little")
  writeBin(rep(1, block_size), con, size = 4L, endian = "little")
  writeBin(as.numeric(diag(block_size)), con, size = 4L, endian = "little")
  close(con)
}
saveRDS(list(beta_per_allele = truth, causal_indices = causal,
              model = "independent Gaussian summary estimates; identity population LD",
              covariance = "diag(1 / (N * 2*f*(1-f)))", phenotype_variance = 1,
              target_genetic_variance = .25, seed = 4101L), file.path(work, "synthetic_truth.rds"))
provenance <- list(
  gwas_source = "Synthetic independent Gaussian summary estimates, seed4101; no individual genotype claims",
  ld_source = "Exact identity correlation, source-verified GCTB binary eigendecomposition; four blocks",
  annotation_source = "Synthetic alternating binary category; sixteen planted effects restricted to category1",
  genome_build = "synthetic coordinates, no human build",
  ancestry = "synthetic, no human ancestry",
  panel_id = "synthetic-chromosome1-256-snp-smoke-v1",
  preprocessing = "Complete common SNP set, exact A/G effect alleles, original per-allele effects, no filtering"
)
fit <- run_sbayesrc_native(file.path(work, "gwas.ma"), ld, file.path(work, "annotations.txt"),
                           file.path(work, "native-fit"), scale = "chromosome",
                           provenance = provenance, gctb_path = binary,
                           n_iter = 300L, burn_in = 150L, seed = 22L, threads = 1L)
stopifnot(length(fit$pip) == p, all(is.finite(fit$pip)),
          all(fit$pip >= 0 & fit$pip <= 1),
          fit$additional$native_executed,
          identical(fit$additional$variant_ids, ids),
          fit$additional$manifest$modeled_variant_count == p,
          identical(fit$additional$manifest$scale, "chromosome"))
saveRDS(fit, file.path(work, "validated_fit.rds"))
summary <- data.frame(implementation = fit$additional$implementation,
                       version = fit$additional$manifest$actual_version,
                       scope = "chromosome", variants = p, blocks = 4L,
                       n_iter = 300L, burn_in = 150L,
                       runtime_seconds = fit$runtime_seconds,
                       pip_min = min(fit$pip), pip_max = max(fit$pip),
                       sum_pip = sum(fit$pip), convergence_established = FALSE)
write.csv(summary, file.path(work, "smoke_summary.csv"), row.names = FALSE)
print(summary)
message("Real official GCTB adapter smoke passed; MCMC convergence and genome-wide scale remain untested.")
