# Official SBayesRC comparator

`run_sbayesrc_native()` invokes official GCTB once for the complete externally
prepared panel. It provides an optional upstream comparator alongside the
legacy `sbayesrc` region surrogate. Joint input includes every supplied LD
block and every supplied SNP; `variant_ids` selects imported outputs only.
The manifest declares `upstream_mode = "sbayes_RC"`: this is standard official
SBayesRC inference. GCTB's newer `--gwfm RC` workflow has extra preparation,
convergence and credible-set machinery; this adapter does not reproduce that
workflow. The [official GCTB tutorials](https://gctbhub.cloud.edu.au/software/gctb/)
describe the modes separately.
The function is deliberately separate from `run_methods()`, whose adapters
dispatch independent regions. A local region simulation is not automatically
converted into a genome-wide SBayesRC experiment.

The command and output schema were checked against [official GCTB source
2.5.5](https://github.com/jianzeng/GCTB/tree/cc7fa7d765c83a89c6375946cf77fe50ba1a317e).
The verified options are `--ldm-eigen`, `--gwas-summary`, `--sbayes RC`,
`--annot`, `--chain-length`, `--burn-in`, `--seed`, `--thread`, and `--out`.
GCTB's no-argument banner supplies the version; the adapter does not assume
unsupported `--help` or `--version` options. Installation and scientific
context are provided by the [official SBayesRC
README](https://github.com/zhilizheng/SBayesRC/blob/2881f9a6f57374b80b4d19a09dbf9387939a2f46/README.md)
and [SBayesRC paper](https://www.nature.com/articles/s41588-024-01704-y).

## Prepared input contract

Prepare and freeze these files using upstream tools before running the adapter.
They must describe the same SNP set and effect alleles, with documented genome
build, ancestry, imputation, allele harmonization and selection decisions.

| Input | Required content |
|---|---|
| GWAS file | Exactly `SNP A1 A2 freq b se p N`, in that order. `b` and `se` are original per-allele estimates, and `freq` is the A1 allele frequency. |
| `LDdir/snp.info` | The official ordered columns `Chrom ID Index GenPos PhysPos A1 A2 A1Freq N Block`. Index runs from zero in SNP row order. |
| `LDdir/ldm.info` | The official GCTB ordered columns `Block Chrom StartSnpIdx StartSnpID EndSnpIdx EndSnpID NumSnps`. Boundaries and SNP memberships must agree with `snp.info`. |
| `LDdir/block<BLOCK>.eigen.bin` | Every eigen-LD block listed in the metadata, in the upstream binary format. |
| Annotation file | `SNP Intercept` followed by finite numeric annotation categories; the intercept is one for every SNP. |

The adapter requires the exact same full SNP set in GWAS, LD and annotations.
Missing variants must be handled explicitly during external preprocessing;
the adapter does not silently insert zero annotations, impute GWAS effects,
or reduce the modeled panel to the evaluation SNPs. It checks both allele
directions and imports effect signs relative to the prepared GWAS A1. Strand
complements and ambiguous allele decisions belong in external harmonization.
GCTB reads LD metadata positionally, so its modern column order is enforced;
legacy metadata used by some R package resources must be converted first.

`scale` is mandatory. `genome_wide` requires all 22 autosomes in the panel;
`chromosome` requires exactly one chromosome; `subset` makes no full-genome
claim. Presence of all autosomes does not establish high-density coverage.
The manifest records the actual SNP count, chromosomes and LD block count.
Report those counts and the panel selection rule alongside the scope label.

## Run and review

```r
provenance <- list(
  gwas_source = "Frozen discovery summary release, accession and date",
  ld_source = "Upstream eigen-LD panel, release and sample ancestry",
  annotation_source = "Frozen annotation release and selection",
  genome_build = "GRCh37",
  ancestry = "Describe discovery and reference ancestries separately",
  panel_id = "Unique frozen panel identifier",
  preprocessing = "Describe QC, allele harmonization, imputation and exclusions"
)

# Validate and inspect the single joint command without invoking GCTB.
plan <- run_sbayesrc_native(
  "prepared.ma", "prepared_ld", "prepared_annotations.txt", "new-fit",
  scale = "genome_wide", provenance = provenance, dry_run = TRUE
)

# Supply the official executable and use a fresh output directory for each fit.
fit <- run_sbayesrc_native(
  "prepared.ma", "prepared_ld", "prepared_annotations.txt", "new-fit",
  scale = "genome_wide", provenance = provenance,
  gctb_path = "/path/to/official/gctb", expected_version = "2.5.5",
  variant_ids = common_evaluation_ids,
  n_iter = 3000, burn_in = 1000, seed = 22, threads = 1
)
```

`pip` and `additional$posterior_mean_beta` are in the requested stable SNP ID
order. The entire native output is checked before any subset is selected:
duplicate/missing/extra IDs, incompatible alleles, nonfinite effects, negative
posterior SDs, and PIPs outside `[0,1]` reject the result. The official `SE`
field is imported as posterior effect SD. No credible sets are manufactured
from summed marginal PIPs. Inspect `native_process.log`, the original
`.snpRes`, and any `.parRes` file retained by upstream inference.

The output directory must be new or empty. Input GWAS, annotations, LD metadata
and all eigen binaries receive full-file MD5 fingerprints. The executable's
banner/version and checksum are recorded, and changed inputs or executable
bytes during a fit cause rejection. These checks support operational freezing
and auditability; MD5 is not an adversarial authenticity mechanism. Obtain the
binary from the official source and retain its download/source provenance.

Validation reads annotation and metadata tables into R memory; plan RAM for the
full panel. `data.table` provides optional efficient parsing and is required
for individual text files above 128 MiB. Hashing reads every LD block even in
a dry run, which has a real I/O cost at genome scale. A dry run checks schemas,
sets, allele compatibility and block presence/boundaries; upstream GCTB still
validates the binary eigendecompositions when it executes.

## Scientific limits and validation evidence

An official executable is a meaningful comparator only when its inputs and
resources match the stated problem. Include complete discovery/reference
ancestry, actual fitted SNP count, MCMC budget, preprocessing exclusions and
seed in the comparison record. A model jointly fitted to a selected GWAS locus
does not establish genome-wide annotation learning. Evaluation on common SNPs
also needs an explicit missing-causal-variant estimand.

Neither 3000 iterations nor successful exit establishes MCMC convergence.
The adapter records `convergence_claim = "not_established"`; diagnostics,
independent chains, effective sample sizes and posterior stability must be
assessed for scientific claims. Upstream GCTB can set effects/PIPs to zero for
SNPs it flags as unconverged; inspect its logs rather than interpreting zeros
as a convergence certificate. The adapter does not infer calibration or method
accuracy from successful execution.

The unit tests exercise schema, scope, exact SNP-set alignment, effect sign
orientation, full-panel result validation, no output reuse, and frozen-file
fingerprints. Fixture `.snpRes` values are hand-written import contracts; their
placeholder eigen files are never submitted as inference inputs.

`scripts/native/sbayesrc_smoke.R` additionally supplies a bounded real-executable
smoke: 256 synthetic SNPs on one chromosome in four identity-LD blocks, with
original per-allele Gaussian summary estimates and source-verified native
binary eigendecompositions. It fits 300 iterations with 150 burn-in and records
the native version, executable/input fingerprints, log and imported output.
Use `scripts/native/imperial_native_job.sh` for the isolated Imperial PBS run.
This smoke tests the upstream executable and adapter, and makes no human
ancestry, individual-genotype, MCMC-convergence or genome-wide scale claim.

As checked on 10 October 2026, the smallest public complete HapMap3 reference
linked by upstream was approximately 3.14 GB compressed. A full released
genome-wide panel is not bundled here, and execution of that resource remains
a separate validation requirement.
