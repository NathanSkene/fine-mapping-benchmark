# Audit of the original result claims

Audit date: 10 October 2026. The comparison baseline is upstream commit
`f053fdf318a3abe2c44eeb0ac3f3c99a03979b81`, not a moving GitHub branch.
The distinction between an archived claim and a recomputed finding is essential
to the comparative paper. This file records the former; completed campaign
artifacts will supply the latter.

## What can be recovered now

The upstream Git tree contains analysis code and narrative reports, but no raw
PIP archives or collected result tables: the tracked results directory contains
only `.gitkeep`. The result-explorer notebook contains no stored output. The
available Nathan checkout and a filename search across Nathan's Projects and
vault did not locate the original Iteration 002/003 archives. This is an
availability statement about the inspected locations, not evidence that Lucas's
copies no longer exist. Server inventory is a separate audit.

The original reports were already critical of their own method configurations
and metrics. In particular, the [variable-importance analysis specification](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/variable-importance-analysis.md#L647-L750)
diagnoses coarse-threshold AP and pooling of locus rankings. The current repairs
implement and extend that validity work. They should not be presented as if
every concern was unknown to the original author.

## Historical quantitative claims

These are **reported observations**, transcribed from the immutable original
report. They have not yet been independently recomputed from historical PIPs.
AP and mass-ratio numbers inherit the original estimator, aggregation and
simulation conventions. A percentage range below spans the strata described
by its source; it is not an uncertainty interval.

| Historical claim | Setting and interpretation | Immutable source |
|---|---|---|
| `fb_xregion` AP 0.539, described as best deployable method | Continuous annotations, in-sample LD, Iteration 003 | [Iteration 003, lines 256–257](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L256-L257) |
| Functional BEATRICE → joint prior mass ratio 11.45 → 5.04 (continuous), 9.02 → 5.41 (binary) | Reference-panel size 200; a reduction in excess marginal PIP mass, not proof of calibration | [Lines 248–250](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L248-L250) |
| Joint methods improve mass-ratio proximity to one in 68–94% of cells | Paired within-cell comparison to single-locus Functional BEATRICE | [Lines 251–252](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L251-L252) |
| `fb_xregion` beats `fb_pooled` on AP in 62–97% of cells | LassoNet versus linear shared-prior head; source aggregates by annotation and LD | [Lines 253–255](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L253-L255) |
| `fb_pooled` improves AP in only 6–23% of binary-annotation cells versus single-locus FB | A historical basis for dropping the linear-head method; requires renewed assessment after scoring/inference repair | [Lines 254–255](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L254-L255) |
| Enrichment AP slope 0.0085/fold (`fb_xregion`) versus 0.0070/fold (FB); both approximately zero with panel LD, oracle 0.017–0.022 | Binary annotations; a fitted response of ranking to simulator enrichment | [Lines 258–259](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L258-L259) |
| `fb_xregion` beats FINEMAP AP in only 2–10% of cells outside in-sample LD, less than 2% on the FDR criterion | Panel-LD strata; does not establish FINEMAP universal robustness | [Lines 260–262](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L260-L262) |
| At S=1: mass ratio 6.37 (`fb_xregion`), 9.96 (FB), 2.23 (FINEMAP) | Low causal count; source motivates causal-count/effect-prior tuning | [Lines 263–265](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L263-L265) |
| Adopted tuning trial 401: AP 0.7702, mass ratio 1.73, reliability 1.00, FDR@0.95 0 | Tuning data, not independent test performance; selected among 962 completed trials | [Iteration 004, lines 143–154](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-004.md#L143-L154) |
| AP-oriented tuning point AP 0.8109, mass ratio 8.62, reliability 0.90 | Same tuning study; illustrates ranking/probability-quality trade-off | [Lines 150–154](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-004.md#L150-L154) |

Several analysis scripts contain additional interpretation, such as an ABF
single-cause chart title and FINEMAP described as the best baseline. Without
their generated tables, those strings are not independently verified numeric
rankings. Ties are also important: the old pattern-mining script assigned
`seq_len` ranks after sorting, allowing method order to choose a winner at equal
AP. The corrected comparison uses explicit ties and uncertainty.

## Experiment identity and historical artifact routes

### Iteration 002/003

The Iteration 003 report describes 135 grid rows × 250 scenarios = 33,750
scenarios, with ten regions per scenario. The latter have two regions at each
of p = 100, 200, 400, 500 and 1,000. The grid uses n = 1,000; S = 1, 2, 3, 5,
10; phi = 0.0075, 0.05, 0.1, 0.2, 0.4; ten replicates; annotations absent,
binary or continuous; enrichment 2.7, 5.4, 8.1, 10.8 where relevant; in-sample,
reference-500 and reference-200 LD; and sparse or sparse-plus-infinitesimal
architectures with sparse genetic-variance proportions 0.5, 0.7, 0.9 and 1.
The final proportion duplicates the sparse architecture distribution.

The same cached simulations were used to add `fb_pooled` and `fb_xregion` to
14 existing comparators. This gives paired *data exposure* across those methods,
but the original PIPs still need a corrected evaluator before an updated
historical ranking can be asserted.

Reported provenance:

- PBS supplemental array `3449477[]`, 1,350 tasks; collection `3470549`.
- Original scratch `$EPHEMERAL/fmbench_iter002`.
- Per-scenario `job_*/scenario_*/results_supp.rds` and `evaluation_supp.rds`.
- Durable original local `results/iter003/combined_scenario_metrics.{rds,csv}`
  (reported 367,025 cells), `combined_pip_calibration.rds`,
  `combined_fdr_curves.rds`, `combined_evaluation.{rds,csv}`, `run_summary.csv`.
- Generated report `results/iteration-003-joint-prior-analysis.pdf` (ten pages).
- Earlier analysis defaults `results/iter002/` and `results/iter002_fixed/`.

The report explicitly warns that ephemeral storage could purge after about
30 days from 23 July. Several scripts contain a Lucas-local temporary path
under a `-Users-ls1020` session directory. Those are path clues, not stable
available archives. See [the original provenance block](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-003.md#L216-L244).

### Iteration 004

The upstream report is labelled **running**, not a completed final result. It
reports an early QC pass of 2,608/11,250 scenarios, all 19 methods producing
outputs, and structural Funmap absence in the no-annotation arm. That is an
execution check, not a published ranking.

Iteration 004 changes n to 5,000 and p to two regions each at 500, 750, 1,000,
1,500 and 2,000; phi to 0.05, 0.1, 0.2, 0.4, 0.6; sparse genetic-variance
proportions to 0.2, 0.4, 0.6, 0.8; and retains in-sample LD only. Other axes
give 45 rows × 250 scenarios = 11,250 scenarios. Reported array `3599396`,
scratch `$EPHEMERAL/fmbench_iter004`, 450 tasks, 72-hour limits.

The nineteen method IDs are `functional_beatrice`, `fb_pooled`, `fb_xregion`,
`beatrice`, `polyfun_oracle`, `funmap`, `paintor`, `polyfun_est`,
`polyfun_ldsc`, `sbayesrc`, `susie`, `susie_inf`, `finemap`, `finemap_inf`,
`abf`, `marginal_z`, `carma`, `finimom`, `sparsepro`. The legacy `sbayesrc`
is a local surrogate, not native GCTB/GWFM. The oracle receives simulator
information and must be separated from deployable participants.

The report itself [forbids comparing Iteration 004 and 001–003 numbers as the
same benchmark](https://github.com/lucassiu14/fine-mapping-benchmark/blob/f053fdf318a3abe2c44eeb0ac3f3c99a03979b81/docs/autoresearch/iteration-004.md#L7-L10).
Changing several axes simultaneously prevents attributing a difference to
one repair. The new Gaussian challenge changes the experiment again and
therefore supplies new stress evidence rather than a direct rescoring.

## Fair attribution in the corrected comparison

The paper should distinguish four contrasts, retaining every attempted fit:

| Contrast | Inputs held fixed | What a difference can mean |
|---|---|---|
| Historical rescoring | Historical PIPs, truth, candidate universe, locus identities | Evaluation correction changes the measured ranking |
| Inference repair | Original simulation objects, method arguments, declared stochastic seeds/budgets; original versus repaired implementation | A code repair changes predictions; use the same corrected scorer |
| Simulation repair | Matched source genotypes and random streams where meaningful; original versus corrected phenotype construction | Changed generative assumptions alter the question; not merely rescoring |
| New challenge | Frozen new inputs, participant contracts and held-out assessment labels | Robustness/coverage under new mechanisms; not a reversal of a historical ranking |

If original PIPs are unavailable, regenerate the original grid with the
immutable upstream code, frozen dependencies and source panels. Such results
are **reproduced historical experiments**, not recovered originals. A same-input
original-versus-repaired refit gives a valid implementation contrast; it does
not prove byte-for-byte recovery of Lucas's historical run.

The historical and repaired model versions need separate labels. In particular,
our corrected benchmark and new exact/ABF controls are not independent
inventions of Lucas's Functional BEATRICE or its shared-prior variants. The
vault's supervision history describes a collaborative programme; attribution
must reflect that rather than frame the paper as one person's methods defeating
another person's methods.

### CARMA set structure and the comparison endpoint

The original wrapper flattened `fit[["Credible set"]][[2]]` into one vector,
but the pinned original CARMA commit returns a list of LD-grouped sets. This
loses native grouping and can undercount recovery of distinct causal targets.
The corrected wrapper preserves those groups and labels them
`ld_clustered_pip_mass_atleastone`. This is an implementation correction,
separate from any measured ranking change.

The official construction scans absolute-LD thresholds from .5 to .99,
accumulates descending-PIP members within a neighbourhood until their summed
PIP exceeds the chosen level, and retains the threshold yielding the most
sets. The intended empirical coverage endpoint is whether each reported set
contains at least one causal variant, rather than whether a global union
contains every cause. Summed marginal PIPs give expected causal count; without
additional restrictions they do not establish the joint probability of at
least one cause. Set hit rates, support sizes and distinct-target matching are
therefore reported with that construction disclosed, without inventing a
nominal joint-probability guarantee. See the [pinned CARMA source](https://github.com/ZikunY/CARMA/blob/51f0fcf6a74ee74c75fd2ce70912fdd6cdd214fe/R/CARMA.R#L146-L200)
and the authors' [assessment definitions](https://www.nature.com/articles/s41588-023-01392-0).

Within each setting, report localisation, calibrated discovery, uncertainty
coverage, completion and resource use together. A method can lead in AP and
lose on high-confidence false discoveries. A method unsupported on a declared
input track has no scientific-quality comparison there; its operational
availability still remains visible. Neither a universal overall winner nor an
unweighted count of track wins is justified by an arbitrary stress-suite mix.
