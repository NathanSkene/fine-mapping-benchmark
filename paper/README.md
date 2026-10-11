# Comparative manuscript

`manuscript.md` is the authored paper source. It is explicitly a working
draft while the full optional/native campaign is prepared and executed. It contains the
mathematical and experimental comparison contract, primary references,
respectful historical attribution and verified development evidence.
Its final-result requirements are in `results-requirements.md`.

`manuscript.tex` is generated standalone LaTeX for the Codex source editor
and live PDF preview. The built-in editor compiler successfully compiled
the first draft. Export a local PDF with `bash paper/render.sh` using Pandoc
and pdfLaTeX; the PDF and build logs are ignored under `results/`.

The historical evidence and archive search are documented in
`docs/original-results-audit.md`. No original PIP archive has been recovered
from inspected Git/local locations, so the planned source-frozen replay
must be labelled as a reproduction rather than an exact recovery.

`plot_results.R` exports setting-specific PDF figures directly from a sealed
collector analysis. It validates the analysis commitment and table hashes,
keeps methods in alphabetical order, preserves event-alignment warnings, and
uses the collector's existing paired method intervals. A separate regional
stage table adds the fixed analysis specified below; it does not replace the
collector's method estimates or intervals. There is no overall leaderboard.
Native points retain independent genotype-draw outcomes. The figure manifest
records every input/output hash, including the stage helper and specification.

Regional exports require the sealed `task-metrics.csv` and its
`task-metrics-schema.rds` in complete publication mode. The schema supplies
column types and the plan-specific missing-value marker; the census reads
17-digit decimal tokens and uses locus scores only to identify endpoint
availability, finiteness and the existing contributing masks. It does not
recompute scores. `regional-target-census.csv` records counts and minimum,
quartile, median and maximum distributions of requested/retained p, candidate
gaps, target prevalence and synthetic experiment N for each endpoint.
The `n_basis` distinguishes actual simulated GWAS/reference individual counts
in the legacy experiment from the Gaussian score design parameter `n`.
Gaussian scores are generated directly, without simulated individual rows.
The `meta_missingness` and `meta_known_overlap` tracks use variant-specific
sample sizes; the locus CSV records the design N, not their distribution.
That distribution is explicitly unavailable here: the committed
`datasets/generation_<generation_id>/public/inputs.rds` contains `n` and
`n_per_variant`, reached through the plan's dataset-to-generation mapping and
checked by the public release manifest. No experiment N is a
participant count from a published GWAS.
The retained prevalence is observed positive candidates / retained candidates.
Generating-null loci and missing causes are counted separately: the legacy
filter acts before causal selection, while deliberate Gaussian candidate
restrictions can remove a cause. A Gaussian size gap is not labeled as the
legacy monomorphism filter.

The helper streams chunks, writes selected rows to temporary partitions and
processes at most 100,000 rows from one selected design stratum at a time.
That bound covers all selected stages and methods together. The optional
`on_stratum(d, stratum)` callback receives the validated, typed source rows
plus private normalized metadata for exactly that selection. It supplies no
excluded stages; downstream stage attribution must explicitly leave a
partially selected stage comparison unassessed.
Duplicate identities, inconsistent dimensions/truth, summary denominators or
paired denominators fail closed. `regional-summary-publication.csv` preserves
the collector's values and adds the all-attempted census;
`regional-paired-publication.csv` preserves every source contrast/interval and
adds each side's census on exactly its finite-difference pairs. AP uses finite
observed-positive scores; recall includes recorded failure zeros; log-loss
summary support includes infinity. The failure-zero AP composite retains only
observed-positive attempts. Brier point alignment uses its contributing locus
mask; unknown alignment is not displayed as established alignment.
Sparse membership and any-nonzero total-effect
censuses are separate, including the all-positive no-data baseline and the
`p_causal=1` exception as determined by actual counts. Dense alignment is not
inferred from the source's sparse-target flag. Missing census evidence in a
development export is explicitly unassessed.

Setting labels retain **requested** p; captions give retained-p ranges.
Reference-panel rows use different genotype draws and can retain different
candidates, so they are reference-panel regime comparisons, not isolated
LD-only paired contrasts. Regional summary Brier is a mean over locus scores;
pooled evaluator Brier weights variants. The existing paired table compares
methods within a stage. `regional-stage-attribution.csv` separately compares
corrected refits with original refits, then corrected phenotype resimulation
with corrected refits. [The separate helper](stage_attribution.R) and its
[registered specification](stage-attribution-spec.json) distinguish the whole
implementation/acceptance change from a changed generating experiment.
The specification was frozen before full-campaign outcomes; development
results had already been seen. These are analytical differences of stored
scores, with no refitting or rescoring. The implementation contrast holds
observed inputs and truth fixed, but attribution to a particular repair needs
a separate ablation. The resimulation contrast holds genotypes and design
fixed; it cannot alone establish an inference improvement or reversal of
Lucas's historical rankings.

The stage table reports availability, AP on finite paired scores, AP with
failure zeros on observed-positive attempts, Brier, attempted recall and
distinct-target credible-set recall. Endpoint denominators, recovered/lost
fits and unmatched attempts remain explicit. A fixed 1,000-draw percentile
bootstrap resamples declared independent panel/region clusters; methods that
share a learned prior use whole panels. A single panel gives no confidence
interval. The [reconstructed historical grid](../scripts/hpc/generate_params_grid.R)
has two loci per requested-p class and reuses
their genotypes across phenotype repetitions. Thus these non-pooled strata
usually have only two independent genotype clusters, and pooled strata one.
With two independent symmetric cluster effects, half of samples have equal
signs and their bootstrap envelope excludes the true zero mean: its nominal
95% coverage is not justified. The export therefore labels intervals as
descriptive envelopes, records the low cluster count, and does not claim
statistical distinguishability from them. These are pointwise envelopes,
without multiplicity control across the grid. Their interpretation conditions
on the supplied genotype design. Brier alignment
is checked on its exact contributing pairs and is an operational diagnostic
when unqualified. AP among successful pairs can favour methods that fail on
difficult loci, so it must be read alongside availability and failure-zero AP.

`regional-stage-selection.csv` records whether each selected stratum/method
contains all three stages. Excluding a stage with `--setting-regex` leaves the
comparison unassessed, rather than treating excluded rows as failures.
Gaussian score tasks have no implementation-stage pair. Complete collections
must have matched attempted locus identities. Incomplete exports retain their
unmatched counts and development labels. Scientific stage analysis is bound
to the exact full-plan SHA in the specification; a different full plan needs
its own specified analysis. This additional export is bounded at 100,000 rows.

Native views stream the sealed `task-metrics.csv` truth counts into bounded
totals for the selected draw groups. A complete publication export requires
that census and checks its locus, successful-fit and completed-fit counts
against the draw summaries. No scores or intervals are recomputed. The derived
`native-target-census.csv` and `native-draw-publication.csv` preserve the source
values and record whether each locus has positive and negative candidates.
`native-paired-publication.csv` preserves the collector's contrasts and adds
each side's probability event, alignment, conditioning and truth census.
Paired annotations use exactly the finite, complete draw pairs included in
that stored metric contrast and verify its paired-draw count. Separately
prefixed `selected_scope_` columns retain the wider selected-draw census,
including draws omitted from a particular endpoint. Neither view changes the
source numerical estimates or intervals.

All-positive AP is 1 regardless of ranking. Recall on those loci is a
threshold-call fraction, with no localisation discrimination. Figures mark
these endpoints with triangles and explain the constant-one baseline: without
data, PIP 1 gives AP 1, recall 1, Brier 0 and log loss 0 for an all-positive
target. Dense-event probability scores and sparse-membership localisation
remain separate views. The software derives these warnings from the observed
sealed census, rather than a model-name assumption. Older development fixtures
without truth counts are labelled unassessed and cannot pass the publication
gate.
Artificial validation fixtures carry an additional explicit synthetic label
on each PDF; ordinary measured integration views keep their development label.
An explicit runtime test-fixture marker cannot qualify as scientific evidence
by changing the plan's profile to `full`.

```sh
Rscript --vanilla paper/plot_results.R /absolute/campaign/analysis /absolute/figures
```

Choose publication panels with `--setting-regex PATTERN` and, where useful,
`--methods abf,susie,finemap`. These filters and source hashes are recorded in
the figure manifest. The default export limit is 100 settings and 2,400 paired
contrasts; exceeding it requires a narrower declared view. Filters select
already computed outcomes and do not retune methods or recompute the existing
method intervals. Stage differences use the same fixed specification for
every eligible selected group.
The CSV reader scans 10,000 rows at a time and retains only selected rows,
with a separate 100,000-row retention limit for setting tables. The helper's
source hash accompanies the main plotting-source hash. Large unselected paired
tables therefore do not need to fit in RAM. Run its parser/bounds checks with
`Rscript --vanilla paper/test_csv_reader.R`.
Run the native truth-census, event-metadata, export and gate checks with
`Rscript --vanilla paper/test_native_publication.R`.
Run the regional size/truth census, endpoint masks, bounded streaming and
publication gate checks with `Rscript --vanilla paper/test_regional_publication.R`.
Run the stage mathematics and selection/provenance checks with
`Rscript --vanilla paper/test_stage_attribution.R` and
`Rscript --vanilla paper/test_stage_publication.R`.

Publication-mode exports require a complete full-profile collection.
`--allow-incomplete` enables watermarked installation/development views, which
cannot replace the paper's final campaign results.

Do not remove the working-draft status, add a universal winning method,
or describe the full campaign as complete until its frozen results and
terminal scheduler evidence satisfy the completion contract. Human
publication authorship and collaborator approval have not been assigned.
