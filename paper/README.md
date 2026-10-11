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
uses the collector's existing paired intervals. It does not recalculate a
bootstrap or invent an overall leaderboard. Native points retain independent
genotype-draw outcomes. The figure manifest records every input/output hash.

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
already computed outcomes and do not retune methods or recompute intervals.
The CSV reader scans 10,000 rows at a time and retains only selected rows,
with a separate 100,000-row retention limit for setting tables. The helper's
source hash accompanies the main plotting-source hash. Large unselected paired
tables therefore do not need to fit in RAM. Run its parser/bounds checks with
`Rscript --vanilla paper/test_csv_reader.R`.
Run the native truth-census, event-metadata, export and gate checks with
`Rscript --vanilla paper/test_native_publication.R`.

Publication-mode exports require a complete full-profile collection.
`--allow-incomplete` enables watermarked installation/development views, which
cannot replace the paper's final campaign results.

Do not remove the working-draft status, add a universal winning method,
or describe the full campaign as complete until its frozen results and
terminal scheduler evidence satisfy the completion contract. Human
publication authorship and collaborator approval have not been assigned.
