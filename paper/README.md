# Comparative manuscript

`manuscript.md` is the authored paper source. It is explicitly a working
draft while the full optional/native campaign runs. It contains the
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

```sh
Rscript --vanilla paper/plot_results.R /absolute/campaign/analysis /absolute/figures
```

Choose publication panels with `--setting-regex PATTERN` and, where useful,
`--methods abf,susie,finemap`. These filters and source hashes are recorded in
the figure manifest. The default export limit is 100 settings and 2,400 paired
contrasts; exceeding it requires a narrower declared view. Filters select
already computed outcomes and do not retune methods or recompute intervals.

Publication-mode exports require a complete full-profile collection.
`--allow-incomplete` enables watermarked installation/development views, which
cannot replace the paper's final campaign results.

Do not remove the working-draft status, add a universal winning method,
or describe the full campaign as complete until its frozen results and
terminal scheduler evidence satisfy the completion contract. Human
publication authorship and collaborator approval have not been assigned.
