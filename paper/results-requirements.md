# Evidence needed to finalise the comparative paper

This is the completion contract for `manuscript.md`. It is not a substitute
for executing the user's full benchmarking campaign. The paper remains a
working draft until the results and provenance below are present and checked.

1. **Method inventory:** every original registered method, its original/native
   identity, intended input regime, arguments, dependency commit/version,
   execution status, warnings, failures, timeouts and nonconvergence. Include
   optional methods, original local surrogates, and native implementations as
   distinct rows. Do not rename one to imply fidelity to the other.
2. **Original-comparison evidence:** either recover historical per-locus PIPs
   and truth, or regenerate the historical experiments at the frozen original
   source. Identify which happened. Original and corrected scorers on identical
   predictions isolate evaluator effects; repaired method code on identical
   inputs isolates the combined implementation/adapter/acceptance change.
   A particular repair needs a separate ablation. New simulation data cannot
   establish either contrast alone.
3. **Full experiment identity:** original-grid family, repaired genotype/OLS
   family, new Gaussian stress family and native genome-wide family each need
   frozen design files, source/input hashes, resource limits, seeds, variant
   universes, truth definitions, probability-event alignment, exposure history
   and independent cluster IDs. The native 128-block comparison is an empirical
   subset, not dense genome-wide coverage.
   Historical-grid p is a requested design size: report realized candidate
   distributions and filter losses alongside it. Ns are synthetic experiment
   parameters: legacy N counts simulated individuals, while Gaussian N is a
   score-design parameter; variant-wise N distributions require the committed
   public inputs. None is a published GWAS participant total.
   Reference-panel regime rows are not isolated LD-only pairs because
   their genotype draws and retained candidates can differ.
4. **Final assessment:** hyperparameters and contrasts fixed before final truth
   inspection; author defaults, controlled ablations and competitive tuning
   reported separately. The exact reference and true-prior controls are labelled
   controls, not ordinary competitors. All attempts are retained.
5. **Completed output archives:** per-attempt predictions and status; variant
   alignment; raw and assessed results; task-level and grouped metrics;
   posterior-set semantics; resource measurements. A submission/job ID alone
   does not prove completion. Record native posterior conditioning and null
   support, including conditional explored-support set events. Every reported
   null forecast must identify its candidate-model or generating-universe
   event; operational proxy scores are separate from same-event calibration.
6. **Primary evidence tables:** per-setting mean/macro AP with ties; attempted
   recovery and completion; PIP raw Brier/log loss; pooled FDR and discovery
   counts at fixed prespecified PIP thresholds; null false-call rates; native
   set event hit rates, distinct-target recovery and sizes; full-configuration
   coverage where declared; CPU/wall time and memory; chain assessability,
   diagnostics and final-summary/trace coherence. Native imported outcomes
   remain separate from diagnostically qualified accuracy; four imported
   chains alone cannot establish mixing or permit a native winner claim.
7. **Paired uncertainty:** use the same independent locus/panel clusters for
   every method in a contrast. Preserve phenotype reuse and shared-prior
   dependence. Report effect differences and confidence intervals, not only
   win counts. State whether intervals condition on existing LD or generalise
   to new genomic loci. Treat all unregistered contrasts as exploratory.
   The separate regional stage specification fixes implementation/acceptance
   and generating-experiment contrasts before full-campaign outcomes, after
   development outcomes were seen. Preserve endpoint-specific finite pairs,
   availability, recovered/lost fits, failure-zero AP and unmatched attempts;
   a selected-out stage is unassessed. One shared-prior panel gives no interval.
   Historical per-p strata usually contain two genotype clusters regardless
   of phenotype repetition count: their percentile envelopes are descriptive,
   have no validated nominal coverage and cannot establish statistical
   distinguishability. Generalising across genomic loci needs more independent
   genotype panels, rather than more repetitions on the same loci.
8. **Historical claim checks:** for each claim in `original-results-audit.md`,
   classify as reproduced, revised by scoring, revised by inference, revised
   by simulation, untestable from available historical data, or outside the new
   experiment's scope. Do not force a binary confirm/refute answer.
9. **Setting-specific guidance:** explain which methods lead on which declared
   endpoint, under which settings, with trade-offs and uncertainty. Native
   genome-wide learning must use a sufficiently broad input universe; a small
   identity-LD smoke run cannot establish the 2026 GWFM paper's advantages.
10. **Paper delivery:** replace the explicit working-draft results section with
    verified tables/figures/narrative; update abstract and conclusions to match
    measured data; retain unsuccessful methods and limiting evidence; validate
    all quantitative claims against the frozen files; render and visually
    inspect the final PDF. Confirm author/contribution attribution with the
    user before publication, not by inventing an author list.

Suggested primary contrasts, chosen before final results are inspected:

- Corrected minus original AP on the same predictions, especially compressed
  score distributions and exact ties.
- Original versus repaired learned-prior posterior calculation on the same
  inputs, with final-prior consistency as a mathematical diagnostic.
- Conditional-union set selection versus the original summed-PIP diagnostic
  ablation on the same explored support and conditioning event. Report whether
  the conditional union target was met; do not claim an unconditional or
  exhaustive-model coverage guarantee.
- Joint `fb_xregion` versus single-locus Functional BEATRICE, preserving the
  controlled likelihood/optimisation budget and documenting prior reuse.
- Functional methods versus their annotation-free counterparts with useful,
  irrelevant and reversed annotations; separate fixed-prior transfer from
  annotation learning.
- Sparse and sparse-plus-infinitesimal inference within in-sample/compatible
  LD, then separately as out-of-scope stress under external/mismatched LD.
- Null-aware versus forced-signal controls for null and selected-null inputs.
- Native GCTB mode versus same-input regional strategies, separated from any
  genome-wide information advantage; report chain diagnostics before accuracy
  claims.
   Regional publication figures require their sealed locus metric CSV/type
   schema and endpoint-specific realized-size/truth census. Retain generating
   nulls, observed-positive eligibility, missing causes, failures and infinite
   losses with their proper denominators. Paired census metadata must use the
   exact source contributing mask, without altering scores or intervals.

For manuscript phrasing, "leads" means the strongest observed result on a
named endpoint; "statistically distinguishable" additionally requires an
appropriate paired uncertainty contrast. Avoid a universal "winner" unless
the target distribution and endpoint weights were declared in advance.
