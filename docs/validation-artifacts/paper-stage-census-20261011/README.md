# Regional publication and stage-analysis validation

The [source-bound proof](proof.json) records 285 passing checks: 11 CSV-reader,
63 unchanged native-publication, 85 regional-census, 96 stage-mathematics and
30 stage-selection/provenance assertions. A prior sealed development collection
also passed the export: its 148 locus rows produced 940 endpoint census rows,
94 summary rows, 540 unchanged method-contrast rows and 24 separate analytical
stage-contrast rows. These are software and interpretation checks, not full
campaign results or a method leaderboard.

The [stage specification](../../../paper/stage-attribution-spec.json) was fixed
after development outcomes had been seen and before full-campaign outcomes.
It separates the implementation/adapter/acceptance bundle on fixed observed
inputs and truth from a changed phenotype/truth experiment on fixed genotypes.
Selected-out stages stay unassessed. Scores, existing method contrasts and
their intervals are preserved; the additional stage table computes only its
specified paired differences and fixed cluster bootstrap.

The reconstructed historical grid has two genotype regions per requested-p
class, reused across phenotype repetitions. Shared-prior methods use one
whole-panel cluster. A two-cluster symmetric-effect counterexample shows why
the resulting percentile envelopes cannot justify nominal 95% coverage.
Exports therefore retain the numbers as descriptive envelopes, record cluster
counts, mark coverage unvalidated and make no statistical-distinguishability
or multiplicity-controlled claim. One cluster gives no interval.

The census preserves requested/retained sizes, endpoint-specific prevalence,
nulls, missing causes, failures and infinite losses. Its synthetic N labels
distinguish simulated individual counts from Gaussian score-design parameters;
variant-wise N distributions require the committed public inputs. Reference
panel regimes do not isolate an LD-only intervention.

The [render record](render-and-review.json) binds the updated working draft,
successful PDF and compiler checks, and independent review. Raw fixtures,
logs, PDFs and the complete development export remain in ignored local
validation directories. Logical labels and hashes are not download URLs.
The server's inference source remains frozen at `0e4a7ead774b`; its fresh
integration, resource-probe, native-convergence and full-campaign gates are
separate from this publication-only change.
