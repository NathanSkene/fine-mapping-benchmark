# Figure export validation

This record documents software and presentation checks. It contains no full
scientific campaign rankings.

The current `plot_results.R` exported collector tables from the restored sealed
integration at `results/campaign-final-integrity-fixture/plan-v2/analysis`.
That run contains 28 terminal tasks and 292 locus rows. Selecting the
`allele_flip|original_refit` settings and `abf,sbayesrc,susie` methods produced
nine setting pages and 36 paired-contrast pages under
`results/paper-plot-regional-final-v5/`.

Every page says **DEVELOPMENT ONLY: profile integration collection complete**.
The figure manifest marks `scientific_full_profile=false` and retains the
explicit filters. Its input, output and plotting-script SHA256 values were
independently checked against the saved bytes. Setting and contrast pages were
visually inspected, including event-misaligned open Brier markers, conditional
successful-fit AP, counted attempts, paired cluster counts and footer spacing.
The intervals are read from the collector rather than reconstructed by the
figure script.

A separately labelled synthetic native-format software fixture tested native
schema handling: one setting page and two paired-contrast pages under
`results/paper-plot-fixtures/native-plot-v7/`. Every page is marked development
and incomplete, and says native accuracy eligibility is unavailable for this
older synthetic schema. Panel captions describe imported fits. These arbitrary fixture values test plotting behaviour and
must never be incorporated as scientific results.

The exporter rejected publication-mode export of development inputs and a
changed sealed CSV, produced no false figure directory, and released its
collector lock after rejection. Current full-campaign scientific plots remain
pending verified full-profile collected evidence.

The current exporter requires explicit native accuracy eligibility for
publication figures and rejects finite ordinary native accuracy values in
unqualified rows. Native completion remains an execution/import measure.
The updated setting panels use compact wrapped headings and a two-line footer;
the native fixture was visually inspected after correcting heading overlap.

The final streamed exporter reads 10,000 CSV rows at a time, applies the
declared selection before retaining rows and enforces the existing setting
and paired-contrast limits. Eleven parser/bounds checks cover quoted commas,
quotes, embedded newlines, missing/infinite outcomes, empty selections and
150,000 rows with a small selection near EOF. All values and strings selected
from the fresh sealed streamed integration at
`results/campaign-bundle-cli-integration-5/campaign/analysis` match the
full-memory reader at zero numeric tolerance (four summary rows and eight
paired rows); integer/double storage types may differ without changing values.
All ten rendered pages match the full-memory version of the same plotting
code pixel for pixel. The final export is at
`results/paper-plot-streamed-final/`, with two settings and eight paired pages.

The final native-format fixture export is at
`results/paper-plot-fixtures/native-streamed-final/`. It remains explicitly
development/incomplete and says native eligibility is unavailable. Figure
input/output/helper hashes and both development labels were verified. Native
eligibility checks run on every chunk before filters, so an invalid
unqualified row cannot be hidden by selecting other methods. These checks
validate presentation and bounded parsing; they add no scientific rankings.
