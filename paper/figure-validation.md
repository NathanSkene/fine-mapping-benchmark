# Figure export validation

This record documents software and presentation checks. It contains no full
scientific campaign rankings.

The current `plot_results.R` exported collector tables from the restored sealed
integration at `results/campaign-final-integrity-fixture/plan-v2/analysis`.
That run contains 28 terminal tasks and 292 locus rows. Selecting the
`allele_flip|original_refit` settings and `abf,sbayesrc,susie` methods produced
nine setting pages and 36 paired-contrast pages under
`results/paper-plot-regional-final-v3/`.

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
`results/paper-plot-fixtures/native-plot-v5/`. Every page is marked development
and incomplete. These arbitrary fixture values test plotting behaviour and
must never be incorporated as scientific results.

The exporter rejected publication-mode export of development inputs and a
changed sealed CSV, produced no false figure directory, and released its
collector lock after rejection. Current full-campaign scientific plots remain
pending verified full-profile collected evidence.
