# Integration repair evidence — 11 October 2026

The small integration campaign exposed six unexpected output-validation
failures, alongside 240 declared input or algorithm absences. This ledger
preserves that outcome and records the separate repair checks. It contains no
method ranking, convergence certificate or repaired full-campaign result.
The machine-readable [evidence ledger](evidence.json) gives logical artifact
labels, SHA-256 hashes, source bindings and the limits of each check.

The integration used immutable source
`eb2d9bc34190c9fc89cfa25b43ecb2f070139fc8` and plan
`f0cfb18aadab117ac4bec9a48d158fb8575ef66c484b05f0a0ba7ef327dd8724`.
All 168 task receipts and 282 dataset records were present. Of 1,404 locus
attempts, **1,158 succeeded and 246 recorded a failure**:

| Recorded outcome | Locus attempts | Interpretation |
|---|---:|---|
| Scalar sample size unsupported | 152 | Nineteen methods across eight non-scalar-N challenge loci |
| Gaussian score input unsupported by INF adapters | 72 | Thirty-six each for SuSiE-inf and FINEMAP-inf |
| Exact/forced-one reference lacks required overlap information | 4 | Declared reference-model absence |
| Funmap input has no annotations | 12 | Four loci in each of three legacy stages |
| FINEMAP serialized-probability coherence rejection | 2 | Rounded native PIPs; diagnosed below |
| SparsePro native credible-set parser rejection | 4 | Three challenge loci and one corrected-resimulation locus |

The completion census passed its receipt and identity checks. Passing that
operational gate does not erase the six unexpected failures or certify the
scientific predictions. The original predictions, errors and sealed analysis
remain unchanged.

## Known-overlap construction

A separate replay of an actual planned `p = 1000, n = 10000` generation case
found 10,506 zero and 10,956 subnormal LD entries. Recovering overlap by
dividing the noise covariance by LD produced 10,506 NaNs and a maximum
subnormal error of 0.2928932188. The repair constructs the overlap matrix
directly as the normalized Gram matrix of the study-membership vectors,
`V_j = (1, a_j) / sqrt(1 + a_j)`. Thus `W = V V^T` has unit diagonal and
rank at most two, and `Cov(z) = LD * W` by entrywise multiplication.

The sealed replay verified that every other public-input field and every truth
field was identical. The input validator and PSD thresholds were unchanged.
This is one exact generation replay; it does not establish completion of the
larger panel or downstream method success. The earlier replay is retained;
the v2 report additionally covers the subnormal distortion.

## FINEMAP serialized probabilities

PBS `4314729.pbs-7` completed with exit status zero. Its two native PIP vectors
were exactly identical to the archived vectors, with sums 0.9999998 and
0.9999995. Inspection of the exact executable
`41e0af31c18e2e92d2a5fd1fbf35c9832d62ae737746ecd0aa4629a5925c9bcc`
established general-format six-significant-digit SNP-probability output.
The active rounding mode was not separately measured. Conservative adjacent
last-place envelopes cover nearest and directed rounding, including the
asymmetric spacing at decimal powers. Both vectors admit a coherent underlying
union probability of one. PIPs, scores and native tokens are preserved exactly.

This allowance is an explicit assessor choice with public evaluation **strict
by default**. Only the controlled FINEMAP adapter route can opt in after its
source, executable, literal tokens, identifiers and excluded-null support have
been checked. A self-consistent descriptor is not proof of custody. Unknown
executables and generic client predictions retain strict point validation.
The zero null probability here is structural: native SNP PIPs condition on
nonempty modeled configurations. It is not a calibrated forecast of the
whole generating locus being null. See the
[native contract](../../finemap-native-contract.md).

The same executable prints `1.4.2` to help/stdout and `1.4.1` in its native
log header. The v2 proof records those channels separately and preserves the
earlier proof. The model-size log parser's existing three-digit half-unit
policy is a separate contract; its lower allowance at a decimal power is now
correctly ten times smaller. Its active log-writer rounding mode remains
unverified.

## SparsePro transport and strict parser

PBS `4314778.pbs-7` completed with exit status zero. Five actual native output
probes covered the four reported cases and a deterministic multi-member-set
case. Paired parsers consumed the same saved outputs; imported PIPs were
identical to native values. Slash-separated set members were supported by the
repaired parser. This v2 proof binds its own wrapper snapshot, not the later
strict wrapper. Its preceding v1 harness failure is retained and is not counted
as a successful native check.

The final parser also validates set-column schema, member/probability/effect
alignment and corrupted text bytes, and rejects ambiguous input identifiers.
The v3 parser-only harness (`4314810.pbs-7`) failed before parsing because its
Linux symlink check treated a missing output path as occupied. Its code,
submission lineage and errors are preserved; it executed no native program.
The replacement v4 proof (`4314827.pbs-7`, exit status zero) passed on all five
sealed outputs. Every case reproduced all three original written-input hashes;
both output-copy stages matched exact bytes and SHA-256 hashes. Independent
readers confirmed PIPs, set members, group probabilities and effects, with
exactly three native set fields. All 40 output artifacts were verified. The
five intercepted calls executed **zero new SparsePro inference runs**. No
benchmark prediction has been overwritten by these probes.

## Remaining gates and historical evidence

The resource probe (`4314711.pbs-7`, exit status one) exposed a fourth,
separate issue: an original-generator draw can lose monomorphic variants after
simulation, so realized locus size can be smaller than the nominal request.
The audit now records requested, realized and dropped counts/fractions,
GWAS/reference sample sizes and every warning. It permits known filtering
only for upstream genotype-source SHA
`651dbb1640b63187028c12406a82801a676f10be1f9a52b635aedad5a49be139`.
Canonical VCF MAF undersupply and bundled-VCF capping warnings remain fatal,
as do unexplained dimension loss and contradictory matrices or identifiers.
Genotype/phenotype payloads and RNG are preserved; there is no resampling.

All **74 generation-audit assertions** passed, as did 147 existing campaign
contract and 81 precision assertions in the source-bound 302-assertion proof.
Independent review repeated the 74 focused checks. These are software checks,
not a regenerated full panel or a new resource estimate. The old failed
attempt is preserved and a fresh generation execution gate is still required.

The final repaired-source v3 package check passed **Status: OK**, with **1,824
passing expectations, 24 explicit optional/repository-only skips and zero
errors, warnings or notes**. `NOT_CRAN=true` enabled the pipeline test. Separate
source-only native orchestration checks passed **218 expectations**, with no
failures, errors, warnings or skips. The exact archive SHA is
`6d3812db885521a290bd03a6625e396e85caadf9e207b75f399ec0d1efe02ffd`;
all 159 archive file hashes and 264 current scientific/source files were
independently checked against their recorded manifests. This certifies the
tested software snapshot, not native accuracy or full campaign execution.

A new immutable source snapshot, fresh plans and their preflight, integration
and resource gates are required before campaign launch.
An intermediate staged-v2 package check passed on its own archived source;
later timer-accounting and generation-audit edits mean it is not the final
repaired-source check. FINEMAP's executable verification remains outside its
native elapsed timer, while task-level PBS costs include that verification.
The older native integration quality check passed its specified execution
contracts with 252 event rows, 59 supported regional cells and one structural
absence. Its 20 retained draws per chain do not establish convergence.

Most valuable for the historical comparison are **Lucas's original raw
simulation objects and per-method prediction/result files**, with seeds and
configuration logs. Collected tables and figures are secondary. The
[original-results audit](../../original-results-audit.md) explains the recovery
request and attribution limits. Until those files are recovered, original-code
replay is a new experiment, not direct rescoring of historical results.

Full raw files are retained in ignored execution archives. Logical labels and
hashes here identify those records; they do not make the underlying data or
server environment publicly downloadable. These compact records accompany the
[validation record](../../validation.md) and
[artifact index](../README.md).
