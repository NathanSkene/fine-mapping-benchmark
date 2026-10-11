# FINEMAP native probability and model-size contract

The complete model-size posterior comes from FINEMAP's native log, not from
the sum of its retained `.config` rows. The latter are a truncated collection
of high-probability configurations. Selected-model credible sets retain their
paired native identifier/probability columns.

The inspected server executable has SHA256
`41e0af31c18e2e92d2a5fd1fbf35c9832d62ae737746ecd0aa4629a5925c9bcc`.
The author download label and actual help/console say **1.4.2**, while the
same executable's native `.log_sss` header says **1.4.1**. These are distinct
observed channels. The adapter selects the new SNP precision contract by
executable bytes, never by either version string. Earlier execution artifacts
remain unchanged.

## SNP probability serialization

Two Gaussian integration outputs had printed PIP sums of 0.9999998 and
0.9999995, although their native prior excluded the null and their reported
model-size posterior selected one cause. The previous point-value union check
rejected them. The native literals reproduce the stored numerical PIPs exactly;
this discrepancy is output serialization, not a reason to renormalize scores.

Inspection of the exact executable verified the call from
`SNPs::writeLine`'s probability field to `floatToString(double, 6)`. That routine
constructs a fresh stringstream, sets precision six, and inserts the double
without fixed or scientific format flags. This agrees with C++'s
[general-format numeric conversion contract](https://eel.is/c++draft/facet.num.put.virtuals).
The author manual supplies no independent probability-rounding guarantee, and
the active floating-point rounding mode was not verified. Therefore the new
envelope conservatively allows nearest and directed rounding rather than
assuming nearest rounding.

For a positive reported value q of decimal order k and precision d = 6, let

    s = 10^(k − d + 1).

The upper allowance is s. The lower allowance is s except at an exact decimal
power q = 10^k, where it is s/10: the preceding decimal grid is ten times finer.
Thus a printed PIP of 1 has lower bound 0.999999, rather than 0.99999. Bounds
are intersected with [0,1]. A literal zero remains zero; general formatting
does not round a nonzero double to a zero token. At subnormal scales the helper
uses at least one representable binary subnormal quantum for a nonzero token's
allowance. Decimal order and power boundaries are derived from the literal
token, avoiding `log10` boundary classification errors.

The actual two vectors have conservative summed envelopes
[0.9999965, 1.0000031] and [0.9999962, 1.0000028], both containing one. The
validator tests the PIP/null union inequalities for feasibility within these
aligned envelopes. It does not claim to verify the full joint posterior,
search convergence, or calibration. PIPs, null point probabilities, proper
scores, ranks, and confidence intervals are never replaced by envelope bounds.

The wrapper retains `snp_probability_tokens` in input-variant order. Its
`native_probability_precision` record declares the exact executable SHA,
writer, precision, unknown rounding mode, lower/upper bounds, bijective native
identifier map, and structural null-zero basis. A record hash detects changes
to its representation. Validation also checks token/PIP equality, the native
diagnostics, both identifier maps, and event/support declarations. Known
metadata cannot excuse a genuine union contradiction. Other executable bytes,
missing provenance, or a nonzero null prior receive no new envelope allowance;
malformed supplied metadata is rejected.

**A consistent record is not authentication or custody.** Public
`run_challenge()` and default `evaluate_challenge()` retain strict point-value
validation, including for a client-created or copied descriptor. The assessor
must explicitly set `allow_native_probability_precision = TRUE` only after
establishing custody of outputs from the audited executable and adapter.
The frozen campaign worker enables this option only for its controlled
FINEMAP route after plan/source/runtime/input checks. Successful prediction
diagnostics retain the required probability fields rather than copying large
configuration tables. Failed native diagnostics remain available for audit.

The source-bound integration repair evidence is recorded in
[the repair evidence ledger](validation-artifacts/integration-repairs-20261011/evidence.json).

## Complete model-size posterior

The existing complete-log parser handles model-size probabilities separately.
The earlier observed executable output uses three significant digits and suppresses
trailing zeros. The two annotated installation loci produced posterior sums
0.9996 and 1.00112. Rejecting either solely because its distance from one
exceeds 0.0001 incorrectly rejects a usable native fit.

For a nonzero printed value q with d significant digits, the adapter uses a
half-unit rounding allowance

    h = 0.5 × 10^(floor(log10(q)) − d + 1).

The existing parser policy uses at least three digits; greater printed precision
narrows the allowance. At an exact decimal power the lower half-unit is h/10,
while the upper half-unit remains h. For example, at q = 1 and d = 3 the lower
allowance is 0.0005, not 0.005. Consequently, under this existing half-unit
policy, the tokens 0, 1, 0.0006 cannot form a normalized posterior, while
0, 1, 0.0004 remain feasible. This narrow
boundary correction preserves every reported value and adds no blanket tolerance.
Each interval is intersected with [0,1]. A complete vector
passes normalisation validation only when the summed intervals contain one.
This check permits native output rounding while rejecting incompatible sums,
missing model sizes, duplicate sizes, nonfinite values and values outside [0,1].
The parser **does not renormalise** printed probabilities. The reported sum and
tokens and intervals remain in `additional$k_posterior_rounding`. `best_k` maximises the
reported values; near ties remain limited by native output precision.
This historical three-digit log policy is distinct from the newly inspected
six-digit SNP writer proof. The log's half-unit policy corresponds to nearest
decimal rounding; historical observations and the SNP proof do not separately
establish the log writer's active rounding mode.

## Conditional PIPs and the null event

The same log declares a zero prior on k=0 in this executable's default run.
Consequently, its PIPs condition on at least one modelled cause. The adapter
records the complete reported prior, `native_null_support`,
`posterior_conditioning` and the reported null posterior. A zero posterior
under this prior is not a data-driven rejection of the null. Null-locus false
calls remain useful operational stress outcomes, with this model restriction
made explicit.

The structural zero names the empty configuration within the non-null
conditioned candidate model support. The challenge scorer's generating-null
truth is a different event: scoring that structural zero against a genuinely
null locus is an operational model-support stress test, not calibration of
the same event. Metadata makes this distinction explicit and marks conditional
PIP scores as operational rather than unconditional calibration.

The [author manual](http://www.christianbenner.com/) states that `--prior-k0`
affects the posterior number of causes but not fine-mapping itself. A nonzero
count-null prior therefore does not establish that `.snp` PIPs are a joint,
null-including marginal posterior. The wrapper keeps this conditioning
distinction and provides no new serialization relaxation for that regime.

## Historical execution evidence and limits

PBS job `4310528.pbs-7` reran the actual executable through the repaired R
adapter on two N=200, p=24 annotated loci. Both completed, kept the native
probabilities and imported valid identifiers and sets. The small
[execution proof](../scripts/campaign/finemap-parser-proof.json) records the
executable, tested wrapper and output hashes. This establishes execution and
parsing; it makes no scientific accuracy or search-convergence claim. A later
source change requires fresh source-bound preflights.

The native help lists no seed option, and this wrapper does not pass an external
random seed. R's random seed alone therefore does not establish reproducible
native stochastic search. Binary identity and observed-input identity are
recorded separately from that unresolved RNG contract.

Focused regression cases are in
`tests/testthat/test-native-parser-events.R` and
`tests/testthat/test-native-probability-precision.R`. They cover actual rounded
vectors, decimal powers, generic client strictness, tampered provenance,
genuine contradictions, unchanged scores, and controlled parser alignment.
The FINEMAP semantics discussed
here are distinct from the [PAINTOR null and enumeration contract](paintor-native-contract.md)
and the [general evaluation contract](evaluation.md).
