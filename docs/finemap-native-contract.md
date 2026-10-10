# FINEMAP model-size posterior contract

The complete model-size posterior comes from FINEMAP's native log, not from
the sum of its retained `.config` rows. The latter are a truncated collection
of high-probability configurations. Selected-model credible sets retain their
paired native identifier/probability columns.

The actual server executable reports FINEMAP **1.4.2**. Its SHA256 is
`41e0af31c18e2e92d2a5fd1fbf35c9832d62ae737746ecd0aa4629a5925c9bcc`.
It prints model-size probabilities with three significant digits and suppresses
trailing zeros. The two annotated installation loci produced posterior sums
0.9996 and 1.00112. Rejecting either solely because its distance from one
exceeds 0.0001 incorrectly rejects a usable native fit.

For a nonzero printed value q with d significant digits, the adapter uses a
half-unit rounding allowance

    h = 0.5 × 10^(floor(log10(q)) − d + 1).

The format contract uses at least three digits; greater printed precision
narrows the allowance. Each interval is clipped to [0,1]. A complete vector
passes normalisation validation only when the summed intervals contain one.
This check permits native output rounding while rejecting incompatible sums,
missing model sizes, duplicate sizes, nonfinite values and values outside [0,1].
The parser **does not renormalise** printed probabilities. The reported sum and
intervals remain in `additional$k_posterior_rounding`. `best_k` maximises the
reported values; near ties remain limited by native output precision.

The same log declares a zero prior on k=0 in this executable's default run.
Consequently, its PIPs condition on at least one modelled cause. The adapter
records the complete reported prior, `native_null_support`,
`posterior_conditioning` and the reported null posterior. A zero posterior
under this prior is not a data-driven rejection of the null. Null-locus false
calls remain useful operational stress outcomes, with this model restriction
made explicit.

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
`tests/testthat/test-native-parser-events.R`. The FINEMAP semantics discussed
here are distinct from the [PAINTOR null and enumeration contract](paintor-native-contract.md)
and the [general evaluation contract](evaluation.md).
