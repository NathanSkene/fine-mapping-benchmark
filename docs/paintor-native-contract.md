# Audited PAINTOR native probability and resource contract

This audit concerns upstream `gkichaev/PAINTOR_V3.0` commit
`ffce82a7114eeb98220db39ba5d6b69e5a419a9e`, inspected on 2026-10-10. The
repository directory says V3.0 while source headers say V3.1; use the commit
and byte hashes to identify this implementation. Source inspection establishes
its contract; a campaign's runtime fingerprint separately identifies its actual
executable. It does not validate arbitrary binaries supplied to the wrapper.

## Source identity and probability support

The inspected server checkout is
`/rds/general/user/nskene/home/fmbenchmark-campaign-20261010/tools/PAINTOR_V3.0`.
Its native modeling and I/O sources match the upstream commit. The checkout
contains installation changes and generated objects under bundled NLopt.

| Artifact | SHA256 |
|---|---|
| `Functions_model.cpp` | `115f3106a12881af3d88744d60664b0b50506c9229f52da7e5e833b2d0c37731` |
| `Functions_model.hpp` | `32382bceb103c97faf94d3a651160abfc3adfe57b6a424b8871b2672ad0e4906` |
| `Functions_IO.cpp` | `8748c2ba366dbc360f3584725b0134a3077407b8c6a04bbb76a39fe6b118207e` |
| `main.cpp` | `8df0ca9e0ea599c79edc672e865ae292b72ca74c08fe7fd672d07b85f3f21a4a` |
| Server executable `PAINTOR` | `f7affb3fde82913177e2d61aacf92ad447694004554a428f3ccc3c935340b67b` |

Enumeration initializes zeros at [Functions_model.cpp:337–343](https://github.com/gkichaev/PAINTOR_V3.0/blob/ffce82a7114eeb98220db39ba5d6b69e5a419a9e/Functions_model.cpp#L337-L343), then
advances the configuration before evaluating it at [lines 124–151](https://github.com/gkichaev/PAINTOR_V3.0/blob/ffce82a7114eeb98220db39ba5d6b69e5a419a9e/Functions_model.cpp#L124-L151). The all-zero state is consequently absent. At a two-effect
limit, support consists of `p` singleton and `choose(p, 2)` pair configurations.
For `p=3`, an independent traversal produces `{1}`, `{2}`, `{3}`, `{1,2}`,
`{1,3}`, `{2,3}`: six configurations, without null. The regression test also
checks `p=2,5,12` against independently listed subsets. This establishes the
conditioning of enumeration PIPs on at least one modeled causal variant;
it is not a numerical calibration test. MCMC initializes at null at
[lines 190–222 and 346–354](https://github.com/gkichaev/PAINTOR_V3.0/blob/ffce82a7114eeb98220db39ba5d6b69e5a419a9e/Functions_model.cpp#L190-L222) and includes null in its search space. Neither output supplies a
locus-null posterior. The wrapper records these distinctions without inventing
such a probability, and does not change either upstream algorithm.

## Pooling and resource implications

PAINTOR shares annotation coefficients across supplied loci, while keeping
separate float32 LD matrices in each `Locus`, as declared at
[Functions_model.hpp:127–140](https://github.com/gkichaev/PAINTOR_V3.0/blob/ffce82a7114eeb98220db39ba5d6b69e5a419a9e/Functions_model.hpp#L127-L140).
The R setup writes separate locus files. There is no dense covariance of all
concatenated variants. For 220 loci of 1,000 variants, persistent native LD is
approximately `4 * 220 * 1000^2 = 0.88 GB`; the R inputs require approximately
1.76 GB for that collection of LD matrices. Reader strings, per-locus
factorization workspaces, annotations and retained assessor matrices add
memory. These are allocation components, not a certified peak-RSS bound.

Runtime is more concerning. Enumeration recomputes a prior term over all `p`
variants for each configuration, as shown at
[Functions_model.cpp:83–95](https://github.com/gkichaev/PAINTOR_V3.0/blob/ffce82a7114eeb98220db39ba5d6b69e5a419a9e/Functions_model.cpp#L83-L95).
At `p=1000`, the two-effect support has 500,500 configurations per locus.
A 220-locus EM sweep therefore visits 110,110,000 configurations and evaluates
approximately 110,110,000,000 variant prior terms. The default maximum is
15 EM sweeps at [Functions_model.hpp:47](https://github.com/gkichaev/PAINTOR_V3.0/blob/ffce82a7114eeb98220db39ba5d6b69e5a419a9e/Functions_model.hpp#L47);
convergence can stop sooner. For the separate 128-block, 200-variant native
panel, one sweep visits 2,572,800 configurations. These counts identify work;
they do not predict elapsed hours or assert successful full-size execution.

The challenge census has 220 tasks, but current single-cohort regional dispatch
rejects four unsupported tracks before fitting; supported pooled input is at
most 180 loci, including slightly shorter missing-causal loci. Keep all 220
planned tasks in the completion denominator. The full native panel pools all
128 blocks for one observed phenotype and exposure. Slicing either supplied
pool into independent fits changes learned coefficients and the declared
transductive experiment.

## Required execution gates and possible improvements

A two-locus installation preflight cannot establish full-panel affordability.
Before launching a scientific campaign, run the frozen resource-probe dataset,
full method arguments and final binary in scheduled compute jobs. Retain
scheduler `resources_used.mem`, walltime, exit status, input dimensions, the
number of supported pooled loci, successful-fit census and receipt hashes.
For the configured 32 GB/72 h worker, a practical launch gate is successful
terminal completion with peak total job memory below 80% of allocated memory
and runtime below 80% of walltime. Those margins are operational policy, not
mathematical guarantees. A failure requires a new declared resource allocation,
a new prespecified algorithm/tuning arm, or an explicit failed attempt; never
silently replace enumeration with MCMC or shrink the training pool.

An architecture-preserving optimization could compute the all-null prior once
per locus/EM sweep, add selected-variant log odds for each nonempty
configuration, and reuse configuration likelihoods while enrichment changes.
This removes redundant work while preserving the current non-null support;
it must reproduce native PIPs, enrichment updates and convergence on small
independent fixtures before use. Introducing a null model would be a separate
scientific algorithm change. This audit implements metadata and tests only.
