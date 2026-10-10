# Executed RL validation

Validated locally on 10 October 2026. The package is a separate Python artifact;
the larger R/server fine-mapping campaign has its own validation record.

## Tests and interfaces

`FINEMAP_RL_DOCKER=1 python -m unittest discover -s rl_envs/tests -v` completed
**30 tests with no failures, errors or skips**. The frozen
[test log](artifacts/tests.log) records every executed check. This includes:

- Private-seed/key separation, disjoint split identities, observation allowlists,
  immutable sparse labels, release integrity and public/private path separation.
- Proper-score algebra, attempted-failure denominators, missing-target event
  semantics, posterior configuration accounting, exact-LD exchangeability,
  variant permutation and consistent allele-sign invariance.
- Assay likelihood updates, budget/repetition limits, deterministic episodes,
  one terminal submission, no terminal truth breakdown, and a value-of-information
  stopping check when the declared cost exceeds available model-based gain.
- Actual contextual-bandit parameter updates from reward; bounded validation
  queries; final client-seed rejection; unknown/non-finite/oversized JSON requests.
- Rejection of final resets before submitting the active task, preventing a
  client from skipping difficult cases through the reset operation.
- Official Gymnasium **1.3.0** `check_env` for all three numeric adapters.
  Its warnings concern intentionally unbounded Gaussian Box observation spaces.
- Actual OpenEnv **0.3.0** persistent WebSocket reset/state/submit lifecycle,
  rejected repeated scoring/final seeds, and simultaneous independent-key sessions.
- Actual Docker predictor execution and held-out scoring, non-root/read-only
  source, absent assessor files/socket/environment secrets, disabled networking,
  timeout, output flood, malformed JSON and failed-task accounting.

The first Docker flood probe found a cleanup timeout. Continued output draining
with bounded buffers and disabled container logging corrected that path; the
complete repaired suite above passed. Current OpenEnv HTTP controls are
stateless per request; the implemented hosted episode route uses `/ws`.

Local interpreter: Python **3.14.7**, NumPy **2.5.3**. The predictor runtime uses
Python **3.11.14** and NumPy **2.2.6**, in Linux ARM containers on Docker **29.2.1**.
The resolved tested predictor image is
`sha256:ec80cbe8ad8612c58e4db36cf7f7af982a059f4c95e97310494fca2df65d743f`.
Its Python base was
`python:3.11.14-slim-bookworm@sha256:65a93d69fa75478d554f4ad27c85c1e69fa184956261b4301ebaf6dbb0a3543d`.
Image identity establishes this local run's provenance; rebuilding an image can
produce a different digest.

## Reward-based training and example trajectories

The [frozen provenance](artifacts/provenance.json) and
[training log](artifacts/training.csv) record **250 numerical contextual-bandit
updates**, followed by **100 distinct validation loci per policy, six policies
(600 episodes)**, on 12-variant tasks from eleven stress/control mechanisms.
The training key is explicitly public and the run is marked `blind: false`.
The learned policy uses input features and training reward only. It updates
regression estimates for eight prespecified inference configurations; it does
not train an LLM or synthesize a new inference algorithm.

| Policy | Mean utility on 100 validation loci |
|---|---:|
| Learned contextual bandit | 0.928781 |
| Active marginal uncertainty | 0.926501 |
| Active model-based myopic value of information | 0.907386 |
| Working exact 0/1/2-effect reference | 0.898580 |
| Independent ABF | 0.884009 |
| Prior only | 0.819194 |

Utility includes each policy's declared acquisition costs. These means compare
a public synthetic development mixture, not a blinded biological competition.
The [paired comparisons](artifacts/paired-comparisons.json) resample the same loci.
The learned-policy difference from the working exact reference was **+0.030201**,
with a 95% percentile paired bootstrap interval **[+0.015384, +0.048274]**.
The uncertainty-acquisition difference was **+0.027921**, interval
**[-0.006574, +0.075107]**. The model-based myopic difference was **+0.008806**,
interval **[-0.004904, +0.023194]**. These intervals describe this simulation
experiment conditional on its specification and trained policy; they do not
account for selection across alternate mixtures, protocols or future experiments.

The two acquisition rules illustrate why a working likelihood matters. A policy
that stops using its posterior's predicted risk reduction can miss useful
measurements when that posterior is overconfident under a stress mechanism.
The aggregate result alone does not identify which mechanism caused the gap;
mechanism-stratified and independently replicated evaluation is needed before
making that attribution.

The actual Docker code-improvement example evaluated a forced-one starter and
an independently written ABF predictor on six development loci. Development
loss fell from **0.149947** to **0.112943**; the second predictor was selected
using those development scores only. A distinct six-locus final panel returned
**0.881110 utility after two query costs**, with **zero failed tasks**. This
demonstrates runnable source improvement and held-out scoring, not trained
coding-agent performance. All held-out labels remained outside the predictor
container. [Public example trajectories](artifacts/example-trajectories.json)
contain observations/actions/rewards rather than hidden labels.

## R challenge import

The export bridge was executed with a trusted eight-task R panel containing
matched-prior, missing-causal, known-overlap and shared-cohort tracks, two loci
each. It exported **four supported tasks** and explicitly recorded the other
four as unsupported multi-cohort/variant-wise-n inputs. Python verified both
SHA-256 manifests and loaded the aligned assessor labels successfully. The
[bridge report](artifacts/bridge-report.json) preserves this support boundary.
Missing generating targets are retained in whole-locus null semantics; observed
membership labels do not impose a false lower union bound on observed PIP mass.

## Packaging and remaining scope

Wheel and source-distribution builds are checked separately from scientific
correctness. Build hashes, a clean Linux wheel installation and its dependency
versions are recorded in `artifacts/build-provenance.json` after the final build.
These files are software provenance, not evidence of a commercial transaction.
The clean Linux wheel test runs the core suite on Python 3.11.14/NumPy 2.2.6;
Docker/Gymnasium/OpenEnv integration cases are explicitly skipped in that minimal
container. All integrations were executed in the full local suite above. The
build-provenance file is excluded from the source distribution to avoid requiring
an archive to contain its own hash.

No marketplace submission, receiving-lab acceptance, LLM fine-tuning, remote
multi-tenant hosting, biological assay validation or certified sandbox audit was
performed. A fixed-key hosted final service needs persistent coordinator replay
protection across restarts and independent custody/rate limits. Optional native
fine-mapping binaries are not bundled with this package. The current real-code
runner is opt-in and requires Docker; Gymnasium/OpenEnv are explicit extras.
