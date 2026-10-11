# Assessor boundaries and red-team coverage

The trusted components are the task factory, private labels, verifier, episode
coordinator and Docker host. The agent controls structured actions, probability
vectors and, in code tasks, up to 64 KiB of Python source. Treat submitted code
as untrusted. Direct Python/Gym access is only a development interface for a
trusted training process; Python attribute names do not enforce confidentiality.

| Failure or attack | Implemented control | Verification / remaining limit |
|---|---|---|
| Seed or task ID reveals labels | HMAC private-key domains, opaque IDs, random variant order; no latent fields in input schema | Repeatability, different-key and disjoint-split tests; private keys still require operator custody |
| Private truth included in observations/state | Exact public allowlist; copied inputs; public lifecycle state; no truth-derived terminal breakdown | Core/JSON/OpenEnv tests search response fields; legitimate inputs can reveal a mechanism statistically |
| Nulls disappear from reward | Separate positively weighted whole-locus-null Brier term | Analytic proper-score and null event tests |
| Duplicate/missing/wrong task identities gain credit | Complete unique IDs, task identity checks and ID alignment | Malformed and cross-task submission tests |
| Failed hard loci improve average score | Loss 1 for each attempted failure; no filtered denominator | Panel and actual failed-code tests |
| Repeated reward queries recover final labels | One terminal submission; later steps rejected; no test seed selection; issued episode cap; reset cannot skip an active final task | JSON and OpenEnv tests; fixed-key restarts need a durable service ledger |
| Assays become free truth | Noisy public channel, finite budget and repetition bounds, fixed negative acquisition cost | Determinism, legal-action and exhausted-budget tests; operating rates are stylised |
| Unlimited development feedback | Fixed query budget, explicit costs, disjoint final loci | Design tests; development feedback intentionally supports learning and can reveal aggregate label information |
| Agent modifies reward or reads key | Predictor source only mounted; public JSON input; verifier/private tasks outside container | Actual Docker hidden-file/source-write probes; no certified host/kernel boundary |
| Network exfiltration | `--network none`, no Docker socket, non-root/read-only runtime, dropped capabilities | Actual network and socket/path probes |
| Compute/memory/process abuse | CPU/memory/swap/PID limits; bounded time; ephemeral container | Actual infinite-loop timeout; no exhaustive cgroup/kernel stress or escape assessment |
| Stdout/stderr flood | Bounded buffers, continued draining during kill, logging disabled | Actual flood probe and cleanup check; Docker daemon health remains an operator concern |
| JSON or arbitrary protocol introspection | Bounded JSON-lines input; non-finite rejection; limited reset/step/state operations; no `eval`/`exec` | Non-finite/oversized/unknown-operation tests; OpenEnv needs gateway request/rate limits |
| Split or prior shift conceals fragility | Explicit generator version and stress mechanisms; no universal exact-ceiling claim | Sign/permutation/perfect-LD/prior/annotation checks; realistic genotypes and biological truth require additional panels |

Intermediate assay and query rewards are declared costs, independent of hidden
final correctness. The verifier uses raw Brier losses without rounding or bins.
Necessary coherence requires each observed PIP to be at most the non-null
probability. We do not enforce `sum observed PIPs >= P(non-null)`: a missing
generating variant makes that inequality invalid. An explicit complete-universe
contract would support additional checks in a future specification.

The first Docker output-flood test caused a cleanup timeout. The runner was
changed to keep draining excess output without retaining it, disable container
logging and report cleanup failure explicitly. The repaired full Docker suite
was rerun successfully. This is evidence of a corrected tested path, not a claim
that arbitrary adversarial code cannot affect an infrastructure worker.

No model-judged scientific narrative enters the reward. Statements of biological
causality, population realism or commercial quality need evidence outside these
numeric verifiers. Passing schemas and tests does not make those claims true.
