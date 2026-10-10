# FineMap RL

Four verifiable environments train agents to make calibrated inferences, choose
measurements, improve inference programs and repair actual predictor code.
They turn the fine-mapping benchmark's scientific contract into observable
actions and assessor-side rewards. No model judge or external API key is needed.

| Environment | Agent action | Feedback | Terminal objective |
|---|---|---|---|
| `InferenceEnv` | Submit PIPs and a null probability | One terminal reward | Honest observed-variant sparse-membership and whole-locus-null probabilities |
| `ExperimentEnv` | Buy a noisy assay, or submit | Assay result and declared cost | Probability accuracy minus acquisition cost |
| `MethodDesignEnv` | Evaluate or submit a declarative inference program | Budgeted development scores | Separate held-out locus accuracy minus query costs |
| `CodeImprovementEnv` | Evaluate or submit Python predictor source | Budgeted development scores in Docker | Separate held-out accuracy; every failed locus is charged |

Install from the repository root:

```sh
python -m pip install './rl_envs[gym]'
python -m unittest discover -s rl_envs/tests -v
finemap-rl smoke --output results/rl-development --episodes 250 --test-episodes 100
```

For real code submissions, build the dependency-only runtime, then run its tests
and example. The runner resolves the image to its immutable ID before execution.

```sh
docker build -t finemap-rl-predictor:0.1 -f rl_envs/Dockerfile.predictor rl_envs
FINEMAP_RL_DOCKER=1 python -m unittest discover -s rl_envs/tests -v
finemap-rl smoke --output results/rl-docker-development --docker
```

The numerical smoke actually updates a contextual bandit from reward feedback.
It produces training/evaluation CSVs, bootstrap summaries, example trajectories
and source hashes. It establishes that training and scoring work; it does not
demonstrate LLM training, a commercial buyer or improved biological discoveries.

## Minimal interaction

```python
from finemap_rl import TaskFactory, InferenceEnv
from finemap_rl.baselines import predict

factory = TaskFactory()  # fresh private key; retain it only in the assessor
env = InferenceEnv(factory)
observation, info = env.reset(seed=1)
prediction = predict(observation["inputs"])
observation, reward, terminated, truncated, info = env.step(
    {"type": "submit", "prediction": prediction}
)
```

Use local Python/Gym objects only when the agent is trusted. Python underscores
do not prevent an untrusted agent from reading private state. For agent training
or hosted evaluation, run the assessor in a separate process/worker and expose
only the bounded JSON protocol:

```sh
finemap-rl serve --environment experiment --split train --max-episodes 1000
```

Send one JSON request per line, such as `{"op":"reset"}`, then
`{"op":"step","action":{"type":"assay","variant_id":"<supplied ID>"}}`.
`{"op":"state"}` returns public lifecycle metadata only. Requests have a 1 MiB
limit. Held-out `--split test` rejects client-selected seeds and caps issued
episodes; final rewards must not become an unlimited tuning service. A secure
remote gateway, authentication, rate limits and tenant separation remain an
operator responsibility.

## Inputs, actions and rewards

Each public task contains exactly `task_id`, `variant_id`, `z`, `LD`,
`annotations`, `n` and `statistic`. The convention is
`z | beta ~ Normal(sqrt(n) R beta, R)` with unit Gaussian residual scale.
Task IDs and variant order do not directly encode generating mechanisms; the assessor retains
sparse labels, whole-locus null status and mechanism. Inputs do not include the
private key, seed, causal count, effects or truth. Two annotations are supplied
without declaring whether they help. Eleven procedural tracks include matched
priors, nulls, multiple/large causal counts, opposing effects, background effects,
LD shift, exact LD equivalence, annotation enrichment/reversal and allele error.
These are synthetic stress tests, not empirical ancestry or assay validation.

A prediction contains exactly `task_id`, `variant_id`, `pip` and
`null_probability`. IDs must cover every supplied variant once. Predictions are
aligned by ID; non-finite, missing, duplicated and incoherent values fail. Each
PIP must not exceed `1 - null_probability`. There is no lower bound on the sum
of observed PIPs: a generating target can lie outside a frozen observed universe.
These event bounds do not establish a unique joint posterior or a credible set.

The terminal score is

```text
1 - [0.5 mean_j (PIP_j - gamma_sparse,j)^2
     + 0.5 (P_null - 1{no sparse target anywhere in locus})^2].
```

Both positive weighted terms are strictly proper for their declared marginal
events in expectation. Weighting the null separately prevents the sparse-label
majority from hiding locus false positives. This score does not evaluate joint
configuration calibration, and cannot justify cumulative-PIP credible sets.
Malformed actions terminate with reward -1; failed predictions in a panel incur
loss 1. The denominator includes every attempted task. There is no per-variant
reward or terminal label breakdown to use as an extra oracle.

`ExperimentEnv` provides a **stylised** noisy sparse-membership assay:
sensitivity .85, specificity .90, budget 6, at most two repetitions per variant,
and cost .003 utility units per measurement by default. The parameters are
public. Each acquisition receives only its noisy result and negative cost;
final scoring occurs once. Intermediate reward never depends directly on hidden
truth. Repetitions produce independent draws. This is a testbed for value of
information, not a claim that MPRA, CRISPR or replication studies directly report
GWAS causality with these operating characteristics.

Two acquisition baselines are supplied. `uncertainty_assay` buys variants with
the highest marginal posterior variance. `value_information_assay` evaluates
both possible assay outcomes under the working joint posterior, computes their
expected weighted Brier-risk reduction across all PIPs and the null event, and
stops when the best legal assay's gain fails to exceed its cost. It is a myopic
model-based policy; a misspecified posterior can underestimate useful information.

`MethodDesignEnv` admits `exact`, `independent_abf` and `prior` programs with
explicit effect SD, count prior, annotation weight and LD ridge. It gives four
development queries by default. The exact working reference enumerates zero,
one and two effects, and is limited to 64 variants. It is exact only for compatible
matched-prior inputs. Stress scenarios intentionally exceed or violate this
working model. No reference is a universal ceiling.

`CodeImprovementEnv` executes up to 64 KiB of source in a fresh non-root,
read-only Docker container without network, assessor mounts or secrets. It has
CPU, memory, PID, output and time limits. The predictor reads a JSON task list
from stdin and writes the prediction list to stdout; NumPy is available. The
starter is a forced-one approximation for the agent to improve. Assessor truth
and scoring execute outside the container. Docker is not a proven boundary
against kernel attacks: production needs dedicated workers and security review.

## Gymnasium

```python
import gymnasium as gym
from finemap_rl.gym_adapter import register_environments
register_environments()
env = gym.make("FineMappingExperiment-v0")
observation, info = env.reset(seed=42)
```

IDs: `FineMappingInference-v0`, `FineMappingExperiment-v0` and
`FineMappingMethodDesign-v0`. `info["public_observation"]` carries the structured
input contract and public assay/query state; it contains no truth. Inference and
experiment action dictionaries have `operation` (0 submit, 1 assay), `variant`
(index) and `probabilities` (PIPs followed by null probability). In numeric
observations, z comes first, then row-major LD, then row-major annotations.
Experiment vectors append attempts, successes and remaining budget. Method design
uses a documented candidate catalog: the first half of discrete actions evaluate,
the second half submit. Its vector appends four fields per query (present,
catalog index / catalog size, development loss, failure fraction) and remaining
queries. The structured observation makes these semantics inspectable.

The adapters follow the current [Gymnasium API](https://gymnasium.farama.org/api/env/),
including separate termination/truncation and seeded resets. Compatibility must
be checked with Gymnasium's checker, not inferred from method names.

The optional `openenv` extra provides typed OpenEnv 0.3.0 server interfaces and
an `openenv.yaml` manifest. Its persistent `/ws` endpoint is tested for complete
episodes and independent simultaneous sessions. Use WebSockets for this adapter:
the current framework's HTTP control endpoints create a fresh environment per
request. See the dossier for commands and the hosted-deployment limits.

## Frozen benchmark bridge

`scripts/rl/export_challenge.R panel.rds public_dir private_dir` exports trusted R
challenge panels. It explicitly reports unsupported multi-cohort/variant-wise-n
tasks; it does not delete them from the full R campaign. Missing-cause tasks retain
their whole-locus null label and observed sparse-membership labels. The bridge is
for single-cohort constant-n known-noise score data, not marginal OLS statistics.
Use `read_release(public_dir)` for inputs only, and
`read_release(public_dir, private_dir)` only in the assessor. JSON files have
SHA-256 manifests; private directories are separate, with mode 0700 when created.
Hashes detect corruption but do not confer blindness or custody.

See [submission dossier](docs/submission-dossier.md),
[threat model and red-team coverage](docs/threat-model.md) and
[validation record](docs/validation.md) for readiness, evidence and remaining gaps.
