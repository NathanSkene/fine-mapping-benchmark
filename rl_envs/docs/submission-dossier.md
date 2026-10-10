# Submission and commercial evaluation dossier

Prepared 10 October 2026. This is a concrete software offer and evaluation plan,
not a claim that a frontier lab has accepted, bought or requested it. No outreach,
account creation, upload or publication to an environment marketplace has been
performed. Repository delivery and an externally hosted service are separate steps.

## The product

FineMap RL offers four connected tasks: numerical probability forecasting, noisy
measurement allocation, budgeted method design, and actual Python predictor repair.
The scientific theme is transferable: distinguish evidence from latent truth,
account for unavailable information, spend a finite evaluation budget, and produce
calibrated forecasts when a model's assumptions fail. Fine mapping supplies a
domain with a tractable generative model, exact small-model checks, correlated
candidates, meaningful uncertainty and realistic pipeline failure modes.

The deliverable includes an installable MIT Python package, deterministic
assessor-side verifiers, a finite candidate language, a Docker code runner,
Gymnasium adapters, a typed OpenEnv server adapter, JSON-only frozen releases,
R benchmark import hooks, baselines, actual reward-based numerical training,
sample trajectories, adversarial tests and provenance. Raw target labels and
secret final keys are not package contents. The compact development example uses
a public key and is explicitly unsuitable as a secret assessment panel.

| Buyer use case | Supplied environment | What evidence it measures |
|---|---|---|
| Calibration after post-training | Inference | Proper forecasts of explicit sparse and null events |
| Scientific tool choice | Experiment | Value of information with noisy assays and finite costs |
| Research optimisation | Method design | Improvements selected using bounded development feedback, then evaluated on distinct loci |
| Scientific coding agents | Code improvement | Executable predictor changes verified against held-out outcomes outside the agent container |

Claims of transfer to biological discovery, robust LLM training or a prospective
community assessment need further evidence. The assay is stylised; it is not an
MPRA, CRISPR or replication simulator validated against experimental error rates.
The current exact baseline has a two-effect cap and cannot represent the entire
fine-mapping research frontier. The code runner is tested locally on Linux ARM
containers; production isolation is not certified by those tests.

## Integration routes checked against current primary documentation

**Gymnasium 1.3.0:** numeric observation/action spaces, seeded reset, separate
termination/truncation and official environment checks are implemented. This
supports conventional policy optimisation without an LLM service. The contract
follows [Farama's environment API](https://gymnasium.farama.org/api/env/).

**OpenEnv 0.3.0:** `FineMappingAction`, `FineMappingObservation` and
`FineMappingState` implement the installed framework's typed server interfaces.
The application factory supports inference, experiment and method design, with
independent private keys per persistent session. Persistent `/ws` interactions
are tested. OpenEnv's current HTTP control endpoints instantiate and close a
fresh environment per request, so they are unsuitable for a multi-step episode
here; use `/ws` rather than chaining HTTP `/reset` and `/step`. These choices
follow the framework's [client-server concepts](https://github.com/huggingface/OpenEnv/blob/main/docs/source/guides/concepts.md)
and [official server implementation](https://github.com/huggingface/OpenEnv/blob/main/src/openenv/core/env_server/http_server.py).
The optional package extra and `openenv.yaml` declare the integration. An
OpenEnv marketplace upload, remote Docker deployment and its platform-specific
validator remain unexecuted.

**Prime Intellect / Verifiers:** current official documentation provides
taskset/harness composition and a community environment hub. The observed
step/reset protocol and quantitative verifier could be wrapped as a hosted
tool or framework harness. We do not ship a Verifiers adapter or claim its
conformance: its current API must be pinned and tested with the receiving
operator. The [official environment guide](https://docs.primeintellect.ai/verifiers/environments)
is an integration route, not an acceptance policy or purchase commitment.

Run the OpenEnv training application locally after installing its extra:

```sh
python -m pip install './rl_envs[openenv]'
python -m uvicorn finemap_rl.openenv_adapter:app --host 127.0.0.1 --port 8000
```

WebSocket requests use `{"type":"reset","data":{}}`,
`{"type":"step","data":{"command":<environment action>}}` and
`{"type":"state"}`. `/health` and `/schema` expose framework diagnostics.
The JSON-lines assessor is an alternative with fewer dependencies and an explicit
request limit. Neither interface supplies authentication or a tenant firewall.

## Licensing and provenance

The package retains the repository's Lucas Siu copyright and adds contributor
attribution. Its MIT terms permit commercial distribution subject to retaining
the notice. This is a description of the supplied license text, not clearance
of institutional IP ownership or third-party data rights. A proposed supplier
must establish ownership/authority to sign a commercial contract independently.
No participant genotypes, 1000 Genomes files, FINEMAP/GCTB binaries, licensed
reference panels or restricted annotations are bundled in this package.

Procedural synthetic inputs are generated from repository-authored NumPy code.
The Docker runtime contains pinned Python/NumPy dependencies; their own license
notices and the base image's distribution licenses remain applicable. Record an
image digest and dependency manifest for each deployment. The local tested image
and package/build hashes appear in the [validation record](validation.md).

## Assessment custody and reward governance

Separate task generation, predictor execution and truth scoring operationally.
For development, public fixed keys enable reproducible debugging. For a closed
assessment, the operator supplies fresh unreleased keys or independently held
frozen labels and does not mount them in an agent container. Public/private
directory separation and hashes alone do not establish security or blindness.

The terminal scalar reward conveys information about labels. A final episode is
consumed after submission; repeated stepping is rejected. Final client-selected
seeds are rejected, issued episodes are capped, and OpenEnv sessions receive
fresh keys. Active final tasks cannot be skipped by requesting another reset.
A service must count all issued tasks, including disconnected/abandoned episodes,
and give an unsubmitted task its declared failure score. A service with
persistent fixed keys must also keep a durable
coordinator ledger across process restarts, maintain participant access rules,
rate-limit new sessions and forbid using closed assessment feedback to tune.
Those service-level controls are not replaced by the local in-memory cap.

Acquisition costs and validation-query penalties are explicit task utilities,
not empirical prices or hardware runtime charges. Publish all changes to their
values as a specification version. Keep task mixture weights visible: improving
utility in one synthetic mixture does not prove universal superiority.

## Concrete pilot acceptance plan

1. A receiving lab installs the wheel in a clean environment and runs the core,
   Gymnasium and hosted API checks; a dedicated worker runs Docker adversarial
   tests. Record its platform, image digest and limits.
2. Freeze a training distribution and a withheld evaluation panel under separate
   custody. Include distinct loci, generator-family shifts and unreleased
   experimental evidence when the intended claim concerns biology.
3. Evaluate the lab's existing agent before training. Compare fixed baselines,
   its trained policy and a blinded ablation on the same loci. Report all attempts,
   costs, failures, calibration and paired uncertainty, not only a mean reward.
4. Red-team label recovery, reward-code tampering, output/compute abuse, final
   replay and query-budget bypass on the receiving infrastructure.
5. Agree a measurable pilot deliverable, deployment support, ownership, hosting
   responsibilities, update cadence and acceptance criteria before pricing.

There is no invented revenue forecast or list of alleged buyers. Commercial
value would come from reliable task renewal, withheld evidence, integration and
support, rather than exclusivity over an already open MIT codebase.
