# Reproducible development tuning

Functional BEATRICE's tuning driver is a development tool. A search optimum on
its fixed panel is an in-sample selection result. It is not a comparative
fine-mapping benchmark result or evidence of calibration on future loci. Freeze
candidate hyperparameters and run every candidate on independently generated
validation tasks before comparing methods. A competitive evaluation also needs
comparable declared tuning budgets for the competing methods.

The revised contract (`attempted-locus-macro-ap-v2`) makes these choices explicit:

* Every cached scenario is evaluated. `--scenarios` must equal the scenario
  count in the RDS, so extra hard cells cannot be silently excluded.
* The same first `--regions` loci are used in every scenario (`0` means all).
  Requested and realised locus sizes are subset with the genotype and truth
  lists. These loci are fixed experimental blocks, not random samples from a
  population of genomic loci.
* `--rotate-regions` is rejected. Coprime cycle lengths alone do not make a
  one-locus-per-scenario schedule balanced within each architecture cell.
  In the former default design, locus size was paired with particular `S` and
  `phi` cells. To reduce cost, generate smaller loci or use a common fixed
  subset across the complete architecture grid.
* `ap` is the evaluator's exact tie-aware **all-attempt locus-macro AP** on
  non-null loci. It cannot fall back to pooled variant AP or success-only AP.
  Null scenarios have undefined AP and contribute their available calibration
  diagnostics. An entirely null panel has no ranking objective.
* Python reads and validates `n_failed`. A partial-success configuration becomes
  an Optuna `FAIL`, instead of competing using the favourable surviving fits.
  A nonzero subprocess exit also overrides any printed `status=ok`.
* The mass objective is the mean of each scenario's absolute count-ratio error.
  The previous absolute error of the mean ratio let over- and underprediction
  cancel. The ratio is undefined when the scenario contains no causal variants.
* Pareto trials run the complete panel. AP-only quantile aborts are rejected,
  including in scalar mode; they can remove candidates that improve another
  objective. Scalar mode retains Optuna's pruner, which receives the same
  scalar formula used at completion. This remains a heuristic based on a
  scenario prefix and can miss candidates whose benefit appears in later cells.

The count ratio is an aggregate diagnostic, and high-PIP ECE is a conditional
binned diagnostic. Neither is a strictly proper score for all submitted PIPs.
FDR estimates and high-PIP reliability are undefined when there are no selected
variants. Existing tuning objectives retain their documented missing-value
penalties; they do not turn abstention into evidence of failed probability
calibration. The challenge suite's raw Brier and log scores, native credible
coverage, null tasks and availability are better tools for final assessment.

## Study identity and random streams

Before optimisation, the driver records SHA-256 hashes of the simulation RDS,
R evaluator/wrapper sources, tuning scripts and BEATRICE Python sources, together
with the objective, search settings, Optuna version and base seed. R loads the
fingerprinted checkout through `pkgload`, so an older installed package cannot
silently replace the evaluator. A journal sidecar is created atomically; workers
with different contracts cannot initialise the same new study concurrently.
Changed data, source, scoring, search settings or seed require a new study.
Historical studies without fingerprints remain readable with `--report` but
cannot be resumed under the revised objective. Keep each `.contract.json` sidecar
with its journal when moving or archiving a study.

Use `--seed` for a declared sampler base seed and `--worker-id` for distinct
worker streams. The default worker ID is `PBS_ARRAY_INDEX` or zero, and the
worker seed is `(base_seed + 1000003 * worker_id) mod 2^32`. This makes initial
single-worker runs with the same environment and schedule repeatable. It does
**not** promise an identical resumed or asynchronous search trajectory: journal
storage preserves trials, not the sampler's complete random state; scheduling,
completion order and worker count affect proposals. Archive environment/package
versions and worker logs as well as the study contract for reproducibility.

## Verification

`tests/testthat/test-tuning-contract.R` checks the common locus subset,
authoritative metadata, invalid rotation and the AP estimand.
`tests/python/test_tuning_contract.py` checks partial-failure rejection,
subprocess status, grid completeness, immutable study identity, data/source hash
changes, and opposing count errors. These tests run without a neural campaign.
The CLI preflight and current Optuna journal backend were also exercised against
a two-scenario synthetic metadata fixture. No full hyperparameter campaign was
rerun, and historical winning configurations have not been revalidated.
