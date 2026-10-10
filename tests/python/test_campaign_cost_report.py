"""Artificial sealed checkpoints test resource arithmetic, not scientific costs."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
RSCRIPT = os.environ.get("FMB_RSCRIPT") or (
    "/opt/homebrew/bin/Rscript" if Path("/opt/homebrew/bin/Rscript").exists() else shutil.which("Rscript"))


@unittest.skipUnless(RSCRIPT, "R is needed for the campaign resource-contract test")
class CostReportTests(unittest.TestCase):
    def test_legacy_and_challenge_rows_preserve_separate_resource_estimates(self):
        self.check_retention("loose")

    def test_sealed_task_bundles_include_all_headers_ledgers_and_receipts(self):
        self.check_retention("bundled")

    def check_retention(self, retention):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root / "source"
            shutil.copytree(REPO / "R", source / "R")
            (source / "scripts/campaign").mkdir(parents=True)
            for name in ("common.R", "cost_report.R"):
                shutil.copy2(REPO / "scripts/campaign" / name, source / "scripts/campaign" / name)
            setup = root / "fixture.R"
            setup.write_text(r'''
args <- commandArgs(TRUE); repo <- normalizePath(args[1L]); root <- normalizePath(args[2L])
storage <- args[3L]
for (p in list.files(file.path(repo, "R"), "[.]R$", full.names = TRUE)) source(p)
tasks <- data.frame(task_id = c("task_1", "task_2"), index = 1:2,
  dataset_ids = c("legacy", "challenge"), stage = c("corrected_refit", "challenge"), method = "abf")
datasets <- data.frame(dataset_id = c("legacy", "challenge"),
  family = c("legacy_reconstruction", "challenge"), row = c(NA_integer_, 1L))
for (profile in c("resource_probe", "full")) {
  directory <- file.path(root, profile); dir.create(directory)
  planned_tasks <- tasks; planned_datasets <- datasets
  if(profile == "full") {
    extra<-datasets[rep(1L,4L),];extra$dataset_id<-paste0("legacy",2:5)
    planned_datasets<-rbind(datasets,extra)
    planned_tasks$dataset_ids[1L]<-paste(c("legacy",extra$dataset_id),collapse="|")
  }
  plan <- list(profile = profile, original_source = "", corrected_source = .campaign_snapshot(repo),
    method_args = list(abf = list()), runtime_fingerprint = .campaign_runtime_fingerprint(list(abf = list())),
    tasks = planned_tasks, datasets = planned_datasets, generation = data.frame(generation_id = 1L, replicates = 10L), seed = 1L)
  plan$plan_sha256 <- .campaign_object_sha(plan)
  .campaign_atomic_rds(plan, file.path(directory, "plan.rds"))
  if (profile == "full") next
  for (i in 1:2) {
    task <- tasks[i, ]; d <- datasets[i, ]; n <- if (i == 1L) 10L else 220L
    output_dir <- file.path(directory, "tasks", task$task_id)
    artifact <- file.path(output_dir, paste0(d$dataset_id, ".rds"))
    # No inference took place: declared durations exercise only the arithmetic.
    output <- list(terminal = TRUE, plan_sha256 = plan$plan_sha256, task = list(task = task, dataset = d),
      method_task = task, metadata = as.list(d), result = list(abf = list(total_runtime_seconds = 10)),
      elapsed_seconds = if (i == 1L) 10 else 66,
      table = data.frame(actual_p = if (i == 1L) rep(c(500L, 750L, 1000L, 1500L, 2000L), each = 2L) else 1000L,
        n_gwas = 5000L, n = if (i == 1L) 5000L else 10000L, success = rep(TRUE, n)))
    .campaign_atomic_rds(output, artifact)
    receipt <- paste0(d$dataset_id, ".terminal.json")
    .campaign_receipt(output_dir, artifact, plan$plan_sha256,
      paste(task$task_id, d$dataset_id, sep = "."), receipt)
    if (storage == "bundled") {
      .campaign_start_attempt(output_dir,plan,task,d)
      .campaign_compact_task(output_dir,plan,task)
    } else {
      completed <- file.path(output_dir, "completed.rds")
      .campaign_atomic_rds(list(terminal = TRUE, plan_sha256 = plan$plan_sha256, task = task,
        artifacts = .campaign_artifact_manifest(c(artifact, file.path(output_dir, receipt)))), completed)
      .campaign_receipt(output_dir, completed, plan$plan_sha256, task$task_id)
    }
    receipt_path<-file.path(output_dir,"terminal.json"); value<-jsonlite::fromJSON(receipt_path)
    value$resource_measurement<-list(schema=1L,n_new=1L,n_datasets=1L,
      script_elapsed_seconds=output$elapsed_seconds+i+1L,data_processing_seconds=output$elapsed_seconds,overhead_seconds=i+1L)
    jsonlite::write_json(value,receipt_path,auto_unbox=TRUE,digits=NA)
  }
}
''')
            env = os.environ.copy()
            if (REPO / ".local-lib").exists():
                env["R_LIBS"] = str(REPO / ".local-lib")
                env["FMB_R_LIB"] = str(REPO / ".local-lib")
            built = subprocess.run([RSCRIPT, "--vanilla", str(setup), str(source), str(root), retention],
                                   env=env, text=True, capture_output=True)
            self.assertEqual(built.returncode, 0, built.stderr)
            output = root / "cost.csv"
            command = [RSCRIPT, "--vanilla", str(source / "scripts/campaign/cost_report.R"),
                       "--canary-plan", str(root / "resource_probe/plan.rds"),
                       "--full-plan", str(root / "full/plan.rds"), "--out", str(output)]
            result = subprocess.run(command, env=env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open() as handle:
                rows = {row["stage"]: row for row in csv.DictReader(handle)}
            self.assertEqual(set(rows), {"corrected_refit", "challenge"})
            self.assertEqual(int(rows["corrected_refit"]["planned_locus_fits"]), 50)
            self.assertEqual(int(rows["challenge"]["planned_locus_fits"]), 220)
            self.assertAlmostEqual(float(rows["corrected_refit"]["worker_hours_median"]), 50 / 3600)
            self.assertAlmostEqual(float(rows["challenge"]["worker_hours_median"]), 66 / 3600)
            self.assertAlmostEqual(float(rows["corrected_refit"]["end_to_end_worker_hours_median"]), 52 / 3600)
            self.assertAlmostEqual(float(rows["challenge"]["end_to_end_worker_hours_median"]), 69 / 3600)
            self.assertTrue(all(float(row["checkpoint_storage_gb_median"]) > 0 for row in rows.values()))
            with (root / "cost-measured.csv").open() as handle:
                measured = list(csv.DictReader(handle))
            self.assertTrue(all(float(row["primary_retention_bytes"]) > float(row["checkpoint_bytes"])
                                for row in measured))
            if retention == "bundled":
                self.assertTrue(all(row["retention"].startswith("sealed_task_bundle") for row in measured))
                for i in (1, 2):
                    self.assertEqual(len([x for x in (root / f"resource_probe/tasks/task_{i}").rglob("*")
                                          if x.is_file()]), 2)
            with (root / "cost-file-counts.csv").open() as handle:
                census = {row["retention"]: row for row in csv.DictReader(handle)}
            self.assertEqual(int(census["sealed_task_bundles"]["primary_terminal_regular_files"]), 4)
            self.assertEqual(int(census["loose_checkpoints_with_one_attempt"]["primary_terminal_regular_files"]), 22)
            checkpoint = root / "resource_probe/tasks/task_2" / (
                "completed.rds" if retention == "bundled" else "challenge.rds")
            checkpoint.write_bytes(checkpoint.read_bytes() + b"changed valid-prefix RDS")
            rejected = subprocess.run(command, env=env, text=True, capture_output=True)
            self.assertNotEqual(rejected.returncode, 0)
            self.assertIn("bytes changed", rejected.stderr)


if __name__ == "__main__":
    unittest.main()
