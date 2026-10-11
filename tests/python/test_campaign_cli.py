"""Real small R fits exercise bundled CLI integration, never scientific ranks."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
RSCRIPT = os.environ.get("FMB_RSCRIPT") or shutil.which("Rscript")


@unittest.skipUnless(RSCRIPT, "Rscript required for actual campaign CLI integration")
class CampaignCliTests(unittest.TestCase):
    def test_actual_fits_paired_reuse_bundle_resume_collection_and_tamper(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = os.environ.get("FMB_CAMPAIGN_TEST_KEEP_DIR")
            root = Path(destination) if destination else Path(temporary)
            root.mkdir(parents=True, exist_ok=True)
            source = root / "source"
            reuse_fixture = os.environ.get("FMB_CAMPAIGN_TEST_REUSE_FIXTURE") == "1"
            if reuse_fixture:
                self.assertTrue(destination and (root / "campaign/plan.rds").exists())
                # A collector-test retry can reuse already executed real fits
                # only when all method R bytes and scientific worker scripts
                # match. Controllers may be newer than the frozen workers.
                for path in (REPO / "R").glob("*.R"):
                    self.assertEqual(path.read_bytes(), (source / "R" / path.name).read_bytes())
                for name in ("common.R", "generate.R", "run_task.R", "collect.R", "collect_streamed.R"):
                    self.assertEqual((REPO / "scripts/campaign" / name).read_bytes(),
                                     (source / "scripts/campaign" / name).read_bytes())
            else:
                shutil.copytree(REPO / "R", source / "R")
                shutil.copytree(REPO / "scripts/campaign", source / "scripts/campaign")
            original = root / "same-source-control"
            if not reuse_fixture:
                shutil.copytree(source / "R", original / "R")
            env = os.environ.copy()
            if (REPO / ".local-lib").exists():
                env["R_LIBS"] = str(REPO / ".local-lib")
                env["FMB_R_LIB"] = str(REPO / ".local-lib")
            fixture = root / "prepare.R"
            fixture.write_text(r'''
a <- commandArgs(TRUE); repo <- normalizePath(a[1L]); root <- normalizePath(a[2L])
if(nzchar(Sys.getenv("FMB_R_LIB"))) .libPaths(c(Sys.getenv("FMB_R_LIB"),.libPaths()))
for(x in list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE)) source(x)
campaign <- file.path(root,"campaign"); dir.create(campaign)
set.seed(121); geno <- lapply(1:2,function(i) {
 X <- scale(matrix(rnorm(80*12),80,12))
 list(X=X,X_raw=X,n=80L,p=12L,LD=cor(X),maf=rep(.3,12),variant_ids=paste0("r",i,":",1:12),
      annotations_matrix=cbind(category=rep(0:1,6),decoy=rnorm(12))) })
scenarios <- lapply(1:5,function(i) list(scenario_id=i,S=1L,phi=.2,iter=i,model="sparse",
 regions=simulate_phenotypes(geno,S=1L,phi=.2,annotations="binary",n_annotations=2L,
                            enrichment=c(2,1),seed=500+i,verbose=FALSE)))
sim <- list(genotypes=geno,scenarios=scenarios,params=list(annotation_type="binary",model="sparse",
              n_annotations=2L,enrichment=c(2,1),annotation_correlation=0))
generation <- data.frame(generation_id=1:2,family=c("legacy_reconstruction","challenge"),
 row=1L,p=c(NA_integer_,12L),n=c(80L,100L),replicates=c(NA_integer_,1L),seed=1001:1002)
sim_path <- file.path(campaign,"datasets/generation_00001/simulation.rds")
ids <- paste0("controlled_sc",1:5)
datasets <- data.frame(dataset_id=c(ids,"challenge"),family=c(rep("legacy_reconstruction",5),"challenge"),
 row=c(rep(1L,5),2L),scenario=c(1:5,NA_integer_),
 sim_file=c(rep(sim_path,5),file.path(campaign,"datasets/generation_00002/panel.rds")),
 sim_sha256="generated",result_files="",result_sha256="",genotypes_sha256="generated",
 original_label="same-source synthetic-genotype integration control")
tasks <- list(); k <- 0L
for(stage in c("original_refit","corrected_refit","corrected_resimulation"))
 for(method in c("abf","marginal_z")) {
  k<-k+1L; tasks[[k]] <- data.frame(task_id=sprintf("task_%07d",k),index=k,
   dataset_ids=paste(ids,collapse="|"),stage=stage,method=method,seed=100L+k) }
for(method in c("abf","marginal_z","exact","forced_one")) {
 k<-k+1L; tasks[[k]] <- data.frame(task_id=sprintf("task_%07d",k),index=k,
   dataset_ids="challenge",stage="challenge",method=method,seed=100L+k) }
tasks <- do.call(rbind,tasks)
method_args <- list(abf=list(prior_variance=.04),marginal_z=list(coverage=.95))
plan <- list(profile="integration",specification="same-source small real-fit CLI control; no historical or full-rank claims",
 repo=repo,original_source=normalizePath(file.path(root,"same-source-control")),
 corrected_source=.campaign_snapshot(repo),
 original_source_hashes=.campaign_snapshot(normalizePath(file.path(root,"same-source-control"))),
 methods=c("abf","marginal_z"),method_args=method_args,
 runtime_fingerprint=.campaign_runtime_fingerprint(method_args),seed=143L,
 generation=generation,datasets=datasets,tasks=tasks)
plan$plan_sha256 <- .campaign_object_sha(plan)
.campaign_atomic_rds(plan,file.path(campaign,"plan.rds"))
writeLines(plan$plan_sha256,file.path(campaign,"plan.sha256"))
for(nm in c("tasks","datasets","generation")) write.csv(plan[[nm]],file.path(campaign,paste0(nm,".csv")),row.names=FALSE)
files<-c("tasks.csv","datasets.csv","generation.csv")
jsonlite::write_json(list(schema=1L,plan_sha256=plan$plan_sha256,
 plan_rds_sha256=.challenge_sha256(file.path(campaign,"plan.rds")),
 files=stats::setNames(as.list(vapply(file.path(campaign,files),.challenge_sha256,"")),files)),
 file.path(campaign,"dispatch-manifest.json"),auto_unbox=TRUE)
.campaign_atomic_rds(sim,sim_path)
checkpoint<-file.path(dirname(sim_path),"generated.rds")
.campaign_atomic_rds(list(terminal=TRUE,plan_sha256=plan$plan_sha256,task=generation[1,],
 genotypes_sha256=.campaign_object_sha(geno),artifacts=.campaign_artifact_manifest(sim_path)),checkpoint)
.campaign_receipt(dirname(sim_path),checkpoint,plan$plan_sha256,"generation_1")
''')
            def run(arguments):
                result = subprocess.run([RSCRIPT, "--vanilla", *arguments], env=env,
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
                return result
            if not reuse_fixture:
                run([str(fixture), str(source), str(root)])
            plan = root / "campaign/plan.rds"
            script = source / "scripts/campaign"
            if not reuse_fixture:
                for index in (1, 2): run([str(script / "preflight.R"), "--plan", str(plan), "--index", str(index)])
                run([str(script / "generate.R"), "--plan", str(plan), "--index", "2"])
                for index in range(1, 11):
                    run([str(script / "run_task.R"), "--plan", str(plan), "--index", str(index)])
            tasks = root / "campaign/tasks"
            self.assertEqual(len(list(tasks.rglob("*.rds"))), 10)
            self.assertEqual(len(list(tasks.rglob("*.json"))), 10)
            self.assertFalse((root / "campaign/fit-cache").exists())
            for receipt in tasks.rglob("terminal.json"):
                timing = json.loads(receipt.read_text())["resource_measurement"]
                self.assertEqual(timing["n_new"], timing["n_datasets"])
                self.assertGreater(timing["script_elapsed_seconds"], 0)
                self.assertAlmostEqual(timing["script_elapsed_seconds"],
                    timing["data_processing_seconds"] + timing["overhead_seconds"], places=5)
                self.assertIn("excludes R process startup", timing["scope"])
            before = {str(p): p.read_bytes() for p in tasks.rglob("*") if p.is_file()}
            result = run([str(script / "run_task.R"), "--plan", str(plan), "--index", "3"])
            self.assertIn("already complete", result.stdout)
            for path, content in before.items(): self.assertEqual(Path(path).read_bytes(), content)
            run([str(script / "collect.R"), "--plan", str(plan), "--bootstrap-draws", "20",
                 "--backend", "memory"])
            analysis = root / "campaign/analysis"
            reference = root / f"memory-reference-{os.getpid()}"
            shutil.copytree(analysis, reference)
            run([str(script / "collect.R"), "--plan", str(plan), "--bootstrap-draws", "20"])
            comparison = root / "compare.R"
            comparison.write_text(r'''
a <- commandArgs(TRUE)
reference <- a[1L]; streamed <- a[2L]
source(file.path(a[3L],"collect_streamed.R"))
schema <- readRDS(file.path(streamed,"task-metrics-schema.rds"))
classes <- unname(ifelse(schema$storage_classes=="double","numeric",schema$storage_classes))
x <- readRDS(file.path(reference,"task-metrics.rds"))
y <- read.csv(file.path(streamed,"task-metrics.csv"),stringsAsFactors=FALSE,check.names=FALSE,
              colClasses=classes,na.strings=schema$missing_marker,fileEncoding="UTF-8")
stopifnot(identical(names(x),setdiff(names(y),"stratum")))
y <- y[names(x)]; rownames(x)<-rownames(y)<-NULL
stopifnot(isTRUE(all.equal(x,y,tolerance=1e-14,check.attributes=TRUE)))
blocks<-list()
.campaign_each_metric_block(file.path(streamed,"task-metrics.blocks"),function(z) {
 z<-.campaign_normalize_metric_block(z,schema$storage_classes)
 blocks[[length(blocks)+1L]]<<-z[names(x)]
})
binary<-do.call(rbind,blocks);rownames(binary)<-NULL
stopifnot(identical(x,binary))
for(nm in c("status.csv","planned-coverage.csv","scoring-attribution.csv",
            "regime-summary.csv","paired-comparisons.csv")) {
 x<-read.csv(file.path(reference,nm),stringsAsFactors=FALSE,check.names=FALSE)
 y<-read.csv(file.path(streamed,nm),stringsAsFactors=FALSE,check.names=FALSE)
 stopifnot(identical(names(x),names(y)))
 keys<-intersect(c("task_id","dataset_id","stage","method","stratum","method_a","method_b","metric"),names(x))
 if(length(keys)) {
  x<-x[do.call(order,c(unname(as.list(x[keys])),list(method="radix"))),,drop=FALSE]
  y<-y[do.call(order,c(unname(as.list(y[keys])),list(method="radix"))),,drop=FALSE]
 }
 rownames(x)<-rownames(y)<-NULL
 stopifnot(isTRUE(all.equal(x,y,tolerance=1e-13,check.attributes=TRUE)))
}
index<-readRDS(file.path(streamed,"task-metrics.rds"))
stopifnot(index$schema==2L,index$format=="metric_block_stream",index$rows==nrow(y) || index$rows>0)
''')
            run([str(comparison), str(reference), str(analysis), str(script)])
            with (analysis / "status.csv").open() as handle: status = list(csv.DictReader(handle))
            self.assertTrue(all(row["complete"] == "TRUE" for row in status))
            self.assertTrue(all(row["storage"] == "sealed_task_bundle" for row in status))
            with (analysis / "task-metrics.csv").open() as handle: metrics = list(csv.DictReader(handle))
            reused = [row for row in metrics if row["stage"] == "corrected_refit"]
            self.assertEqual(len(reused), 20)
            self.assertTrue(all(row["fit_reused"] == "TRUE" for row in reused))
            manifest = json.loads((analysis / "analysis-manifest.json").read_text())
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["manifest_sha256"], hashlib.sha256(
                (analysis / "analysis-manifest.rds").read_bytes()).hexdigest())
            checkpoint = tasks / "task_0000003/completed.rds"
            original_bytes = checkpoint.read_bytes()
            checkpoint.write_bytes(original_bytes + b"changed valid-prefix checkpoint")
            failed = subprocess.run([RSCRIPT, "--vanilla", str(script / "collect.R"),
                "--plan", str(plan), "--bootstrap-draws", "20"], env=env, text=True, capture_output=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("bytes changed", failed.stderr)
            self.assertFalse((analysis / "task-metrics.csv").exists())
            self.assertEqual(json.loads((analysis / "analysis-manifest.json").read_text())["status"],
                             "failed_or_interrupted")
            checkpoint.write_bytes(original_bytes)
            run([str(script / "collect.R"), "--plan", str(plan), "--bootstrap-draws", "20"])
            bounded_failure = subprocess.run([RSCRIPT, "--vanilla", str(script / "collect.R"),
                "--plan", str(plan), "--bootstrap-draws", "20", "--max-stratum-rows", "1"],
                env=env, text=True, capture_output=True)
            self.assertNotEqual(bounded_failure.returncode, 0)
            self.assertIn("explicit collection memory bound", bounded_failure.stderr)
            self.assertFalse((analysis / "regime-summary.csv").exists())
            self.assertEqual(json.loads((analysis / "analysis-manifest.json").read_text())["status"],
                             "failed_or_interrupted")
            missing = tasks / "task_0000003"
            hidden = root / "temporarily-unavailable-task"
            missing.rename(hidden)
            run([str(script / "collect.R"), "--plan", str(plan), "--bootstrap-draws", "20"])
            self.assertEqual(json.loads((analysis / "analysis-manifest.json").read_text())["status"],
                             "incomplete")
            hidden.rename(missing)
            all_hidden = root / "temporarily-unavailable-tasks"
            tasks.rename(all_hidden)
            run([str(script / "collect.R"), "--plan", str(plan), "--bootstrap-draws", "20"])
            self.assertFalse((analysis / "task-metrics.csv").exists())
            self.assertEqual(json.loads((analysis / "analysis-manifest.json").read_text())["status"],
                             "no_scientific_output")
            all_hidden.rename(tasks)
            run([str(script / "collect.R"), "--plan", str(plan), "--bootstrap-draws", "20"])


if __name__ == "__main__": unittest.main()
