"""Design-only scheduling tests; no simulation or inference is performed."""
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
BUILDER = REPO / "scripts/campaign/build_priority_schedule.R"


@unittest.skipUnless(shutil.which("Rscript"), "Rscript required for design projection")
class PriorityBuilderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.storage = tempfile.TemporaryDirectory()
        root = Path(cls.storage.name)
        original = root / "upstream/R"; original.mkdir(parents=True)
        (original / "run_simulation.R").write_text("""run_simulation <- function(S,phi,n_iter,p_causal,model) {
if (model == 'sparse') expand.grid(S=S,phi=phi,iter=seq_len(n_iter))
else expand.grid(S=S,phi=phi,p_causal=p_causal,iter=seq_len(n_iter)) }
""")
        args = root / "args.rds"
        prepare = root / "args.R"
        prepare.write_text('''a <- commandArgs(TRUE)
for(x in list.files(file.path(a[1L],"R"),"[.]R$",full.names=TRUE)) source(x)
settings <- stats::setNames(rep(list(list()),length(.FM_REGISTRY)),names(.FM_REGISTRY))
fields <- c(finemap="finemap_path",paintor="paintor_path",caviar="caviar_path",
 sparsepro="sparsepro_dir",finemap_inf="finemap_inf_dir",beatrice="beatrice_dir")
for(m in names(fields)) settings[[m]][[fields[[m]]]] <- paste0("/nonexistent-priority-fixture/",m)
saveRDS(settings,a[2L])
''')
        ready = subprocess.run(["Rscript", "--vanilla", str(prepare), str(REPO), str(args)],
                               text=True, capture_output=True)
        if ready.returncode: raise RuntimeError(ready.stderr)
        cls.full_plan = root / "full/plan.rds"
        inputs = root / "planning-inputs"; inputs.mkdir()
        (inputs / "source.txt").write_text("planning-only provenance fixture\n")
        made = subprocess.run(["Rscript", "--vanilla", str(REPO / "scripts/campaign/create_plan.R"),
            "--out", str(cls.full_plan.parent), "--original-source", str(original.parent),
            "--method-args", str(args), "--vcf-dir", str(inputs), "--genetic-map-dir", str(inputs)],
            text=True, capture_output=True)
        if made.returncode: raise RuntimeError(made.stderr)

    @classmethod
    def tearDownClass(cls): cls.storage.cleanup()

    def test_full_design_counts_and_whole_boundary_chunks(self):
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / "design-test.R"
            script.write_text('''args <- commandArgs(TRUE)
source(args[1L]); p <- readRDS(args[2L])
s <- .campaign_priority_design(p)
stopifnot(identical(s$summary$priority_method_tasks, 32324L),
          identical(s$summary$priority_generation_units, 132L),
          identical(s$summary$audit_actual_scenarios, 2160L),
          identical(s$summary$audit_boundary_extra_scenarios, 0L),
          identical(as.numeric(s$counts$locus_method_attempts), c(268400,1296000,25704000)))
stopifnot(all(p$tasks$stage[match(head(s$method_task_indices,1220L),p$tasks$index)] == "challenge"))
# Deliberately move the first five-scenario chunk across a fixed audit boundary.
# All three stages/methods sharing it remain whole; the extra scenario is visible.
old <- p$tasks$dataset_ids[1L]
ids <- strsplit(old,"|",fixed=TRUE)[[1L]]
replacement <- paste(c(ids[1:4],"reconstructed_row001_sc006"),collapse="|")
p$tasks$dataset_ids[p$tasks$dataset_ids == old] <- replacement
b <- .campaign_priority_design(p)
stopifnot(b$summary$audit_boundary_extra_scenarios == 1L,
          "reconstructed_row001_sc006" %in% b$audit_datasets$dataset_id)
selected <- p$tasks[p$tasks$index %in% b$method_task_indices,]
stopifnot(all(p$tasks$index[p$tasks$dataset_ids == replacement] %in% selected$index),
          all(lengths(strsplit(selected$dataset_ids[selected$dataset_ids==replacement],"|",fixed=TRUE)) == 5L))
cat("Full design and intact boundary chunks verified\\n")
''')
            result = subprocess.run(["Rscript", "--vanilla", str(script), str(BUILDER),
                str(self.full_plan)], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_new_builder_validates_distinct_unchanged_worker_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); fixture = root / "fixture.R"
            fixture.write_text('''args <- commandArgs(TRUE); repo <- args[1L]; root <- args[2L]
.libPaths(c(file.path(repo,".local-lib"),.libPaths()))
worker <- file.path(root,"frozen-worker"); dir.create(file.path(worker,"R"),recursive=TRUE)
for (name in c("campaign.R","challenge_io.R"))
  file.copy(file.path(repo,"R",name),file.path(worker,"R",name))
original <- file.path(root,"upstream"); dir.create(file.path(original,"R"),recursive=TRUE)
writeLines(c("run_simulation <- function(S,phi,n_iter,p_causal,model) {",
 "if (model == 'sparse') expand.grid(S=S,phi=phi,iter=seq_len(n_iter))",
 "else expand.grid(S=S,phi=phi,p_causal=p_causal,iter=seq_len(n_iter)) }"),
 file.path(original,"R/run_simulation.R"))
for(x in list.files(file.path(worker,"R"),"[.]R$",full.names=TRUE)) source(x)
p <- readRDS(args[3L])
p$repo <- normalizePath(worker); p$original_source <- normalizePath(original)
p$corrected_source <- .campaign_snapshot(p$repo)
p$original_source_hashes <- .campaign_snapshot(p$original_source)
p$method_args <- list(); p$runtime_fingerprint <- .campaign_runtime_fingerprint(p$method_args)
p$plan_sha256 <- NULL; p$plan_sha256 <- .campaign_object_sha(p)
campaign <- file.path(root,"campaign"); dir.create(campaign)
saveRDS(p,file.path(campaign,"plan.rds"),version=2)
writeLines(p$plan_sha256,file.path(campaign,"plan.sha256"))
for(name in c("tasks","datasets","generation"))
 write.csv(p[[name]],file.path(campaign,paste0(name,".csv")),row.names=FALSE)
names <- c("tasks.csv","datasets.csv","generation.csv")
jsonlite::write_json(list(schema=1L,plan_sha256=p$plan_sha256,
 plan_rds_sha256=.challenge_sha256(file.path(campaign,"plan.rds")),
 files=stats::setNames(as.list(vapply(file.path(campaign,names),.challenge_sha256,"")),names)),
 file.path(campaign,"dispatch-manifest.json"),auto_unbox=TRUE)
''')
            prepared = subprocess.run(["Rscript", "--vanilla", str(fixture), str(REPO), str(root), str(self.full_plan)],
                                      text=True, capture_output=True)
            self.assertEqual(prepared.returncode, 0, prepared.stderr)
            plan = root / "campaign/plan.rds"; worker = root / "frozen-worker"
            before = {str(p): p.read_bytes() for p in worker.rglob("*") if p.is_file()}
            before[str(plan)] = plan.read_bytes()
            output = root / "external-priority.json"
            command = ["Rscript", "--vanilla", str(BUILDER), "--plan", str(plan),
                       "--worker-repo", str(worker), "--out", str(output)]
            result = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = json.loads(output.read_text())
            self.assertEqual(manifest["schema"], 1)
            self.assertEqual(len(manifest["method_task_indices"]), 32324)
            self.assertEqual(len(manifest["generation_indices"]), 132)
            self.assertFalse(manifest["summary"]["schedule_uses_results"])
            self.assertTrue(manifest["summary"]["full_remainder_retained"])
            for path, content in before.items(): self.assertEqual(Path(path).read_bytes(), content)
            duplicate = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(duplicate.returncode, 0)
            self.assertIn("immutable", duplicate.stderr)
            (root / "campaign/tasks.csv").write_text("tampered projection\n")
            command[-1] = str(root / "different-priority.json")
            changed = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(changed.returncode, 0)
            self.assertIn("Changed frozen projection", changed.stderr)


if __name__ == "__main__": unittest.main()
