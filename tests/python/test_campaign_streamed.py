"""Bounded collector serialization contracts; no benchmark measurements."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
RSCRIPT = os.environ.get("FMB_RSCRIPT") or shutil.which("Rscript")


@unittest.skipUnless(RSCRIPT, "Rscript required for collector contract tests")
class StreamedCollectionTests(unittest.TestCase):
    def run_r(self, code):
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / "test.R"
            script.write_text(code)
            result = subprocess.run([RSCRIPT, "--vanilla", str(script), str(REPO), directory],
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)

    def test_double_logical_missing_text_and_block_round_trip(self):
        self.run_r(r'''
a<-commandArgs(TRUE); source(file.path(a[1L],"scripts/campaign/collect_streamed.R"))
x<-data.frame(real=c(pi,1+.Machine$double.eps,1e-300,Inf,-Inf,NaN,NA_real_,0),
 integer=c(1L,2L,NA_integer_,4L,5L,6L,7L,8L),
 flag=c(TRUE,FALSE,NA,TRUE,FALSE,TRUE,FALSE,TRUE),
 text=c("NA","NaN","comma, quote\"","new\nline","β","",NA_character_,"FALSE"),
 stringsAsFactors=FALSE)
missing<-"<unique-test-missing>"; csv<-file.path(a[2L],"metrics.csv")
.campaign_csv_block(x[1:4,],csv,missing=missing)
.campaign_csv_block(x[5:8,],csv,missing=missing)
classes<-.campaign_metric_classes(x)
y<-read.csv(csv,stringsAsFactors=FALSE,check.names=FALSE,
 colClasses=unname(ifelse(classes=="double","numeric",classes)),na.strings=missing)
# R's decimal parser can differ for subnormal/extreme exponents. The retained
# binary archive and hexadecimal partition path establish exact preservation.
stopifnot(identical(x[-3,],y[-3,]),identical(is.nan(x$real),is.nan(y$real)))
hex<-file.path(a[2L],"partition.csv")
.campaign_csv_block(x,hex,missing=missing,double_format="hexadecimal")
y<-read.csv(hex,stringsAsFactors=FALSE,check.names=FALSE,
 colClasses=unname(ifelse(classes=="double","numeric",classes)),na.strings=missing)
stopifnot(identical(x,y))
blocks<-file.path(a[2L],"metrics.blocks"); con<-file(blocks,"wb")
.campaign_write_metric_block(x[1:4,],con); .campaign_write_metric_block(x[5:8,],con);close(con)
ys<-list(); census<-.campaign_each_metric_block(blocks,function(z) ys[[length(ys)+1L]]<<-z)
stopifnot(census$blocks==2L,census$rows==8L,identical(x,do.call(rbind,ys)))
con<-file(blocks,"ab");writeBin(as.raw(c(1,2,3)),con);close(con)
err<-tryCatch(.campaign_each_metric_block(blocks,function(z) NULL),error=function(e)e)
stopifnot(inherits(err,"error"),grepl("Truncated metric block header",conditionMessage(err)))
err<-tryCatch(.campaign_csv_block(data.frame(x=missing),csv,missing=missing),error=function(e)e)
stopifnot(inherits(err,"error"),grepl("collides",conditionMessage(err)))
writer<-.campaign_partition_writer(2L)
for(i in 1:12) writer$write(x,file.path(a[2L],paste0("part",(i%%3)+1L,".csv")),missing)
stopifnot(writer$peak()==2L);writer$close();writer$close()
for(i in 1:3) {
 z<-read.csv(file.path(a[2L],paste0("part",i,".csv")),stringsAsFactors=FALSE,check.names=FALSE,
  colClasses=unname(ifelse(classes=="double","numeric",classes)),na.strings=missing)
 stopifnot(nrow(z)==4L*nrow(x),identical(x,z[1:nrow(x),]))
}
''')

    def test_infinite_loss_and_fixed_eligibility_summary(self):
        self.run_r(r'''
a<-commandArgs(TRUE); source(file.path(a[1L],"scripts/campaign/collect_streamed.R"))
.campaign_ap_failure_zero<-function(ap,success,n_positive) ifelse(n_positive>0,ifelse(success,ap,0),NA_real_)
n<-3L
x<-as.data.frame(setNames(rep(list(rep(NA_real_,n)),22L),
 c("ap","brier","recall","log_loss","brier_any_nonzero_total_effect","log_loss_any_nonzero_total_effect",
   "cs_distinct_recall","native_all_causal_coverage_reported","native_all_causal_recovery_all_attempts",
   "cs_causal_variant_recall_all_attempts","cs_support_size","runtime_seconds","tp","fp",
   "n_observed_causal","null_model_in_support","null_probability_reported","probability_target_aligned",
   "pip_event","posterior_conditioning","native_null_support","probability_score_interpretation")),
 stringsAsFactors=FALSE)
x$stage<-"challenge";x$method<-"test";x$stratum<-"null-and-positive-control"
x$success<-c(TRUE,TRUE,FALSE);x$ap<-c(.5,NA_real_,NA_real_);x$n_observed_causal<-c(1L,0L,0L)
x$brier<-c(.2,.1,NA_real_);x$log_loss<-c(Inf,.25,NA_real_);x$recall<-c(.5,NA_real_,NA_real_)
x$tp<-c(1L,0L,0L);x$fp<-c(1L,0L,0L);x$runtime_seconds<-c(1,2,3)
x$null_model_in_support<-c(TRUE,TRUE,NA);x$null_probability_reported<-c(FALSE,FALSE,FALSE)
x$probability_target_aligned<-c(TRUE,FALSE,NA)
x$pip_event<-"sparse_component_membership";x$score_target<-"sparse_component_membership"
x$posterior_conditioning<-"unconditional";x$native_null_support<-"included"
x$probability_score_interpretation<-"operational";x$cluster_id<-c("locus1","locus2","locus2")
x$panel_cluster_id<-"panel1";x$source_locus_cluster<-c("source1","source2",NA_character_)
y<-.campaign_streamed_group(x)
stopifnot(is.infinite(y$log_loss_successful),y$n_log_loss_infinite==1L,y$n_log_loss_observed==2L,
 y$n_attempted==3L,y$n_ap_eligible_attempts==1L,y$n_ap_eligible_successful==1L,
 y$ap_observed_positive_attempts_zero_for_failure==.5,y$n_genotype_clusters==2L,
 y$n_panel_clusters==1L,y$n_probability_target_aligned_successful==1L,
 y$n_probability_target_misaligned_successful==1L,y$n_probability_target_unknown_successful==0L)
''')


if __name__ == "__main__": unittest.main()
