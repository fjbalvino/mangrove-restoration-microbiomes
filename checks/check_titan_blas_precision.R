#!/usr/bin/env Rscript
# Replay the first five bootstrap draws from an external archived TITAN2 fit.
# Run in separate processes with different BLAS backends to diagnose tie sensitivity.
# Usage: Rscript check_titan_blas_precision.R reference.rds seed backend_label output_prefix
args <- commandArgs(TRUE)
stopifnot(length(args)==4L)
suppressPackageStartupMessages(library(TITAN2))
stopifnot(as.character(packageVersion('TITAN2'))=='2.4.4')
reference <- readRDS(args[1]); seed <- as.integer(args[2]); label <- args[3]
out <- args[4]; dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
par <- reference$arguments
stopifnot(par['nBoot']>=5,par['ncpus']==1,par['memory']==0)
set.seed(seed)
invisible(capture.output(fit <- titan(reference$env,reference$taxa,
 minSplt=par['minSplt'],numPerm=par['numPerm'],nBoot=5,ncpus=1,memory=FALSE,
 imax=as.logical(par['imax']),ivTot=as.logical(par['ivTot']),
 pur.cut=par['pur.cut'],rel.cut=par['rel.cut'],messaging=FALSE)))
result <- do.call(rbind,lapply(seq_len(4),function(j) {
 x<-fit$metricArray[,j,];y<-reference$metricArray[,j,seq_len(5)]
 data.frame(backend=label,metric=j,n_different=sum(x!=y,na.rm=TRUE),
   max_absolute_difference=max(abs(x-y),na.rm=TRUE),
   n_missingness_mismatches=sum(is.na(x)!=is.na(y)))
}))
write.csv(result,paste0(out,'.csv'),row.names=FALSE)
writeLines(c(paste0('LD_PRELOAD=',Sys.getenv('LD_PRELOAD')),
             capture.output(sessionInfo())),paste0(out,'_session.txt'))
print(result)
