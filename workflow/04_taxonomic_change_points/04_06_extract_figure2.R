#!/usr/bin/env Rscript
# Extract one consistent set of Figure 2 tables directly from complete TITAN2 fits.
# Usage: Rscript 04_06_extract_figure2.R run_manifest.csv output_directory
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
manifest_path <- normalizePath(args[1], mustWork=TRUE)
manifest <- read.csv(manifest_path, check.names=FALSE, stringsAsFactors=FALSE)
out <- args[2]; dir.create(out, recursive=TRUE, showWarnings=FALSE)
gradients <- c('MHI_local','moisture_stress_PC1','nutrients_redox_PC1',
               'physicochemical_PC1','vegetation_landscape_PC1','water_inundation_PC1')
strata <- c('Global','5','20','40')
expected <- as.vector(outer(strata,gradients,paste,sep='__'))
keys <- paste(manifest$stratum,manifest$gradient,sep='__')
stopifnot(nrow(manifest)==24L,!anyDuplicated(keys),setequal(keys,expected))
thresholds <- indicators <- checks <- list()
for(i in seq_len(nrow(manifest))) {
  r <- manifest[i,]
  object_path <- if (grepl("^/", r$rds)) r$rds else file.path(dirname(manifest_path), r$rds)
  fit <- readRDS(object_path)
  par <- fit$arguments
  wanted <- c(minSplt=5,numPerm=999,boot=1,nBoot=500,imax=0,ivTot=0,
              pur.cut=.95,rel.cut=.95,ncpus=1,memory=0)
  stopifnot(all(names(wanted) %in% names(par)),all(par[names(wanted)]==wanted))
  spp <- as.data.frame(fit$sppmax,check.names=FALSE)
  ntax <- c(Global=1031L,'5'=1015L,'20'=978L,'40'=997L)[[r$stratum]]
  stopifnot(nrow(fit$taxa)==17L,ncol(fit$taxa)==ntax,nrow(spp)==ntax,
            identical(dim(fit$metricArray),c(ntax,4L,500L)),
            identical(dim(fit$maxSumz),c(500L,2L)),
            identical(dim(fit$maxFsumz),c(500L,2L)),
            !anyDuplicated(rownames(spp)),all(is.finite(fit$taxa)))
  retained <- spp$purity>=.95 & spp$reliability>=.95
  stopifnot(identical(as.numeric(ifelse(retained,spp$maxgrp,0)),as.numeric(spp$filter)))
  # Bootstrap resampling must reproduce the stored purity/reliability exactly.
  pure <- ifelse(spp$maxgrp==1,
                 rowMeans(fit$metricArray[,1,,drop=FALSE]==1),
                 rowMeans(fit$metricArray[,1,,drop=FALSE]==2))
  reli <- rowMeans(fit$metricArray[,4,,drop=FALSE]<.05,na.rm=TRUE)
  stopifnot(max(abs(pure-spp$purity))<1e-12,max(abs(reli-spp$reliability))<1e-12)
  indicators[[i]] <- data.frame(stratum=r$stratum,gradient=r$gradient,
    taxon_id=rownames(spp),direction=ifelse(spp$maxgrp==1,'z-','z+'),
    cp=spp$zenv.cp,z_score=spp$zscore,purity=spp$purity,reliability=spp$reliability,
    frequency=spp$freq,retained=retained,
    q05=spp[['5%']],q10=spp[['10%']],q50=spp[['50%']],q90=spp[['90%']],q95=spp[['95%']])
  cp <- fit$sumz.cp
  for(d in 1:2) {
    sign <- c('-','+')[d]; n <- sum(spp$filter==d)
    vals <- as.numeric(cp[paste0('fsumz',sign),c('cp','0.05','0.10','0.50','0.90','0.95')])
    # TITAN can return placeholders for a filtered direction without indicators.
    if(n==0L) vals[]<-NA_real_
    if(n>0L) stopifnot(all(is.finite(vals)),all(diff(vals[-1])>=-1e-12))
    thresholds[[length(thresholds)+1L]] <- data.frame(stratum=r$stratum,
      gradient=r$gradient,direction=paste0('z',sign),n_profiles=17L,n_taxa=ntax,
      n_indicators=n,cp=vals[1],q05=vals[2],q10=vals[3],q50=vals[4],q90=vals[5],q95=vals[6])
  }
  checks[[i]] <- data.frame(stratum=r$stratum,gradient=r$gradient,n_profiles=17L,
    n_taxa=ntax,nBoot=500L,numPerm=999L,retained=sum(retained),
    object_md5=unname(tools::md5sum(object_path)),status='PASS')
}
all_ind <- do.call(rbind,indicators)
write.csv(do.call(rbind,thresholds),file.path(out,'Figure2_source_data.csv'),row.names=FALSE,na='')
write.csv(all_ind,file.path(out,'TITAN_all_taxon_gradient_associations.csv'),row.names=FALSE,na='')
write.csv(all_ind[all_ind$retained,],file.path(out,'TITAN_retained_indicators.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,checks),file.path(out,'TITAN_extraction_checks.csv'),row.names=FALSE)
cat('PASS: all 24 complete fits, all bootstrap metrics, all counts and directions verified.\n')
