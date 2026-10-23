#!/usr/bin/env Rscript
# Check exact numerical/RNG equivalence of execution aids with unmodified TITAN2.
# Optional first argument: RDS list(env=<vector>, taxa=<matrix>) with real inputs.
args <- commandArgs(TRUE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
root <- normalizePath(file.path(dirname(script),'..'))
lib <- file.path(root,'workflow/04_taxonomic_change_points')
suppressPackageStartupMessages(library(TITAN2))
stopifnot(as.character(packageVersion('TITAN2'))=='2.4.4')
source(file.path(lib,'stream_permutations.R'))
source(file.path(lib,'checkpoint_boot.R'))
source(file.path(lib,'prepare_global_rng.R'))
source(file.path(lib,'validated_metadata.R'))
official_getivz <- get('getivz',asNamespace('TITAN2'))
official_boot <- get('boot.titan',asNamespace('TITAN2'))
check_case <- function(env,taxa,seed,nperm,nboot) {
  fit <- function() titan(env,taxa,minSplt=5,numPerm=nperm,nBoot=nboot,
    ncpus=1,memory=FALSE,imax=FALSE,ivTot=FALSE,pur.cut=.95,rel.cut=.95,messaging=FALSE)
  assignInNamespace('getivz',official_getivz,ns='TITAN2')
  assignInNamespace('boot.titan',official_boot,ns='TITAN2')
  set.seed(seed); invisible(capture.output(reference <- fit())); final_rng <- .Random.seed
  set.seed(seed); advance_titan_rng(env,numPerm=nperm,nBoot=nboot)
  stopifnot(identical(.Random.seed,final_rng))
  assignInNamespace('getivz',stream_getivz,ns='TITAN2')
  set.seed(seed); invisible(capture.output(streamed <- fit()))
  stopifnot(identical(reference,streamed),identical(.Random.seed,final_rng))
  assignInNamespace('boot.titan',checkpoint_boot_titan,ns='TITAN2')
  checkpoint <- tempfile(fileext='.rds')
  options(titan.checkpoint.path=checkpoint,titan.checkpoint.stop_after=2L)
  set.seed(seed)
  stopped <- tryCatch({invisible(capture.output(fit()));FALSE},error=function(e) {
    if(!grepl('INTENTIONAL_CHECKPOINT_TEST_STOP',conditionMessage(e))) stop(e)
    TRUE
  })
  stopifnot(stopped,readRDS(checkpoint)$completed==2L)
  options(titan.checkpoint.stop_after=Inf)
  set.seed(seed); invisible(capture.output(resumed <- fit()))
  stopifnot(identical(reference,resumed),identical(.Random.seed,final_rng))
  unlink(c(checkpoint,paste0(checkpoint,'.progress')))
  cat('PASS: ',nrow(taxa),' units, ',ncol(taxa),' taxa, ',nperm,' permutations, ',
      nboot,' bootstraps; all 13 object components and RNG states identical.\n',sep='')
}
taxa <- matrix(((seq_len(17*12)*13) %% 41)+1,nrow=17)
taxa[seq_len(nrow(taxa)) %% 5 == 0,seq(2,12,2)] <- 0
colnames(taxa) <- paste0('taxon_',seq_len(ncol(taxa)))
check_case(sort(rep(seq_len(9),length.out=17)),taxa,seed=123,nperm=29,nboot=5)
if(length(args)>0) {
  actual <- readRDS(args[1]); check_case(actual$env,actual$taxa,seed=123,nperm=999,nboot=5)
}
meta_path <- file.path(root,'source_data/indices/metadata_003b.csv')
meta <- read.csv(meta_path,check.names=FALSE,stringsAsFactors=FALSE)
validate_titan_metadata(meta[rev(seq_len(nrow(meta))),],meta_path)
bad <- meta;bad$MHI_local[1] <- bad$MHI_local[1]+.1
stopifnot(inherits(try(validate_titan_metadata(bad,meta_path),silent=TRUE),'try-error'))
bad <- meta;bad$sample_id[1] <- bad$sample_id[2]
stopifnot(inherits(try(validate_titan_metadata(bad,meta_path),silent=TRUE),'try-error'))
cat('PASS: sample-ID joins accepted; obsolete MHI values and duplicate IDs rejected.\n')
