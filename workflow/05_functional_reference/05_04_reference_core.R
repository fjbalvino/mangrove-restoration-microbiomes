#!/usr/bin/env Rscript
# Core analyses specified in the current manuscript, reconstructed because
# the earlier repository did not contain their original producer script.
# Read-based length-adjusted fragments; full-layer denominators; site-level tests.
options(stringsAsFactors=FALSE, warn=1)
args <- commandArgs(trailingOnly=TRUE)
arg <- function(flag, default) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (length(i) != 1L || i == length(args)) stop("Invalid argument: ", flag)
  args[i+1L]
}
base <- arg("--base",getwd())
results <- arg("--out_root", file.path(base,"local_runs","functional_reference"))
matrix_dir <- arg("--matrix_dir", file.path(base,"source_data","functional","read_based"))
out <- file.path(results,"605_reference_core")
if (dir.exists(out)) stop("Output exists; use a fresh --out_root")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
prep <- readRDS(file.path(results,"602_composition_reference/rds/602_primary_preprocessed_data.rds"))
md <- prep$metadata
stages <- levels(md$restoration4)
threshold <- 0.0001
summary <- list(); coverage <- list(); membership <- list(); profile_stats<-list()
sample_distances<-list(); profile_distances<-list(); tests<-list(); overlaps<-list(); stage_sets<-list()
write_tsv <- function(x,file) write.table(x,file.path(out,file),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
clr <- function(x) {z<-log(x+1);sweep(z,1,rowMeans(z),"-")}
partial_f <- function(y,dat) {
  f0<-lm(y~locality,data=dat);f1<-lm(y~locality+restoration4,data=dat)
  c(F=((sum(resid(f0)^2)-sum(resid(f1)^2))/3)/(sum(resid(f1)^2)/8),df1=3,df2=8)
}
for (layer in c("KO","PFAM")) {
  cat("Core analysis:",layer,"\n")
  ontology<-if(layer=="KO") "KEGG_ko" else "PFAM"
  x<-read.delim(file.path(matrix_dir,paste0(layer,"_complete_RPK.tsv")),check.names=FALSE,row.names=1)
  m<-t(as.matrix(x[,md$sample_id,drop=FALSE]))
  stopifnot(identical(rownames(m),md$sample_id),all(is.finite(m)),all(m>=0))
  rel<-sweep(m,1,rowSums(m),"/")
  detect<-rel>=threshold
  cs<-md$restoration4=="Conserved";ds<-md$restoration4=="Degraded"
  stopifnot(sum(cs)==9,sum(ds)==9)
  cores<-sapply(stages,function(st) colSums(detect[md$restoration4==st,,drop=FALSE])==sum(md$restoration4==st))
  core<-cores[,"Conserved"]
  cov<-rowSums(detect[,core,drop=FALSE])/sum(core)*100
  coverage[[layer]]<-data.frame(layer=layer,md[,c("sample_id","profile_id","locality","depth_cm","restoration4")],detected_core_features=rowSums(detect[,core,drop=FALSE]),core_features=sum(core),coverage_percent=cov)
  summary[[layer]]<-data.frame(layer=layer,total_features=ncol(m),core_features=sum(core),universally_detected=sum(colSums(detect)==51),degraded_median_percent=median(cov[ds]),degraded_min_percent=min(cov[ds]),degraded_max_percent=max(cov[ds]))
  membership[[layer]]<-data.frame(layer=layer,feature_id=colnames(m),core,universally_detected=colSums(detect)==51,detected_samples=colSums(detect),detected_degraded=colSums(detect[ds,,drop=FALSE]))
  for (st in stages) stage_sets[[paste(layer,st)]]<-data.frame(layer=layer,stage=st,n_samples=sum(md$restoration4==st),n_profiles=length(unique(md$profile_id[md$restoration4==st])),core_features=sum(cores[,st]))
  keys<-apply(cores,1,function(v) paste(as.integer(v),collapse=""))
  tab<-table(keys[keys!="00000"])
  overlaps[[layer]]<-data.frame(layer=layer,stage_order=paste(stages,collapse="|"),membership_pattern=names(tab),features=as.integer(tab))
  z<-clr(m[,core,drop=FALSE])
  d<-as.matrix(dist(z));refidx<-vapply(seq_len(nrow(md)),function(i) {
    k<-which(md$locality==md$locality[i]&md$depth_cm==md$depth_cm[i]&cs)
    stopifnot(length(k)==1);k
  },integer(1))
  dist_ref<-d[cbind(seq_len(nrow(md)),refidx)]
  baseline_idx<-vapply(seq_len(nrow(md)),function(i) {
    k<-which(md$locality==md$locality[i]&md$depth_cm==md$depth_cm[i]&ds)
    stopifnot(length(k)==1);k
  },integer(1))
  sample_tbl<-data.frame(layer=layer,md[,c("sample_id","profile_id","locality","depth_cm","restoration4")],conserved_sample_id=md$sample_id[refidx],degraded_sample_id=md$sample_id[baseline_idx],distance_to_conserved=dist_ref,degraded_distance_to_conserved=dist_ref[baseline_idx],distance_reduction=dist_ref[baseline_idx]-dist_ref,used_for_inference=!cs)
  sample_tbl$distance_to_conserved[cs]<-NA_real_
  sample_tbl$distance_reduction[cs]<-NA_real_
  sample_distances[[layer]]<-sample_tbl
  pt<-aggregate(cbind(distance_to_conserved,distance_reduction)~profile_id+locality+restoration4,data=sample_tbl[!cs,],FUN=mean)
  pt<-pt[order(pt$locality,pt$profile_id),];pt$restoration4<-droplevels(pt$restoration4)
  stopifnot(nrow(pt)==14,qr(model.matrix(~locality+restoration4,pt))$rank==6)
  perms<-readRDS(file.path(results,"602_composition_reference/rds",paste0("602_",ontology,"_reference_permutations.rds")))
  check_order<-read.delim(file.path(results,"602_composition_reference/tables/602_reference_distance_by_profile.tsv"))
  stopifnot(identical(as.character(pt$profile_id),check_order$profile_id[check_order$ontology==ontology]))
  reduced<-lm(distance_to_conserved~locality,data=pt)
  obs<-partial_f(pt$distance_to_conserved,pt)
  sim<-apply(perms,1,function(p)partial_f(fitted(reduced)+resid(reduced)[p],pt)[["F"]])
  tests[[layer]]<-data.frame(layer=layer,n_profiles=nrow(pt),df_stage=3,df_residual=8,F=obs[["F"]],p_value=(1+sum(sim>=obs[["F"]]-sqrt(.Machine$double.eps)))/(nrow(perms)+1),permutations=nrow(perms),response="profile mean core distance to matched Conserved",randomization="locality-only residuals within locality")
  profile_distances[[layer]]<-data.frame(layer=layer,pt)
  profile_stats[[layer]]<-aggregate(distance_reduction~locality+restoration4,data=pt,FUN=mean)
  profile_stats[[layer]]$layer<-layer
  write_tsv(data.frame(sample_id=md$sample_id,z,check.names=FALSE),paste0(layer,"_core_CLR.tsv"))
}
write_tsv(do.call(rbind,summary),"core_summary.tsv")
write_tsv(do.call(rbind,coverage),"core_coverage_by_sample.tsv")
write_tsv(do.call(rbind,membership),"core_membership.tsv")
write_tsv(do.call(rbind,stage_sets),"stage_core_sizes.tsv")
write_tsv(do.call(rbind,overlaps),"stage_core_exact_intersections.tsv")
write_tsv(do.call(rbind,sample_distances),"core_distances_by_sample.tsv")
write_tsv(do.call(rbind,profile_distances),"core_distances_by_profile.tsv")
write_tsv(do.call(rbind,profile_stats),"core_reduction_by_locality_stage.tsv")
write_tsv(do.call(rbind,tests),"core_reference_stage_tests.tsv")
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
print(do.call(rbind,summary));print(do.call(rbind,tests))
cat("CORE_SECTION_COMPLETE\n")
