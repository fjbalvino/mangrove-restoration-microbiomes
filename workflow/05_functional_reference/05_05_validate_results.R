#!/usr/bin/env Rscript
# Compare the run with reviewed results and independently reconstruct Figure 3c,d.
options(stringsAsFactors=FALSE)
args <- commandArgs(trailingOnly=TRUE)
arg <- function(flag, default) {
  i <- which(args == flag)
  if (!length(i)) return(default)
  if (length(i) != 1L || i == length(args)) stop("Invalid argument: ", flag)
  args[i+1L]
}
base <- arg("--base", getwd())
out <- arg("--out_root", file.path(base,"local_runs","functional_reference"))
read_tsv <- function(path) read.delim(path,check.names=FALSE)
prep <- readRDS(file.path(out,"602_composition_reference/rds/602_primary_preprocessed_data.rds"))
res <- readRDS(file.path(out,"602_composition_reference/rds/602_analysis_results.rds"))
md <- prep$metadata; p <- res$permutations
stopifnot(nrow(md)==51L, length(unique(md$profile_id))==17L,
          identical(dim(p),c(9999L,51L)), nrow(unique(p))==9999L,
          qr(model.matrix(~locality+depth_f+restoration4,md))$rank==9L)
idx <- split(seq_len(nrow(md)),md$profile_id)
for (i in seq_len(nrow(p))) {
  z <- p[i,]
  stopifnot(setequal(z,1:51),!identical(as.integer(z),1:51),
            all(md$locality[z]==md$locality), all(md$depth_cm[z]==md$depth_cm),
            all(vapply(idx,function(j)length(unique(md$profile_id[z[j]]))==1L,logical(1))))
}
fixture <- read_tsv(file.path(base,"checks/functional_reference_expected.tsv"))
for (rel in unique(fixture$table)) {
  actual <- read_tsv(file.path(out,rel))
  tests <- fixture[fixture$table==rel,]
  for (j in seq_len(nrow(tests))) {
    t <- tests[j,]; cols <- strsplit(t$key_columns,"|",fixed=TRUE)[[1]]
    key <- apply(actual[,cols,drop=FALSE],1,paste,collapse="|")
    at <- which(key==t$key_values)
    if (length(at)!=1L) stop("Missing or duplicated result: ",rel," ",t$key_values)
    value <- actual[[t$column]][at]
    if (!isTRUE(all.equal(as.numeric(value),as.numeric(t$expected),tolerance=1e-9))) {
      stop("Reviewed result differs: ",rel," ",t$key_values," ",t$column,
           " expected=",t$expected," actual=",value)
    }
  }
}
contrasts <- read_tsv(file.path(out,"604_depth_reference_contrasts/tables/604_distances_and_local_contrasts.tsv"))
stopifnot(nrow(contrasts)==84L)
reconstructed <- list()
for (layer in c("KO","PFAM")) {
  ontology <- if(layer=="KO") "KEGG_ko" else "PFAM"
  files <- file.path(base,"source_data/functional/read_based",paste0(layer,"_complete_",c("RPK","counts"),".tsv"))
  matrices <- lapply(files,function(f) {
    x <- read_tsv(f); ids <- x[[1]]
    stopifnot(!anyDuplicated(ids),setequal(names(x)[-1],md$sample_id))
    m <- t(as.matrix(x[,md$sample_id,drop=FALSE]));colnames(m)<-ids;m
  })
  x <- matrices[[1]]; counts <- matrices[[2]]
  stopifnot(identical(dimnames(x),dimnames(counts)))
  keep <- colSums(x>0)>=26 & colSums(counts)>=10 & apply(log1p(x),2,var)>0
  logx <- log(x[,keep,drop=FALSE]+1); clr <- sweep(logx,1,rowMeans(logx),"-")
  stopifnot(isTRUE(all.equal(unname(clr),unname(prep[[ontology]]$clr),tolerance=1e-10)))
  z <- contrasts[contrasts$ontology==ontology,]
  for (i in seq_len(nrow(z))) {
    row <- z[i,]; a <- match(row$sample_id,md$sample_id)
    c <- which(md$locality==row$locality & md$depth_cm==row$depth_cm & md$restoration4=="Conserved")
    d <- which(md$locality==row$locality & md$depth_cm==row$depth_cm & md$restoration4=="Degraded")
    stopifnot(length(c)==1L,length(d)==1L,md$sample_id[c]==row$conserved_sample_id,
              md$sample_id[d]==row$degraded_sample_id)
    da <- sqrt(sum((clr[a,]-clr[c,])^2)); dd <- sqrt(sum((clr[d,]-clr[c,])^2))
    stopifnot(abs(da-row$distance)<1e-9,abs(dd-row$degraded_distance)<1e-9,
              abs(da-dd-row$delta_distance)<1e-9)
  }
  reconstructed[[layer]] <- z[order(z$sample_id),c("sample_id","delta_distance")]
}
stopifnot(identical(reconstructed$KO$sample_id,reconstructed$PFAM$sample_id),
          !isTRUE(all.equal(reconstructed$KO$delta_distance,reconstructed$PFAM$delta_distance)))
report <- data.frame(check=c("51 samples and 17 complete profiles","9999 valid distinct profile permutations",
                            "full-rank composition model","reviewed numerical expectations",
                            "84 contrasts independently reconstructed","KO and PFAM are distinct matrices"),pass=TRUE)
write.table(report,file.path(out,"functional_reference_validation.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
print(report)
cat("FUNCTIONAL_REFERENCE_VALIDATED; numerical expectations:",nrow(fixture),"\n")
