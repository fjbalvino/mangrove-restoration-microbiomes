suppressPackageStartupMessages(library(lmerTest))
root <- if(length(commandArgs(TRUE))) commandArgs(TRUE)[1] else '.'
m <- read.csv(file.path(root,'prepared/metadata.csv'))
selected <- read.csv(file.path(root,'validation/selected_validation_models.csv'))
Y <- read.csv(file.path(root,'validation/responses.csv'),check.names=FALSE)
m$depth <- factor(m$depth_cm,levels=c(5,20,40))
m$locality <- factor(m$locality,levels=c('Carmen','Cozumel','Tuxpan'))
m$profile_id <- factor(m$profile_id)
rows <- list();k<-0L
for(ax in unique(selected$axis)) {
 m$x <- as.numeric(scale(m[[ax]]))
 for(idx in selected$feature_index[selected$axis==ax]) {
  m$y <- as.numeric(Y[Y$feature_index==idx,-1])
  for(kind in c('common','interaction')) {
   notes <- character()
   f <- withCallingHandlers(lmer(if(kind=='common') y~x+depth+locality+(1|profile_id) else y~x*depth+locality+(1|profile_id),data=m,REML=TRUE,control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))),warning=function(w){notes<<-c(notes,conditionMessage(w));invokeRestart('muffleWarning')},message=function(w){invokeRestart('muffleMessage')})
   if(kind=='common'){L<-matrix(0,1,6);L[1,2]<-1} else {L<-matrix(0,3,8);L[,2]<-1;L[2,7]<-1;L[3,8]<-1}
   for(j in seq_len(nrow(L))){
    tt<-contest1D(f,L[j,],confint=TRUE);k<-k+1L
    rows[[k]]<-data.frame(axis=ax,feature_index=idx,model=kind,contrast=if(kind=='common')'common' else paste0('slope',c(5,20,40)[j]),estimate=tt$Estimate,se=tt$'Std. Error',df=tt$df,p=tt$'Pr(>|t|)',ci_low=tt$lower,ci_high=tt$upper,theta=lme4::getME(f,'theta'),sigma=sigma(f),singular=lme4::isSingular(f,tol=1e-4),warnings=paste(notes,collapse=' | '))
   }
   if(kind=='interaction'){
    tt<-contestMD(f,rbind(L[2,]-L[1,],L[3,]-L[1,]));k<-k+1L
    rows[[k]]<-data.frame(axis=ax,feature_index=idx,model=kind,contrast='joint',estimate=NA,se=NA,df=tt$DenDF,p=tt$'Pr(>F)',ci_low=NA,ci_high=NA,theta=lme4::getME(f,'theta'),sigma=sigma(f),singular=lme4::isSingular(f,tol=1e-4),warnings=paste(notes,collapse=' | '))
   }
  }
 }
 cat(ax,'complete\n');flush.console()
}
write.csv(do.call(rbind,rows),file.path(root,'validation/lmerTest_results.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'validation/R_session_info.txt'))
