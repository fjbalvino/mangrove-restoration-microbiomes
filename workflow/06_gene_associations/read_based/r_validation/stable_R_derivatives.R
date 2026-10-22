suppressPackageStartupMessages(library(lmerTest))
root <- if(length(commandArgs(TRUE))) commandArgs(TRUE)[1] else '.'
m<-read.csv(file.path(root,'prepared/metadata.csv'));sel<-read.csv(file.path(root,'validation/stable_refits.csv'));Y<-read.csv(file.path(root,'validation/responses.csv'),check.names=FALSE)
m$depth<-factor(m$depth_cm,levels=c(5,20,40));m$locality<-factor(m$locality,levels=c('Carmen','Cozumel','Tuxpan'));m$profile_id<-factor(m$profile_id)
Z<-model.matrix(~0+profile_id,m);ZZ<-tcrossprod(Z);rows<-list();k<-0L
for(ii in seq_len(nrow(sel))){
 ax<-sel$axis[ii];idx<-sel$feature_index[ii];kind<-sel$model[ii];m$x<-as.numeric(scale(m[[ax]]));m$y<-as.numeric(Y[Y$feature_index==idx,-1])
 f<-suppressMessages(suppressWarnings(lmer(if(kind=='common') y~x+depth+locality+(1|profile_id) else y~x*depth+locality+(1|profile_id),data=m,REML=TRUE,control=lme4::lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000)))))
 X<-model.matrix(f);y<-m$y;n<-nrow(X);p<-ncol(X)
 dev<-function(par){V<-diag(n)+par[1]*ZZ;W<-solve(V);AA<-crossprod(X,W%*%X);bb<-solve(AA,crossprod(X,W%*%y));e<-y-X%*%bb;as.numeric(determinant(V,logarithm=TRUE)$modulus+determinant(AA,logarithm=TRUE)$modulus+crossprod(e,W%*%e)/par[2]^2+(n-p)*log(2*pi*par[2]^2))}
 covb<-function(par){V<-diag(n)+par[1]*ZZ;par[2]^2*solve(crossprod(X,solve(V,X)))}
 par<-c(unname(lme4::getME(f,'theta'))^2,sigma(f))
 H<-numDeriv::hessian(dev,par,method.args=list(eps=1e-3,d=1e-3,zero.tol=1,r=4));J<-numDeriv::jacobian(covb,par,method.args=list(eps=1e-3,d=1e-3,zero.tol=1,r=4))
 f@vcov_varpar<-2*solve(H);f@Jac_list<-lapply(1:2,function(j)matrix(J[,j],p,p))
 if(kind=='common'){L<-matrix(0,1,6);L[1,2]<-1}else{L<-matrix(0,3,8);L[,2]<-1;L[2,7]<-1;L[3,8]<-1}
 for(j in seq_len(nrow(L))){tt<-contest1D(f,L[j,]);k<-k+1L;rows[[k]]<-data.frame(axis=ax,feature_index=idx,model=kind,contrast=if(kind=='common')'common' else paste0('slope',c(5,20,40)[j]),df=tt$df,p=tt$'Pr(>|t|)')}
 if(kind=='interaction'){tt<-contestMD(f,rbind(L[2,]-L[1,],L[3,]-L[1,]));k<-k+1L;rows[[k]]<-data.frame(axis=ax,feature_index=idx,model=kind,contrast='joint',df=tt$DenDF,p=tt$'Pr(>F)')}
}
write.csv(do.call(rbind,rows),file.path(root,'validation/R_stable_derivatives.csv'),row.names=FALSE)
