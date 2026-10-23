# Memory-bounded equivalent of TITAN2 2.4.4 getivz (GPL-2).
# Original authors: Matthew E. Baker, Ryan S. King and David Kahle.
# Process the same permutation matrices one partition at a time instead of
# retaining the full permutation x taxon x partition array. All native
# indval/permiv calculations, sd, colMeans and RNG calls remain unchanged.
# Optional execution aid; --stream_permutations false uses the unmodified package.
stream_getivz <- function(clss, spp, ivTot=ivTot, nPerm=nPerm, numClass=numClass,
                         imax=imax, messaging=FALSE) {
  ns <- asNamespace('TITAN2'); n <- ncol(spp)
  obs <- apply(clss,2,get(if(ivTot) 'indvals' else 'indval', ns),spp=spp)
  obsiv <- obs[(n+1):(2*n),,drop=FALSE]; obsmx <- obs[1:n,,drop=FALSE]
  ivz <- iv.p <- matrix(NA_real_,n,ncol(clss))
  permiv <- get('permiv',ns)
  for(j in seq_len(ncol(clss))) {
    piv <- permiv(clss[,j],spp,ivTot=ivTot,nPerm=nPerm)
    mn <- colMeans(piv,na.rm=TRUE); st <- apply(piv,2,sd); st[which(st<.1)]<-1
    ivz[,j] <- (obsiv[,j]-mn)/st
    for(i in seq_len(n)) iv.p[i,j] <- (sum(piv[,i]>=obsiv[i,j],na.rm=TRUE)+1)/nPerm
  }
  rbind(obsmx,ivz,obsiv,iv.p)
}
