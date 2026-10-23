# Execution-only wrapper: leaves TITAN2's tboot calculation and RNG calls intact.
checkpoint_boot_titan <- function(env, taxa, ivTot=ivTot, boot=boot,
 ncpus=ncpus, nBoot=nBoot, minSplt=minSplt, nPerm=250, memory=memory,
 imax=imax, numUnit=numUnit) {
  stopifnot(ncpus == 1, identical(memory, FALSE), nBoot > 0)
  path <- getOption("titan.checkpoint.path")
  stopifnot(is.character(path), length(path)==1)
  identity <- list(env=env,taxa=taxa,ivTot=ivTot,boot=boot,ncpus=ncpus,
    nBoot=nBoot,minSplt=minSplt,nPerm=nPerm,memory=memory,imax=imax,
    numUnit=numUnit,R=R.version.string,TITAN2=as.character(packageVersion("TITAN2")),
    RNGkind=RNGkind(),initial_rng=get(".Random.seed",envir=.GlobalEnv))
  completed <- 0L
  values <- vector("list",nBoot)
  if(file.exists(path)) {
    state <- readRDS(path)
    stopifnot(identical(identity,state$identity),length(state$values)==nBoot,
      state$completed>=0,state$completed<=nBoot,
      all(!vapply(state$values[seq_len(state$completed)],is.null,logical(1))))
    values <- state$values
    completed <- state$completed
    assign(".Random.seed",state$rng,envir=.GlobalEnv)
    cat("RESUME",completed,"/",nBoot,"\n")
    flush.console()
  }
  tboot <- get("tboot",envir=asNamespace("TITAN2"))
  if(completed<nBoot) for(i in seq.int(completed+1L,nBoot)) {
    values[[i]] <- tboot(i,env=env,taxa=taxa,ivTot=ivTot,
      minSplt=minSplt,nPerm=nPerm,memory=memory,imax=imax)
    stop_after <- getOption("titan.checkpoint.stop_after",Inf)
    if(i %% 10L == 0L || i==nBoot || i==stop_after) {
      state <- list(identity=identity,completed=i,values=values,
                    rng=get(".Random.seed",envir=.GlobalEnv))
      tmp <- paste0(path,".tmp")
      saveRDS(state,tmp,compress=FALSE)
      if(!file.rename(tmp,path)) stop("Atomic checkpoint rename failed")
      writeLines(as.character(i),paste0(path,".progress"))
      cat("CHECKPOINT",i,"/",nBoot,as.character(Sys.time()),"\n")
      flush.console()
    }
    if(i==stop_after) stop("INTENTIONAL_CHECKPOINT_TEST_STOP")
  }
  list(0,values)
}
