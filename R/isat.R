isat <-
function(y, mc=TRUE, ar=NULL, ewma=NULL, mxreg=NULL,
  iis=FALSE, sis=TRUE, tis=FALSE, uis=FALSE, blocks=NULL,
  ratio.threshold=0.8, max.block.size=30, t.pval=0.001,
  wald.pval=t.pval, vcov.type=c("ordinary", "white", "newey-west"),
  do.pet=FALSE, ar.LjungB=NULL, arch.LjungB=NULL,
  normality.JarqueB=NULL, user.diagnostics=NULL,
  info.method=c("sc", "aic", "hq"), include.gum=FALSE,
  include.empty=FALSE, max.paths=NULL, parallel.options=NULL,
  tol=1e-07, LAPACK=FALSE, max.regs=NULL, print.searchinfo=TRUE,
  plot=NULL, alarm=FALSE)
{

  ##arguments:
  isat.call <- sys.call()
  vcov.type <- match.arg(vcov.type)
  info.method <- match.arg(info.method)

  ##max paths argument:
  if( !is.null(max.paths) && max.paths < 1 ){
    stop("'max.paths' cannot be smaller than 1")
  }

  ##parallel.options argument:
  if(!is.null(parallel.options)){
    if(is.numeric(parallel.options)){
      clusterSpec <- parallel.options
    }
    OScores <- detectCores()
    if(parallel.options > OScores){
      stop("parallel.options > number of cores/threads")
    }
  }

  ##estimate gum:
  mod <- arx(y, mc=mc, ar=ar, ewma=ewma, mxreg=mxreg,
    vcov.type=vcov.type, qstat.options=NULL,
    user.diagnostics=user.diagnostics, tol=tol, LAPACK=LAPACK,
    plot=FALSE)
  y <- mod$aux$y
  y.n <- mod$aux$y.n
  y.index <- mod$aux$y.index
  y.index.as.char <- as.character(y.index)
  y.name <- mod$aux$y.name
  mX <- mod$aux$mX #NULL if is.null(mX)
  mXnames <- mod$aux$mXnames #NULL if is.null(mX)
  colnames(mX) <- mXnames
  mXncol <- mod$aux$mXncol
  vcov.type <- mod$aux$vcov.type
  qstat.options <- mod$aux$qstat.options
  if(is.null(mX)){ mxkeep <- NULL }else{ mxkeep <- 1:mXncol }

  ##indicator saturation matrices:
  ISmatrices <- list()

  if(iis){ #impulse indicators
    mIIS <- matrix(0,y.n,y.n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste("iis", y.index.as.char, sep="")
    ISmatrices <- c(ISmatrices,list(IIS=mIIS))
  }

  if(sis){ #step-shift indicators
    mSIS <-matrix(0,y.n,y.n) #replace by , y.n, y.n-1 ?
    loop.indx <- 1:y.n #replace by 2:y.n ?
    tmp <- function(i){ mSIS[i,1:i] <<- 1 }
    tmp <- sapply(loop.indx,tmp)
    colnames(mSIS) <- paste("sis", y.index.as.char, sep="")
    mSIS <- mSIS[,-1]
    ISmatrices <- c(ISmatrices,list(SIS=mSIS))
  }

  if(tis){ #trend indicators
    mTIS <- matrix(0,y.n,y.n)
    v1n <- seq(1,y.n)
    loop.indx <- 1:y.n
    tmp <- function(i){
      mTIS[c(i:y.n),i] <<- v1n[1:c(y.n-i+1)]
    }
    tmp <- sapply(loop.indx,tmp)
    colnames(mTIS) <- paste("tis", y.index.as.char, sep="")
    mTIS <- mTIS[,-1]
    ISmatrices <- c(ISmatrices,list(TIS=mTIS))
  }

  ##user-defined indicators/variables:
  ##----------------------------------

  #if uis is a matrix:
  if(!is.list(uis) && !identical(as.numeric(uis),0)){

    ##handle colnames:
    uis <- as.zoo(cbind(uis))
    uis.names <- colnames(uis)
    if(is.null(uis.names)){
      uis.names <- paste("uisxreg", 1:NCOL(uis), sep="")
    }
    if(any(uis.names == "")){
      missing.colnames <- which(uis.names == "")
      for(i in 1:length(missing.colnames)){
       uis.names[i] <- paste("uisxreg", missing.colnames[i], sep="")
      }
    }
    #for the future?: uis.names <- make.names(uis.names)

    ##select sample:
    uis <- na.trim(uis, sides="both", is.na="any")
    uis.index.as.char <- as.character(index(uis))
    t1 <- which(uis.index.as.char==y.index.as.char[1])
    t2 <- which(uis.index.as.char
      ==y.index.as.char[length(y.index.as.char)])
    uis <- coredata(uis)
    uis <- window(uis, start=t1, end=t2)
    uis <- cbind(coredata(as.zoo(uis)))
    colnames(uis) <- uis.names

    #check nrow(uis):
    if(nrow(uis) != y.n) stop("nrow(uis) is unequal to no. of observations")
    ISmatrices <- c(ISmatrices,list(UIS=uis))
  }

  ##if uis is a list of matrices:
  if(is.list(uis)){
    #check nrow(uis[[i]]):
    for(i in 1:length(uis)){
      uis[[i]] <- as.matrix(coredata(as.zoo(uis[[i]])))
      if(nrow(uis[[i]]) != y.n){
        stop(paste("nrow(uis[[",i,"]]) is unequal to no. of observations",
          sep=""))
      }
    } #end check nrow
    uis.names <- paste("UIS", 1:length(uis), sep="")
    if(is.null(names(uis))){
      names(uis) <- uis.names
    }else{
      for(i in 1:length(uis)){
        if(names(uis)[i]==""){
          names(uis)[i] <- uis.names[i]
        }else{
          names(uis)[i] <- paste(uis.names[i], ".", names(uis)[i],
            sep="")
        } #close if..else
      } #close for..loop
    }
    ISmatrices <- c(ISmatrices,uis)
  }

  ##check blocks:
  if(is.list(blocks)){
    if(length(ISmatrices)!=length(blocks)){
      stop("No. of IS matrices is unequal to length(blocks)")
    }
    blocks.is.list <- TRUE
    ISblocks <- blocks
  }else{
    blocks.is.list <- FALSE
    ISblocks <- list()
  }

  ##loop on ISmatrices:
  ISfinalmodels <- list()
  for(i in 1:length(ISmatrices)){

    ##blocks:
    if(!blocks.is.list){

      ncol.adj <- NCOL(ISmatrices[[i]])

      if(is.null(blocks)){
        blockratio.value <- ncol.adj/(ratio.threshold*ncol.adj - mXncol)
        blocksize.value <- ncol.adj/min(y.n*ratio.threshold, max.block.size)
#OLD:
#        blocksize.value <- ncol.adj/max.block.size
        no.of.blocks <- max(2,blockratio.value,blocksize.value)
        no.of.blocks <- ceiling(no.of.blocks)
        no.of.blocks <- min(ncol.adj, no.of.blocks) #ensure blocks < NCOL
      }else{
        no.of.blocks <- blocks
      }

      blocksize <- ceiling(ncol.adj/no.of.blocks)
      partitions.t2 <- blocksize
      for(j in 1:no.of.blocks){
        if( blocksize*j <= ncol.adj ){
          partitions.t2[j] <- blocksize*j
        }
      }
      #check if last block contains last indicator:
      if(partitions.t2[length(partitions.t2)] < ncol.adj){
        partitions.t2 <- c(partitions.t2, ncol.adj)
      }
      blocksadj <- length(partitions.t2)
      partitions.t1 <- partitions.t2 + 1
      partitions.t1 <- c(1,partitions.t1[-blocksadj])

      tmp <- list()
      for(j in 1:blocksadj){
        tmp[[j]] <- partitions.t1[j]:partitions.t2[j]
      }
      ISblocks[[i]] <- tmp

    } #end if(!blocks.is.list)

    ##make blocks function:
    ISblocksFun <- function(j, i, ISmatrices, ISblocks, mX,
      parallel.options, tol, LAPACK, print.searchinfo, y,
      vcov.type, qstat.options, user.diagnostics, mxkeep,
      t.pval, wald.pval, do.pet, ar.LjungB, arch.LjungB,
      normality.JarqueB, info.method, include.empty,
      max.paths, max.regs){

      ##check if block contains 1 regressor:
      if( length(ISblocks[[i]][[j]])==1 ){
        tmp <- colnames(ISmatrices[[i]])[ ISblocks[[i]][[j]] ]
        mXis <- cbind(ISmatrices[[i]][, ISblocks[[i]][[j]] ])
        colnames(mXis) <- tmp
        mXis <- cbind(mX, mXis)
      }else{
        mXis <- cbind(mX,ISmatrices[[i]][, ISblocks[[i]][[j]] ])
      }

      ##apply dropvar:
      mXis <- gets:::dropvar(mXis, tol=tol, LAPACK=LAPACK,
        silent=print.searchinfo)

      ##print info:
      if(is.null(parallel.options)){
        if(print.searchinfo){
          message("\n", appendLF=FALSE)
          message(names(ISmatrices)[i],
            " block ", j, " of ", length(ISblocks[[i]]), ":",
            appendLF=TRUE)
          #message("\n", appendLF=FALSE)
        }
      }

      ##gum:
      mod <- gets:::arx(y, mxreg=mXis, vcov.type=vcov.type,
        qstat.options=qstat.options, user.diagnostics=user.diagnostics,
        tol=tol, LAPACK=LAPACK, plot=FALSE)
#future?: keep=NULL instead of mxkeep?
#      getsis <- getsm(mod, keep=NULL, t.pval=t.pval,
      getsis <- gets:::getsm(mod, keep=mxkeep, t.pval=t.pval,
        wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, info.method=info.method,
        include.empty=include.empty, max.regs=max.regs,
        estimate.specific=FALSE,
        print.searchinfo=print.searchinfo, plot=FALSE)

      if(is.null(getsis$specific.spec)){
        ISspecific.models <- NULL
#        ISspecific.models[[j]] <- NULL
      }else{
        ISspecific.models <- names(getsis$specific.spec)
#        ISspecific.models[[j]] <- names(getsis$specific.spec)
#For the future?:
#        ISgums[[j]] <- getsis$gum.mean
#        ISpaths[[j]] <- getsis$paths
#        ISterminals.results[[j]] <- getsis$terminals.results
      }

      ##return
      return(ISspecific.models)

    } #close ISblocksFun

#OLD:
#    ##print info:
#    if(print.searchinfo){
#      message("\n", appendLF=FALSE)
#      message(names(ISmatrices)[i],
#        " block ", j, " of ", length(ISblocks[[i]]), ":",
#        appendLF=TRUE)
#      message("\n", appendLF=FALSE)
#    }

    ##do gets on each block:
    if(is.null(parallel.options)){
      ISspecific.models <- lapply(1:length(ISblocks[[i]]),
        ISblocksFun, i, ISmatrices, ISblocks, mX,
        parallel.options, tol, LAPACK, print.searchinfo, y,
        vcov.type, qstat.options, user.diagnostics, mxkeep,
        t.pval, wald.pval, do.pet, ar.LjungB, arch.LjungB,
        normality.JarqueB, info.method, include.empty,
        max.paths, max.regs)
    }else{

      ##print info:
      if(print.searchinfo){
        message("\n", appendLF=FALSE)
        message("Preparing parallel computing...",
          appendLF=TRUE)
        message(names(ISmatrices)[i],
          " blocks to search in parallel: ", length(ISblocks[[i]]),
          appendLF=TRUE)
        message("Searching...", appendLF=TRUE)
        #message("\n", appendLF=FALSE)
      }

      #to do:
      # clusterArgs <- list(spec=numberOfCores, type, ...)
      # do.call("makeCluster", clusterArgs)
      blocksClust <- makeCluster(clusterSpec,
        outfile="")
      ISspecific.models <- parLapply(blocksClust,
        1:length(ISblocks[[i]]), ISblocksFun, i, ISmatrices,
        ISblocks, mX, parallel.options, tol, LAPACK,
        print.searchinfo, y, vcov.type, qstat.options,
        user.diagnostics, mxkeep, t.pval, wald.pval, do.pet,
        ar.LjungB, arch.LjungB, normality.JarqueB, info.method,
        include.empty, max.paths, max.regs)
      stopCluster(blocksClust)

    } #end if..

###START OLD:
#
#    ##gets on each block:
#    ISspecific.models <- list()
#    #for the future?: ISgums <- list(); ISpaths <- list(); ISterminals.results <- list()
#    for( j in 1:length(ISblocks[[i]]) ){
#
#      ##print info:
#      if(print.searchinfo){
#        message("\n", appendLF=FALSE)
#        message(names(ISmatrices)[i],
#          " block ", j, " of ", length(ISblocks[[i]]), ":",
#          appendLF=TRUE)
#        message("\n", appendLF=FALSE)
#      }
#
#      ##check if block contains 1 regressor:
#      if( length(ISblocks[[i]][[j]])==1 ){
#        tmp <- colnames(ISmatrices[[i]])[ ISblocks[[i]][[j]] ]
#        mXis <- cbind(ISmatrices[[i]][, ISblocks[[i]][[j]] ])
#        colnames(mXis) <- tmp
#        mXis <- cbind(mX, mXis)
#      }else{
#        mXis <- cbind(mX,ISmatrices[[i]][, ISblocks[[i]][[j]] ])
#      }
#
#      ##apply dropvar:
#      mXis <- dropvar(mXis, tol=tol, LAPACK=LAPACK,
#        silent=print.searchinfo)
#
#      ##gum:
#      mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
#        qstat.options=qstat.options, user.diagnostics=user.diagnostics,
#        tol=tol, LAPACK=LAPACK, plot=FALSE)
##future?: keep=NULL instead of mxkeep?
##      getsis <- getsm(mod, keep=NULL, t.pval=t.pval,
#      getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
#        wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
#        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
#        user.diagnostics=user.diagnostics, info.method=info.method,
#        include.empty=include.empty, paths.method=paths.method,
#        max.regs=max.regs, estimate.specific=FALSE,
#        print.searchinfo=print.searchinfo, plot=FALSE)
#
#      if(is.null(getsis$specific.spec)){
#        ISspecific.models[[j]] <- NULL
#      }else{
#        ISspecific.models[[j]] <- names(getsis$specific.spec)
##For the future?:
##        ISgums[[j]] <- getsis$gum.mean
##        ISpaths[[j]] <- getsis$paths
##        ISterminals.results[[j]] <- getsis$terminals.results
#      }
#
#    } #end for(j in 1:length(ISblocks[[i]]))
#
###END OLD

    ##print info:
    if(print.searchinfo){
      message("\n", appendLF=FALSE)
      message("GETS of union of retained ",
        names(ISmatrices)[i], " variables... ",
        appendLF=TRUE)
    }

    ##if no indicators retained from the blocks:
    if(length(ISspecific.models)==0){
      isNames <- NULL
      ISfinalmodels[[i]] <- NULL
    }

    ##when indicators/variables(uis) retained from the blocks:
    if(length(ISspecific.models)>0){

      isNames <- NULL

      #which indicators/variables(uis) retained?:
      for(j in 1:length(ISspecific.models)){
        #check if mean is non-empty:
        if(!is.null(ISspecific.models[[j]])){
          isNames <- union(isNames, ISspecific.models[[j]])
        }
      } #end for(j) loop
      isNames <- setdiff(isNames, mXnames)
      #isNamesAll[[i]] <- isNames

##BEGIN MODIFIED (panel):

      #redo gets with union of retained indicators:
      if(length(isNames)==0){
        ISfinalmodels[[i]] <- mXnames
      }else{
        mXisNames <- c(mXnames,isNames)
        mXis <- cbind(mX,ISmatrices[[i]][,isNames])
        colnames(mXis) <- mXisNames
        mXis <- dropvar(mXis, tol=tol, LAPACK=LAPACK,
          silent=print.searchinfo)
        mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
          qstat.options=NULL, tol=tol, LAPACK=LAPACK,
          plot=FALSE)
        getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
          do.pet=do.pet, wald.pval=wald.pval, ar.LjungB=ar.LjungB,
          arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
          user.diagnostics=user.diagnostics, info.method=info.method,
          include.gum=include.gum, include.empty=include.empty,
          max.paths=max.paths, max.regs=max.regs,
          print.searchinfo=print.searchinfo,
          estimate.specific=FALSE, plot=FALSE)
        ISfinalmodels[[i]] <- names(getsis$specific.spec)
      }

##END MODIFIED (panel)

    } #end if(length(ISspecific.models > 0)

  } #end for(i) loop (on ISmatrices)

  ##add names to ISblocks:
  names(ISblocks) <- names(ISmatrices)

  ##gets of union of retained impulses:
  if(print.searchinfo){
    message("\n", appendLF=FALSE)
    message("GETS of union of ALL retained variables...",
      appendLF=TRUE)
    message("\n", appendLF=FALSE)
  }

  ##no final models estimated:
  if(length(ISfinalmodels)==0){
    ISfinalmodels <- NULL
    if(is.null(mX)){ mXis <- NULL }else{
      mXis <- zoo(cbind(mX), order.by=y.index)
      colnames(mXis) <- mXnames
    }
  }

  ##final models estimated:
  if(length(ISfinalmodels)>0){

    mIS <- NULL #matrix

    #which indicators were retained?
    for(i in 1:length(ISfinalmodels)){
      isNames <- NULL
      #check if non-empty:
      if(!is.null(ISfinalmodels[[i]])){
        isNames <- setdiff(ISfinalmodels[[i]], mXnames)
      }
      if(length(isNames)>0){
        tmp <- cbind(ISmatrices[[i]][, isNames ])
        colnames(tmp) <- isNames
        mIS <- cbind(mIS, tmp)
      }
    } #end for loop

    mXis <- dropvar(cbind(mX,mIS), tol=tol, LAPACK=LAPACK,
      silent=print.searchinfo)
    mXis <- zoo(mXis, order.by=y.index)
  } #end if(length(ISfinalmodels)>0)

  ##gum and gets:
  y <- zoo(y, order.by=y.index)
  mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
    qstat.options=qstat.options, user.diagnostics=user.diagnostics,
    tol=tol, LAPACK=LAPACK, plot=FALSE)
  getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
    do.pet=do.pet, wald.pval=wald.pval, ar.LjungB=ar.LjungB,
    arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
    user.diagnostics=user.diagnostics, info.method=info.method,
    include.empty=include.empty, max.paths=max.paths,
    max.regs=max.regs, print.searchinfo=print.searchinfo,
    plot=FALSE)

  ##names of retained impulses, mX colnames:
  ISnames <- setdiff(getsis$aux$mXnames, mXnames)
  if(length(ISnames)==0){ ISnames <- NULL }
  colnames(getsis$aux$mX) <- getsis$aux$mXnames

  ##return:
  getsis$gets.type <- "isat"
  getsis$call <- isat.call
  getsis <- c(list(ISfinalmodels=ISfinalmodels,
    ISnames=ISnames), getsis)
  getsis$aux$t.pval <- t.pval #needed for biascorr
  class(getsis) <- "isat"
  if(alarm){ alarm() }
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.isat(getsis, coef.path=TRUE) }
  return(getsis)
}
