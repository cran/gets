isat <-
function(y, mc=TRUE, ar=NULL, ewma=NULL, mxreg=NULL,
  iis=TRUE, sis=FALSE, blocks=NULL, ratio.threshold=0.8,
  max.block.size=30, vcov.type=c("ordinary", "white"),
  t.pval=0.001, do.pet=FALSE, wald.pval=0.001, ar.LjungB=NULL,
  arch.LjungB=NULL, normality.JarqueB=NULL,
  info.method=c("sc", "aic", "hq"), include.gum=FALSE,
  include.empty=FALSE, tol=1e-07, LAPACK=FALSE, max.regs=NULL,
  verbose=TRUE, print.searchinfo=TRUE, alarm=FALSE, plot=TRUE)
{

  ##arguments:
  tis=FALSE #for the future!
  isat.call <- sys.call()
  vcov.type <- match.arg(vcov.type)
  info.method <- match.arg(info.method)
  mod <- arx(y, mc=mc, ar=ar, ewma=ewma, mxreg=mxreg,
    vcov.type=vcov.type, qstat.options=NULL,
    tol=tol, LAPACK=LAPACK, verbose=FALSE, plot=FALSE)
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

  #indicator saturation matrices:
  listIS <- list()
  if(is.null(mX)){ mxkeep <- NULL}else{ mxkeep <- 1:mXncol }

  if(iis){
    mIIS <- matrix(0,y.n,y.n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste("iis", y.index.as.char, sep="")
    listIS <- c(listIS,list(mIIS=mIIS))
  }else{
    mIIS <- NULL
  }

  if(sis){
    mSIS <-matrix(0,y.n,y.n)
    loop.indx <- 1:y.n
    tmp <- function(i){ mSIS[i,1:i] <<- 1 }
    tmp <- sapply(loop.indx,tmp)
    colnames(mSIS) <- paste("sis", y.index.as.char, sep="")
    mSIS <- mSIS[,-1]
    listIS <- c(listIS,list(mSIS=mSIS))
  }else{
    mSIS <- NULL
  }

  if(tis){
#NOTE: This part has to be changed and adapted for TIS labels!
    mTIS <- NULL
    v1n <- seq(1,y.n)
    zeros <- rep(0,y.n)
    loop.indx <- 1:c(y.n-1)
    tmp <- function(i){
      mTIS <<- cbind(mTIS,c(zeros[1:i],v1n[1:c(y.n-i)]))
    }
    tmp <- sapply(loop.indx,tmp)
    colnames(mTIS) <- paste("tis", 2:c(NCOL(mTIS)+1), sep="")
    listIS <- c(listIS,list(mTIS=mTIS))
  }else{
    mTIS <- NULL
  }

  #determine no. of blocks:
  if(is.null(blocks)){
    blockratio.value <- y.n/(ratio.threshold*y.n - NCOL(mX))
    blocksize.value <- y.n/max.block.size
    blocks <- max(2,blockratio.value,blocksize.value)
    blocks <- ceiling(blocks)
  }

  #loop on listIS (i.e. mIIS, mSIS, mTIS):
  ISfinalmodels <- list()
  for(i in 1:length(listIS)){

    #make partitions:
    smplsize <- floor(y.n/blocks)
    partitions.t2 <- rep(NA,blocks)
    partitions.t2[blocks] <- NCOL(listIS[[i]]) #y.n
    for(j in 1:c(blocks-1)){ partitions.t2[j] <- smplsize*j }
    partitions.t1 <- partitions.t2 + 1
    partitions.t1 <- c(1,partitions.t1[-blocks])

    #gets on each block:
    ISspecific.models <- list()
#    ISgums <- list()
#    ISpaths <- list()
#    ISterminals.results <- list()
    for(j in 1:blocks){

      #print info:
      if(print.searchinfo){
        cat("\n")
        cat(substr(names(listIS)[i],2,4),
          " block ", j, " of ", blocks, "...\n", sep="")
        cat("\n")
      }

      mXis <- cbind(mX,listIS[[i]][,partitions.t1[j]:partitions.t2[j]])
      mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
        qstat.options=qstat.options, tol=tol, LAPACK=LAPACK,
        verbose=TRUE, plot=FALSE)
      getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
        do.pet=do.pet, wald.pval=wald.pval, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        info.method=info.method, include.empty=include.empty,
        max.regs=max.regs, estimate.specific=FALSE,
        print.searchinfo=print.searchinfo, plot=FALSE)

      if(is.null(getsis$specific.spec)){
        ISspecific.models[[j]] <- NULL
      }else{
        ISspecific.models[[j]] <- names(getsis$specific.spec)
#        ISgums[[j]] <- getsis$gum.mean
#        ISpaths[[j]] <- getsis$paths
#        ISterminals.results[[j]] <- getsis$terminals.results
      }

    } #end for(j in 1:blocks)

    #print info:
    if(print.searchinfo){
      cat("\n")
      cat("GETS of union of retained ",
        substr(names(listIS)[i],2,4), " indicators... \n",
        sep="")
      cat("\n")
    }

    #if no indicators retained from the blocks:
    if(length(ISspecific.models)==0){
      isNames <- NULL
      ISfinalmodels[[i]] <- NULL
    }

    #when indicators retained from the blocks:
    if(length(ISspecific.models)>0){

      isNames <- NULL

      #which indicators retained?:
      for(j in 1:length(ISspecific.models)){
        #check if mean is non-empty:
        if(!is.null(ISspecific.models[[j]])){
          isNames <- union(isNames, ISspecific.models[[j]])
        }
      } #end for(j) loop
      isNames <- setdiff(isNames, mXnames)
      #isNamesAll[[i]] <- isNames

      #redo gets with union of retained indicators:
      mXisNames <- c(mXnames,isNames)
      mXis <- cbind(mX,listIS[[i]][,isNames])
      colnames(mXis) <- mXisNames
      mXis <- dropvar(mXis)
      mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
        qstat.options=NULL, tol=tol, LAPACK=LAPACK,
        verbose=TRUE, plot=FALSE)
      getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
        do.pet=do.pet, wald.pval=wald.pval, ar.LjungB=ar.LjungB,
        arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
        info.method=info.method, include.gum=include.gum,
        include.empty=include.empty, max.regs=max.regs,
        print.searchinfo=print.searchinfo, estimate.specific=FALSE,
        plot=FALSE)
      ISfinalmodels[[i]] <- names(getsis$specific.spec)

    } #end if(length(IS..)>0

  } #end for(i) loop (on listIS)

  #gets of union of retained impulses:
  if(print.searchinfo){
    cat("\n")
    cat("GETS of union of ALL retained indicators...\n")
    cat("\n")
  }

  #no final models estimated:
  if(length(ISfinalmodels)==0){
    ISfinalmodels <- NULL
    if(is.null(mX)){ mXis <- NULL }else{
      mXis <- zoo(cbind(mX), order.by=y.index)
      colnames(mXis) <- mXnames
    }
  }

  #final models estimated:
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
        tmp <- cbind(listIS[[i]][, isNames ])
        colnames(tmp) <- isNames
        mIS <- cbind(mIS, tmp)
      }
    } #end for loop

    mXis <- dropvar(cbind(mX,mIS))
    mXis <- zoo(mXis, order.by=y.index)
  } #end if(length(ISfinalmodels)>0)

  #gum and gets:
  y <- zoo(y, order.by=y.index)
  mod <- arx(y, mxreg=mXis, vcov.type=vcov.type,
    qstat.options=qstat.options, tol=tol,
    LAPACK=LAPACK, verbose=TRUE, plot=FALSE)
  getsis <- getsm(mod, keep=mxkeep, t.pval=t.pval,
    do.pet=do.pet, wald.pval=wald.pval, ar.LjungB=ar.LjungB,
    arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
    info.method=info.method, include.empty=include.empty,
    max.regs=max.regs, print.searchinfo=print.searchinfo,
    plot=FALSE)

  #names of retained impulses:
  ISnames <- setdiff(getsis$aux$mXnames, mXnames)
  if(length(ISnames)==0){ ISnames <- NULL }

  #return:
  getsis$gets.type <- "isat"
  getsis$call <- isat.call
  getsis <- c(list(ISfinalmodels=ISfinalmodels,
    ISnames=ISnames), getsis)
  class(getsis) <- "gets"
  if(alarm){ alarm() }
  if(plot){ plot.gets(getsis, coef.path=TRUE) }
  return(getsis)
}
