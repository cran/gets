getsm <-
function(object, t.pval=0.05, wald.pval=t.pval, vcov.type=NULL,
  do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025), normality.JarqueB=NULL,
  user.diagnostics=NULL, info.method=c("sc","aic","aicc","hq"),
  gof.function=NULL, gof.method=NULL, keep=NULL, include.gum=FALSE,
  include.1cut=TRUE, include.empty=FALSE, max.paths=NULL,
  turbo=FALSE, print.searchinfo=TRUE, plot=NULL, alarm=FALSE)
{
  ## contents:
  ## 1 arguments
  ## 2 gets modelling
  ## 3 estimate specific
  ## 4 output
  
  ##------------------
  ## 1 arguments
  ##------------------
  
  ##check if mean equation:
  if( is.null(object$aux$mX) ){ stop("Mean equation empty") }

  ##check max.paths:
  if( !is.null(max.paths) && max.paths < 1){
    stop("'max.paths' cannot be smaller than 1")
  }

  ##diagnostics: determine ar and arch lags:
  if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])

  ##user-defined diagnostics?:
  if( is.null(user.diagnostics) ){
    user.diagnostics <- object$call$user.diagnostics  
  }
  
  ##if( user-defined estimator ):
  if( !is.null(object$aux$user.estimator) ){
    user.estimator <- object$call$user.estimator
    if( is.null(plot) || identical(plot,TRUE) ){
      plot <- FALSE
      message("  user-defined estimator: 'plot' set to FALSE")
    }
  } #close if user estimator

  ##if( default estimator ):
  if( is.null(object$call$user.estimator) ){

    ##determine ols method:
    if( is.null(vcov.type) ){ vcov.type <- object$aux$vcov.type }
    vcovTypes <- c("a", "b", "ordinary", "white", "newey-west")
    olsMethod <- charmatch(vcov.type, vcovTypes)
    if( (olsMethod%in%c(3,4,5))==FALSE ){ stop("'vcov.type' invalid") }
    
    ##ols arguments:
    user.estimator <- list()
    user.estimator$name <- "ols"
    user.estimator$tol <- object$aux$tol 
    user.estimator$LAPACK <- object$aux$LAPACK
    user.estimator$method <- olsMethod
    
    ##variance specification:
    if( is.null(object$variance.result) ){
      user.estimator$variance.spec <- NULL
    }else{
      user.estimator$variance.spec <- list(vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg)
    }
    
  } #close if( default estimator )

  ##gof arguments:
  if( is.null(gof.function) ){

    ##determine info method:
    infoTypes <- c("sc","aic","aicc","hq")
    whichMethod <- charmatch(info.method[1], infoTypes)
    info.method <- infoTypes[ whichMethod ]
    
    ##make gof arguments:
    gof.function <- list(name="infocrit", method=info.method)
    gof.method <- "min"
    
  }
  
  ##------------------
  ## 2 gets modelling
  ##------------------

  ##out list:
  out <- list()
  out$time.started <- date()
  out$time.finished <- NA ##added below, towards the end
  out$call <- sys.call() #used by coef.arx

  ##add gum results and diagnostics to out:
  tmp <- matrix(0, NROW(object$mean.results), 2)
  colnames(tmp) <- c("reg.no.", "keep")
  tmp[,1] <- 1:NROW(tmp) #fill reg.no. column
  tmp[keep,2] <- 1 #fill keep column
  out$gum.mean <- cbind(tmp, object$mean.results)
  out$gum.variance <- object$variance.results
  out$gum.diagnostics <- object$diagnostics

  ##do the gets:
  est <- getsFun(object$aux$y, object$aux$mX,
    user.estimator=user.estimator, gum.result=NULL, t.pval=t.pval,
    wald.pval=wald.pval, do.pet=do.pet, ar.LjungB=ar.LjungB,
    arch.LjungB=arch.LjungB, normality.JarqueB=normality.JarqueB,
    user.diagnostics=user.diagnostics, gof.function=gof.function,
    gof.method=gof.method, keep=keep, include.gum=include.gum,
    include.1cut=include.1cut, include.empty=include.empty,
    max.paths=max.paths, turbo=turbo, max.regs=NULL,
    print.searchinfo=print.searchinfo, alarm=alarm)
  est$time.started <- NULL
  est$time.finished <- NULL
  out$time.finished <- date()
  est$call <- NULL
  out <- c(out, est)

  ##---------------------
  ## 3 estimate specific
  ##---------------------

  ## if no search has been undertaken:
  if( is.null(out$terminals.results) ){
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }

  ##if search has been undertaken:
  if( !is.null(out$terminals.results) ){

    if( length(out$specific.spec)>0 ){
      out$specific.spec <- sort(out$specific.spec)
    }

    ##prepare estimation:
    yadj <- zoo(object$aux$y, order.by=object$aux$y.index)
    if( length(out$specific.spec)==0 ){
      mXadj <- NULL
    }else{
      mXadj <- cbind(object$aux$mX[, out$specific.spec ])
      colnames(mXadj) <- object$aux$mXnames[ out$specific.spec ]
      mXadj <- zoo(mXadj, order.by=object$aux$y.index)
    }
    if(is.null(ar.LjungB)){ ar.LjungB <- object$aux$qstat.options[1] }
    if(is.null(arch.LjungB)){ arch.LjungB <- object$aux$qstat.options[2] }
    if( is.null(normality.JarqueB) ){
      normality.JarqueB <- FALSE
    }else{
      normality.JarqueB <- TRUE
    }
        
    ##if( default estimator ):
    if( is.null(object$call$user.estimator) ){
      ##estimate specific model:
      est <- arx(yadj, mxreg=mXadj, vc=object$aux$vc,
        arch=object$aux$arch, asym=object$aux$asym,
        log.ewma=object$aux$log.ewma, vxreg=object$aux$vxreg,
        zero.adj=object$aux$zero.adj,
        vc.adj=object$aux$vc.adj, vcov.type=vcov.type,
        qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
        normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, plot=FALSE)
    } #end if( default estimator )

    ##if( user-defined estimator ):
    if( !is.null(object$call$user.estimator) ){
      ##estimate specific:
      est <- arx(yadj, mxreg=mXadj, user.estimator=user.estimator,
        qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
        normality.JarqueB=normality.JarqueB,
        user.diagnostics=user.diagnostics, tol=object$aux$tol,
        LAPACK=object$aux$LAPACK, plot=FALSE)
    } #end if( user-defined estimator )

    ##delete, rename, add:
    est$call <- est$date <- NULL
    where.diagnostics <- which(names(est)=="diagnostics")
    if(length(where.diagnostics)>0){
      names(est)[where.diagnostics] <- "specific.diagnostics"
    }
    est$aux$y.name <- object$aux$y.name
    est$aux$call.gum <- object$call #used by predict.gets()
    est <- unclass(est)
    out <- c(out,est)

  } #end if( !is.null(out$terminals.results) )

  ##------------------
  ## 4 output
  ##------------------

  ##finalise and return:
  out <- c(list(date=date(), gets.type="getsm"), out)
  class(out) <- "gets"
  if(alarm){ alarm() }
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.gets(out) }
  return(out)

}
