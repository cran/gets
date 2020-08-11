getsv <-
function(object, t.pval=0.05, wald.pval=t.pval,
  do.pet=TRUE, ar.LjungB=list(lag=NULL, pval=0.025),
  arch.LjungB=list(lag=NULL, pval=0.025),
  normality.JarqueB=NULL, user.diagnostics=NULL,
  info.method=c("sc", "aic", "aicc", "hq"),
  gof.function=NULL, gof.method=NULL, keep=c(1),
  include.gum=FALSE, include.1cut=TRUE, include.empty=FALSE,
  max.paths=NULL, tol=1e-07, turbo=FALSE, print.searchinfo=TRUE,
  plot=NULL, alarm=FALSE)
{
  ### ARGUMENTS ###########

  ##obligatory:
  vc=TRUE
  vcov.type <- "ordinary"

  ##zoo and NA related:
  e <- object$residuals #should not contain NAs
  e.index <- index(e) #use object$aux$y.index instead?
  e <- coredata(e)
  e.n <- length(e) #use object$aux$y.n instead?
  eadj <- e[c(e.n-object$aux$loge2.n+1):e.n] #Note: log(eadj^2)=loge2
  eadj.n <- length(eadj)
  eadj.index <- e.index[c(e.n-object$aux$loge2.n+1):e.n]

  ##diagnostics options, max.regs:
  if(!is.null(ar.LjungB) && is.null(ar.LjungB$lag)){
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if(!is.null(arch.LjungB) && is.null(arch.LjungB$lag)){
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])
  #if(is.null(max.regs)){ max.regs <- 10*object$aux$y.n }

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


  ### INITIALISE ##########

  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  loge2 <- object$aux$loge2
  mX <- object$aux$vX
  colnames(mX) <- object$aux$vXnames
  if( !(1 %in% keep) ){
    keep <- union(1,keep)
    warning("Regressor 1 included into 'keep'")
  }

  ##add gum results and diagnostics to out:
  out$gum.mean <- object$mean.results
  tmp <- matrix(0, NROW(object$variance.results), 2)
  colnames(tmp) <- c("reg.no.", "keep")
  tmp[,1] <- 1:NROW(tmp) #fill reg.no. column
  tmp[keep,2] <- 1 #fill keep column
  out$gum.variance <- cbind(tmp, object$variance.results)
  out$gum.diagnostics <- object$diagnostics


  ### DO MULTI-PATH GETS ##########

  ##do the gets:
  est <- getsFun(loge2, mX,
    user.estimator=list(name="ols", untransformed.residuals=eadj,
    tol=object$aux$tol, LAPACK=object$aux$LAPACK, method=6),
    gum.result=NULL, t.pval=t.pval, wald.pval=wald.pval, do.pet=do.pet,
    ar.LjungB=ar.LjungB, arch.LjungB=arch.LjungB,
    normality.JarqueB=normality.JarqueB, user.diagnostics=user.diagnostics,
    gof.function=gof.function, gof.method=gof.method, keep=keep,
    include.gum=include.gum, include.1cut=include.1cut,
    include.empty=include.empty, max.paths=max.paths, turbo=turbo,
    tol=tol, max.regs=NULL, print.searchinfo=print.searchinfo,
    alarm=alarm)
  est$time.started <- NULL
  est$time.finished <- NULL
  est$call <- NULL
  out <- c(out, est)

  ## if no search has been undertaken:
  if(is.null(est$terminals.results)){
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }


  ### ESTIMATE SPECIFIC ################

  ## prepare estimation:
  e <- zoo(cbind(eadj), order.by=eadj.index)
  colnames(e) <- "e"
  specificadj <- setdiff(out$specific.spec, 1)
  if(length(specificadj)==0){
    vXadj <- NULL
  }else{
    vXadj <- cbind(object$aux$vX[,specificadj])
    colnames(vXadj) <- object$aux$vXnames[specificadj]
    vXadj <- zoo(vXadj, order.by=eadj.index)
  }
  if( is.null(ar.LjungB) ){ ar.LjungB <- object$aux$qstat.options[1] }
  if( is.null(arch.LjungB) ){ arch.LjungB <- object$aux$qstat.options[2] }
  if( is.null(normality.JarqueB) ){
    normality.JarqueB <- FALSE
  }else{
    normality.JarqueB <- TRUE
  }

  ## estimate model:
  est <- arx(e, vc=TRUE, vxreg=vXadj,
    zero.adj=object$aux$zero.adj, vc.adj=object$aux$vc.adj,
    qstat.options=c(ar.LjungB[1],arch.LjungB[1]),
    normality.JarqueB=normality.JarqueB,
    user.diagnostics=user.diagnostics, tol=object$aux$tol,
    LAPACK=object$aux$LAPACK, plot=FALSE)

  ## delete, rename and change various stuff:
  est$call <- est$date <- NULL
  where.diagnostics <- which(names(est)=="diagnostics")
  if(length(where.diagnostics)>0){
    names(est)[where.diagnostics] <- "specific.diagnostics"
  }
  est$mean.fit <- object$mean.fit[ index(object$mean.fit) %in% eadj.index ]
  #est$mean.fit <- object$mean.fit[ eadj.index ] #should work, but doesn't!
  est$vcov.mean <- NULL
  est$aux$vxreg <- est$aux$vxreg.index <- NULL
  est$aux$y.name <- "e"

  ## finalise:
  est <- unclass(est)
  out <- c(out,est)

  ### OUTPUT ########

  out$aux$vXnames.gum <- object$aux$vXnames
  out$aux$call.gum <- object$call
  if(is.null(out$aux$vcov.type)){ out$aux$vcov.type <- vcov.type }
  out <- c(list(date=date(), gets.type="getsv"), out)
  out$time.finished <- date()
  class(out) <- "gets"

  if(alarm){ alarm() }
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.gets(out) }
  return(out)
}
