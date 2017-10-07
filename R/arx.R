arx <-
function(y, mc=FALSE, ar=NULL, ewma=NULL, mxreg=NULL,
  vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
  zero.adj=0.1, vc.adj=TRUE,
  vcov.type=c("ordinary", "white", "newey-west"),
  qstat.options=NULL, user.estimator=NULL, user.diagnostics=NULL,
  tol=1e-07, LAPACK=FALSE, plot=NULL)
{
  ### ARGUMENTS: ###########

  vcov.type <- match.arg(vcov.type)

  ##regressand:
  y.name <- deparse(substitute(y))
  #if(is.ts(y)){ y <- as.zooreg(y) }
  if(is.zoo(y)){ y <- cbind(y) }else{ y <- as.zoo(cbind(y)) }
  #OLD: y <- as.zoo(cbind(y))
  y <- cbind(y)
  if(NCOL(y) > 1) stop("Dependent variable not 1-dimensional")
  if( is.null(y.name)){ y.name <- colnames(y)[1] }
  if( y.name[1] =="" ){ y.name <- "y" }
  y.n <- NROW(y)
  y.index <- index(y)
  y <- coredata(y)

  ##regressors:
  mX <- NULL
  mXnames <- NULL

  ##mean intercept:
  if(identical(as.numeric(mc),1)){
    mX <- cbind(rep(1,y.n))
    mXnames  <- "mconst"
  }

  ##ar terms:
  if(!is.null(ar) && !identical(as.numeric(ar),0) ){
    tmp <- NULL
    nas <- rep(NA, max(ar))
    tmpfun <- function(i){
      tmp <<- cbind(tmp, c(nas[1:i],y[1:c(y.n-i)]))
    }
    tmpfun <- sapply(ar,tmpfun)
    mX <- cbind(mX, tmp)
    mXnames <- c(mXnames, paste("ar", ar, sep=""))
  }

  ##ewma term:
  if(!is.null(ewma)){
    tmp <- do.call(eqwma, c(list(y),ewma) )
    mXnames <- c(mXnames, colnames(tmp))
    colnames(tmp) <- NULL
    mX <- cbind(mX, tmp)
  }

  ##adjust for NAs:
  tmp <- zoo(cbind(y,mX), order.by=y.index)
  tmp <- na.trim(tmp, sides="both", is.na="any")
  y <- tmp[,1]
  y.n <- NROW(y) #re-define y.n
  y.index <- index(y) #re-define y.index
  t1 <- y.index[1]
  t2 <- y.index[y.n]
  y <- coredata(y)
  if(!is.null(mX)){
    mX <- tmp[,2:NCOL(tmp)]
    mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL
  }

  ##mxreg:
  if(!is.null(mxreg) && !identical(as.numeric(mxreg),0) ){
    #if(is.ts(mxreg)){ mxreg <- as.zooreg(mxreg) }
    mxreg <- as.zoo(cbind(mxreg))
    mxreg.names <- colnames(mxreg)
    if(is.null(mxreg.names)){
      mxreg.names <- paste("mxreg", 1:NCOL(mxreg), sep="")
    }
    if(any(mxreg.names == "")){
      missing.colnames <- which(mxreg.names == "")
      for(i in 1:length(missing.colnames)){
        mxreg.names[i] <- paste("mxreg", i, sep="")
      }
    }
    #mxreg.names <- make.names(mxreg.names)
    mXnames <- c(mXnames, mxreg.names)
    mxreg <- window(mxreg, start=t1, end=t2)
    mxreg <- cbind(coredata(mxreg))
    mX <- cbind(mX, mxreg)

    ##re-adjust for NAs:
    tmp <- zoo(cbind(y,mX), order.by=y.index)
    tmp <- na.trim(tmp, sides="both", is.na="any")
    y <- tmp[,1]
    y.n <- NROW(y) #re-define y.n
    y.index <- index(y) #re-define y.index
    t1 <- y.index[1]
    t2 <- y.index[y.n]
    y <- coredata(y)
    mX <- tmp[,2:NCOL(tmp)]
    mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL

  } #end if(!is.null(mxreg))

  ##vxreg:
  if(!is.null(vxreg) && !identical(as.numeric(vxreg),0) ){
    vxreg <- as.zoo(cbind(vxreg))
    vxreg.names <- colnames(vxreg)
    if(is.null(vxreg.names)){
      vxreg.names <- paste("vxreg", 1:NCOL(vxreg), sep="")
    }
    if(any(vxreg.names == "")){
      missing.colnames <- which(vxreg.names == "")
      for(i in 1:length(missing.colnames)){
        vxreg.names[i] <- paste("vxreg", i, sep="")
      }
    }
    vxreg <- window(vxreg, start=t1, end=t2)
#OLD:
#    vxreg <- cbind(coredata(vxreg))
    colnames(vxreg) <- NULL
  } #end if(!is.null(vxreg))

  ##determine qstat.options:
  if(is.null(qstat.options)){
    if(is.null(ar)){ar.lag <- 1}else{ar.lag <- max(ar)+1}
    if(is.null(arch)){arch.lag <- 1}else{arch.lag <- max(arch)+1}
    qstat.options <- c(ar.lag, arch.lag)
  }

  ##aux: info for getsm/getsv functions
  aux <- list()
  aux$y <- y
  aux$y.index <- y.index
  aux$y.name <- y.name
  aux$y.n <- y.n
  if(!is.null(mX)){
    colnames(mX) <- NULL
    aux$mX <- mX
    aux$mXnames <- mXnames
    aux$mXncol <- NCOL(mX)
  }
  aux$vc <- vc
  aux$zero.adj <- zero.adj
  aux$vc.adj <- vc.adj
  aux$vcov.type <- vcov.type
  aux$qstat.options <- qstat.options
  aux$user.estimator <- user.estimator
  aux$user.diagnostics <- user.diagnostics
  aux$tol <- tol
  aux$LAPACK <- LAPACK

  ### INITIALISE ##########

  sysCall <- sys.call()
  #for the future: make the following objects part of the out-list?
  vcov.var <- NULL #make sure this object exists
  variance.results <- NULL #make sure this object exists

  #### for the future regarding user.estimator: check if
  #### user.estimator$spec is NULL, "mean", "variance" or "both"
  #### in order to determine what kind of estimator it is

  ##check if mean and log-garch spec:
  meanSpec <- !is.null(mX)
  noVarianceSpec <- if( vc==FALSE && is.null(arch)
    && is.null(asym) && is.null(log.ewma)
    && is.null(vxreg) ){ TRUE }else{ FALSE }

  #### MEAN ###############

  ##estimate:
  if( is.null(user.estimator) ){

    estMethod <- which(vcov.type==c("none", "none", "ordinary",
      "white", "newey-west"))
    out <- ols(y, mX, tol=tol, LAPACK=LAPACK,
      method=estMethod, user.fun=NULL, user.options=NULL)

    ##delete unneeded entries:
    #out$n <- NULL #this might have to be changed in order to enable gum.result in getsFun
    #out$k <- NULL ##this might have to be changed in order to enable gum.result in getsFun
    #out$df <- NULL: Do not delete!
    out$qr <- NULL
    out$rank <- NULL
    out$qraux <- NULL
    out$pivot <- NULL
    out$xtxinv <- NULL
    out$residuals2 <- NULL
    #out$rss <- NULL

    ##add entries if no variance spec:
    if(noVarianceSpec){
      out$var.fit <- rep(out$sigma2, aux$y.n)
      out$std.residuals <- out$residuals/sqrt(out$sigma2)
      aux$loge2.n <- aux$y.n #change to out$n?
    }

  }else{

    ##user-defined estimator:
    out <- do.call(user.estimator$name, list(y, mX),
      envir=.GlobalEnv)
    #delete?:
    if( is.null(out$vcov) && !is.null(out$vcov.mean) ){
      out$vcov <- out$vcov.mean
    }

  } #end if(..)else(..)

  ##make estimation results (a data frame):
  if(meanSpec){
    stderrs <- sqrt(diag(out$vcov))
    colnames(out$vcov) <- mXnames
    rownames(out$vcov) <- mXnames
    t.stat <- out$coefficients/stderrs
    p.val <- pt(abs(t.stat), out$df, lower.tail=FALSE)*2
    out$mean.results <- as.data.frame(cbind(out$coefficients,
      stderrs, t.stat, p.val))
    colnames(out$mean.results) <- c("coef", "std.error",
      "t-stat", "p-value")
    rownames(out$mean.results) <- mXnames
  } #end if(meanSpec)

  ##rename some stuff:
  if( is.null(user.estimator) ){
    outNames <- names(out)
    whereIs <- which(outNames=="vcov")
    if( length(whereIs) > 0 ){ names(out)[whereIs] <- "vcov.mean" }
    whereIs <- which(outNames=="fit")
    names(out)[whereIs] <- "mean.fit"
  }


  #### VARIANCE #############

  ##if log-arch spec:
  if( noVarianceSpec==FALSE ){

    ##regressand
    zero.where <- which(out$residuals==0)
    eabs <- abs(out$residuals)
    if(length(zero.where) > 0){
      eabs[zero.where] <- quantile(eabs[-zero.where], zero.adj)
    }
    loge2 <- log(eabs^2)

    ##regressor matrix:
    vX <- cbind(rep(1,y.n))
    vXnames <- "vconst"

    ##arch terms:
    if(!is.null(arch) && !identical(as.numeric(arch),0) ){
      tmp <- NULL
      nas <- rep(NA, max(arch))
      tmpfun <- function(i){
        tmp <<- cbind(tmp, c(nas[1:i],loge2[1:c(y.n-i)]))
      }
      tmpfun <- sapply(arch,tmpfun)
      vX <- cbind(vX, tmp)
      vXnames <- c(vXnames, paste("arch", arch, sep=""))
    }

    ##asym terms:
    if(!is.null(asym) && !identical(as.numeric(asym),0) ){
      tmp <- NULL
      nas <- rep(NA, max(asym))
      tmpfun <- function(i){
        tmp <<- cbind(tmp, c(nas[1:i],
          loge2[1:c(y.n-i)]*as.numeric(out$residuals[1:c(y.n-i)]<0)))
      }
      tmpfun <- sapply(asym,tmpfun)
      vX <- cbind(vX, tmp)
      vXnames <- c(vXnames, paste("asym", asym, sep=""))
    }

    ##log.ewma term:
    if(!is.null(log.ewma)){
      if(is.list(log.ewma)){
        log.ewma$lag <- 1
      }else{
        log.ewma <- list(length=log.ewma)
      }
      tmp <- do.call(leqwma, c(list(out$residuals),log.ewma) )
      vXnames <- c(vXnames, colnames(tmp))
      colnames(tmp) <- NULL
      vX <- cbind(vX, tmp)
    }

    ##adjust for NAs:
    tmp <- zoo(cbind(loge2,vX), order.by=y.index)
    tmp <- na.trim(tmp, sides="left", is.na="any")
    loge2 <- tmp[,1]
    loge2.n <- NROW(loge2)
    loge2.index <- index(loge2) #re-define y.index
    loge2 <- coredata(loge2)
    vX <- tmp[,2:NCOL(tmp)]
    vX <- coredata(vX)
    vX <- cbind(vX)
    colnames(vX) <- NULL

    ##vxreg:
    if( !is.null(vxreg) && !identical(as.numeric(vxreg),0) ){
      vxreg <- window(vxreg, start=loge2.index[1],
        end=loge2.index[loge2.n])
      vxreg <- cbind(vxreg)
      vX <- cbind(vX, coredata(vxreg))
      vXnames <- c(vXnames, vxreg.names)
      colnames(vxreg) <- vxreg.names

      ##re-adjust for NAs:
      tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
      tmp <- na.trim(tmp, sides="left", is.na="any")
      loge2 <- tmp[,1]
      loge2.n <- NROW(loge2)
      loge2.index <- index(loge2) #re-define index
      loge2 <- coredata(loge2)
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }

    ##aux: more info for getsm/getsv functions
    aux$loge2 <- loge2
    aux$loge2.n <- loge2.n
    aux$vX <- vX
    aux$vXnames <- vXnames
    aux$vXncol <- NCOL(vX)
    aux$arch <- arch
    aux$asym <- asym
    aux$log.ewma <- log.ewma
    aux$vxreg <- vxreg #note: NROW(vxreg)!=NROW(vX) is possible

    ##estimate, prepare results:
    estVar <- ols(loge2, vX, tol=tol, LAPACK=LAPACK,
      method=3)
    s.e. <- sqrt(as.vector(diag(estVar$vcov)))
    estVar$vcov <- as.matrix(estVar$vcov[-1,-1])
    colnames(estVar$vcov) <- vXnames[-1]
    rownames(estVar$vcov) <- vXnames[-1]
    Elnz2 <- -log(mean(exp(estVar$residuals)))
    t.stat <- estVar$coefficients/s.e.
    p.val <- pt(abs(t.stat), estVar$df, lower.tail=FALSE)*2
    if(vc.adj){
      t.stat[1] <- ((estVar$coefficients[1]-Elnz2)^2)/s.e.[1]^2
      p.val[1] <- pchisq(t.stat[1], 1, lower.tail=FALSE)
      estVar$coefficients[1] <- estVar$coefficients[1] - Elnz2
    }

    ##add/modify entries in out:
    out$n <- estVar$n
    out$vcov.var <- estVar$vcov
    out$var.fit <- exp(estVar$fit - Elnz2)
    out$ustar.residuals <- estVar$residuals
    out$std.residuals <- out$residuals[c(y.n-loge2.n+1):y.n]/sqrt(out$var.fit)
    out$logl <- -loge2.n*log(2*pi)/2 - sum(log(out$var.fit))/2 - sum(out$std.residuals^2)/2
    out$Elnz2 <- Elnz2
    out$variance.results <- as.data.frame(cbind(estVar$coefficients, s.e., t.stat, p.val))
    colnames(out$variance.results) <- c("coef", "std.error", "t-stat", "p-value")
    rownames(out$variance.results) <- vXnames

  } #end if( noVarianceSpec==FALSE )


  ### OUTPUT: ######################

  ##diagnostics:
  out$diagnostics <- diagnostics(out,
    ar.LjungB=c(qstat.options[1],0), arch.LjungB=c(qstat.options[2],0),
      normality.JarqueB=0, user.fun=user.diagnostics, verbose=TRUE)

  ##add NAs to variance series:
  if( noVarianceSpec==FALSE ){
    NAs2add <- rep(NA,y.n-loge2.n)
    out$var.fit <- c(NAs2add, out$var.fit)
    out$ustar.residuals <- c(NAs2add, out$ustar.residuals)
    out$std.residuals <- c(NAs2add, out$std.residuals)
  }

  ##add zoo-indices:
  out$mean.fit <- zoo(out$mean.fit, order.by=y.index)
  out$residuals <- zoo(out$residuals, order.by=y.index)
  if(!is.null(out$var.fit)){
    out$var.fit <- zoo(out$var.fit, order.by=y.index)
  }
  if(!is.null(out$ustar.residuals)){
    out$ustar.residuals <- zoo(out$ustar.residuals, order.by=y.index)
  }
  if(!is.null(out$std.residuals)){
    out$std.residuals <- zoo(out$std.residuals, order.by=y.index)
  }

  ##result:
  out <- c(list(call=sysCall, date=date(), aux=aux), out)
  class(out) <- "arx"

  ##plot:
  if( is.null(plot) ){
    plot <- getOption("plot")
    if( is.null(plot) ){ plot <- FALSE }
  }
  if(plot){ plot.arx(out) }

  ##return result:
  return(out)
}
