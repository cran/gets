arx <-
function(y, mc=FALSE, ar=NULL, ewma=NULL, mxreg=NULL,
  vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
  zero.adj=0.1, vc.adj=TRUE,
  vcov.type=c("ordinary", "white", "newey-west"),
  qstat.options=NULL, user.estimator=NULL, user.diagnostics=NULL,
  tol=1e-07, LAPACK=FALSE, plot=NULL)
{
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

  out <- list()
  out$call <- sys.call()
#for the future: make the following objects part of the out-list?
  vcov.var <- NULL #make sure this object exists
  variance.results <- NULL #make sure this object exists

  #### for the future regarding user.estimator: check if
  #### user.estimator$spec is NULL, "mean", "variance" or "both"
  #### in order to determine what kind of estimator it is

  ##check if mean and/or log-garch spec:
  noMeanSpec <- is.null(mX)
  noVarianceSpec <- if( vc==FALSE && is.null(arch)
    && is.null(asym) && is.null(log.ewma)
    && is.null(vxreg) ){ TRUE }else{ FALSE }

  #### MEAN ###############

  if( noMeanSpec ){

    estMean <- list()
    estMean$resids <- aux$y
    estMean$mean.fit <- rep(0, aux$y.n)
    if(noVarianceSpec){
      estMean$sigma2 <- var(estMean$resids)
      estMean$var.fit <- rep(estMean$sigma2, aux$y.n)
      estMean$resids.std <- estMean$resids/sqrt( estMean$sigma2 )
      estMean$logl <- -aux$y.n*log(estMean$sigma2*2*pi)/2 - sum(estMean$resids.std^2)/2
    }

  }else{

    ##estimate:
    if( is.null(user.estimator) ){

      estMethod <- which(vcov.type==c("none", "none", "ordinary",
        "white", "newey-west"))
      estMean <- ols(y, mX, tol=tol, LAPACK=LAPACK,
        method=estMethod, user.fun=NULL, user.options=NULL)

      ##delete unneeded entries:
      estMean$qr <- NULL
      estMean$rank <- NULL
      estMean$qraux <- NULL
      estMean$pivot <- NULL
      estMean$xtxinv <- NULL
      estMean$resids2 <- NULL
      estMean$rss <- NULL
      estMean$n <- NULL

      ##add entries if no variance spec:
      if(noVarianceSpec){
        estMean$var.fit <- rep(estMean$sigma2, aux$y.n)
        estMean$resids.std <- estMean$resids/sqrt(estMean$sigma2)
        estMean$logl <- -aux$y.n*log(estMean$sigma2*2*pi)/2 - sum(estMean$resids.std^2)/2
        aux$loge2.n <- aux$y.n
      }

    }else{

      ##user-defined estimator:
      estMean <- do.call(user.estimator$name, list(y, mX),
        envir=.GlobalEnv)
      if( is.null(estMean$vcov) && !is.null(estMean$vcov.mean) ){
        estMean$vcov <- estMean$vcov.mean
      }

    }

    stderrs <- sqrt(diag(estMean$vcov))
    colnames(estMean$vcov) <- mXnames
    rownames(estMean$vcov) <- mXnames
    t.stat <- estMean$coefficients/stderrs
    p.val <- pt(abs(t.stat), estMean$df, lower.tail=FALSE)*2

    estMean$mean.results <- as.data.frame(cbind(estMean$coefficients,
      stderrs, t.stat, p.val))
    colnames(estMean$mean.results) <- c("coef", "std.error",
      "t-stat", "p-value")
    rownames(estMean$mean.results) <- mXnames

    ##rename some stuff:
    if( is.null(user.estimator) ){
      estMeanNames <- names(estMean)
      whereIs <- which(estMeanNames=="vcov")
      names(estMean)[whereIs] <- "vcov.mean"
      whereIs <- which(estMeanNames=="fit")
      names(estMean)[whereIs] <- "mean.fit"
    }

  } #end if(is.null(mX))else(..)


  #### VARIANCE #############

  ##if log-arch spec:
  if( noVarianceSpec==FALSE ){

    ##regressand
    zero.where <- which(estMean$resids==0)
    eabs <- abs(estMean$resids)
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
          loge2[1:c(y.n-i)]*as.numeric(estMean$resids[1:c(y.n-i)]<0)))
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
      tmp <- do.call(leqwma, c(list(estMean$resids),log.ewma) )
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

    ##estimate:
    est.v <- ols(loge2, vX, tol=tol, LAPACK=LAPACK,
      method=2)
    fit.v <- as.vector(vX%*%cbind(est.v$coefficients))
    ustar <- loge2-fit.v
    ustar2 <- ustar^2
    d.f.v. <- loge2.n - aux$vXncol
    sigma2.v <- sum(ustar2)/d.f.v.

    ##covariance coefficient matrix:
    varcovmat.v <- sigma2.v*est.v$xtxinv
    s.e. <- sqrt(as.vector(diag(varcovmat.v)))
    vcov.var <- as.matrix(varcovmat.v[-1,-1])
    colnames(vcov.var) <- vXnames[-1]
    rownames(vcov.var) <- vXnames[-1]
    t.stat <- est.v$coefficients/s.e.
    p.val <- pt(abs(t.stat), d.f.v., lower.tail=FALSE)*2

    Elnz2 <- -log(mean(exp(ustar)))
    if(vc.adj){
      t.stat[1] <- ((est.v$coefficients[1]-Elnz2)^2)/s.e.[1]^2
      p.val[1] <- pchisq(t.stat[1], 1, lower.tail=FALSE)
      est.v$coefficients[1] <- est.v$coefficients[1] - Elnz2
    }
    fit.v <- exp(fit.v - Elnz2)
    resids.std <- estMean$resids[c(y.n-loge2.n+1):y.n]/sqrt(fit.v)

    variance.results <- as.data.frame(cbind(est.v$coefficients, s.e., t.stat, p.val))
    colnames(variance.results) <- c("coef", "std.error", "t-stat", "p-value")
    rownames(variance.results) <- vXnames
  } #end if(noVarianceSpec)

  ### OUTPUT: ######################

  ##add zoo-indices:
  if(!is.null(estMean$resids)){
    estMean$resids <- zoo(estMean$resids, order.by=y.index)
  }
  if(!is.null(estMean$resids.std)){
    estMean$resids.std <- zoo(estMean$resids.std, order.by=y.index)
  }
  if(!is.null(estMean$mean.fit)){
    estMean$mean.fit <- zoo(estMean$mean.fit, order.by=y.index)
  }
  if(!is.null(estMean$var.fit)){
    estMean$var.fit <- zoo(estMean$var.fit, order.by=y.index)
  }

  out <- c(out, estMean)
  #rm(estMean)? Clean up?

  if( noVarianceSpec==FALSE ){

    ##log-variance
    add.nas2var <- rep(NA,y.n-loge2.n)
    out$var.fit <- zoo(c(add.nas2var,fit.v),
      order.by=y.index)
    out$resids.ustar <- zoo(c(add.nas2var,ustar),
        order.by=y.index)
    out$resids.std <- zoo(c(add.nas2var,resids.std),
      order.by=y.index)
    out$Elnz2 <- Elnz2
#check?, since it is made up of the standardised residuals?
    out$logl <- -loge2.n*log(2*pi)/2 - sum(log(fit.v))/2 - sum(na.trim(out$resids.std)^2)/2

  } #end if( varianceSpec )

  ##diagnostics:
  if( !is.null(out$resids.std) ){
    tmpY <- tmpXreg <- NULL
    if(!is.null(user.diagnostics)){
      tmpY <- out$aux$y
      tmpXreg <- out$aux$mX
    }
    out$diagnostics <- diagnostics( coredata(na.trim(out$resids.std, sides="both", is.na="any")),
      s2=1, y=tmpY, xreg=tmpXreg, ar.LjungB=c(qstat.options[1],0),
      arch.LjungB=c(qstat.options[2],0), normality.JarqueB=0,
      user.fun=user.diagnostics, verbose=TRUE)
  }

  ##result:
  out$vcov.var <- vcov.var
  out$variance.results <- variance.results
  out <- c(list(date=date(),aux=aux), out)
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
