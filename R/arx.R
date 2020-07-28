arx <-
function(y, mc=FALSE, ar=NULL, ewma=NULL, mxreg=NULL,
  vc=FALSE, arch=NULL, asym=NULL, log.ewma=NULL, vxreg=NULL,
  zero.adj=0.1, vc.adj=TRUE,
  vcov.type=c("ordinary", "white", "newey-west"),
  qstat.options=NULL, normality.JarqueB=FALSE, user.estimator=NULL,
  user.diagnostics=NULL, tol=1e-07, LAPACK=FALSE, plot=NULL)
{
  ### ARGUMENTS: ###########

  vcov.type <- match.arg(vcov.type)

  ##regressand, regressors:
  tmp <- regressorsMean(y, mc=mc, ar=ar, ewma=ewma, mxreg=mxreg,
    return.regressand=TRUE, return.as.zoo=TRUE,
    na.trim=TRUE,
    na.omit=FALSE)
  
  ##aux: auxiliary list, also used by getsm/getsv
  aux <- list()
  aux$y <- coredata(tmp[,1])
  aux$y.n <- length(aux$y)
  aux$y.name <- colnames(tmp)[1]
  aux$y.index <- index(tmp)
  if( NCOL(tmp)>1 ){
    aux$mX <- cbind(coredata(tmp[,-1]))
    aux$mXnames <- colnames(tmp)[-1]
    colnames(aux$mX) <- NULL
    aux$mXncol <- NCOL(aux$mX)
  }

  ##modify vxreg:
  if( !is.null(vxreg) ){
    ##time:
    vxreg <- as.zoo(cbind(vxreg))
    vxreg <- window(vxreg, start=aux$y.index[1],
      end=aux$y.index[length(aux$y.index)])
    ##colnames:
    vxreg.names <- colnames(vxreg)
    if(is.null(vxreg.names)){
      vxreg.names <- paste0("vxreg", 1:NCOL(vxreg))
    }
    if( any(vxreg.names == "") ){
      missing.colnames <- which(vxreg.names == "")
      for(i in 1:length(missing.colnames)){
        vxreg.names[missing.colnames[i]] <- paste0("vxreg", i)
      }
    }
    colnames(vxreg) <- vxreg.names
    ##add to aux (necessary for getsm/getsv):
    aux$vxreg <- vxreg #note: NROW(vxreg)!=NROW(vX) is possible
  }

  ##determine qstat.options:
  if(is.null(qstat.options)){
    if(is.null(ar)){ar.lag <- 1}else{ar.lag <- max(ar)+1}
    if(is.null(arch)){arch.lag <- 1}else{arch.lag <- max(arch)+1}
    qstat.options <- c(ar.lag, arch.lag)
  }
  
  ##info for getsm/getsv functions
  aux$vcov.type <- vcov.type
  aux$qstat.options <- qstat.options
  aux$user.estimator <- user.estimator
  aux$user.diagnostics <- user.diagnostics
  aux$tol <- tol
  aux$LAPACK <- LAPACK

  ### INITIALISE ##########

  sysCall <- sys.call()
  #for the future: make sure the following objects are part of the out-list?
  vcov.var <- NULL #make sure this object exists
  variance.results <- NULL #make sure this object exists

  #### for the future regarding user.estimator: check if
  #### user.estimator$spec is NULL, "mean", "variance" or "both"
  #### in order to determine what kind of estimator it is

  ##check if mean and log-garch spec:
  meanSpec <- !is.null(aux$mX)
  varianceSpec <- if( vc==FALSE && is.null(arch)
    && is.null(asym) && is.null(log.ewma)
    && is.null(vxreg) ){ FALSE }else{ TRUE }

  #### DEFAULT ESTIMATOR ###############

  if( is.null(user.estimator) ){

    ##estimate:
    estMethod <- which(vcov.type==c("none", "none", "ordinary",
      "white", "newey-west"))
    varianceSpecArg <- NULL
    if( varianceSpec ){
      ##note: vc must be TRUE
      varianceSpecArg <- list(vc=TRUE, arch=arch, asym=asym,
        log.ewma=log.ewma, vxreg=vxreg)
    }
    out <- ols(aux$y, aux$mX, tol=tol, LAPACK=LAPACK, method=estMethod,
      variance.spec=varianceSpecArg)

    ##delete some unneeded entries:
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

    ##re-organise stuff related to mean spec:
    colnames(out$vcov) <- aux$mXnames
    rownames(out$vcov) <- aux$mXnames
    outNames <- names(out)
    whereIs <- which(outNames=="vcov")
    if( length(whereIs) > 0 ){ names(out)[whereIs] <- "vcov.mean" }
    whereIs <- which(outNames=="fit")
    names(out)[whereIs] <- "mean.fit"

    ##if no variance spec:
    if( varianceSpec==FALSE ){
      out$var.fit <- rep(out$sigma2, aux$y.n)
      out$std.residuals <- out$residuals/sqrt(out$sigma2)
      aux$loge2.n <- aux$y.n #same as out$n, change?
      aux$vc <- FALSE #needed for specific in getsm()
    }

    ##if variance spec:
    if( varianceSpec ){

      ##aux: info for getsm and getsv:
      aux$vc <- TRUE #obligatory if varianceSpec
      aux$zero.adj <- zero.adj
      aux$vc.adj <- vc.adj
      aux$loge2 <- out$regressorsVariance[,1]
      aux$loge2.n <- length(aux$loge2)
      aux$vX <- cbind(out$regressorsVariance[,-1])
      aux$vXnames <- colnames(out$regressorsVariance)[-1]
      colnames(aux$vX) <- NULL
      aux$vXncol <- NCOL(aux$vX)
      aux$arch <- arch
      aux$asym <- asym
      aux$log.ewma <- log.ewma
      out$regressorsVariance <- NULL #delete, not needed anymore

      ##re-organise stuff related to variance spec:
      s.e. <- sqrt(as.vector(diag(out$vcov.var)))
      tmpdf <- aux$loge2.n - length(out$var.coefficients)
      tmpvcov <- as.matrix(out$vcov.var[-1,-1])
      colnames(tmpvcov) <- aux$vXnames[-1]
      rownames(tmpvcov) <- aux$vXnames[-1]
      t.stat <- out$var.coefficients/s.e.
      p.val <- pt(abs(t.stat), tmpdf, lower.tail=FALSE)*2
      t.stat[1] <- ((out$var.coefficients[1]-out$Elnz2)^2)/s.e.[1]^2
      p.val[1] <- pchisq(t.stat[1], 1, lower.tail=FALSE)
      out$var.coefficients[1] <- out$var.coefficients[1] - out$Elnz2
      out$n <- aux$loge2.n
      out$vcov.var <- tmpvcov
      out$variance.results <-
        as.data.frame(cbind(out$var.coefficients, s.e., t.stat, p.val))
      colnames(out$variance.results) <- c("coef", "std.error", "t-stat", "p-value")
      rownames(out$variance.results) <- aux$vXnames
      out$var.coefficients <- NULL

    } #close if( varianceSpec )
        
  } #close if( is.null(user.estimator) )
  

  #### USER-DEFINED ESTIMATOR ###############

  if( !is.null(user.estimator) ){

    ##make user-estimator argument:
    if( is.null(user.estimator$envir) ){ user.estimator$envir <- .GlobalEnv }
    userEstArg <- user.estimator
    userEstArg$name <- NULL
    userEstArg$envir <- NULL
    if( length(userEstArg)==0 ){ userEstArg <- NULL }

    ##user-defined estimator:
    if( is.null(user.estimator$envir) ){
      out <- do.call(user.estimator$name,
        c(list(aux$y,aux$mX), userEstArg))
    }else{
      out <- do.call(user.estimator$name,
        c(list(aux$y,aux$mX), userEstArg), envir=user.estimator$envir)
#        c(list(y=aux$y,x=aux$mX), userEstArg), envir=user.estimator$envir)
    }
    
#delete?:
    ##just in case...:
    if( is.null(out$vcov) && !is.null(out$vcov.mean) ){
      out$vcov <- out$vcov.mean
    }

  } #end if( user.estimator )


  ### OUTPUT: ######################

  ##mean estimation results (a data frame):
  if( meanSpec ){
    if( !is.null(out$vcov) ){
      coefvar <- out$vcov
    }else{
      coefvar <- out$vcov.mean
    }
    stderrs <- sqrt(diag(coefvar))
    t.stat <- out$coefficients/stderrs
    p.val <- pt(abs(t.stat), out$df, lower.tail=FALSE)*2
    out$mean.results <- as.data.frame(cbind(out$coefficients,
      stderrs, t.stat, p.val))
    colnames(out$mean.results) <- c("coef", "std.error",
      "t-stat", "p-value")
    rownames(out$mean.results) <- aux$mXnames
  } #end if(meanSpec)

  ##diagnostics:
  if( any( names(out) %in% c("residuals", "std.residuals") ) ){
    ar.LjungBarg <- c(qstat.options[1],0)
    arch.LjungBarg <- c(qstat.options[2],0)
    if( identical(normality.JarqueB,FALSE) ){
      normality.JarqueBarg <- NULL
    }else{
      normality.JarqueBarg <- as.numeric(normality.JarqueB)
    }
  }else{
    ar.LjungBarg <- NULL
    arch.LjungBarg <- NULL
    normality.JarqueBarg <- NULL
  }
  out$diagnostics <- diagnostics(out,
    ar.LjungB=ar.LjungBarg, arch.LjungB=arch.LjungBarg,
    normality.JarqueB=normality.JarqueBarg,
    user.fun=user.diagnostics, verbose=TRUE)
    
  ##add zoo-indices:
  if(!is.null(out$mean.fit)){
    out$mean.fit <- zoo(out$mean.fit, order.by=aux$y.index)
  }
  if(!is.null(out$residuals)){
    out$residuals <- zoo(out$residuals, order.by=aux$y.index)
  }
  if(!is.null(out$var.fit)){
    out$var.fit <- zoo(out$var.fit, order.by=aux$y.index)
  }
  if(!is.null(out$ustar.residuals)){
    out$ustar.residuals <- zoo(out$ustar.residuals, order.by=aux$y.index)
  }
  if(!is.null(out$std.residuals)){
    out$std.residuals <- zoo(out$std.residuals, order.by=aux$y.index)
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
