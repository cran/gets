predict.arx <-
function(object, spec=NULL, n.ahead=5,
  newmxreg=NULL, newvxreg=NULL, n.sim=1000, innov=NULL,
  plot=TRUE, ...)
{
  ##spec:
  if(is.null(spec)){
    if(!is.null(object$mean.results)) spec <- "mean"
    if(is.null(object$mean.results)
        && !is.null(object$variance.results)) spec <- "variance"
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  } #end if(..)else(..)
  if(is.null(spec)){ stop("No estimated model") }

  ##n.ahead:
  if(n.ahead < 1){ stop("n.ahead must be 1 or greater") }

  ##newindex:
#  if(!is.null(newindex) &&
#    n.ahead!=length(newindex)){ stop("length(newindex) must equal n.ahead") }


  ##-------------
  ##if mean spec:
  ##-------------

  outMean <- NULL
  if(spec=="mean" || spec=="both"){

    coefs <- coef.arx(object, spec="mean")

    ##mc:
    if(!is.null(object$call$mc)){
      mconst <- as.numeric(coefs[1])
      mconstIndx <- 1
    }else{
      mconst <- 0
      mconstIndx <- 0
    }

    ##ar:
    arMax <- 0
    arIndx <- max(mconstIndx)
    if(!is.null(object$call$ar)){
      arEval <- eval(object$call$ar)
      arIndx <- 1:length(arEval) + max(mconstIndx)
      arMax <- max(arEval)
      arCoefs <- rep(0,arMax)
      arCoefs[arEval] <- as.numeric(coefs[arIndx])
    }

    ##ewma:
    ewmaMax <- 0
    ewmaIndx <- max(arIndx)
    if(!is.null(object$call$ewma)) stop("Sorry, 'ewma' not implemented yet")
#    if(!is.null(object$call$ewma)){
#      modify arEval
#      modify arIndx
#      modify arMax
#      modify arCoefs
#    }

    ##backcast length:
    backcastMax <- max(arMax,ewmaMax)

    ##mxreg:
    mxreghat <- rep(0, n.ahead + backcastMax)
    if(!is.null(object$call$mxreg)){

      ##check newmxreg:
      if(is.null(newmxreg)){ stop("'newmxreg' is NULL") }
      if(NROW(newmxreg)!=n.ahead){ stop("NROW(newmxreg) must equal n.ahead") }

      ##newmxreg:
      newmxreg <- coredata(cbind(as.zoo(newmxreg)))
      colnames(newmxreg) <- NULL

      ##mxreghat:
      mxregIndx <- c(max(ewmaIndx)+1):length(coefs)
      mxreghat <-  newmxreg %*% as.numeric(coefs[mxregIndx])
      mxreghat <- c(rep(0,backcastMax),mxreghat)

    } ##end mxreg

    ##prepare prediction:
    yhat <- rep(NA, n.ahead + backcastMax)
    yhat.n <- length(yhat)
    meanFit <- coredata(fitted.arx(object, spec="mean"))
    yhat[1:backcastMax] <- meanFit[c(length(meanFit)-backcastMax+1):length(meanFit)]

    ##predict:
    for(i in c(backcastMax+1):yhat.n){
      yhat[i] <- mconst + sum(arCoefs*yhat[c(i-1):c(i-arMax)]) + mxreghat[i]
    }

    ##out:
    outMean <- yhat[c(yhat.n-n.ahead+1):yhat.n]
    outMean <- as.zoo(outMean)

  } #end mean spec


  ##-----------------
  ##if variance spec:
  ##-----------------

  outVariance <- NULL
  if(spec=="variance" || spec=="both"){

    ##record coef estimates:
    coefs <- as.numeric(coef.arx(object, spec="variance"))
    Elnz2hat <- coefs[length(coefs)]
    coefs <- coefs[-length(coefs)]

    ##vc:
    vconst <- as.numeric(coefs[1])

    ##arch:
    archMax <- 0
    archIndx <- 1
    if(!is.null(object$call$arch)){
      archEval <- eval(object$call$arch)
      archIndx <- 1:length(archEval) + 1
      archMax <- max(archEval)
      archCoefs <- rep(0,archMax)
      archCoefs[archEval] <- as.numeric(coefs[archIndx])
    }

    ##asym:
    asymMax <- 0
    asymIndx <- max(archIndx)
    if(!is.null(object$call$asym)){
      asymEval <- eval(object$call$asym)
      asymIndx <- 1:length(asymEval) + max(archIndx)
      asymMax <- max(asymEval)
      asymCoefs <- rep(0,asymMax)
      asymCoefs[asymEval] <- as.numeric(coefs[asymIndx])
    }

    ##log.ewma:
    logewmaMax <- 0
    logewmaIndx <- max(asymIndx)
    if(!is.null(object$call$log.ewma)) stop("Sorry, 'log.ewma' not implemented yet")

    ##backcast length:
    backcastMax <- max(archMax,asymMax,logewmaMax)

    ##vxreg:
    vxreghat <- rep(0, n.ahead + backcastMax)
    if(!is.null(object$call$vxreg)){

      ##check newvxreg:
      if(is.null(newvxreg)){ stop("'newvxreg' is NULL") }
      if(NROW(newvxreg)!=n.ahead){ stop("NROW(newvxreg) must equal n.ahead") }

      ##newmxreg:
      newvxreg <- coredata(cbind(as.zoo(newvxreg)))
      colnames(newvxreg) <- NULL

      ##vxreghat:
      vxregIndx <- c(max(logewmaIndx)+1):length(coefs)
      vxreghat <-  newvxreg %*% coefs[vxregIndx]
      vxreghat <- c(rep(0,backcastMax),vxreghat)

    } #end vxreg

    ##prepare lnsd2:
    lnsd2hat <- rep(NA, n.ahead + backcastMax)
    lnsd2hat.n <- length(lnsd2hat)
    lnsd2Fit <- log(coredata(fitted.arx(object, spec="variance")))
    lnsd2hat[1:backcastMax] <- lnsd2Fit[c(length(lnsd2Fit)-backcastMax+1):length(lnsd2Fit)]
    mLnsd2Hat <- matrix(NA, lnsd2hat.n, n.sim)
    mLnsd2Hat[,1:NCOL(mLnsd2Hat)] <- lnsd2hat

    ##prepare lnz2:
    if(backcastMax>0){
      lnz2hat <- rep(NA, n.ahead + backcastMax)
      lnz2hat.n <- length(lnz2hat)
      lnz2Fit <- coredata(object$resids.ustar + Elnz2hat)
      lnz2hat[1:backcastMax] <- lnz2Fit[c(length(lnz2Fit)-backcastMax+1):length(lnz2Fit)]
      mLnz2Hat <- matrix(NA, lnz2hat.n, n.sim)
      mLnz2Hat[,1:NCOL(mLnz2Hat)] <- lnz2hat
    }

    ##zhat:
    if(backcastMax>0){

      ##bootstrap:
      if(is.null(innov)){
        zhat <- coredata(na.trim(object$resids.std))
        where.zeros <- which(zhat==0)
        if(length(where.zeros)>0){ zhat <- zhat[-where.zeros] }
        draws <- runif(n.ahead*n.sim, min=0.5+.Machine$double.eps,
          max=length(zhat)+0.5+.Machine$double.eps)
        draws <- round(draws, digits=0)
        zhat <- zhat[draws]
      }

      ##user-provided:
      if(!is.null(innov)){
        if(length(innov)!=n.ahead*n.sim){ stop("length(innov) must equal n.ahead*n.sim") }
        zhat <- as.numeric(innov)
      }

      zhat <- matrix(zhat,n.ahead,n.sim)
      mLnz2Hat[c(backcastMax+1):NROW(mLnz2Hat),] <- log(zhat^2)
    }

    ##prepare asym:
    if(asymMax>0){
      zhatIneg <- rep(NA, n.ahead + backcastMax)
      zhatIneg.n <- length(zhatIneg)
      zhatFit <- coredata(object$resids.std)
      zhatIneg[1:backcastMax] <- zhatFit[c(length(zhatFit)-backcastMax+1):length(zhatFit)]
      zhatIneg <- as.numeric(zhatIneg<0)
      mZhatIneg <- matrix(NA, zhatIneg.n, n.sim)
      mZhatIneg[,1:NCOL(mZhatIneg)] <- zhatIneg
      mZhatIneg[c(backcastMax+1):NROW(mZhatIneg),] <- matrix(as.numeric(zhat<0),NROW(zhat),NCOL(zhat))
    }

    ##predict:
    asymTerm <- 0
    for(j in 1:NCOL(mLnsd2Hat)){
      for(i in c(backcastMax+1):NROW(mLnsd2Hat)){
        archTerm <- sum( archCoefs*mLnsd2Hat[c(i-1):c(i-archMax),j] )
        lnz2Term <- sum( archCoefs*mLnz2Hat[c(i-1):c(i-archMax),j] )
        if(asymMax>0){
          asymTermSd2 <- sum( asymCoefs*mLnsd2Hat[c(i-1):c(i-asymMax),j]*mZhatIneg[c(i-1):c(i-asymMax),j] )
          asymTermLnz2 <- sum( asymCoefs*mLnz2Hat[c(i-1):c(i-asymMax),j]*mZhatIneg[c(i-1):c(i-asymMax),j] )
          asymTerm <- asymTermSd2 + asymTermLnz2
        }
        mLnsd2Hat[i,j] <- vconst + archTerm + lnz2Term + asymTerm + vxreghat[i]
      } ##end for(i)
    } ##end for(j)

    ##out:
    outVariance <- mLnsd2Hat[c(lnsd2hat.n-n.ahead+1):lnsd2hat.n,]
    outVariance <- rowMeans(exp(outVariance))
    outVariance <- as.zoo(outVariance)

  } #end variance spec

  ##out:
  if(spec=="mean"){ out <- outMean }
  if(spec=="variance"){ out <- outVariance }
  if(spec=="both"){
    out <- cbind(outMean,outVariance)
    colnames(out) <- c("mean","variance")
  }
  #if(!is.null(newindex)){
  #  out <- zoo(coredata(out), order.by="index")
  #}

  ##plot:
  if(plot){
#    if(is.null(newindex)){
      xlabArg <- "Step ahead"
#    }else{
#      xlabArg <- ""
#    }
    if(spec=="both"){
      ylabArg <- c("Mean","Variance")
    }else{
      if(spec=="mean"){ ylabArg <- "Mean" }
      if(spec=="variance"){ ylabArg <- "Variance" }
    }
    plot(out, xlab=xlabArg, ylab=ylabArg, main="Forecasts",
      col="red")
  } #end if(plot)

  ##out:
  return(out)

}
