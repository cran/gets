predict.arx <-
function(object, spec=NULL, n.ahead=12,
  newmxreg=NULL, newvxreg=NULL, newindex=NULL,
  n.sim=1000, innov=NULL, return=TRUE, plot=NULL,
  plot.options=list(), ...)
{

  ## contents:
  ## 0 initialise
  ## 1 if mean spec
  ## 2 if variance spec
  ## 3 if plot=TRUE
  ## 4 if return=TRUE

  ##n.ahead:
  if(n.ahead < 1){ stop("n.ahead must be 1 or greater") }

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

  ##newindex:
  yInSample <- zoo(object$aux$y, order.by=object$aux$y.index)

  if(!is.null(newindex)){
    yAsRegular <- FALSE
    if( n.ahead!=length(newindex) ){
      stop("length(newindex) must equal n.ahead")
    }
  }

  if( is.null(newindex) && is.regular(yInSample, strict=TRUE) ){
    endCycle <- cycle(yInSample)
    endCycle <- as.numeric(endCycle[length(endCycle)])
    endYear <- floor(as.numeric(object$aux$y.index[object$aux$y.n]))
    yFreq <- frequency(yInSample)
    yhataux <- rep(NA, n.ahead+1)
    yDeltat <- deltat(yInSample)
    if( yDeltat==1 && yFreq==1 ){
      yhataux <- zoo(yhataux,
        order.by=seq(endYear, endYear+n.ahead, by=1))
      yAsRegular <- FALSE
    }else{
      yhataux <- zooreg(yhataux, start=c(endYear, endCycle),
                frequency=yFreq)
      yAsRegular <- TRUE
    }
    yhataux <- yhataux[-1]
    newindex <- index(yhataux)
  }

  if(is.null(newindex)){ newindex <- 1:n.ahead }

  ##----------------
  ## 1 if mean spec
  ##----------------

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
    arCoefs <- NULL ##J-dog change
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
      #      if(NROW(newmxreg)!=n.ahead
      #        && NROW(newmxreg)!=1){ stop("NROW(newmxreg) must equal 1 or n.ahead") }
      #      if(NROW(newmxreg)!=n.ahead
      #        && NROW(newmxreg)==1){
      #          newmxreg <- coredata(rbind(as.zoo(newmxreg)))
      #          tmp <- matrix(NA,n.ahead,NCOL(newmxreg))
      #          tmp[1:NROW(tmp),] <- newmxreg
      #          newmxreg <- tmp
      #      }

      ##newmxreg:
      newmxreg <- coredata(cbind(as.zoo(newmxreg)))
      colnames(newmxreg) <- NULL

      ##mxreghat:
      mxregIndx <- c(max(ewmaIndx)+1):length(coefs) ##J-dog change

      mxreghat <-  newmxreg %*% as.numeric(coefs[mxregIndx])
      mxreghat <- c(rep(0,backcastMax),mxreghat)

    } ##end mxreg

    ##prepare prediction:
    yhat <- rep(NA, n.ahead + backcastMax)
    yhat.n <- length(yhat)
#OLD:
#    meanFit <- coredata(fitted.arx(object, spec="mean"))
    if(backcastMax>0) {
      ##actual y-values:
      yhat[1:backcastMax] <- tail(object$aux$y, n=backcastMax)
#        OR:
#        object$aux$y[c(length(object$aux$y)-backcastMax+1):length(object$aux$y)]
#        OLD (fitted y-values); this was erroneous!:
#        meanFit[c(length(meanFit)-backcastMax+1):length(meanFit)]
    }

    ##predict:
    for(i in c(backcastMax+1):yhat.n){
      yhat[i] <- mconst + sum(arCoefs*yhat[c(i-1):c(i-arMax)]) + mxreghat[i]
    }

    ##out:
    outMean <- yhat[c(yhat.n-n.ahead+1):yhat.n]
    #outMean <- as.zoo(outMean)

  } #end mean spec


  ##--------------------
  ## 2 if variance spec
  ##--------------------

  outVariance <- NULL
  if(spec=="variance" || spec=="both"){ # || !is.null(plot.options$errors.only)

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
    if(!is.null(object$call$log.ewma)){
      logewmaEval <- eval(object$call$log.ewma)
      if(is.list(logewmaEval)){ logewmaEval <- logewmaEval$length }
      logewmaIndx <- 1:length(logewmaEval) + max(asymIndx)
      logewmaMax <- max(logewmaEval)
      logewmaCoefs <- as.numeric(coefs[logewmaIndx])
    }

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
      lnz2Fit <- coredata(object$ustar.residuals + Elnz2hat)
      lnz2hat[1:backcastMax] <- lnz2Fit[c(length(lnz2Fit)-backcastMax+1):length(lnz2Fit)]
      mLnz2Hat <- matrix(NA, lnz2hat.n, n.sim)
      mLnz2Hat[,1:NCOL(mLnz2Hat)] <- lnz2hat
    }

    ##zhat:
    ##-----

    ##bootstrap innov (not user-provided):
    if(is.null(innov)){
      zhat <- coredata(na.trim(object$std.residuals))
      where.zeros <- which(zhat==0)
      if(length(where.zeros)>0){ zhat <- zhat[-where.zeros] }
      draws <- runif(n.ahead*n.sim, min=0.5+.Machine$double.eps,
                     max=length(zhat)+0.5+.Machine$double.eps)
      draws <- round(draws, digits=0)
      zhat <- zhat[draws]
    }

    ##user-provided innov:
    if(!is.null(innov)){
      if(length(innov)!=n.ahead*n.sim){ stop("length(innov) must equal n.ahead*n.sim") }
      if(any(innov==0)){ stop("innov cannot contain zeros") }
      zhat <- as.numeric(innov)
    }

    zhat <- matrix(zhat,n.ahead,n.sim)
    mZhat2 <- zhat^2
    mLnz2Hat[c(backcastMax+1):NROW(mLnz2Hat),] <- log(mZhat2)

    vEpsilon2 <- rep(NA, n.ahead+backcastMax)
    vEpsilon2[1:backcastMax] <- as.numeric(object$residuals[c(length(object$residuals)-backcastMax+1):length(object$residuals)]^2)
    mZhat2 <- rbind(matrix(NA,backcastMax,NCOL(mZhat2)),mZhat2)


    ##prepare asym:
    if(asymMax>0){
      zhatIneg <- rep(NA, n.ahead + backcastMax)
      zhatIneg.n <- length(zhatIneg)
      zhatFit <- coredata(object$std.residuals)
      zhatIneg[1:backcastMax] <- zhatFit[c(length(zhatFit)-backcastMax+1):length(zhatFit)]
      zhatIneg <- as.numeric(zhatIneg<0)
      mZhatIneg <- matrix(NA, zhatIneg.n, n.sim)
      mZhatIneg[,1:NCOL(mZhatIneg)] <- zhatIneg
      mZhatIneg[c(backcastMax+1):NROW(mZhatIneg),] <- matrix(as.numeric(zhat<0),NROW(zhat),NCOL(zhat))
    }

    ##prepare log.ewma:
    if(logewmaMax>0){
      mLogEwmaHat <- matrix(NA, n.ahead+backcastMax, length(logewmaCoefs))
      colnames(mLogEwmaHat) <- object$aux$vXnames[logewmaIndx]
      mLogEwmaHat[1:backcastMax,] <- object$aux$vX[c(NROW(object$aux$vX)-backcastMax+1):NROW(object$aux$vX),logewmaIndx]
      mLogEwmaHat <- as.matrix(mLogEwmaHat)
    }

    ##predict:
    archTerm <- 0
    lnz2Term <- 0
    asymTerm <- 0
    logewmaTerm <- 0
    for(j in 1:NCOL(mLnsd2Hat)){
      for(i in c(backcastMax+1):NROW(mLnsd2Hat)){
        if(archMax>0){
          archTerm <- sum( archCoefs*mLnsd2Hat[c(i-1):c(i-archMax),j] )
          lnz2Term <- sum( archCoefs*mLnz2Hat[c(i-1):c(i-archMax),j] )
        }
        if(asymMax>0){
          asymTermSd2 <- sum( asymCoefs*mLnsd2Hat[c(i-1):c(i-asymMax),j]*mZhatIneg[c(i-1):c(i-asymMax),j] )
          asymTermLnz2 <- sum( asymCoefs*mLnz2Hat[c(i-1):c(i-asymMax),j]*mZhatIneg[c(i-1):c(i-asymMax),j] )
          asymTerm <- asymTermSd2 + asymTermLnz2
        }
        if(logewmaMax>0){
          for(k in 1:NCOL(mLogEwmaHat)){
            mLogEwmaHat[i,k] <- log( mean(vEpsilon2[c(i-logewmaEval[k]):c(i-1)]) )
          }
          logewmaTerm <- sum( coefs[logewmaIndx] * mLogEwmaHat[i,] )
        }
        mLnsd2Hat[i,j] <- vconst + archTerm + lnz2Term + asymTerm + logewmaTerm + vxreghat[i]
        vEpsilon2[i] <- exp(mLnsd2Hat[i,j])*mZhat2[i,j]
      } ##end for(i)
    } ##end for(j)

    ##out:
    outVariance <- mLnsd2Hat[c(lnsd2hat.n-n.ahead+1):lnsd2hat.n,]
    outVariance <- rowMeans(rbind(exp(outVariance)))

  } #end variance spec


  ##--------------------
  ## 3 out object
  ##--------------------

  ##out:
  if(spec=="mean"){ out <- outMean }
  if(spec=="variance"){ out <- outVariance }
  if(spec=="both"){
    out <- cbind(outMean,outVariance)
    colnames(out) <- c("mean","variance")
  }
  if(yAsRegular){
    startYear <- floor(as.numeric(index(yhataux))[1])
    startCycle <- cycle(yhataux)[1]
    out <- zooreg(out, start=c(startYear, startCycle),
      frequency=yFreq)
  }else{
    out <- zoo(out, order.by=newindex)
  }


  ##----------------
  ## 4 if plot=TRUE
  ##----------------

  ##determine plot-argument:
  plotArg <- plot
  if( is.null(plotArg) ){
    plotArg <- getOption("plot")
    if( is.null(plotArg) ){ plotArg <- FALSE }
  }

  ##if plot=TRUE:
  if(plotArg){

    ##how many observations to retain?:
    if(is.null(plot.options$keep)) {
      plot.options$keep=12L
    }

    ##include retained fitted?
    if(is.null(plot.options$fitted)) {
      plot.options$fitted <- FALSE
    }

#    ##WHAT DOES THIS DO???:
#    ##errors only?
#    if(!is.null(plot.options$errors.only)) {
##      plot.options$errors.only <- TRUE
#      if((spec=="mean")){
#        print("Cannot plot estimated error bars when mean specified")
#      }
#    }

    ##location of legend:
    if(is.null(plot.options$legend.loc)) {
      plot.options$legend.loc <- "topleft"
    }

    ##xlab:
    xlabArg <- ifelse(is.null(newindex), "Step ahead", "")

    ##make dataForPlot:
    dataForPlot <- matrix(NA, NROW(out), 6)
    colnames(dataForPlot) <- c("MeanActual", "MeanFitted",
      "MeanPrediction", "VarianceActual", "VarianceFitted",
      "VariancePrediction")
    if(!is.null(plot.options$newmactual)){
      dataForPlot[1:length(plot.options$newmactual),"MeanActual"] <-
        plot.options$newmactual
    }
    if(!is.null(outMean)){ dataForPlot[,"MeanPrediction"] <- outMean }
    if(!is.null(plot.options$newvactual)){
      dataForPlot[1:length(plot.options$newvactual),"VarianceActual"] <-
        plot.options$newvactual
    }
    if(is.null(outVariance)){
      dataForPlot[,"VariancePrediction"] <- sigma.arx(object)^2
    }else{
      dataForPlot[,"VariancePrediction"] <- outVariance
    }
    if(plot.options$keep > 0){
      retainedData <- matrix(NA, plot.options$keep, NCOL(dataForPlot))
      colnames(retainedData) <- c("MeanActual", "MeanFitted",
        "MeanPrediction", "VarianceActual", "VarianceFitted",
        "VariancePrediction")
      retainedData[,"MeanActual"] <-
        tail(coredata(yInSample), n=plot.options$keep)
      retainedData[,"MeanFitted"] <-
        tail(coredata(object$mean.fit), n=plot.options$keep)
      retainedData[,"VarianceActual"] <-
        tail(coredata(object$residuals^2), n=plot.options$keep)
      retainedData[,"VarianceFitted"] <-
        tail(coredata(object$var.fit), n=plot.options$keep)
      dataForPlot <- rbind(retainedData, dataForPlot)
      tmpIndx <- c(tail(index(yInSample), n=plot.options$keep),
        index(out))
      dataForPlot <- zoo(dataForPlot, order.by=tmpIndx)
    }else{
      dataForPlot <- zoo(dataForPlot, order.by=index(out))
    }

    ##get current par-values:
    def.par <- par(no.readonly=TRUE)

    ##prepare plotting:
    if(xlabArg=="") { #plot margins?
      par(mar = c(2,3,0.5,0.5) + 0.1) #b,l,t,r
    }else{
      par(mar = c(3,3,0.5,0.5) + 0.1) #b,l,t,r
    }
    plotTypeForecast <- ifelse(n.ahead==1, "p", "l")
    plotTypeRetained <- ifelse(plot.options$keep > 1, "l", "p")

    ##if spec="mean" or "both":
    ##-------------------------

    if( spec=="mean" || spec=="both" ){

      ##y-axis (limits):
      ylimVals <- range( na.trim(dataForPlot[,"MeanActual"]),
        na.trim(dataForPlot[,"MeanFitted"]),
        na.trim(dataForPlot[,"MeanPrediction"]
          + 2*sqrt(dataForPlot[,"VariancePrediction"])),
        na.trim(dataForPlot[,"MeanPrediction"]
          - 2*sqrt(dataForPlot[,"VariancePrediction"])) )

      ##plot the mean-predictions:
      plot.zoo(dataForPlot[,"MeanPrediction"], xlab="", ylab="",
        col="red", type=plotTypeForecast, ylim=ylimVals)

      ##add 2 x SEs:
      SE <- sqrt(dataForPlot[,"VariancePrediction"])
      lines(dataForPlot[,"MeanPrediction"] + 2*SE,
        col="darkgreen", lty=2, type=plotTypeForecast)
      lines(dataForPlot[,"MeanPrediction"] - 2*SE,
        col="darkgreen", lty=2, type=plotTypeForecast)

      ##add actuals (pre- and post-prediction):
      if(plot.options$keep > 0){
        lines(dataForPlot[,"MeanActual"], col="blue",
          type=plotTypeRetained)
      }

      ##add fitted (pre-prediction):
      if( plot.options$keep > 0 && plot.options$fitted ){
        lines(dataForPlot[,"MeanFitted"], col="red",
          type=plotTypeRetained)
      }

      ##add text closer to plot than xlab or ylab would do
      if(xlabArg!="") {
        mtext(xlabArg, side=1, line=2)
      }
      mtext("Mean", side=2, line=2)

      ##add legend:
      legend(plot.options$legend.loc, col=c("red","darkgreen","blue"),
        lty=c(1,2,1),
        legend=c("Forecast","+/- 2 x SE of regression","Actual"),
        bty="n")

    } #close if spec="mean" or "both"


    ##if spec="variance":
    ##-------------------

    if(spec=="variance"){

      ##y-axis (limits):
      ylimVals <- range( na.trim(dataForPlot[,"VarianceActual"]),
        na.trim(dataForPlot[,"VarianceFitted"]),
        na.trim(dataForPlot[,"VariancePrediction"]) )

      ##plot the variance-predictions:
      plot.zoo(dataForPlot[,"VariancePrediction"], xlab="", ylab="",
        col="red", type=plotTypeForecast, ylim=ylimVals)

      ##add fitted (pre-prediction):
      if( plot.options$keep > 0 && plot.options$fitted ){
        lines(dataForPlot[,"VarianceFitted"], col="red",
          type=plotTypeRetained)
      }

      ##add actuals (pre- and post-prediction):
      if(plot.options$keep > 0){
        lines(dataForPlot[,"VarianceActual"], col="blue",
          type=plotTypeRetained)
      }

      ##add text closer to plot than xlab or ylab would do
      if(xlabArg!="") {
        mtext(xlabArg, side=1, line=2)
      }
      mtext("Variance", side=2, line=2)

      ##add legend:
      legend(plot.options$legend.loc, col=c("red","blue"),
        lty=c(1,1), legend=c("Forecast (conditional variance)", "Actual (residuals^2)"), bty="n")

    } #close if spec="variance"


    ##add vertical green line between last observation
    ##and first forecast:
    if(plot.options$keep > 0){
      abline( v=( as.numeric(index(dataForPlot))[plot.options$keep] + as.numeric(index(dataForPlot))[plot.options$keep+1] )/2,
        col="darkgreen", lty=3)
    }

    ##return to old par-values:
    par(def.par)

  } #end if(plot)


  ##------------------
  ## 5 if return=TRUE
  ##------------------

  if(return){ return(out) }

}
