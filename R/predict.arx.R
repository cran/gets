predict.arx <-
function(object, spec=NULL, n.ahead=12,
  newmxreg=NULL, newvxreg=NULL, newindex=NULL,
  n.sim=5000, innov=NULL, probs=NULL, ci.levels=NULL,
  quantile.type=7, return=TRUE, verbose=FALSE, plot=NULL,
  plot.options=list(), ...)
{

  ## contents:
  ## 0 initialise
  ## 1 simulate innov
  ## 2 variance predictions
  ## 3 mean predictions
  ## 4 probs (quantiles)
  ## 5 newindex
  ## 6 if plot=TRUE
  ## 7 if return=TRUE

  ##-----------------------
  ## 0 initialise
  ##-----------------------

  ##name of object:
  objectName <- deparse(substitute(object))

  ##check n.ahead:
  if(n.ahead < 1){ stop("n.ahead must be 1 or greater") }

  ##determine spec argument:
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
  
  ##what needs to be predicted?
  predictMean <- switch(spec, "mean"=TRUE, "both"=TRUE,
    "variance"=FALSE)
  predictVariance <- switch(spec, "mean"=FALSE, "both"=TRUE,
    "variance"=TRUE)
  
  ##is there a mean spec?
  coefs <- as.numeric(coef.arx(object, spec="mean"))
  if(length(coefs)>0){ specMean <- TRUE }else{ specMean <- FALSE }

  ##is there a variance spec?
  coefs <- as.numeric(coef.arx(object, spec="variance"))
  if(length(coefs)>0){ specVar <- TRUE }else{ specVar <- FALSE }

  ##determine plot-argument:
  plotArg <- plot
  if( is.null(plotArg) ){
    plotArg <- getOption("plot")
    if( is.null(plotArg) ){ plotArg <- FALSE }
  }

  ##probs argument:
  if( !is.null(probs) ){
    if( any(probs <= 0) || any(probs >= 1) ){
      stop("the values of 'probs' must be between 0 and 1")
    }
    probs <- union(probs,probs) #ensure values are distinct/not repeated
    probs <- probs[order(probs, decreasing=FALSE)] #re-order to increasing
  }
  probsArg <- probs

  ##ci.levels argument:
  if( is.null(ci.levels) && plotArg==TRUE ){
    if( class(object)=="isat" ){ ##if isat:
      ciLevelsArg <- c(0.68,0.95)
    }else{ ##if not isat:
      ciLevelsArg <- c(0.5,0.9)    
    }
  }else{
    ciLevelsArg <- ci.levels
  }
  if( !is.null(ciLevelsArg) ){
    if( any(ciLevelsArg <= 0) || any(ciLevelsArg >= 1) ){
      stop("'ci.levels' must be between 0 and 1")
    }
    ciLevelsArg <- union(ciLevelsArg, ciLevelsArg) #ensure levels are distinct/not repeated
    ciLevelsArg <- ciLevelsArg[order(ciLevelsArg, decreasing=FALSE)] #ensure levels are increasing
    ciLower <- (1-ciLevelsArg)/2
    ciLower <- ciLower[order(ciLower, decreasing=FALSE)]    
    ciUpper <- ciLevelsArg + (1-ciLevelsArg)/2
    ciUpper <- ciUpper[order(ciUpper, decreasing=TRUE)]    
    probsArg <- union(probsArg, c(ciLower,ciUpper) ) #add to probsArg
    probsArg <- probsArg[order(probsArg, decreasing=FALSE)] #ascending
  }
  
  ##are simulations of innov needed?
  doSimulations <- FALSE
  if( !is.null(probsArg) ){ doSimulations <- predictVariance <- TRUE }
  if( predictVariance && specVar ){ doSimulations <- TRUE }


  ##-----------------------
  ## 1 simulate innov
  ##-----------------------

  ##simulate:
  mZhat <- NULL
  if(doSimulations){

    ##bootstrap (innov not user-provided):
    if(is.null(innov)){
      zhat <- coredata(na.trim(object$std.residuals))
      if(specVar){
        where.zeros <- which(zhat==0)
        if(length(where.zeros)>0){ zhat <- zhat[-where.zeros] }
      }
      draws <- runif(n.ahead*n.sim, min=0.5+.Machine$double.eps,
                     max=length(zhat)+0.5+.Machine$double.eps)
      draws <- round(draws, digits=0)
      zhat <- zhat[draws]
    }
    
    ##user-provided innov:
    if(!is.null(innov)){
      if(length(innov)!=n.ahead*n.sim){ stop("length(innov) must equal n.ahead*n.sim") }
      if(specVar){
        if(any(innov==0)){ stop("'innov' cannot contain zeros") }
      }
      zhat <- as.numeric(innov)
    }

    ##matrix of innovations:
    mZhat <- matrix(zhat,n.ahead,n.sim)
    colnames(mZhat) <- paste0("mZhat.", seq(1,n.sim))
  
  } #end if(doSimulations)
  
  
  ##-----------------------
  ## 2 variance predictions
  ##-----------------------

  sd2hat <- NULL #variance predictions
  mEpsilon <- NULL #matrix of simulated errors
  if(predictVariance){
    
    ##there is no variance specification:
    if(specVar==FALSE){
      sigmahat <- sigma.arx(object)
      sd2hat <- rep(sigmahat^2, n.ahead)
    }

    ##there is a variance specification:
    if(specVar==TRUE){

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
      if(backcastMax>0){
        lnsd2hat[1:backcastMax] <- 
          lnsd2Fit[c(length(lnsd2Fit)-backcastMax+1):length(lnsd2Fit)]
      }
      mLnsd2Hat <- matrix(NA, lnsd2hat.n, n.sim) #matrix of lnsd2 predictions
      mLnsd2Hat[,1:NCOL(mLnsd2Hat)] <- lnsd2hat  #fill with backcast values
  
      ##prepare lnz2:
      lnz2hat <- rep(NA, n.ahead + backcastMax)
      lnz2hat.n <- length(lnz2hat)
      lnz2Fit <- coredata(object$ustar.residuals) + Elnz2hat
      if(backcastMax>0){
        lnz2hat[1:backcastMax] <- lnz2Fit[c(length(lnz2Fit)-backcastMax+1):length(lnz2Fit)]
      }
      mLnz2Hat <- matrix(NA, lnz2hat.n, n.sim)
      mLnz2Hat[,1:NCOL(mLnz2Hat)] <- lnz2hat
      mZhat2 <- mZhat^2 #mZhat from section 1
      mLnz2Hat[c(backcastMax+1):NROW(mLnz2Hat),] <- log(mZhat2)
      vEpsilon2 <- rep(NA, n.ahead+backcastMax) #needed for log(ewma) term(s)
      if(backcastMax>0){
        vEpsilon2[1:backcastMax] <- as.numeric(object$residuals[c(length(object$residuals)-backcastMax+1):length(object$residuals)]^2)
        mZhat2 <- rbind(matrix(NA,backcastMax,NCOL(mZhat2)),mZhat2)
      }
   
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
      mSd2Hat <- exp( mLnsd2Hat[c(lnsd2hat.n-n.ahead+1):lnsd2hat.n,] )
      if(n.ahead==1){ mSd2Hat <- rbind(mSd2Hat) } #rbind() needed when n.ahead=1
      sd2hat <- as.vector(rowMeans(mSd2Hat))
  
    } #end if(specVar==TRUE)

    ##matrix of errors:
    if(!is.null(mZhat)){
      mEpsilon <- sqrt(sd2hat)*mZhat
      colnames(mEpsilon) <- paste0("mEpsilon.", seq(1,n.sim))
    }

  } #end if(predictVariance)


  ##-----------------------
  ## 3 mean predictions
  ##-----------------------

  yhat <- NULL #mean predictions
  mY <- NULL #matrix of simulated y's
  if(predictMean){
    
    ##there is no mean specification:
    if(specMean==FALSE){
      yhat <- rep(0, n.ahead)
      if(!is.null(mEpsilon)){
        mY <- yhat + mEpsilon
        colnames(mY) <- paste0("mY.", seq(1,n.sim))
      }
    }

    ##there is a mean specification:
    if(specMean==TRUE){

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
      if( !is.null(object$call$ewma) ){
        ewmaEval <- eval(object$call$ewma)
        if(is.list(ewmaEval)){ ewmaEval <- ewmaEval$length }
        ewmaIndx <- 1:length(ewmaEval) + max(arIndx)
        ewmaMax <- max(ewmaEval)
        ewmaCoefs <- as.numeric(coefs[ewmaIndx])
      }
  
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
      if(backcastMax>0) {
        ##actual y-values:
        yhat[1:backcastMax] <- tail(object$aux$y, n=backcastMax)
      }

      ##prepare ewma:
      if(ewmaMax>0){
        mEwmaHat <- matrix(NA, n.ahead+backcastMax, length(ewmaCoefs))
        colnames(mEwmaHat) <- object$aux$mXnames[ewmaIndx]
        mEwmaHat[1:backcastMax,] <- object$aux$mX[c(NROW(object$aux$mX)-backcastMax+1):NROW(object$aux$mX),ewmaIndx]
        mEwmaHat <- as.matrix(mEwmaHat)
      }

      ##predict yhat:
      arTerm <- 0
      ewmaTerm <- 0
      for(i in c(backcastMax+1):yhat.n){
        if( arMax>0 ){ arTerm <- sum(arCoefs*yhat[c(i-1):c(i-arMax)]) }
        if( ewmaMax>0 ){
          for(k in 1:NCOL(mEwmaHat)){
            mEwmaHat[i,k] <- mean( yhat[c(i-ewmaEval[k]):c(i-1)] )
          }
          ewmaTerm <- sum( coefs[ewmaIndx] * mEwmaHat[i,] )
        }
        yhat[i] <- mconst + arTerm + ewmaTerm + mxreghat[i]
      } #end loop
  
      ##out:
      yhat <- yhat[c(yhat.n-n.ahead+1):yhat.n]
  
      ##simulate mY?:
      if( !is.null(mEpsilon) ){
        
        ##loop on j:
        mY <- matrix(NA, NROW(mEpsilon), NCOL(mEpsilon))
        for(j in 1:NCOL(mEpsilon)){

          ##prepare prediction no. j:
          yhatadj <- rep(NA, n.ahead + backcastMax)
          if(backcastMax>0) {
            ##actual y-values:
            yhatadj[1:backcastMax] <- tail(object$aux$y, n=backcastMax)
          }

          ##prepare ewma:
          if(ewmaMax>0){
            mEwmaHat <- matrix(NA, n.ahead+backcastMax, length(ewmaCoefs))
            colnames(mEwmaHat) <- object$aux$mXnames[ewmaIndx]
            mEwmaHat[1:backcastMax,] <- object$aux$mX[c(NROW(object$aux$mX)-backcastMax+1):NROW(object$aux$mX),ewmaIndx]
            mEwmaHat <- as.matrix(mEwmaHat)
          }

          ##predict yhatadj:
          arTerm <- 0
          ewmaTerm <- 0
          for(i in c(backcastMax+1):yhat.n){
            if( arMax>0 ){ arTerm <- sum(arCoefs*yhatadj[c(i-1):c(i-arMax)]) }
            if( ewmaMax>0 ){
              for(k in 1:NCOL(mEwmaHat)){
                mEwmaHat[i,k] <- mean( yhatadj[c(i-ewmaEval[k]):c(i-1)] )
              }
              ewmaTerm <- sum( coefs[ewmaIndx] * mEwmaHat[i,] )
            }
            yhatadj[i] <- mconst + arTerm + ewmaTerm +
              mxreghat[i] + mEpsilon[c(i-backcastMax),j]
          } #end loop
  
          ##store the simulation of yhatadj:
          mY[,j] <- yhatadj[c(yhat.n-n.ahead+1):yhat.n]

        } #end for(j)

        ##add colnames:
        colnames(mY) <- paste0("mY.", seq(1,n.sim))

      } #end if( not null(mEpsilon) )
    
    } #end if(specMean==TRUE)

  } #end if(predictMean)


  ##-----------------------
  ## 4 probs (quantiles)
  ##-----------------------

  ##mean:
  mMeanQs <- NULL
  if( predictMean && !is.null(probsArg) ){
    mMeanQs <- matrix(NA, n.ahead, length(probsArg))
    for(i in 1:NROW(mY)){
      mMeanQs[i,] <- quantile(mY[i,], probs=probsArg, type=quantile.type)
    }
    colnames(mMeanQs) <- paste0(probsArg)
  }
   
  ##variance:
  mVarianceQs <- NULL
  if( predictVariance && specVar && !is.null(probsArg) ){
    mVarianceQs <- matrix(NA, n.ahead, length(probsArg))
    for(i in 1:NROW(mSd2Hat)){
      mVarianceQs[i,] <- quantile(mSd2Hat[i,], probs=probsArg, type=quantile.type)
    }
    colnames(mVarianceQs) <- paste0(probsArg)
  }


  ##-----------------------
  ## 5 newindex
  ##-----------------------
  
  ##in-sample:
  yInSample <- zoo(object$aux$y, order.by=object$aux$y.index)

  #newindex user-provided:
  if( !is.null(newindex) ){
    yAsRegular <- FALSE
    if( n.ahead!=length(newindex) ){
      stop("length(newindex) must equal 'n.ahead'")
    }
    newindexInSample <- any( newindex %in% object$aux$y.index )
  }else{ newindexInSample <- FALSE }

  #in-sample index regular:
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

  ##neither user-provided nor regular:
  if( is.null(newindex) ){ newindex <- 1:n.ahead }

  ##add index to results:
  if(!is.null(mZhat)){ mZhat <- zoo(mZhat, order.by=newindex) }
  if(!is.null(sd2hat)){ sd2hat <- zoo(sd2hat, order.by=newindex) }
  if(!is.null(mEpsilon)){ mEpsilon <- zoo(mEpsilon, order.by=newindex) }
  if(!is.null(yhat)){ yhat <- zoo(yhat, order.by=newindex) }
  if(!is.null(mY)){ mY <- zoo(mY, order.by=newindex) }
  if(!is.null(mMeanQs)){ mMeanQs <- zoo(mMeanQs, order.by=newindex) }
  if(!is.null(mVarianceQs)){ mVarianceQs <- zoo(mVarianceQs, order.by=newindex) }

  
  ##-----------------------
  ## 6 if plot=TRUE
  ##-----------------------

  ##check special case:
  if( plotArg && spec=="variance" && specVar==FALSE ){
    message("Set 'vc = TRUE' to enable a plot of the variance predictions")
    plotArg <- FALSE #change argument
  }

  ##check another special case (plot.zoo does not work if
  ##the index is not unique):
  if( plotArg && newindexInSample==TRUE ){
    message("'newindex' not entirely out-of-sample, so no plot produced")
    plotArg <- FALSE #change argument
  }
  
  ##idea?: the special case where the out-of-sample index
  ##is not of the same type as the in-sample index. for 
  ##regular zoo-series, this can possibly be checked by checking
  ##whether the numeric delta is the same both in-sample and
  ##out-of-sample.
  
  ##way to check for this is 
  ##plot?:
  if( plotArg ){

    ##some of the plot.options:
    ##-------------------------

    ##remember: both ylab and hlines argument can be specified,
    ##even though they do not appear below
                
    ##how many in-sample observations to include in plot?:
    if( is.null(plot.options$keep) ){ plot.options$keep <- 12L }
    if( plot.options$keep < 1 ){
      plot.options$keep <- 1L
      message("'plot.options$keep' changed to 1")
    }
    
    ##if "main" argument:
    if( is.null(plot.options$main) ){
      parMarVals <- c(2.1,3.1,0.6,0.6) #bottom,left,top,right
    }else{
      parMarVals <- c(2.1,3.1,1.5,0.6) #bottom,left,top,right
    }
        
    ##linetype (solid=1, dashed=2, etc.). Order: lty=c(Forecast,Actual)
    if( is.null(plot.options$lty )){ plot.options$lty <- c(1,1) }

    ##linewidth. Order: lwd=c(Forecast,Actual)
    if( is.null(plot.options$lwd) ){ plot.options$lwd <- c(1,1) }
    
    ##colours. Order: col=c(Forecast,Actual)
    if( is.null(plot.options$col) ){ plot.options$col <- c("red","blue") }
    
    ##text for legend:
    if( is.null(plot.options$legend.text) ){
      if( spec == "variance" ){
        plot.options$legend.text <- c("Forecast", "Squared residuals")
      }else{
        plot.options$legend.text <- c("Forecast", "Actual")
      }
    }
    
    ##whether to include retained fitted or not:
    if( is.null(plot.options$fitted) ){ plot.options$fitted <- FALSE }

    ##should predictions start at origin?:
    if( is.null(plot.options$start.at.origin) ){
      plot.options$start.at.origin <- TRUE
    }

    ##add dot at forecast origin?:
    if( is.null(plot.options$dot.at.origin) ){    
      plot.options$dot.at.origin <- TRUE
    }
    
    ##add vertical line at forecast origin?:
    if( is.null(plot.options$line.at.origin) ){
      plot.options$line.at.origin <- FALSE
    }
      
    ##check if( !is.null(shades.of.grey) ):
    if( !is.null(plot.options[["shades.of.grey"]]) ){
      message(
        "argument 'shades.of.grey' has changed name to 'shades',\n",
        "'shades.of.grey' will be deprecated in future versions"
      )
      if( is.null(plot.options[["shades"]]) ){
        plot.options[["shades"]] <- plot.options[["shades.of.grey"]]
      }
    }


    ##start preparing:
    ##----------------
    
    ##select the shades of grey for the ci's:
    if( is.null(plot.options$shades) ){
      shadesOfGrey <- 40:90 #1 to 100 is possible
      shadesOfGrey <- quantile(shadesOfGrey, probs=ciLevelsArg) 
      shadesOfGrey <- shadesOfGrey[length(shadesOfGrey):1] #invert, i.e. last to first
      plot.options$shades <- round(as.numeric(shadesOfGrey))
    }
    greySelection <- paste0("grey", plot.options$shades)
    
    ##make dataForPlot:
    dataForPlot <- matrix(NA, n.ahead, 6)
    colnames(dataForPlot) <- c("MeanActual", "MeanFitted",
      "MeanPrediction", "ResidualsSquared", "VarianceFitted",
      "VariancePrediction")
    if(!is.null(plot.options$newmactual)){
      dataForPlot[1:length(plot.options$newmactual),"MeanActual"] <-
        plot.options$newmactual
    }
    if(!is.null(plot.options$newvactual)){
      dataForPlot[1:length(plot.options$newvactual),"ResidualsSquared"] <-
        plot.options$newvactual
    }
    if(!is.null(yhat)){
      dataForPlot[,"MeanPrediction"] <- coredata(yhat)
    }
    if(!is.null(sd2hat)){
      dataForPlot[,"VariancePrediction"] <- coredata(sd2hat)
    }
    retainedData <- matrix(NA, plot.options$keep, NCOL(dataForPlot))
    colnames(retainedData) <- colnames(dataForPlot)
    retainedData[,"MeanActual"] <-
      tail(coredata(yInSample), n=plot.options$keep)
    retainedData[,"ResidualsSquared"] <-
      tail(coredata(object$residuals)^2, n=plot.options$keep)
    retainedData[,"MeanFitted"] <-
      tail(coredata(object$mean.fit), n=plot.options$keep)
    retainedData[,"VarianceFitted"] <-
      tail(coredata(object$var.fit), n=plot.options$keep)
    if( plot.options$start.at.origin ){ ##let predictions start.at.origin:      
      retainedData[NROW(retainedData),"MeanPrediction"] <- 
        retainedData[NROW(retainedData),"MeanActual"]
      retainedData[NROW(retainedData),"VariancePrediction"] <- 
        retainedData[NROW(retainedData),"VarianceFitted"]
    }
    if( !plot.options$start.at.origin && plot.options$fitted ){
      retainedData[,"MeanPrediction"] <- retainedData[,"MeanFitted"]
    }
    dataForPlot <- rbind(retainedData, dataForPlot)
    tmpIndx <- c(tail(index(yInSample), n=plot.options$keep), newindex)
    dataForPlot <- zoo(dataForPlot, order.by=tmpIndx)
        
    ##if spec="mean" or "both":
    ##-------------------------

    if( spec %in% c("mean","both") ){

      ##create polygon index:
      i1 <- ifelse(plot.options$start.at.origin, 0, 1)   
      polygonIndx <- c(NROW(dataForPlot)-n.ahead+i1):NROW(dataForPlot)
      polygonIndx <- c(polygonIndx, polygonIndx[c(length(polygonIndx):1)])
      polygonIndx <- index(dataForPlot)[polygonIndx]

      ##matrices with the ci's:
      mCiLowerValsMean <- cbind(coredata(mMeanQs[,as.character(ciLower)]))
      colnames(mCiLowerValsMean) <- as.character(ciLower)
      mCiUpperValsMean <- cbind(coredata(mMeanQs[,as.character(ciUpper)]))
      colnames(mCiUpperValsMean) <- as.character(ciUpper)
      mCiUpperValsMean <-
        mCiUpperValsMean[NROW(mCiUpperValsMean):1,] #invert (first to last, last to first)
      if(n.ahead==1){ ##ensure they are still matrices:
        mCiLowerValsMean <- rbind(mCiLowerValsMean)
        mCiUpperValsMean <- rbind(mCiUpperValsMean)
      }
            
      ##add actual value at forecast origin to ci matrices?:
      if( plot.options$start.at.origin ){
        actualValue <- retainedData[NROW(retainedData),"MeanActual"]  
        mCiLowerValsMean <- rbind(actualValue,mCiLowerValsMean)
        mCiUpperValsMean <- rbind(mCiUpperValsMean,actualValue)
      }
      
      ##y-axis (limits):
      if( is.null(plot.options$ylim) ){

        ylimArg <- c(coredata(mCiLowerValsMean[,1]),
          coredata(mCiUpperValsMean[,1]))
        if( plot.options$keep > 0 ){
          ylimArg <- c(ylimArg,
            tail(coredata(yInSample), n=plot.options$keep) )
          if(plot.options$fitted){
            ylimArg <- c(ylimArg,
              tail(coredata(object$mean.fit), n=plot.options$keep))
          }
        }
        if(!is.null(plot.options$newmactual)){
          ylimArg <- c(ylimArg, coredata(plot.options$newmactual))
        }
        ylimArg <- range(ylimArg)
        eps <- abs(ylimArg[2]-ylimArg[1])
        ylimArg[2] <- ylimArg[2] + eps*0.15 #add more space at the top
        ylimArg[1] <- ylimArg[1] - eps*0.05 #add more space at the bottom

      }else{ ylimArg <- plot.options$ylim }
           
      ##get current par-values:
      def.par <- par(no.readonly=TRUE)
  
      ##margins:
      par(mar=parMarVals) 
  
      ##plot actual values in white (i.e. create plot):
      plot.zoo(dataForPlot[,"MeanActual"], xlab="", ylab="",
        main=plot.options$main, lty=plot.options$lty[2],
        col="white", lwd=plot.options$lwd[2], ylim=ylimArg)

      ##add start line?:
      if( plot.options$line.at.origin ){
        startlineIndx <- rep( index(dataForPlot)[plot.options$keep], 2)
        eps <- abs(ylimArg[2]-ylimArg[1])
        startlineVals <- c(ylimArg[1]-eps*0.05, ylimArg[2]/1.2)
        polygon(startlineIndx, startlineVals, col="grey",
          border="grey", lwd=plot.options$lwd)
      }
      
      ##add ci's:
      for(i in 1:length(ciLevelsArg)){
        polygon( polygonIndx,
          c(mCiLowerValsMean[,i],mCiUpperValsMean[,i]),
          col=greySelection[i], border=greySelection[i] )  
      }
                  
      ##add horisontal lines?:
      if(!is.null(plot.options$hlines)){
        abline(h=plot.options$hlines, col="grey", lty=3)
      }

      ##add prediction:
      lines(dataForPlot[,"MeanPrediction"], lty=plot.options$lty[1],
        col=plot.options$col[1], lwd=plot.options$lwd[1])
    
      ##add actual:
      lines(dataForPlot[,"MeanActual"], lty=plot.options$lty[2],
        col=plot.options$col[2], lwd=plot.options$lwd[2], type="l")
        
      ##add fitted (pre-prediction):
      if( plot.options$keep > 0 && plot.options$fitted ){
        lines(dataForPlot[,"MeanFitted"], lty=plot.options$lty[2],
          col=plot.options$col[1], lwd=plot.options$lwd[1],
          type="l")
      }

      ##add point at forecast origin?:
      if( plot.options$dot.at.origin ){
        points(index(dataForPlot)[NROW(retainedData)],
          retainedData[NROW(retainedData),"MeanActual"],
          pch=19, col=plot.options$col[2], lwd=plot.options$lwd[2])
      }
      
      ##add actual values out-of-sample:
      if( !is.null(plot.options$newmactual) ){
        lines(dataForPlot[,"MeanActual"], lty=plot.options$lty[2],
          col=plot.options$col[2], lwd=plot.options$lwd[2],
          type="l")
      }

      ##add text closer to plot than xlab or ylab would do
      mtextValue <- ifelse(is.null(plot.options$ylab),
        "Mean", plot.options$ylab)
      mtext(mtextValue, side=2, line=2)
  
      ##add plot-legend:
      legend("top", lty=plot.options$lty, col=plot.options$col,
        lwd=plot.options$lwd, legend=plot.options$legend.text,
        bg="white", bty="n")
  
      ##add ci-legend:
      legendArg <- ciLevelsArg[length(ciLevelsArg):1]*100
      legendArg <- paste0(legendArg, "%")      
      legend("topright", lty=c(1,1), lwd=13, bty="n",
        col=greySelection[c(1,length(ciLevelsArg))],
        legend=legendArg[c(1,length(ciLevelsArg))])
      
      ##return to old par-values:
      par(def.par)

    } #end if(spec %in% c("mean","both"))


    ##if spec="variance":
    ##-------------------

    if( spec=="variance" ){

      ##create polygon index:
      i1 <- ifelse(plot.options$start.at.origin, 0, 1)   
      polygonIndx <- c(NROW(dataForPlot)-n.ahead+i1):NROW(dataForPlot)
      polygonIndx <- c(polygonIndx, polygonIndx[c(length(polygonIndx):1)])
      polygonIndx <- index(dataForPlot)[polygonIndx]

      ##matrices with the ci's:
      mCiLowerValsVar <- cbind(coredata(mVarianceQs[,as.character(ciLower)]))
      colnames(mCiLowerValsVar) <- as.character(ciLower)
      mCiUpperValsVar <- cbind(coredata(mVarianceQs[,as.character(ciUpper)]))
      colnames(mCiUpperValsVar) <- as.character(ciUpper)
      mCiUpperValsVar <-
        mCiUpperValsVar[NROW(mCiUpperValsVar):1,] #invert (first to last, last to first)

      ##add fitted value to forecast origin?:
      if( plot.options$start.at.origin ){
        fittedValue <- retainedData[NROW(retainedData),"VarianceFitted"]  
        mCiLowerValsVar <- rbind(fittedValue,mCiLowerValsVar)
        mCiUpperValsVar <- rbind(mCiUpperValsVar,fittedValue)
      }

      ##y-axis (limits):
      if(is.null(plot.options$ylim)){

        ylimArg <- c(coredata(mCiLowerValsVar[,1]),
          coredata(mCiUpperValsVar[,1]))
        if( plot.options$keep > 0 ){
          ylimArg <- c(ylimArg,
            tail(coredata(object$residuals)^2, n=plot.options$keep) )
          if(plot.options$fitted){
            ylimArg <- c(ylimArg,
              tail(coredata(object$var.fit), n=plot.options$keep))
          }
        }
        if(!is.null(plot.options$newvactual)){
          ylimArg <- c(ylimArg, coredata(plot.options$newvactual))
        }
        ylimArg <- range(ylimArg)
        eps <- abs(ylimArg[2]-ylimArg[1])
        ylimArg[2] <- ylimArg[2] + eps*0.15 #add more space at the top
        ylimArg[1] <- ylimArg[1] - eps*0.05 #add more space at the bottom
     
      }else{ ylimArg <- plot.options$ylim }
      
      ##get current par-values:
      def.par <- par(no.readonly=TRUE)
  
      ##margins:
      par(mar=parMarVals) 
  
      ##plot the actual values:
      plot.zoo(dataForPlot[,"ResidualsSquared"], xlab="", ylab="",
        main=plot.options$main, lty=plot.options$lty[2],
        col=plot.options$col[2], lwd=plot.options$lwd,
        ylim=ylimArg)
    
      ##add start line?:
      if( plot.options$line.at.origin ){
        startlineIndx <- rep( index(dataForPlot)[plot.options$keep], 2)
        eps <- abs(ylimArg[2]-ylimArg[1])
        startlineVals <- c(ylimArg[1]-eps*0.05, ylimArg[2]/1.2)
        polygon(startlineIndx, startlineVals, col="grey",
          border="grey", lwd=plot.options$lwd)
      }
  
      ##add ci's:
      for(i in 1:length(ciLevelsArg)){
        polygon( polygonIndx,
          c(mCiLowerValsVar[,i],mCiUpperValsVar[,i]),
          col=greySelection[i], border=greySelection[i] )  
      }

      ##add horisontal lines?:
      if(!is.null(plot.options$hlines)){
        abline(h=plot.options$hlines, col="grey", lty=3)
      }
  
      ##add prediction:
      lines(dataForPlot[,"VariancePrediction"], lty=plot.options$lty[1],
        col=plot.options$col[1], lwd=plot.options$lwd,
        type="l")
    
      ##add fitted (in-sample):
      if( plot.options$keep > 0 && plot.options$fitted ){
        lines(dataForPlot[,"VarianceFitted"], lty=plot.options$lty[2],
          lwd=plot.options$lwd, col=plot.options$col[1],
          type="l")
      }
  
      ##add point at forecast origin?:
      if(plot.options$dot.at.origin){
        points(index(dataForPlot)[NROW(retainedData)],
          retainedData[NROW(retainedData),"VarianceFitted"], 
          #OLD: fittedValue,
          pch=19, col=plot.options$col[1], lwd=plot.options$lwd)
      }

      ##add actual values of residuals squared out-of-sample:
      if( !is.null(plot.options$newvactual) ){
        lines(dataForPlot[,"ResidualsSquared"], lty=plot.options$lty[2],
          col=plot.options$col[2], lwd=plot.options$lwd,
          type="l")
      }

      ##add text closer to plot than xlab or ylab would do
      mtextValue <- ifelse(is.null(plot.options$ylab),
        "Variance", plot.options$ylab)
      mtext(mtextValue, side=2, line=2)

      ##add plot-legend:
      legend("top", lty=plot.options$lty, col=plot.options$col,
        lwd=plot.options$lwd, legend=plot.options$legend.text,
        bty="n")
      
      ##add ci-legend:
      legendArg <- ciLevelsArg[length(ciLevelsArg):1]*100
      legendArg <- paste0(legendArg, "%")      
      legend("topright", lty=c(1,1), lwd=13, bty="n",
        col=greySelection[c(1,length(ciLevelsArg))],
        legend=legendArg[c(1,length(ciLevelsArg))])
      
      ##return to old par-values:
      par(def.par)

    } #end if(spec %in% c("mean","both"))

  } #end if(plotArg)

      
  ##-----------------------
  ## 7 if return=TRUE
  ##-----------------------

  if(return){

    ##change colnames on quantiles:
    if(!is.null(mMeanQs)){ colnames(mMeanQs) <- paste0("y",probsArg) }
    if(!is.null(mVarianceQs)){ colnames(mVarianceQs) <- paste0("sd2",probsArg) }
    
    ##return everything:
    if(verbose){
      result <- NULL
      if(!is.null(yhat)){ result <- cbind(yhat) }
      if(!is.null(mMeanQs)){
        if(is.null(result)){ result <- mMeanQs }else{ result <- cbind(result,mMeanQs) }
      }
      if(!is.null(mY)){
        if(is.null(result)){ result <- mY }else{ result <- cbind(result,mY) }
      }
      if(!is.null(sd2hat)){
        if(is.null(result)){ result <- sd2hat }else{ result <- cbind(result,sd2hat) }
      }
      if(!is.null(mVarianceQs)){
        if(is.null(result)){ result <- mVarianceQs }else{ result <- cbind(result,mVarianceQs) }
      }
      if(!is.null(mEpsilon)){
        if(is.null(result)){ result <- mEpsilon }else{ result <- cbind(result,mEpsilon) }
      }
      if(!is.null(mZhat)){
        if(is.null(result)){ result <- mZhat }else{ result <- cbind(result, mZhat) }
      }
    } #end if(verbose)
    
    ##do not return everything:
    if(!verbose){

      resultMean <- NULL
      resultVariance <- NULL
      
      ##mean specification:
      if( spec %in% c("mean","both" ) ){
        resultMean <- yhat
        if( !is.null(probs) || !is.null(ci.levels) ){
          resultMean <- cbind(yhat,mMeanQs)
        }
      }

      ##mean specification:
      if( spec %in% c("variance","both" ) ){
        resultVariance <- sd2hat
        if( !is.null(probs) || !is.null(ci.levels) ){
          resultVariance <- cbind(sd2hat,mVarianceQs)
        }
      }

      ##combine:
      if(is.null(resultMean)){ result <- resultVariance }
      if(is.null(resultVariance)){ result <- resultMean }
      if(!is.null(resultMean) && !is.null(resultVariance) ){
        result <- cbind(resultMean,resultVariance)
        colnames(result) <- c("yhat", "sd2hat")
      }
          
    } #end if(!verbose)

    ##return the result:
    return(result)

  } #end if(return)
  
}
