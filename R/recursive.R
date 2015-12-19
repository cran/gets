recursive <-
function(object, spec="mean",
  std.errors=TRUE, from=40, tol=1e-07, LAPACK=FALSE,
  plot=TRUE, return=TRUE)
{
  ##which specification:
  specType <- c("mean", "variance", "both")
  whichType <- charmatch(spec, specType)
  if(whichType==3){ stop("Sorry, 'both' not possible.") }
  spec <- specType[whichType]

  ##if mean-specification:
  if(spec=="mean"){
    if(is.null(object$mean.results)){
      stop("No mean-equation")
    }
    vY <- object$aux$y
    yNrow <- NROW(vY)
    mX <- object$aux$mX
    mXncol <- object$aux$mXncol
    mXnames <- object$aux$mXnames
    mainlab <- "Recursive estimates: mean equation"
  }

  ##if variance-specification:
  if(spec=="variance"){
    if(is.null(object$variance.results)){
      stop("No variance-equation")
    }
    vY <- object$aux$loge2
    yNrow <- NROW(vY)
    mX <- object$aux$vX
    mXncol <- object$aux$vXncol
    mXnames <- object$aux$vXnames
    mainlab <- "Recursive estimates: log-variance equation"
  }

  ##determine ols method:
  if(spec=="mean"){
    if(std.errors){
      if(object$aux$vcov.type=="ordinary"){ olsMethod=3 }
      if(object$aux$vcov.type=="white"){ olsMethod=4 }
      if(object$aux$vcov.type=="newey-west"){ olsMethod=5 }
    }else{
      olsMethod=1
    }
  } #close if(mean)

  if(spec=="variance"){
    if(std.errors){
      olsMethod=3
    }else{
      olsMethod=2
    }
  } #close if(variance)

  ##initialise:
  colnames(mX) <- mXnames
  recursiveEstimates <- matrix(NA, yNrow, mXncol)
  colnames(recursiveEstimates) <- mXnames
  if(std.errors){
    recursiveStdErrs <- recursiveEstimates
  }
  startIndx <- max(mXncol, min(from, yNrow))
  compute.at <- seq.int(from=yNrow, to=startIndx, by=-1)

  ##recursion:
  for(i in 1:length(compute.at)){

    ##estimate:
    vY <- vY[1:compute.at[i]]
    mXnames <- colnames(mX)
    NCOLmX <- NCOL(mX)
    mX <- dropvar(as.matrix(mX[1:compute.at[i], ]), tol=tol,
      LAPACK=LAPACK, silent=TRUE)
    if(NCOLmX==1){ colnames(mX) <- mXnames }
    tmpEst <- ols(vY, mX, tol=tol, LAPACK=LAPACK,
      method=olsMethod)
    recursiveEstimates[compute.at[i],colnames(mX)] <- tmpEst$coefficients
    if(std.errors){
      recursiveStdErrs[compute.at[i],colnames(mX)] <- sqrt(diag(tmpEst$vcov))
    }

    ##if variance-specification:
    if(spec=="variance"){
      Elnz2est <- -log(mean(exp(tmpEst$residuals)))
      recursiveEstimates[compute.at[i], "vconst"] <- Elnz2est + recursiveEstimates[compute.at[i], "vconst"]
    }

  } #close for loop

  ##rename std.errors columns:
  if(std.errors){
    colnames(recursiveStdErrs) <- paste(colnames(recursiveStdErrs),
      "SE", sep="")
  }

  ##set vconstSE to NA:
  if(std.errors==TRUE && spec=="variance"){
    recursiveStdErrs[,1] <- NA
  }

  ##handle zoo-index:
  naDiff <- object$aux$y.n - yNrow
  if(naDiff==0){
    zooIndx <- object$aux$y.index
  }else{
    zooIndx <- object$aux$y.index[-c(1:naDiff)]
  }
  recursiveEstimates <- zoo(recursiveEstimates,
    order.by=zooIndx)
  if(std.errors){
    recursiveStdErrs <- zoo(recursiveStdErrs,
      order.by=zooIndx)
  }

  ##if return=TRUE:
  if(return){
    out <- list()
    out$estimates <- recursiveEstimates
    if(std.errors){
      out$standard.errors <- recursiveStdErrs
    }
  }

  ##if plot=TRUE:
  if(plot){
    recursiveEstimates <- na.trim(recursiveEstimates,
      is.na="all")
    plot(recursiveEstimates, main=mainlab, xlab="", col="blue")
  }

  ##out:
  if(return){
    return(out)
  }
}
