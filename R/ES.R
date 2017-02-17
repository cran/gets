ES <-
function(object, level=0.99, type=7, ...)
{
  ##check whether class is valid:
  classType <- class(object)
  if( !classType %in% c("arx", "gets", "isat") ){
    stop("object not of class 'arx', 'gets' nor 'isat'")
  }

  ##check the risk-levels:
  riskLevel <- 1-level
  if( any(riskLevel > 1) || any(riskLevel < 0) ){
    stop("risk-level(s) must be in the 0 to 1 interval")
  }

  ##fitted sd, standardised residuals, quantile:
  meanFit <- fitted(object, spec="mean")
  sdFit <- sqrt( fitted(object, spec="variance") )
  residsStd <- residuals(object, std=TRUE)
  qValue <- quantile(residsStd, probs=riskLevel, type=type,
    names=FALSE, na.rm=TRUE)
  colNames <- paste("ES", level, sep="")
  mExpShortF <- matrix(NA, length(sdFit), length(colNames))
  for(i in 1:length(colNames)){
    whereExceeds <- which( residsStd < qValue[i] )
    if( length(whereExceeds) == 0 ){
      stop("no standardised residual smaller than ", qValue[i])
    }else{
      ExpShortF <- mean( residsStd[whereExceeds] )
    }
    mExpShortF[,i] <- sdFit*ExpShortF
  }
  colnames(mExpShortF) <- colNames
  if(NCOL(mExpShortF)==1){ mExpShortF <- as.vector(mExpShortF) }
  mExpShortF <- zoo(mExpShortF, order.by=index(sdFit))
  mExpShortF <- meanFit + mExpShortF

  ##return
  return(-mExpShortF)

}
