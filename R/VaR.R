VaR <-
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

  ##fitted mean and sd, standardised residuals, quantile:
  meanFit <- fitted(object, spec="mean")
  sdFit <- sqrt( fitted(object, spec="variance") )
  residsStd <- residuals(object, std=TRUE)
  qValue <- quantile(residsStd, probs=riskLevel, type=type,
    names=FALSE, na.rm=TRUE)
  if( length(riskLevel)==1 ){
    VaR <- meanFit + sdFit*qValue
  }else{
    colNames <- paste("VaR", level, sep="")
    VaR <- matrix(NA, length(sdFit), length(colNames))
    for(i in 1:length(colNames)){
      VaR[,i] <- meanFit + sdFit*qValue[i]
    }
    colnames(VaR) <- colNames
    VaR <- zoo(VaR, order.by=index(sdFit))
  }

  ##return
  return(-VaR)
}
