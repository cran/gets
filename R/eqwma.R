eqwma <-
function(x, length=5, k=1, p=1, abs=FALSE, log=FALSE,
  as.vector=FALSE, lag=NULL, start=NULL)
{
  ##deprecated arguments:
  if(!is.null(start)){ stop("'start' has been deprecated") }
  if(!is.null(lag)){ stop("'lag' has been deprecated, use 'k' instead") }

  ##zoo related:
  if(is.zoo(x)){
    isZoo <- TRUE
    xindex <- index(x)
    x <- coredata(x)
  }else{
    isZoo <- FALSE
  }

  ##is x a matrix?
  if(NCOL(x)>1){
    stop("'x' must be one-dimensional")
  }else{
    x <- as.vector(x)
  }

  ##na's, power p, absolute value:
  xn <- length(x)
  x <- na.trim(x, sides="left")
  xnadj <- length(x)
  if(xn > xnadj){nachk <- TRUE}else{nachk <- FALSE}
  if(abs){xabs <- abs(x)}else{xabs <- x}
  if(p!=1){xabsp <- xabs^p}else{xabsp <- xabs}

  ##compute eqwma:
  result <- NULL
  for(i in 1:length(length)){
    movavg <- rollmean(xabsp, length[i], fill=NA, align="r")
    result <- cbind(result,movavg)
  }

  ##lag?:
  if(k>0){
    result <- cbind(result[ 1:c(NROW(result)-k), ]) #lag
    result <- rbind( matrix(NA,k,NCOL(result)), result) #add NAs
  }
  
  #apply log?:
  if(log){
    result <- log(result)
    colnames(result) <- paste0("logEqWMA(", length, ")")
  }else{
    colnames(result) <- paste0("EqWMA(", length, ")")
  }
  if(nachk){ result <- rbind( matrix(NA,c(xn-xnadj),NCOL(result)), result) }
  if(as.vector && NCOL(result)==1 ){ result <- as.vector(result) }
  if(isZoo){ result <- zoo(result, order.by=xindex) }

  return(result)
}
