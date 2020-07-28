regressorsVariance <-
function(e, vc=TRUE, arch=NULL, asym=NULL,
  log.ewma=NULL, vxreg=NULL, zero.adj=0.1, vc.adj=TRUE,
  return.regressand=TRUE, return.as.zoo=TRUE, na.trim=TRUE,
  na.omit=FALSE)
{

  ##regressand:
  if(is.zoo(e)){ e <- cbind(e) }else{ e <- as.zoo(cbind(e)) }
  if(NCOL(e) > 1) stop("Dependent variable not 1-dimensional")
  e.n <- NROW(e)
  loge2.index <- index(e)
  e <- coredata(e)
  t1 <- loge2.index[1]
  t2 <- loge2.index[e.n]
  zero.where <- which(e==0)
  eabs <- abs(e)
  if( length(zero.where)>0 ){
    eabs[zero.where] <- quantile(eabs[-zero.where], zero.adj, na.rm=TRUE)
  }
  loge2 <- log(eabs^2)

  ##create regressor matrix:
  vX <- NULL
  vXnames <- NULL

  ##variance intercept:
  if( identical(as.numeric(vc),1) ){
    vX <- cbind(rep(1,e.n))
    vXnames <- "vconst"
  }

  ##arch terms:
  if(!is.null(arch) && !identical(as.numeric(arch),0) ){
    tmp <- NULL
    nas <- rep(NA, max(arch))
    tmpfun <- function(i){
      tmp <<- cbind(tmp, c(nas[1:i],loge2[1:c(e.n-i)]))
    }
    tmpfun <- sapply(arch,tmpfun)
    vX <- cbind(vX, tmp)
    vXnames <- c(vXnames, paste0("arch", arch))
  }

  ##asym terms:
  if(!is.null(asym) && !identical(as.numeric(asym),0) ){
    tmp <- NULL
    nas <- rep(NA, max(asym))
    tmpfun <- function(i){
      tmp <<- cbind(tmp, c(nas[1:i],
        loge2[1:c(e.n-i)]*as.numeric(e[1:c(e.n-i)]<0)))
    }
    tmpfun <- sapply(asym,tmpfun)
    vX <- cbind(vX, tmp)
    vXnames <- c(vXnames, paste0("asym", asym))
  }

  ##log.ewma term:
  if(!is.null(log.ewma)){
    if(is.list(log.ewma)){
      log.ewma$k <- 1
    }else{
      log.ewma <- list(length=log.ewma)
    }
    tmp <- do.call(leqwma, c(list(e),log.ewma) )
    vXnames <- c(vXnames, colnames(tmp))
    colnames(tmp) <- NULL
    vX <- cbind(vX, tmp)
  }

  ##trim for NAs:
  if(na.trim){
    tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
    tmp <- na.trim(tmp, sides="both", is.na="any")
    loge2.n <- NROW(tmp)
    loge2.index <- index(tmp) #re-define index
    t1 <- loge2.index[1] #re-define t1
    t2 <- loge2.index[loge2.n] #re-define t2
    loge2 <- tmp[,1]
    loge2 <- coredata(loge2)
    if(!is.null(vX)){ #re-define vX
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }
  }

  ##vxreg:
  if(!is.null(vxreg)){
    vxreg <- as.zoo(cbind(vxreg))
    vxreg.names <- colnames(vxreg)
    if(is.null(vxreg.names)){
      vxreg.names <- paste0("vxreg", 1:NCOL(vxreg))
    }
    if(any(vxreg.names == "")){
      missing.colnames <- which(vxreg.names == "")
      for(i in 1:length(missing.colnames)){
        vxreg.names[missing.colnames[i]] <- paste0("vxreg", i)
      }
    }
    vXnames <- c(vXnames, vxreg.names)
    vxreg <- window(vxreg, start=t1, end=t2)
    vxreg <- cbind(coredata(vxreg))
    vX <- cbind(vX, vxreg)
    colnames(vxreg) <- NULL

    ##re-trim for NAs:
    if(na.trim){
      tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
      tmp <- na.trim(tmp, sides="both", is.na="any")
      loge2.n <- NROW(tmp)
      loge2.index <- index(tmp) #re-define index
      t1 <- loge2.index[1] #re-define t1
      t2 <- loge2.index[loge2.n] #re-define t2
      loge2 <- coredata(tmp[,1])
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }

  } #end if(!is.null(vxreg))

  ##remove rows with NAs:
  if(na.omit){
    tmp <- zoo(cbind(loge2,vX), order.by=loge2.index)
    tmp <- na.omit(tmp)
    loge2.n <- NROW(tmp) #re-define
    loge2.index <- index(tmp) #re-define
    t1 <- loge2.index[1] #re-define t1
    t2 <- loge2.index[loge2.n] #re-define t2
    loge2 <- coredata(tmp[,1]) #re-define
    if(!is.null(vX)){ #re-define vX
      vX <- tmp[,2:NCOL(tmp)]
      vX <- coredata(vX)
      vX <- cbind(vX)
      colnames(vX) <- NULL
    }
  }

  ### OUTPUT: ######################

  if(return.regressand){
    result <- cbind(loge2, vX)
    colnames(result) <- c("loge2", vXnames)
  }else{
    result <- vX
    if(!is.null(result)){ colnames(result) <- vXnames }
  }
  if(return.as.zoo && !is.null(result) ){ result <- zoo(result, order.by=loge2.index) }
  return(result)

}
