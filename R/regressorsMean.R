regressorsMean <-
function(y, mc=FALSE, ar=NULL, ewma=NULL, mxreg=NULL,
  return.regressand=TRUE, return.as.zoo=TRUE, na.trim=TRUE,
  na.omit=FALSE)
{

  ##regressand:
  y.name <- deparse(substitute(y))
  if(is.zoo(y)){ y <- cbind(y) }else{ y <- as.zoo(cbind(y)) }
  if(NCOL(y) > 1) stop("Dependent variable not 1-dimensional")
  if( is.null(y.name)){ y.name <- colnames(y)[1] }
  if( y.name[1] =="" ){ y.name <- "y" }
  y.n <- NROW(y)
  y.index <- index(y)
  y <- coredata(y)
  t1 <- y.index[1]
  t2 <- y.index[y.n]

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
    ewma$as.vector <- FALSE #force result to be a matrix
    tmp <- do.call(eqwma, c(list(y),ewma) )
    mXnames <- c(mXnames, colnames(tmp))
    colnames(tmp) <- NULL
    mX <- cbind(mX, tmp)
  }

  ##trim for NAs:
  if(na.trim){
    tmp <- zoo(cbind(y,mX), order.by=y.index)
    tmp <- na.trim(tmp, sides="both", is.na="any")
    y.n <- NROW(tmp) #re-define y.n
    y.index <- index(tmp) #re-define y.index
    t1 <- y.index[1] #re-define t1
    t2 <- y.index[y.n] #re-define t2
    y <- coredata(tmp[,1]) #re-define y
    if(!is.null(mX)){ #re-define mX
      mX <- tmp[,2:NCOL(tmp)]
      mX <- coredata(mX)
      mX <- cbind(mX)
      colnames(mX) <- NULL
    }
  }
  
  ##mxreg:
  if( !is.null(mxreg) ){
#OLD (until version 0.23):
#  if( !is.null(mxreg) && !identical(as.numeric(mxreg),0) ){
    mxreg <- as.zoo(cbind(mxreg))
    mxreg.names <- colnames(mxreg)
    if(is.null(mxreg.names)){
      mxreg.names <- paste("mxreg", 1:NCOL(mxreg), sep="")
    }
    if(any(mxreg.names == "")){
      missing.colnames <- which(mxreg.names == "")
      for(i in 1:length(missing.colnames)){
        mxreg.names[missing.colnames[i]] <- paste0("mxreg", i)
      }
    }
    #mxreg.names <- make.names(mxreg.names)
    mXnames <- c(mXnames, mxreg.names)
    mxreg <- window(mxreg, start=t1, end=t2)
    mxreg <- cbind(coredata(mxreg))
    mX <- cbind(mX, mxreg)

    ##re-trim for NAs:
    if(na.trim){
      tmp <- zoo(cbind(y,mX), order.by=y.index)
      tmp <- na.trim(tmp, sides="both", is.na="any")
      y.n <- NROW(tmp) #re-define y.n
      y.index <- index(tmp) #re-define y.index
      t1 <- y.index[1] #re-define t1
      t2 <- y.index[y.n] #re-define t2
      y <- coredata(tmp[,1])
      mX <- tmp[,2:NCOL(tmp)]
      mX <- coredata(mX)
      mX <- cbind(mX)
      colnames(mX) <- NULL
    }

  } #end if(!is.null(mxreg))

  ##remove rows with NAs:
  if(na.omit){
    tmp <- zoo(cbind(y,mX), order.by=y.index)
    tmp <- na.omit(tmp)
    y.n <- NROW(tmp) #re-define y.n
    y.index <- index(tmp) #re-define y.index
    t1 <- y.index[1] #re-define t1
    t2 <- y.index[y.n] #re-define t2
    y <- coredata(tmp[,1]) #re-define y
    if(!is.null(mX)){ #re-define mX
      mX <- tmp[,2:NCOL(tmp)]
      mX <- coredata(mX)
      mX <- cbind(mX)
      colnames(mX) <- NULL
    }
  }

  ### OUTPUT: ######################

  if(return.regressand){
    result <- cbind(y, mX)
    colnames(result) <- c(y.name, mXnames)
  }else{
    result <- mX
    if(!is.null(result)){ colnames(result) <- mXnames }
  }
  if(return.as.zoo && !is.null(result) ){ result <- zoo(result, order.by=y.index) }
  return(result)
  
}
