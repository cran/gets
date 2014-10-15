regs.mean <-
function(y, mc=NULL, ar=NULL, ewma=NULL, mxreg=NULL)
{
y.n <- length(y)
mX <- NULL

#constant
if(!is.null(mc) && as.numeric(mc)==1){
  mX <- cbind(mX,rep(1,y.n)); colnames(mX) <- "mconst"
}

#ar terms:
arlags <- NULL
if(!is.null(ar)){
  for(i in 1:length(ar)){
    arlags <- cbind(arlags, glag(y, k=ar[i]))
  }
  colnames(arlags) <- paste("ar", ar, sep="")
}
mX <- cbind(mX, arlags)

#EqWMA term:
if(is.null(ewma)){EqWMA <- NULL}else{
  EqWMA <- do.call(eqwma, c(list(y),ewma) )
}
mX <- cbind(mX, EqWMA)

#create matrix of mean regressors mxreg:
if(!is.null(mxreg)){
  mxregnames <- colnames(mxreg)
  if(is.null(mxregnames)){
    mxregnames <- paste("mxreg", 1:ncol(mxreg), sep="")
  }
  if(any(mxregnames == "")){
    missing.colnames <- which(mxregnames == "")
    for(i in 1:length(missing.colnames)){
      mxregnames[i] <- paste("mxreg", i, sep="")
    }
  }
  colnames(mxreg) <- mxregnames
}
mX <- cbind(mX, mxreg)

#out-matrix:
out <- mX

return(out)

}
