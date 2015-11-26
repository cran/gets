iim <-
function(x, which.ones=NULL)
{
  if(NROW(x)==1){
    n <- x
    mIIS <- matrix(0,n,n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste("iis", 1:n, sep="")
    if(!is.null(which.ones)){ mIIS <- mIIS[,which.ones] }
    mIIS <- as.zoo(mIIS)
  }else{
    n <- NROW(x)
    x <- as.zoo(x)
    x.index <- index(x)
    mIIS <- matrix(0,n,n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste("iis",
      as.character(x.index), sep="")
    mIIS <- zoo(mIIS, order.by=x.index)
    if(!is.null(which.ones)){
      where.indicators <- which(index(mIIS) %in% which.ones)
      if(length(where.indicators > 0)){
        mIIS <- cbind(mIIS[,where.indicators])
      }else{
        print("'which.ones' not in index")
      }
    }
  } #end if(NROW(x)==1)else..
  return(mIIS)
}
