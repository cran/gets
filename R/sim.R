sim <-
function(x, which.ones=NULL)
{
  if(NROW(x)==1){
    x.is.scalar <- TRUE
    n <- x
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which(1:n %in% which.ones)
    }
  }else{
    x.is.scalar <- FALSE
    n <- NROW(x)
    x <- as.zoo(x)
    x.index <- index(x)
    if(is.null(which.ones)){
      where.indicators <- 2:n
    }else{
      where.indicators <- which(x.index %in% which.ones)
    }
  }
  n.where.indicators <- length(where.indicators)
  loop.indx <- 1:n.where.indicators
  mSIS <-matrix(0,n,n.where.indicators)
  tmp <- function(i){ mSIS[ c(where.indicators[i]:n) ,i] <<- 1 }
  tmp <- sapply(loop.indx,tmp)
  if(x.is.scalar){
    colnames(mSIS) <- paste("sis", where.indicators, sep="")
    mSIS <- as.zoo(mSIS)
  }else{
    colnames(mSIS) <- paste("sis",
      as.character(x.index)[where.indicators], sep="")
    mSIS <- zoo(mSIS, order.by=x.index)
  }
  return(mSIS)
}
