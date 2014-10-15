vcov.gets <-
function(object, spec=NULL,  ...)
{
  #determine type:
  if(is.null(spec)){
    if(as.character(object$call)[1]=="getsm"){ spec <- "mean" }
    if(as.character(object$call)[1]=="getsv"){ spec <- "variance" }
  }else{
    spec.type <- c("mean", "variance")
#    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  }

  #vcov:
  if(spec=="mean"){
    result <- object$vcov.mean
  }
  if(spec=="variance"){
    result <- object$vcov.var
  }
#  if(spec=="both"){
#    result <- list(vcov.mean=object$vcov.mean,
#      vcov.var=object$vcov.var)
#  }
  return(result)
}
