fitted.gets <-
function(object, spec=NULL, ...)
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

  #obtain fitted:
  if(spec=="mean"){
    result <- object$mean.fit
  }
  if(spec=="variance"){
    result <- object$var.fit
  }
#  if(spec=="both"){
#    result <- cbind(object$mean.fit, object$var.fit)
#  }

  return(result)
}
