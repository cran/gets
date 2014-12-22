fitted.gets <-
function(object, spec=NULL, ...)
{
  #determine type:
  if(is.null(spec)){
    spec <- switch(object$gets.type, getsm="mean",
      getsv="variance", isat="mean")
#OLD:
#    spec <- switch(as.character(object$call)[1],
#      getsm="mean", getsv="variance")
  }else{
    spec.type <- c("mean", "variance")
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
