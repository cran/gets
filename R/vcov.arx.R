vcov.arx <-
function(object, spec=c("mean","variance"),
  ...)
{
  spec <- match.arg(spec)
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
