fitted.arx <-
function(object, spec=c("mean","variance"),
  ...)
{
  spec <- match.arg(spec)
  #mean:
  if(spec=="mean"){
    result <- object$mean.fit
  }
  #variance:
  if(spec=="variance"){
    result <- object$var.fit
  }
  return(result)
}
