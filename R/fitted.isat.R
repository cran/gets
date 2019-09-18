fitted.isat <-
function(object, ...)
{
  result <- object$mean.fit
  if(is.null(result)){ result <- object$fit }
  return(result)
}
