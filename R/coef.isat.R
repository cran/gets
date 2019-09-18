coef.isat <-
function(object, ...)
{
  result <- object$coefficients
  if(!is.null(result)){ names(result) <- names(object$specific.spec) }
  return(result)
}
