sigma.isat <-
function(object, ...)
{
  if(is.null(object$residuals)){
    result <- NULL
  }else{
    RSS <- sum(object$residuals^2)
    result <- sqrt(RSS/(object$n - object$k))
  }
  return(result)
}
