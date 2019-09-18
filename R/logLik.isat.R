logLik.isat <-
function(object, ...)
{
  result <- object$logl
  if(is.null(result)){
    result <- numeric(0)
    warning("'object$logl' is NULL")
  }else{
    attr(result, "df") <- length(object$coefficients)
    attr(result, "nobs") <- object$n
  }
  class(result) <- "logLik"
  return(result)
}
