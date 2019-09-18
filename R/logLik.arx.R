logLik.arx <-
function(object, ...)
{
  ## in the future: add a df.method argument with
  ## optional values "mean-coefficients" (default),
  ## "variance-coefficients" and "both"?
  
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
