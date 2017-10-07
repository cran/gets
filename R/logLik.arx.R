logLik.arx <-
function(object, ...)
{
  ## in the future: add a df.method argument with
  ## optional values "mean-coefficients" (default),
  ## "variance-coefficients" and "both"?
  
  result <- object$logl
  attr(result, "df") <- length(object$coefficients)
  attr(result, "nobs") <- object$n
  class(result) <- "logLik"
  return(result)

#OLD:
#  if( is.null(object$aux$user.estimator) ){
#
#    resids <- residuals.arx(object, std=FALSE)
#    sdhat <- sqrt(fitted.arx(object, spec="variance"))
#    sdhat <- na.trim(sdhat)
#    resids <- window(resids, start=index(sdhat)[1],
#      end=index(sdhat)[length(sdhat)])
#    result <- sum(dnorm(resids, sd=sdhat, log=TRUE))
#    attr(result, "df") <- length(coef.arx(object, spec="mean"))
#    attr(result, "nobs") <- length(sdhat)
#    class(result) <- "logLik"
#
#  }else{
#
#    result <- object$logl
#    if( !is.null(result) ){
#      attr(result, "df") <- length(object$coefficients)
#      if( "residuals" %in% names(object) ){
##OLD:
##      if( !is.null(objects$residuals) ){
#        attr(result, "nobs") <- length(na.trim(object$residuals))
#      }
#      class(result) <- "logLik"
#    }
#
#  }
#
#  return(result)

}
