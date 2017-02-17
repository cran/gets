logLik.arx <-
function(object, ...)
{
  if( is.null(object$aux$user.estimator) ){

    resids <- residuals.arx(object, std=FALSE)
    sdhat <- sqrt(fitted.arx(object, spec="variance"))
    sdhat <- na.trim(sdhat)
    resids <- window(resids, start=index(sdhat)[1],
      end=index(sdhat)[length(sdhat)])
    result <- sum(dnorm(resids, sd=sdhat, log=TRUE))
    attr(result, "df") <- length(coef.arx(object, spec="mean"))
    attr(result, "nobs") <- length(sdhat)
    class(result) <- "logLik"

  }else{

    result <- object$logl
    if( !is.null(result) ){
      attr(result, "df") <- length(object$coefficients)
      if( !is.null(objects$resids) ){
        attr(result, "nobs") <- length(na.trim(objects$resids))
      }
      class(result) <- "logLik"
    }
    
  }
  
  return(result)
}
