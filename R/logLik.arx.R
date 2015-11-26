logLik.arx <-
function(object, ...)
{
  resids <- residuals.arx(object, std=FALSE)
  sdhat <- fitted.arx(object, spec="variance")
  sdhat <- na.trim(sdhat)
  resids <- window(resids, start=index(sdhat)[1],
    end=index(sdhat)[length(sdhat)])
  result <- sum(dnorm(resids, sd=sdhat, log=TRUE))
  attr(result, "df") <- length(coef.arx(object))
  attr(result, "nobs") <- length(sdhat)
  class(result) <- "logLik"
  return(result)
}
