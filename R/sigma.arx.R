sigma.arx <-
function(object, ...)
{
  residsTrimmed <- na.trim(object$residuals)
  RSS <- sum(residsTrimmed^2)
  nobs <- length(residsTrimmed)
  DFs <- length(coef.arx(object, spec="mean"))
  return( sqrt(RSS/(nobs-DFs)) )
}
