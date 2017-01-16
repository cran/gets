rsquared <-
function(object, ...)
{
  classOK <- class(object) %in% c("arx", "gets", "isat")
  if(!classOK){ stop("Object not of class 'arx', 'gets' nor 'isat'") }
  TSS <- sum( (object$aux$y - mean(object$aux$y))^2 )
  residsTrimmed <- na.trim(object$resids)
  RSS <- sum(residsTrimmed^2)
  result <- 1 - RSS/TSS
  ##to do: adjusted R-squared
  return(result)
}
