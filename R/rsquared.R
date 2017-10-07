rsquared <-
function(object, ...)
{
  classOK <- class(object) %in% c("arx", "gets", "isat")
  if(!classOK){ stop("Object not of class 'arx', 'gets' nor 'isat'") }
  if( class(object) == "gets" ){
    specType <- switch(as.character(object$call)[1],
      getsm="mean", getsv="variance")
  }
  if( class(object) == "gets" && specType=="variance" ){
    Rsquared <- NA
  }else{
    TSS <- sum( (object$aux$y - mean(object$aux$y))^2 )
    residsTrimmed <- na.trim(object$residuals)
    RSS <- sum(residsTrimmed^2)
    Rsquared <- 1 - RSS/TSS
    ##to do: adjusted R-squared
  }
  return(Rsquared)
}
