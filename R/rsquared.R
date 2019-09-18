rsquared <-
function(object, adjusted=FALSE, ...)
{
  classOK <- class(object) %in% c("arx", "gets", "isat")
  if(!classOK){ message("object not of class 'arx', 'gets' or 'isat'") }
  if( class(object) == "gets" ){
    specType <- switch(as.character(object$call)[1],
      getsm="mean", getsv="variance")
  }
  if( class(object) == "gets" && specType=="variance" ){
    result <- NA
#OLD:
#    Rsquared <- NA
  }else{
    TSS <- sum( (object$aux$y - mean(object$aux$y))^2 )
    residsTrimmed <- na.trim(object$residuals)
    RSS <- sum(residsTrimmed^2)
    Rsquared <- 1 - RSS/TSS
    if(adjusted){
      result <- 1 - (1-Rsquared)*(object$n-1)/(object$n-object$k)
    }else{
      result <- Rsquared
    }
  }
  return(result)
}
