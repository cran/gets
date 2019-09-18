residuals.isat <-
function(object, std=FALSE, ...)
{
  if(is.null(object)){
    result <- NULL
  }else{
    if(std){
      result <- object$residuals/sigma.isat(object)
    }else{
      result <- object$residuals
    }
  }
  return(result)
}
