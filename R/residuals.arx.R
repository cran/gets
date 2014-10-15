residuals.arx <-
function(object, std=FALSE, ...)
{
  if(std){
    result <- object$resids.std
  }else{
    result <- object$resids
  }
  return(result)
}
