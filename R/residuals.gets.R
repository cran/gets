residuals.gets <-
function(object, std=NULL, ...)
{
  #determine type:
  if(is.null(std)){
    std <- switch(object$gets.type, getsm=FALSE, getsv=TRUE,
      isat=FALSE)
#OLD:
#    std <- switch(as.character(object$call)[1],
#      getsm=FALSE, getsv=TRUE)
  }

  if(std){
    result <- object$resids.std
  }else{
    result <- object$resids
  }
  return(result)
}
