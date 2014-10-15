residuals.gets <-
function(object, std=NULL, ...)
{
  #determine type:
  if(is.null(std)){
    if(as.character(object$call)[1]=="getsm"){ std <- FALSE }
    if(as.character(object$call)[1]=="getsv"){ std <- TRUE }
  }#else{
#    std.type <- c(FALSE, TRUE)
#    which.type <- charmatch(spec, spec.type)
#    spec <- spec.type[which.type]
#  }

  if(std){
    result <- object$resids.std
  }else{
    result <- object$resids
  }
  return(result)
}
