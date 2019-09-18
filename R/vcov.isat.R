vcov.isat <-
function(object, ...)
{
  result <- object$vcov #also works if object$vcov.mean exists
  if(!is.null(result)){
    if(is.null(colnames(result))){
      colnames(result) <- names(object$specific.spec)
    }
    if(is.null(rownames(result))){
      rownames(result) <- names(object$specific.spec)
    }
  }
  return(result)  
}
