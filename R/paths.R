paths <-
function(object, ...)
{
  if(class(object)=="gets"){
    return(object$paths)
  }else{
    cat("The object does not belong to the 'gets' class\n")
  }
}
