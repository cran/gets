terminals <-
function(object, ...)
{
  if(class(object)=="gets"){
    return(object$terminals)
  }else{
    cat("The object does not belong to the 'gets' class\n")
  }
}
