terminals <-
function(object, ...)
{
  if(class(object)=="gets" || class(object)=="isat"){
    return(object$terminals)
  }else{
    stop("object not of class 'gets' or 'isat'")
  }
}
