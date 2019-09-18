paths <-
function(object, ...)
{
  if(class(object)=="gets" || class(object)=="isat"){
    return(object$paths)
  }else{
    stop("object not of class 'gets' or 'isat' class")
  }
}
