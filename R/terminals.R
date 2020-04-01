terminals <-
function(object, ...)
{
  if(class(object)=="gets" || class(object)=="isat"){
    return(object$terminals)
  }else{
    cat("object not of class 'gets' or 'isat' class\\n")
  }
}
