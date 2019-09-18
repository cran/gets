vcov.arx <-
function(object, spec=NULL, ...)
{

  ##spec argument:
  specOriginal <- spec
  if(is.null(spec)){
    if(!is.null(object$mean.results)){
      spec <- "mean"
    }
    if(is.null(object$mean.results)
      && !is.null(object$variance.results) ){
      spec <- "variance"
    }
  }else{
    spec.type <- c("mean", "variance")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  }

  ##create result:
  result <- NULL
  
  ##if mean:
  if(spec=="mean"){
    result <- object$vcov.mean
  }

  ##if variance:
  if(spec=="variance"){
    result <- object$vcov.var
  }

#  ##check and change if 0 x 0?:
#  if(all(dim(result)==0)){ result <- NULL }

  ##if user-specified estimator:
  if( !is.null(object$aux$user.estimator) && is.null(specOriginal) ){
    result <- object$vcov
  }
  
  return(result)
}
