coef.arx <-
function(object, spec=c("mean","variance","both"),
  ...)
{
  spec <- match.arg(spec)
  if(spec=="mean"){
    result <- object$mean.results[,1]
    names(result) <- rownames(object$mean.results)
  }
  if(spec=="variance"){
    result <- object$variance.results[,1]
    names(result) <- rownames(object$variance.results)
    result <- c(result,object$Elnz2)
    names(result)[length(result)] <- "Elnz2"
  }
  if(spec=="both"){
    result1 <- object$mean.results[,1]
    names(result1) <- rownames(object$mean.results)
    result2 <- object$variance.results[,1]
    names(result2) <- rownames(object$variance.results)
    result2 <- c(result2,object$Elnz2)
    names(result2)[length(result2)] <- "Elnz2"
    result <- c(result1,result2)
  }
  return(result)
}
