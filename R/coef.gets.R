coef.gets <-
function(object, spec=NULL,
  ...)
{
  gets.type <- object$gets.type
  if(is.null(spec)){
    spec <- switch(gets.type, getsm="mean", getsv="variance",
      isat="mean")
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  }

  #if(getsm):
  if(gets.type=="getsm" || gets.type=="isat"){

    #mean:
    result1 <- NULL
    if(is.null(spec) || spec=="mean" || spec=="both"){
      if( !is.null(object$specific.mean)
        && object$specific.mean!="empty" ){
        result1 <- object$specific.mean[,"coef"]
        names(result1) <- rownames(object$specific.mean)
      }
    } #end is.null(spec)

    #variance:
    result2 <- NULL
    if(spec=="variance" || spec=="both"){
      if(!is.null(object$specific.variance)){
        result2 <- c(object$specific.variance[,"coef"],
          object$Elnz2)
        names(result2) <- c(rownames(object$specific.variance),
          "Elnz2")
      } #end if(!is.null(..))
    }

    result <- c(result1,result2)
  } #end if("getsm")

  #if(getsv):
  if(gets.type=="getsv"){

    #mean:
    result1 <- NULL

    #variance:
    result2 <- NULL
    if(is.null(spec) || spec=="variance" || spec=="both"){
      if(!is.null(object$specific.variance)){
        result2 <- c(object$specific.variance[,"coef"],
          object$Elnz2)
        names(result2) <- c(rownames(object$specific.variance),
          "Elnz2")
      }
    } #end is.null(spec)

    result <- c(result1,result2)
  } #end if("getsv")

  return(result)
}
