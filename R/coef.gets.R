coef.gets <-
function(object, spec=NULL,
  ...)
{
  gets.type <- as.character(object$call)[1]
  if(is.null(spec)){
    if(gets.type=="getsm"){ spec <- "mean" }
    if(gets.type=="getsv"){ spec <- "variance" }
  }else{
    spec.type <- c("mean", "variance", "both")
    which.type <- charmatch(spec, spec.type)
    spec <- spec.type[which.type]
  }

  #if(getsm):
  if(gets.type=="getsm"){

    #mean:
    result1 <- NULL
    if(is.null(spec) || spec=="mean" || spec=="both"){
      if( !is.null(object$specific.mean)
        && object$specific.mean!="empty" ){
        result1 <- object$specific.mean[,2]
        names(result1) <- rownames(object$specific.mean)
      }
    } #end is.null(spec)

    #variance:
    result2 <- NULL
    if(spec=="variance" || spec=="both"){
      if(!is.null(object$specific.variance)){
        result2 <- c(object$specific.variance[,1],
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
        result2 <- c(object$specific.variance[,2],
          object$Elnz2)
        names(result2) <- c(rownames(object$specific.variance),
          "Elnz2")
      }
    } #end is.null(spec)

    result <- c(result1,result2)
  } #end if("getsv")

  return(result)
}
