isatvarcorrect <-
function(x,   mcor = 1){
  
  if (class(x)=="isat"){
    if (x$call$iis==TRUE) 
    {
      x$vcov.mean <- x$vcov.mean * as.numeric(isvarcor(x$aux$t.pval, 1)[2]^2)
      x$vcov.mean[x$keep, x$keep] <- x$vcov.mean[x$keep, x$keep] * as.numeric(isvareffcor(x$aux$t.pval, 1, mcor)[2]^2)
      
      x$mean.results$std.error <- sqrt(diag(x$vcov.mean))
      x$mean.results$`t-stat` <- x$mean.results$coef/x$mean.results$std.error
      x$mean.results$`p-value` <-  pt(abs(x$mean.results$`t-stat`), x$df, lower.tail=FALSE)*2   
      
      x$sigma2 <- x$sigma2*as.numeric(isvarcor(x$aux$t.pval, 1)[2]^2)
      x$logl <- -x$n*log(2*x$sigma2*pi)/2 - x$rss/(2*x$sigma2)
      
      return(x)
    } else {
      stop("iis not TRUE")
    }
  } else {
    stop("x must be an isat object")
  }  
  
}
