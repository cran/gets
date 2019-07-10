isatloop <-
function(num=c(seq(from=20, to=1, by=-1)), t.pval.spec = FALSE, print=FALSE, y, ar=NULL, iis=TRUE,  sis=FALSE,  ...){
  
  
  initialm <- arx(y=y, ar=ar) 
  n <- initialm$n
  num <- num[rev(order(num))]
  
  if (t.pval.spec == FALSE){
    p.num <- num/n # scaling significance levels
  } else { #if p-values are pre-specified
    p.num <- num
  }
  
  K <- length(p.num)
  obs.gauge <- data.frame(matrix(NA, nrow=K, ncol=4)) # record significance p value, sample gauge, expected number and sample number of outliers
  names(obs.gauge) <- c("t.pval", "obs.prop", "p.num", "obs.num")
  for (k in 1:K) 
  {
    pval <- p.num[k]
    if (print==TRUE)
    {
      print(paste("k:", k, "/", K, ", p:", pval, sep=""))
    }
    x <- isat(y=y, iis=iis, sis=sis, t.pval=pval, ar=ar, ...) 
    
    #mxreg=mxreg, mc=mc, 
    if (!is.null(x$ISnames)) {
      ISnames <- c(x$ISnames[grep("iis", x$ISnames)])
      noutl= length(ISnames) 
      is.gauge <- noutl/n
      is.num <- noutl
    }  else   {
      is.gauge <- 0
      is.num <- 0
    }
    if (print==TRUE)
    {
      print(paste("gauge k:", is.gauge, ", number k:", "is.num", sep=""))     
    }
    obs.gauge$t.pval[k] <- pval 
    obs.gauge$obs.prop[k] <- is.gauge
    obs.gauge$p.num[k] <- pval*n
    obs.gauge$obs.num[k] <- is.num    
  }
  
  out <- list(n, obs.gauge)
  names(out) <- c("n", "obs.gauge")
  return(out)
}
