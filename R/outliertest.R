outliertest <-
function(x=NULL, noutl=NULL, t.pval=NULL, T=NULL,  m=1, infty=FALSE, alternative="two.sided"){
  
  # noutl=x$obs.gauge$obs.num[i]
  # t.pval=x$obs.gauge$t.pval[i]
  # T=x$n
  # 
  
  if (!is.null(x)){
    
    if (class(x)=="isat"){
      
      if (any(x$call$sis, x$call$tis, x$call$uis)==TRUE){
        stop("Test only valid for iis")
      } else {
        if (  x$call$iis == TRUE){
          ISnames <- c(x$ISnames[grep("iis", x$ISnames)])
          noutl= length(ISnames) 
          t.pval <- x$aux$t.pval
          T <- x$n
        } else {
          stop("iis must be TRUE")
        }
      }
      
      
    } else {
      stop("x must be an isat object")
    }
  } #x null closed
  
  gauge.null <- t.pval
  gauge <- noutl/T
  
  #### Standard Normal Test
  gauge.sd <- vargaugeiis(t.pval=gauge.null, T=T, infty=infty, m=m)$sd_iisgauge 
  gaugetestval <- (gauge - gauge.null)/gauge.sd
  
  if (alternative=="two.sided"){
    pval <- 2*(pnorm(-abs(gaugetestval)))
  }
  if (alternative=="less"){
    pval <- pnorm(gaugetestval)  
  }
  if (alternative=="greater"){
    pval <- 1-  pnorm(gaugetestval) 
  }
  
  ### Poisson Test
  gauge_n <- noutl
  gauge_null_n <- gauge.null*T
  poistestval <- poisson.test(gauge_n, r=gauge_null_n, alternative=alternative)
  
  rval_norm <- list(statistic = gaugetestval, p.value = pval, estimate=gauge, null.value = gauge.null, alternative = alternative, method="Jiao-Pretis Outlier Proportion Test", data.name="Proportion of detected outliers")
  attr(rval_norm, "class") <- "htest"
  rval_pois <- list(statistic = gauge_n, p.value = poistestval$p.value, estimate=gauge_n, null.value = gauge_null_n, alternative = alternative, method="Jiao-Pretis Outlier Count Test", data.name="Number of detected outliers")
  attr(rval_pois, "class") <- "htest"
  out <- list(rval_norm, rval_pois)
  names(out) <- c("proportion", "count")
  
  return(out)
}
