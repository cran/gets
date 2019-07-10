outlierscaletest <-
function(x, nsim = 10000){
  ###################
  
  obs.gauge <- x$obs.gauge
  n <- x$n ##need to extract n
  p.num <- obs.gauge$t.pval
  c <- qnorm(1 - p.num/2) # corresponding cut-offs # assume standardised error follows standard normal
  fc <- dnorm(c) 
  psi <- 1 - p.num
  tau_2_c <- psi - 2*c*fc
  tau_4 <- 3
  K <- length(p.num)
  
  #######################
  ##### Scale Sum Test
  ######################
  
  stand.obs.gauge <- n^(1/2)*(obs.gauge$obs.prop - obs.gauge$t.pval) # standardise gauge
  scalesum.stat <- sum(stand.obs.gauge) # scaling sum statistic
  
  if (K == 1) # covariance structure between test statistics of different significance levels
  {
    covar <- 0
  }  else  {
    covar <- 0 
    for (k in 1:(K - 1)) 
    {
      for (l in (k + 1):K)
      {
        covar <- covar + obs.gauge$t.pval[l] - obs.gauge$t.pval[k]*obs.gauge$t.pval[l]
      }
    }
  }
  var.scalesum.stat <- sum(obs.gauge$t.pval*(1 - obs.gauge$t.pval)) + 2*covar + sum(c*fc)^(2)*(tau_4 - 1) + 2*sum(c*fc)*sum(tau_2_c + obs.gauge$t.pval - 1) # variance for scaling sum test statistic
  sd.scalesum.stat <- sqrt(var.scalesum.stat)
  
  #### Sum Test Output
  stand.scalesum.stat <- scalesum.stat/sd.scalesum.stat
  scalesum.pval <- 2*(pnorm(-abs(stand.scalesum.stat)))
  
  #######################
  ##### Scale Sup Test
  ######################
  
  scalesup.stat <- max(abs(stand.obs.gauge)) # scaling sup statistic
  
  N <- nsim # sample size used to simulate the limiting distribution (could also be added as the argument of function)
  mu <- matrix(0, K, 1) # mean of GP
  Sigma <- matrix(NA, K, K) # covariance of GP
  for (s in 1:K)
  {
    for (t in 1:K)
    {
      if (s <= t)
      {
        Sigma[s, t] <- obs.gauge$t.pval[t]*(1 - obs.gauge$t.pval[s]) + c[s]*fc[s]*(tau_2_c[t] + obs.gauge$t.pval[t] - 1) + c[t]*fc[t]*(tau_2_c[s] + obs.gauge$t.pval[s] - 1) + c[s]*c[t]*fc[s]*fc[t]*(tau_4 - 1)
      }
      else
      {
        Sigma[s, t] <- Sigma[t, s]
      }
    }
  }
  GPsample <- mvrnormsim(N, mu, Sigma) # generate multivariate normal and dim of object is N by K
  Limitsample <- apply(abs(GPsample), 1, max) # find largest value over rows of absolute of GPsample
  
  pvalsample <- Limitsample # use empirical distribution to draw p value
  pvalsample[Limitsample <= scalesup.stat] <- 0
  pvalsample[Limitsample > scalesup.stat] <- 1
  scalesup.pval <- mean(pvalsample)
  

  ##################################
  ##### Output of Sum and Sup Tests
  
  rval_sum <- list(statistic = stand.scalesum.stat, p.value = scalesum.pval, estimate=NULL, null.value = NULL, alternative = NULL, method="Jiao-Pretis Outlier Scaling Sum Test", data.name="Scaling proportion of detected outliers (Sum)")
  attr(rval_sum, "class") <- "htest"
  rval_sup <- list(statistic = scalesup.stat, p.value = scalesup.pval, estimate=NULL, null.value = NULL, alternative = NULL, method="Jiao-Pretis Outlier Scaling Sup Test", data.name="Scaling proportion of detected outliers (Sup)")
  attr(rval_sup, "class") <- "htest"
  
  out <- list(rval_sum, rval_sup)
  names(out) <- c("sum", "sup")
  
  return(out)
  
}
