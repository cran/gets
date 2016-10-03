isvarcor <-
function(t.pval, sigma) {
  
  alpha <- t.pval
  sigmals <- sigma
  
  c <- abs(qnorm(alpha/2))
  
  psi <- pnorm(c) - pnorm(-c) 
  tau <- psi - 2*c*dnorm(c)
  
  xi_sq <- tau/psi
  
  xi <- sqrt(xi_sq)
  corrxi <- 1/xi
  
  sigmacorr <- sigmals * corrxi
  
  object <-  data.frame(cbind(sigmacorr, corrxi))
  names(object) <- c("sigma.cor", "corxi")
  return(object)
}
