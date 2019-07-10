vargaugeiis <-
function(t.pval, T, infty=FALSE, m=1){
  
  alpha <- t.pval
  c <- abs(qnorm(alpha/2))
  fc <- dnorm(c)
  psi <- pnorm(c) - pnorm(-c) 
  tau <- psi - 2*c*dnorm(c)
  chi <- 3*psi - 2*c*(c^2+3)*fc
  rho_sig <- (c^2 - tau/psi)*c*fc/tau
  k <- 3
  
  ###m iterations
  m_min1 <- m-1
  eta_sigma_m_min1 <- ( ((1-rho_sig^m_min1)/((1-rho_sig)*tau))^2 + 2* ((1-rho_sig^m_min1)/((1-rho_sig)*tau)  )*rho_sig^m_min1  ) *(chi-tau^2/psi)/(k-1)+rho_sig^(2*m_min1)
  eta_gamma_m <- (c*fc)^2*eta_sigma_m_min1*(k-1) + 2*c*fc*rho_sig^(m-1)*(tau-psi)
  
  ###infinite iterations
  eta_infty <- (chi-tau^2/psi)*(c*fc)^2/((1-rho_sig)*tau)^2
  
  if (infty==TRUE){
    viis <- psi*(1-psi) + eta_infty  
  } else {
    viis <- psi*(1-psi) + eta_gamma_m
  }
  
  sdiis <- sqrt(viis)
  viis_T <- viis/T
  sdiis_T <- sqrt(viis/T)
  
  out <- data.frame(cbind(viis_T, sdiis_T, viis, sdiis))
  names(out) <- c("var_iisgauge", "sd_iisgauge", "asy_var_iisgauge", "asy_sd_iisgauge")
  return(out)
}
