isvareffcor <-
function(t.pval, se, m=1) {

  alpha <- t.pval

  c <- abs(qnorm(alpha/2))

  psi <- pnorm(c) - pnorm(-c)
  tau <- psi - 2*c*dnorm(c)

  rhobeta <- 2*c*dnorm(c)/psi
  etam <- (((1-rhobeta^m)/((1-rhobeta)*psi))^2 + 2 * (1-rhobeta^m)/((1-rhobeta)*psi)*rhobeta^m)*tau+rhobeta^(2*m)

  se_cor <- se*sqrt(etam)

  output <-  data.frame(cbind(se_cor, sqrt(etam) ))
  names(output) <- c("se.cor", "eta.m")
  return(output)
}
