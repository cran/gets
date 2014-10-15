ols <-
function(y, x, tol = 1e-07 , LAPACK = FALSE,
  method=1)
{
if(method==1){
  out <- list()
  qx <- qr(x, tol, LAPACK = LAPACK)
  out <- c(out, qx)
  out$coefficients <- as.vector(solve.qr(qx, y, tol = tol))
}
if(method==2){
  tmp <- crossprod(x)
  out <- qr(tmp)
  out$xtxinv <- solve.qr(out, tol=tol, LAPACK=LAPACK)
  out$xtx <- tmp
  out$xty <- crossprod(x,y)
  out$coefficients <- as.vector(out$xtxinv%*%out$xty)
}
return(out)
}
