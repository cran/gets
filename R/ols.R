ols <-
function(y, x, tol=1e-07 , LAPACK=FALSE,
  method=1, user.fun=NULL, user.options=NULL)
{

  ##to do number 1:
  ## - rename ols to estFun
  ## - merge user.fun and user.options into a single argument,
  ## user.estimator, which is a list containing at least one
  ## entry, name, i.e. the name of the estimator-function
  ##
  ##to do number 2:
  ## - split estFun into two functions, estFun and vcovFun

  ##user-specified:
  if(method==0){
    user.options <- c(list(y=y, x=x), user.options)
    out <- do.call(user.fun, user.options, envir=.GlobalEnv)
  }

  ##fastest (usually only for estimates):
  if(method==1){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK)
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
  }

  ##second fastest:
  if(method==2){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
    out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
    out$fit <- as.vector(x %*% out$coefficients)
    out$resids <- y - out$fit
#OLD:
#    out$residuals <- y - out$fit
  }

  ##ordinary vcov:
  if(method==3){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
    out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
    out$fit <- as.vector(x %*% out$coefficients)
    out$resids <- y - out$fit
#OLD:
#    out$residuals <- y - out$fit
    out$resids2 <- out$resids^2
    out$rss <- sum(out$resids2)
    out$n <- length(y)
    out$df <- out$n - NCOL(x)
    out$sigma2 <- out$rss/out$df
    out$vcov <- out$sigma2 * out$xtxinv
  }

  ##white vcov:
  if(method==4){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
    out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
    out$fit <- as.vector(x %*% out$coefficients)
    out$resids <- y - out$fit
#OLD:
#    out$residuals <- y - out$fit
    out$resids2 <- out$resids^2
    out$rss <- sum(out$resids2)
    out$n <- length(y)
    out$df <- out$n - NCOL(x)
    out$sigma2 <- out$rss/out$df
    out$omegahat <- crossprod(x, x*out$resids2)
    out$vcov <- out$xtxinv %*% out$omegahat %*% out$xtxinv
  }

  ##newey-west vcov:
  if(method==5){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
    out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
    out$fit <- as.vector(x %*% out$coefficients)
    out$resids <- y - out$fit
#OLD:
#    out$residuals <- y - out$fit
    out$resids2 <- out$resids^2
    out$rss <- sum(out$resids2)
    out$n <- length(y)
    out$df <- out$n - NCOL(x)
    out$sigma2 <- out$rss/out$df

    y.n <- length(y)
    iL <- round(y.n^(1/4), digits=0)
    vW <- 1 - 1:iL/(iL+1)
    vWsqrt <- sqrt(vW)
    mXadj <- out$resids*x
    mS0 <- crossprod(mXadj)

    mSum <- 0
    for(l in 1:iL){
      mXadjw <- mXadj*vWsqrt[l]
      mXadjwNo1 <- mXadjw[-c(1:l),]
      mXadjwNo2 <- mXadjw[-c(c(y.n-l+1):y.n),]
      mSum <- mSum + crossprod(mXadjwNo1, mXadjwNo2) + crossprod(mXadjwNo2, mXadjwNo1)
    }

    out$omegahat <- mS0 + mSum
    out$vcov <- out$xtxinv %*% out$omegahat %*% out$xtxinv
  }

  ##result:
  return(out)

}
