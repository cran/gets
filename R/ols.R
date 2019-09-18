ols <-
function(y, x, untransformed.residuals=NULL, tol=1e-07,
  LAPACK=FALSE, method=3, ...)
{

  ##for the future:
  ## - rename ols to estFun? Split estFun into two functions,
  ## estFun and vcovFun?

  ##user-specified:
  if(method==0){
    stop("method=0 has been deprecated")
  }

  ##fastest (usually only for estimates):
  if(method==1){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK)
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
  }

  ##second fastest (slightly more output):
  if(method==2){
    out <- list()
    qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
    out <- c(out, qx)
    out$coefficients <- solve.qr(qx, y, tol=tol)
    out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
    out$fit <- as.vector(x %*% out$coefficients)
    out$residuals <- y - out$fit
  }

  ##ordinary vcov:
  if(method==3){
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
      out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df
    if(out$k>0){
      out$vcov <- out$sigma2 * out$xtxinv
    }
    out$logl <- -out$n*log(2*out$sigma2*pi)/2 - out$rss/(2*out$sigma2)
  }

  ##white (1980) vcov:
  if(method==4){
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
      out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df
    if(out$k>0){
      out$omegahat <- crossprod(x, x*out$residuals2)
      out$vcov <- out$xtxinv %*% out$omegahat %*% out$xtxinv
    }
    out$logl <- -out$n*log(2*out$sigma2*pi)/2 - out$rss/(2*out$sigma2)
  }

  ##newey-west(1987) vcov:
  if(method==5){
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k>0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
      out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df

    if(out$k>0){
      iL <- round(out$n^(1/4), digits=0)
      vW <- 1 - 1:iL/(iL+1)
      vWsqrt <- sqrt(vW)
      mXadj <- out$residuals*x
      mS0 <- crossprod(mXadj)

      mSum <- 0
      for(l in 1:iL){
        mXadjw <- mXadj*vWsqrt[l]
        mXadjwNo1 <- mXadjw[-c(1:l),]
        mXadjwNo2 <- mXadjw[-c(c(out$n-l+1):out$n),]
        mSum <- mSum + crossprod(mXadjwNo1, mXadjwNo2) + crossprod(mXadjwNo2, mXadjwNo1)
      }

      out$omegahat <- mS0 + mSum
      out$vcov <- out$xtxinv %*% out$omegahat %*% out$xtxinv
    } #end if(out$k>0)
    
    out$logl <- -out$n*log(2*out$sigma2*pi)/2 - out$rss/(2*out$sigma2)
  }

  ##log-variance w/ordinary vcov (note: y = log(e^2)):
  if(method==6){
    out <- list()
    out$n <- length(y)
    if(is.null(x)){ out$k <- 0 }else{ out$k <- NCOL(x) }
    out$df <- out$n - out$k
    if(out$k > 0){
      qx <- qr(x, tol, LAPACK=LAPACK) ## compute qr-decomposition of x
      out <- c(out, qx)
      out$coefficients <- solve.qr(qx, y, tol=tol)
      out$xtxinv <- chol2inv(qx$qr, LINPACK=FALSE) #(x'x)^-1
      out$fit <- as.vector(x %*% out$coefficients)
    }else{
      out$fit <- rep(0, out$n)
    }
    out$residuals <- y - out$fit #residuals of AR-X representation
    out$residuals2 <- out$residuals^2
    out$rss <- sum(out$residuals2)
    out$sigma2 <- out$rss/out$df
    if(out$k>0){
      out$vcov <- out$sigma2 * out$xtxinv
    }
    ##log-variance part:
    out$Elnz2 <- -log(mean(exp(out$residuals)))
    out$var.fit <- exp(out$fit - out$Elnz2)
    out$std.residuals <- untransformed.residuals/sqrt(out$var.fit)
    out$logl <- -out$n*log(2*pi)/2 - sum(log(out$var.fit))/2 - sum(untransformed.residuals^2/out$var.fit)/2
  }

  ##result:
  return(out)

}
