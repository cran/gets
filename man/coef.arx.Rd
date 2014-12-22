\name{coef.arx}
\alias{coef.arx}
\alias{fitted.arx}
\alias{plot.arx}
\alias{print.arx}
\alias{residuals.arx}
\alias{summary.arx}
\alias{vcov.arx}

\title{Extraction functions for 'arx' objects}

\description{Extraction functions for objects of class 'arx'}

\usage{
\method{coef}{arx}(object, spec=NULL, ...)
\method{fitted}{arx}(object, spec=NULL, ...)
\method{plot}{arx}(x, y=NULL, col=c("red","blue"), lty=c("solid","solid"),
  lwd=c(1,1), spec=NULL, ...)
\method{print}{arx}(x, ...)
\method{residuals}{arx}(object, std=FALSE, ...)
\method{summary}{arx}(object, ...)
\method{vcov}{arx}(object, spec=NULL, ...)
}

\arguments{
  \item{object}{an object of class 'arx'}
  \item{x}{an object of class 'arx'}
  \item{spec}{NULL, "mean", "variance" or, except \code{vcov.arx}, "both". If NULL, then it is automatically determined whether information pertaining to the mean or variance specification should be returned. This can be overridden by either "mean", "variance" or, except for \code{vcov.arx}, "both"}
  \item{std}{logical. If FALSE (default), then the mean residuals are returned. If TRUE, then the standardised residuals are returned}
  \item{y}{unused}
  \item{col}{colours of fitted (default=red) and actual (default=blue) lines}
  \item{lty}{types of fitted (default=solid) and actual (default=solid) lines}
  \item{lwd}{widths of fitted (default=1) and actual (default=1) lines}
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{a numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with fitted values}
  \item{plot:}{a plot of the fitted values and the residuals}
  \item{print:}{a print of the estimation results}
  \item{residuals:}{a \code{\link{zoo}} object with the residuals}
  \item{summary:}{a print of the items in the \code{\link{arx}} object}
  \item{vcov:}{variance-covariance matrix}
}

\author{
James Reade,  \url{https://sites.google.com/site/jjamesreade/},\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{arx}}}

\examples{
##simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 100)

##simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*100), 100, 4)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean, and log-ARCH(3) in the variance:
mymod <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:3)

##print results:
print(mymod)

##plot the fitted vs. actual values, and the residuals:
plot(mymod)

##print the entries of object 'mymod':
summary(mymod)

##extract coefficient estimates (automatically determined):
coef(mymod)

##extract mean coefficients only:
coef(mymod, spec="mean")

##extract log-variance coefficients only:
coef(mymod, spec="variance")

##extract all coefficient estimates:
coef(mymod, spec="both")

##extract variance-covariance matrix of mean equation:
vcov(mymod)

##extract variance-covariance matrix of log-variance equation:
vcov(mymod, spec="variance")

##extract and plot the fitted mean values (automatically determined):
mfit <- fitted(mymod)
plot(mfit)

##extract and plot the fitted variance values:
vfit <- fitted(mymod, spec="variance")
plot(vfit)

##extract and plot both the fitted mean and variance values:
vfit <- fitted(mymod, spec="both")
plot(vfit)

##extract and plot the fitted mean values:
vfit <- fitted(mymod, spec="mean")
plot(vfit)

##extract and plot residuals:
epshat <- residuals(mymod)
plot(epshat)

##extract and plot standardised residuals:
zhat <- residuals(mymod, std=TRUE)
plot(zhat)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
