\name{coef.gets}
\alias{coef.gets}
\alias{fitted.gets}
\alias{logLik.gets}
\alias{plot.gets}
\alias{print.gets}
\alias{residuals.gets}
\alias{summary.gets}
\alias{vcov.gets}

\title{Extraction functions for 'gets' objects}

\description{Extraction functions for objects of class 'gets'}

\usage{
\method{coef}{gets}(object, spec=NULL, ...)
\method{fitted}{gets}(object, spec=NULL, ...)
\method{logLik}{gets}(object, ...)
\method{plot}{gets}(x, spec=NULL, col=c("red","blue"),
  lty=c("solid","solid"), lwd=c(1,1), ...)
\method{print}{gets}(x, ...)
\method{residuals}{gets}(object, std=NULL, ...)
\method{summary}{gets}(object, ...)
\method{vcov}{gets}(object, spec=NULL, ...)
}

\arguments{
  \item{object}{an object of class 'gets'}
  \item{x}{an object of class 'gets'}
  \item{spec}{NULL, "mean", "variance" or, in some instances, "both". When \code{NULL} is a valid value, then it is automatically determined whether information pertaining to the mean or variance specification should be returned}
  \item{std}{logical. If FALSE (default), then the mean residuals are returned. If TRUE, then the standardised residuals are returned}
  \item{col}{colours of fitted (default=red) and actual (default=blue) lines}
  \item{lty}{types of fitted (default=solid) and actual (default=solid) lines}
  \item{lwd}{widths of fitted (default=1) and actual (default=1) lines}
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{a numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with fitted values}
  \item{logLik:}{a numeric, the log-likelihood (normal density)}
  \item{plot:}{a plot of the fitted values and the residuals}
  \item{print:}{a print of the estimation results}
  \item{residuals:}{a \code{\link{zoo}} object with the residuals}
  \item{summary:}{a print of the items in the \code{\link{gets}} object}
  \item{vcov:}{a variance-covariance matrix}
}

\author{
James Reade, \url{https://sites.google.com/site/jjamesreade/},\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{ \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}} }

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 100)

##Simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*100), 100, 4)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean, and a log-ARCH(3) in the variance:
mymod <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:3)

##General-to-Specific (GETS) model selection of the mean:
meanmod <- getsm(mymod)

##General-to-Specific (GETS) model selection of the variance:
varmod <- getsv(mymod)

##print results:
print(meanmod)
print(varmod)

##plot the fitted vs. actual values, and the residuals:
plot(meanmod)
plot(varmod)

##print the entries of object 'gets':
summary(meanmod)
summary(varmod)

##extract coefficients of the simplified (specific) model:
coef(meanmod) #mean spec
coef(varmod) #variance spec

##extract log-likelihood:
logLik(mymod)

##extract variance-covariance matrix of simplified
##(specific) model:
vcov(meanmod) #mean spec
vcov(varmod) #variance spec

##extract and plot the fitted values:
mfit <- fitted(meanmod) #mean fit
plot(mfit)
vfit <- fitted(varmod) #variance fit
plot(vfit)

##extract and plot residuals:
epshat <- residuals(meanmod)
plot(epshat)

##extract and plot standardised residuals:
zhat <- residuals(varmod)
plot(zhat)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
