\name{coef.arx}
\alias{coef.arx}
\alias{fitted.arx}
\alias{print.arx}
\alias{residuals.arx}
\alias{summary.arx}
\alias{vcov.arx}

\title{ Extraction functions for 'arx' objects }

\description{ Extraction functions for objects of class 'arx' }

\usage{
\method{coef}{arx}(object, spec = c("mean", "variance", "both"), ...)
\method{fitted}{arx}(object, spec = c("mean", "variance"), ...)
\method{print}{arx}(x, ...)
\method{residuals}{arx}(object, std = FALSE, ...)
\method{summary}{arx}(object, ...)
\method{vcov}{arx}(object, spec = c("mean", "variance"), ...)
}

\arguments{
  \item{object}{an object of class 'arx'}
  \item{x}{an object of class 'arx'}
  \item{spec}{extract values for the 'mean' or 'variance' specification, or 'bothl' (where applicable)}
  \item{std}{logical. If FALSE (default), then the mean residuals are returned. If TRUE, then the standardised residuals are returned}
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{A numeric vector containing the parameter estimates}
  \item{fitted:}{A \code{\link{zoo}} object}
  \item{print:}{Prints the most important parts of the estimation results}
  \item{residuals:}{A \code{\link{zoo}} object with the residuals}
  \item{summary:}{A print of the items in the \code{\link{arx}} object}
  \item{vcov:}{The variance-covariance matrix}
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}
\seealso{ \code{\link{arx}} }

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 100)

##Simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*100), 100, 4)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean, and log-ARCH(3) in the variance:
mymod <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:3)

##print results:
print(mymod)

##print the entries of object 'mymod':
summary(mymod)

##extract mean coefficients:
coef(mymod)

##extract log-variance coefficients:
coef(mymod, spec="variance")

##extract all coefficients:
coef(mymod, spec="both")

##extract variance-covariance matrix of mean equation:
vcov(mymod)

##extract variance-covariance matrix of log-variance equation:
vcov(mymod, spec="variance")

##extract and plot the fitted mean values:
mfit <- fitted(mymod)
plot(mfit)

##extract and plot the fitted variance:
vfit <- fitted(mymod, spec="variance")
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
