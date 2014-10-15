\name{coef.gets}
\alias{coef.gets}
\alias{fitted.gets}
\alias{paths}
\alias{print.gets}
\alias{residuals.gets}
\alias{summary.gets}
\alias{terminals}
\alias{vcov.gets}

\title{ Extraction functions for 'gets' objects }

\description{ Extraction functions for objects of class 'gets' }

\usage{
\method{coef}{gets}(object, spec = NULL, ...)
\method{fitted}{gets}(object, spec = NULL, ...)
## extraction function for class 'gets'
paths(object, ...)
\method{print}{gets}(x, ...)
\method{residuals}{gets}(object, std = NULL, ...)
\method{summary}{gets}(object, ...)
## extraction function for class 'gets'
terminals(object, ...)
\method{vcov}{gets}(object, spec = NULL, ...)
}

\arguments{
  \item{object}{an object of class 'gets'}
  \item{x}{an object of class 'gets'}
  \item{spec}{NULL (default), "mean" or "variance". Determines whether values associated with the mean or variance specification should be extracted. If NULL, then it is automatically determined whether \code{\link{getsm}} or \code{\link{getsv}} has been used (getsm returns values pertaining to the mean specification, getsv returns values pertaining to the variance specification)}
  \item{std}{NULL (default), FALSE or TRUE. If std = FALSE, then the mean residuals are returned. If std = TRUE, then the standardised residuals are returned. If std = NULL, then it is automatically determined whether \code{\link{getsm}} or \code{\link{getsv}} has been used (getsm returns the residuals, getsv returns the standardised residuals) }
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{A numeric vector containing the parameter estimates}
  \item{fitted:}{A \code{\link{zoo}} object}
  \item{paths:}{A \code{\link{list}} with the paths searched (each number refers to the regressor ordering in the GUM) }
  \item{print:}{Prints the most important parts of the search and estimation results}
  \item{residuals:}{A \code{\link{zoo}} object with the residuals}
  \item{summary:}{A print of the items in the \code{\link{gets}} object}
  \item{terminals:}{A \code{\link{list}} with the terminal models (each number refers to the regressor ordering in the GUM)}
  \item{vcov:}{The variance-covariance matrix}
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}
\seealso{ \code{\link{getsm}}, \code{\link{getsv}} }

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 100)

##Simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*100), 100, 4)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean, and log-ARCH(3) in the variance:
mymod <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:3)

##undertake General-to-Specific (GETS) model selection of
##the mean:
meanmod <- getsm(mymod)

##undertake General-to-Specific (GETS) model selection of
##the log-variance:
varmod01 <- getsv(mymod)

##undertake General-to-Specific (GETS) model selection of
##the log-variance (simplified model):
varmod02 <- getsv(meanmod)

##print results:
print(meanmod)
print(varmod01)
print(varmod02)

##print the entries of object 'gets':
summary(meanmod)
summary(varmod01)
summary(varmod02)

##extract coefficients of the simplified (specific) model:
coef(meanmod) #mean spec
coef(varmod01) #log-variance spec 1
coef(varmod02) #log-variance spec 2

##extract the paths searched:
paths(meanmod) #mean
paths(varmod01) #log-variance spec 1
paths(varmod02) #log-variance spec 2

##extract the terminal models:
terminals(meanmod) #mean
terminals(varmod01) #log-variance spec 1
terminals(varmod02) #log-variance spec 2

##extract variance-covariance matrix of simplified
##(specific) model:
vcov(meanmod) #mean spec
vcov(varmod01) #log-variance spec 1
vcov(varmod02) #log-variance spec 2

##extract and plot the fitted values:
mfit <- fitted(meanmod) #mean fit
plot(mfit)
vfit01 <- fitted(varmod01) #variance fit
plot(vfit01)
vfit02 <- fitted(varmod02) #variance fit
plot(vfit02)

##extract and plot residuals:
epshat <- residuals(meanmod) #mean residuals
plot(epshat)
zhat01 <- residuals(varmod01) #standardised residuals
plot(zhat01)
zhat02 <- residuals(varmod02) #standardised residuals
plot(zhat02)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
