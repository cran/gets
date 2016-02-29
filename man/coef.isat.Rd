\name{coef.isat}
\alias{coef.isat}
\alias{fitted.isat}
\alias{logLik.isat}
\alias{plot.isat}
\alias{print.isat}
\alias{residuals.isat}
\alias{summary.isat}
\alias{vcov.isat}

\title{Extraction functions for 'isat' objects}

\description{Extraction functions for objects of class 'isat'}

\usage{
\method{coef}{isat}(object, ...)
\method{fitted}{isat}(object, ...)
\method{logLik}{isat}(object, ...)
\method{plot}{isat}(x, col=c("red","blue"), lty=c("solid","solid"),
  lwd=c(1,1), coef.path=TRUE, ...)
\method{print}{isat}(x, ...)
\method{residuals}{isat}(object, std=FALSE, ...)
\method{summary}{isat}(object, ...)
\method{vcov}{isat}(object, ...)
}

\arguments{
  \item{object}{an object of class 'isat'}
  \item{x}{an object of class 'isat'}
  \item{std}{logical. If FALSE (default), then the mean residuals are returned. If TRUE, then the standardised residuals are returned}
  \item{col}{colours of fitted (default=red) and actual (default=blue) lines}
  \item{lty}{types of fitted (default=solid) and actual (default=solid) lines}
  \item{lwd}{widths of fitted (default=1) and actual (default=1) lines}
  \item{coef.path}{logical. Only applicable if there are retained indicators after the application of \code{isat}}
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with fitted values}
  \item{logLik:}{a numeric, the log-likelihood (normal density)}
  \item{plot:}{plot of the fitted values and the residuals}
  \item{print:}{a print of the estimation results}
  \item{residuals:}{a \code{\link{zoo}} object with the residuals}
  \item{summary:}{a print of the items in the \code{\link{isat}} object}
  \item{vcov:}{variance-covariance matrix}
}

\author{
James Reade, \url{https://sites.google.com/site/jjamesreade/},\cr
Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{ \code{\link{paths}}, \code{\link{terminals}}, \code{\link{coef.gets}}, \code{\link{getsm}}, \code{\link{arx}} }

\examples{
##step indicator saturation:
set.seed(123)
y <- rnorm(30)
isatmod <- isat(y)

##print results:
print(isatmod)

##plot the fitted vs. actual values, and the residuals:
plot(isatmod)

##print the entries of object 'isatmod':
summary(isatmod)

##extract coefficients of the simplified (specific) model:
coef(isatmod)

##extract log-likelihood:
logLik(isatmod)

##extract variance-covariance matrix of simplified
##(specific) model:
vcov(isatmod)

##extract and plot the fitted values:
mfit <- fitted(isatmod)
plot(mfit)

##extract and plot (mean) residuals:
epshat <- residuals(isatmod)
plot(epshat)

##extract and plot standardised residuals:
zhat <- residuals(isatmod, std=TRUE)
plot(zhat)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
