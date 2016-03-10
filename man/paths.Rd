\name{paths}
\alias{paths}
\alias{terminals}

\title{Extraction functions for 'gets' and 'isat' objects}

\description{Extraction functions for objects of class 'gets' and 'isat'}

\usage{
paths(object, ...)
terminals(object, ...)
}

\arguments{
  \item{object}{an object of class 'gets' or 'isat'}
  \item{\dots}{additional arguments}
}

\value{
  \item{paths:}{a \code{\link{list}} with the paths searched (each number refers to a regressor in the GUM)}
  \item{terminals:}{a \code{\link{list}} with the terminal models (each number refers to a regressor in the GUM)}
}

\author{
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
##regressors in the mean:
mymod <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs)

##General-to-Specific (GETS) modelling of the mean:
meanmod <- getsm(mymod)

##extract the paths searched:
paths(meanmod)

##extract the terminal models:
terminals(meanmod)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
