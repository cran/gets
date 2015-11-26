\name{iim}
\alias{iim}
\alias{sim}
\alias{tim}

\title{Make Indicator Matrices (Impulses, Steps, Trends)}

\description{
Auxiliary functions to make, respectively, matrices of impulse indicators (iim), step indicators (sim) and trend indicators (tim)
}

\usage{
##make matrix of impulse indicators:
iim(x, which.ones = NULL)

##make matrix of step indicators:
sim(x, which.ones = NULL)

##make matrix of trend indicators:
tim(x, which.ones = NULL, log.trend = FALSE)
}

\arguments{
  \item{x}{either an integer (the length of the series in question) or the series itself}
  \item{which.ones}{the locations of the impulses. If NULL (the default), then all impulses are returned}
  \item{log.trend}{logical. If TRUE, then the natural log is applied on the trends}
}

\details{
If x is a series or vector of observations, then the index of x will be used as for the \code{\link{zoo}}, and for the labelling of the impulses.

Note: For \code{sim} and \code{tim} the first indicator is removed, since it is exactly colinear with the others.
}

\value{
A \code{\link{zoo}} matrix containing the impulses
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/} }

\seealso{
\code{\link{isat}}, \code{\link{zoo}}
}
\examples{
##generate series:
y <- rnorm(50)

##make matrix of impulse indicators:
mIIM <- iim(50)

##make matrix of step indicators:
mSIM <- sim(y)

##make matrix of trend indicators, but only every third:
mTIM <- tim(y, which.ones=seq(1,50,3))
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
