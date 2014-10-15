\name{glog.ep}
\alias{glog.ep}

\title{ Adjust for zero values and compute log(abs(e)^p) }

\description{
Adjusts a series called e - typically a series of residuals or
financial returns - for zero values, so that the logarithm can be
applied on the absolute pth. exponentiated values. Next,
log(abs(e)^p) is computed
}

\usage{
glog.ep(e, zero.adj = 0.1, p = 2, na.replace = NA)
}

\arguments{
  \item{e}{numeric vector, time series or \code{\link{zoo}} object}
  \item{zero.adj}{
numeric value between 0 and 1 (the quantile adjustment for zero values).
The default 0.1 means zeros are replaced by the 10 percent quantile
of abs(e) before taking the logarithm
}
  \item{p}{numeric value greater than zero. The power of the log-volatility specification}
  \item{na.replace}{the value to replace NA values with. Default: na.replace=NA}
}

\value{log(abs(e)^p), a numeric, where the zeros in e have been adjusted}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
