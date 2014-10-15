\name{regs.var}
\alias{regs.var}

\title{Create the regressors of a log-ARCH-X model}

\description{Creates the regressors of a log-ARCH-X model, see \code{\link{arx}} and \code{\link{getsv}}}

\usage{
regs.var(e, vc = TRUE, arch = NULL, asym = NULL, log.ewma = NULL, vxreg = NULL,
  p = 2, zero.adj = 0.1)
}

\arguments{
  \item{e}{
numeric vector, time-series or \code{\link{zoo}} object. Missing values in the
beginning or at the end of the series is allowed, as they are removed with the
\code{\link{na.trim}} command from the \code{\link{zoo}} package
}
  \item{vc}{
logical, TRUE (default) or FALSE. TRUE creates an intercept, FALSE does not
}
  \item{arch}{
integer vector, say, c(1,3) or 2:5. The ARCH-lags to include in the
log-variance specification
}
  \item{asym}{
integer vector, say, c(1) or 1:3. The asymmetry or leverage terms to
include in the log-variance specification
}
  \item{log.ewma}{
NULL (default) or a list. If NULL then log(EWMA) is not included as
volatility proxy. If a list, then log(EWMA) is included as a volatility
proxy.
}
  \item{vxreg}{
numeric matrix, time-series or \code{\link{zoo}} object of conditioning variables.
Note that missing values in the beginning or at the end of the series
is allowed, as they are removed with the \code{\link{na.trim}} command from the \code{\link{zoo}}
package
}
  \item{p}{
numeric value greater than zero. The power in e^p
}
  \item{zero.adj}{
numeric value between 0 and 1. The quantile adjustment for zero values.
The default 0.1 means zeros on e are replaced by the 10 percent quantile of
the absolute residuals before taking the logarithm
}
}

\value{Matrix with regressors}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{regs.mean}}, \code{\link{arx}} and \code{\link{getsv}}
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}