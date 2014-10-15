\name{regs.mean}
\alias{regs.mean}

\title{Create the regressors of an AR-X model}

\description{Creates the regressors of an AR-X model, see \code{\link{arx}}. }
\usage{
regs.mean(y, mc = NULL, ar = NULL, ewma = NULL, mxreg = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
numeric vector, time-series or zoo object. Note that missing values
in the beginning or at the end of the series is allowed, as they are
removed with the na.trim command from the zoo package
}
  \item{mc}{
logical, TRUE or FALSE (default). TRUE includes intercept in the
specification, FALSE does not
}
  \item{ar}{
integer vector, say, c(2,4) or 1:4. The AR-lags to include in the
specification
}
  \item{ewma}{
NULL or a list of arguments sent to the eqwma function
}
  \item{mxreg}{
numeric matrix, time-series or \code{\link{zoo}} object of conditioning variables.
Note that missing values in the beginning or at the end of the series
is allowed, as they are removed with the \code{\link{na.trim}} command from the \code{zoo}
package
  }
}

\value{Matrix with regressors}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{regs.var}}, \code{\link{arx}} and \code{\link{getsm}}
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}