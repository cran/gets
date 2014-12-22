\name{diagnostics}
\alias{diagnostics}

\title{Diagnostics tests of residuals}

\description{Auxiliary function (i.e. not intended for the average user) called by the \code{\link{getsm}}, \code{\link{getsv}} and \code{\link{isat}} functions. The \code{diagnostics} function tests for autocorrelation, ARCH and non-normality in a residual series. The autocorrelation and ARCH tests are conducted as Ljung and Box (1979) tests of autocorrelation in the residuals and squared residuals, respectively, whereas the test for non-normality is that of Jarque and Bera (1980)}

\usage{
diagnostics(x, s2=1, ar.LjungB=c(1,0.025), arch.LjungB=c(1,0.025),
  normality.JarqueB=NULL, verbose=FALSE)
}

\arguments{
  \item{x}{numeric vector, typically the residuals from a regression}
  \item{s2}{the standard deviation of x}
  \item{ar.LjungB}{a two element vector or NULL. In the former case, the first element contains the AR-order, the second element the p-value. If NULL, then a test for autocorrelation is not conducted}
  \item{arch.LjungB}{a two element vector or NULL. In the former case, the first element contains the ARCH-order, the second element the p-value. If NULL, then a test for ARCH is not conducted}
  \item{normality.JarqueB}{a value between 0 and 1 or NULL. In the former case, a test for non-normality is conducted using a significance level equal to the numeric value. If NULL, then no test for non-normality is conducted}
  \item{verbose}{unused}
}

\value{logical. If TRUE, then the residuals series passes the diagnostics tests}

\references{
G. Ljung and G. Box (1979): 'On a Measure of Lack of Fit in Time
  Series Models'. Biometrika 66, pp. 265-270

C. Jarque and A. Bera (1980): 'Efficient Tests for Normality,
  Homoscedasticity and Serial Independence'. Economics Letters 6,
  pp. 255-259
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/} }

\seealso{\code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}} }

\examples{
##check for autocorrelation, ARCH and non-normality:
diagnostics(rnorm(40))
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
