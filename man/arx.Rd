\name{arx}
\alias{arx}
\title{Estimate an AR-X model with log-ARCH-X errors}
\description{
Estimation is by OLS in two steps. In the first the AR-X mean specification is estimated, whereas in the second step the residuals from the first are used to fit a log-ARCH-X model to the log-variance. The natural logarithm of the squared residuals constitutes the regressand in the second step.

The AR-X mean specification can contain an intercept, AR-terms, lagged moving averages of the regressand and other conditioning covariates ('X'). The log-variance specification can contain log-ARCH terms, asymmetry or 'leverage' terms, log(EqWMA) where EqWMA is a lagged equally weighted moving average of past squared residuals (a volatility proxy) and other conditioning covariates ('X').
}
\usage{
arx(y, mc = NULL, ar = NULL, ewma = NULL, mxreg = NULL, arch = NULL,
  asym = NULL, log.ewma = NULL, vxreg = NULL, zero.adj = 0.1, vc.adj = TRUE,
  vcov.type = c("ordinary", "white"), qstat.options = NULL, tol = 1e-07,
  LAPACK = FALSE, verbose = TRUE)
}
\arguments{
  \item{y}{
numeric vector, time-series or \code{\link{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command from the \code{\link{zoo}} package}
  \item{mc}{
logical or NULL (default). TRUE includes intercept in the mean specification, whereas any other value does not
}
  \item{ar}{
integer vector, say, c(2,4) or 1:4. The AR-lags to include in the mean specification
}
  \item{ewma}{
either NULL (default) or a list with arguments sent to the \code{\link{eqwma}} function. In the latter case a lagged moving average of y is included as a regressor
}
  \item{mxreg}{
numeric vector or matrix, say, a \code{\link{zoo}} object, of conditioning variables. Note that missing values in the beginning or at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command from the \code{\link{zoo}}. Note also that, if both y and mxreg are zoo objects, then their samples are chose to match
}
  \item{arch}{
integer vector, say, c(1,3) or 2:5. The ARCH-lags to include in the log-variance specification
}
  \item{asym}{
integer vector, say, c(1) or 1:3. The asymmetry or leverage terms to include in the log-variance specification
}
  \item{log.ewma}{
either NULL (default) or a list with arguments sent to the \code{\link{leqwma}} function. In the latter case (essentially the log of a volatility proxy) the log of a lagged moving average of the residuals e is included as a regressor
}
  \item{vxreg}{
numeric vector or matrix, say, a \code{\link{zoo}} object, of conditioning variables. Note that missing values in the beginning or at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command from the \code{\link{zoo}}. Note also that, if both y and mxreg are zoo objects, then their samples are chose to match
}
  \item{zero.adj}{
numeric value between 0 and 1. The quantile adjustment for zero values. The default 0.1 means the zero residuals are replaced by the 10 percent quantile of the absolute residuals before taking the logarithm
}
  \item{vc.adj}{
logical, TRUE (default) or FALSE. If true then the log-variance intercept is adjusted by the estimate of E[ln(z^2)]. This adjustment is needed for the conditional scale of e to be equal to the conditional standard deviation. If FALSE then the log-variance intercept is not adjusted
}
  \item{vcov.type}{
character vector, "ordinary" (default) or "white". If "ordinary", then the ordinary variance-covariance matrix is used for inference. Otherwise the White (1980) heteroscedasticity robust matrix is used
}
  \item{qstat.options}{
NULL (default) or an integer vector of length two, say, c(2,5). The first value sets the order of the AR diagnostic test, whereas the second value sets the order of the ARCH diagnostic test. NULL sets the vector automatically
}
  \item{tol}{
numeric value (default = 1e-07). The tolerance for detecting linear dependencies in the columns of the regressors (see \code{\link{qr}}) function). Only used if LAPACK is FALSE (default)
}
  \item{LAPACK}{
logical, TRUE or FALSE (default). If true use LAPACK otherwise use LINPACK (see \code{\link{qr}} function)
}
  \item{verbose}{
logical, TRUE (default) or FALSE. FALSE returns less output and is therefore (slighthly) faster
}
}
\details{ See Sucarrat and Escribano (2012) }

\value{A list of class 'arx'}

\references{
Genaro Sucarrat and Alvaro Escribano (2012): 'Automated Financial Model Selection: General-to-Specific Modelling of the Mean and Volatility Specifications', Oxford Bulletin of Economics and Statistics 74, Issue no. 5 (October), pp. 716-735

Halbert White (1980): 'A Heteroskedasticity-Consistent Covariance Matrix Estimator and a Direct Test for Heteroskedasticity', Econometrica 48, pp. 817-838
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
Extraction functions (S3 methods): \code{\link{coef.arx}}, \code{\link{fitted.arx}}, \code{\link{print.arx}}, \code{\link{residuals.arx}}, \code{\link{summary.arx}} and \code{\link{vcov.arx}}
Related functions: \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{eqwma}}, \code{\link{leqwma}}
}

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 100)

##Simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*100), 100, 4)

##estimate an AR(2) with intercept:
arx(y, mc=TRUE, ar=1:2)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean:
arx(y, mc=TRUE, ar=1:2, mxreg=xregs)

##estimate a log-variance specification with a log-ARCH(4)
##structure:
arx(y, arch=1:4)

##estimate a log-variance specification with a log-ARCH(4)
##structure and an asymmetry or leverage term:
arx(y, arch=1:4, asym=1)

##estimate a log-variance specification with a log-ARCH(4)
##structure, an asymmetry or leverage term, a 10-period log(EWMA) as
##volatility proxy, and the log of the squareds of the conditioning
##regressors in the log-variance specification:
arx(y, arch=1:4, asym=1, log.ewma=list(length=10), vxreg=log(xregs^2))

##estimate an AR(2) with intercept and four conditioning regressors
##in the mean, and a log-variance specification with a log-ARCH(4)
##structure, an asymmetry or leverage term, a 10-period log(EWMA) as
##volatility proxy, and the log of the squareds of the conditioning
##regressors in the log-variance specification:
arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:4, asym=1,
  log.ewma=list(length=10), vxreg=log(xregs^2))
}
\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
