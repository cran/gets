\name{getsm}
\alias{getsm}
\alias{getsv}

\title{General-to-Specific (GETS) Model Selection of an AR-X model with log-ARCH-X errors}

\description{
The starting model, typically an object of the 'arx' class, is referred to as the General Unrestricted Model (GUM). The \code{getsm} function undertakes multi-path GETS model selection of the mean specification, whereas \code{getsv} does the same for the log-variance specification. The diagnostic tests are optionally undertaken on the standardised residuals.
}
\usage{
getsm(object, keep = NULL, vcov.type = NULL, t.pval = 0.05, do.pet = TRUE,
  wald.pval = 0.05, ar.LjungB = list(lag = NULL, pval = 0.025),
  arch.LjungB = list(lag = NULL, pval = 0.025),
  info.method = c("sc", "aic", "hq"), include.empty = FALSE, zero.adj = NULL,
  vc.adj = NULL, tol = NULL, LAPACK = NULL, max.regs = 1e+05, verbose = TRUE,
  alarm = FALSE)
getsv(object, keep = c(1), t.pval = 0.05, do.pet = TRUE, wald.pval = 0.05,
  ar.LjungB = list(lag = NULL, pval = 0.025),
  arch.LjungB = list(lag = NULL, pval = 0.025),
  info.method = c("sc", "aic", "hq"), include.empty = FALSE, zero.adj = NULL,
  vc.adj = NULL, tol = NULL, LAPACK = NULL, max.regs = 1e+05, verbose = TRUE,
  alarm = FALSE)

}

\arguments{
  \item{object}{an object of class 'arx'. For \code{getsv} the object can instead be an object of class 'gets' obtained with \code{getsm} }
  \item{keep}{the regressors excluded from removal in the specification search. keep=c(1), obligatory when using \code{getsv}, excludes the constant from removal. The regressor numbering is contained in the reg.no column of the GUM}
  \item{vcov.type}{type of variance-covariance matrix used. If NULL (default), then the type used in the estimation of the 'arx' object is used. This can be overridden by "ordinary" (i.e. the ordinary variance-covariance matrix) or "white" (i.e. the White (1980) heteroscedasticity robust variance-covariance matrix)}
  \item{t.pval}{numeric value between 0 and 1. The significance level used for the two-sided regressor significance t-tests}
  \item{do.pet}{logical, TRUE (default) or FALSE. If TRUE then a Parsimonious Encompassing Test (PET) against the GUM is undertaken at each regressor removal for the joint significance of all the deleted regressors along the current path. If FALSE, then a PET is not undertaken at each regressor removal}
  \item{wald.pval}{numeric value between 0 and 1. The significance level used for the PETs}
  \item{ar.LjungB}{a list where the \code{lag} determines the order of a Ljung and Box (1979) test for serial correlation in the standardised residuals, and where \code{pval} contains the significance level. If lag=NULL (default), then the order used is that of the estimated 'arx' object }
  \item{arch.LjungB}{a list where the \code{lag} determines the order of a Ljung and Box (1979) test for serial correlation in the squared standardised residuals, and where \code{pval} contains the significance level. If lag=NULL (default), then the order used is that of the estimated 'arx' object }
  \item{info.method}{character string, "sc" (default), "aic" or "hq", which determines the information criterion used to select among terminal models. The abbreviations are short for the Schwarz or Bayesian information criterion (sc), the Akaike information criterion (aic) and the Hannan-Quinn (hq) information criterion}
  \item{include.empty}{logical. If TRUE, then an empty model is included among the terminal models, if it passes the diagnostic tests, even if it is not equal to one of the terminals. If FALSE (default), then the empty model is not included (unless it is one of the terminals)}
  \item{zero.adj}{numeric value between 0 and 1 or NULL (default). If NULL, then the value used in the estimation of the 'arx' object is used (see \code{arx}) }
  \item{vc.adj}{logical or NULL (default). If TRUE, then the log-variance constant is adjusted by the estimate of E[log(z^2)]. This adjustment is needed for the standardised residuals to have unit variance. If FALSE, then the log-variance constant is not adjusted. If NULL, the default, then the option used in estimating the 'arx' object is used }
  \item{tol}{numeric value or NULL (default). The tolerance for detecting linear dependencies in the columns of the regressors (see \code{qr} function). Only used if LAPACK is FALSE. If NULL, the default, then the value of \code{tol} from the estimation of the 'arx' object is used}
  \item{LAPACK}{logical or NULL (default). If TRUE, then LAPACK is used. Otherwise LINPACK is used (see \code{qr} function). If NULL, the default, the the value of LAPACK from the estimation of the 'arx' object is used }
  \item{max.regs}{integer, the maximum number of regressions along a deletion path}
  \item{verbose}{logical. FALSE returns less output and is therefore faster}
  \item{alarm}{logical. If TRUE, then a sound or beep is emitted when the specification search terminates (in order to alert the user)}
}

\details{ See Sucarrat and Escribano (2012) }

\value{A list of class 'gets'}

\references{
Genaro Sucarrat and Alvaro Escribano (2012): 'Automated Financial
  Model Selection: General-to-Specific Modelling of the Mean and
  Volatility Specifications', Oxford Bulletin of Economics and
  Statistics 74, Issue no. 5 (October), pp. 716-735

G. Ljung and G. Box (1979): 'On a Measure of Lack of Fit in Time
  Series Models'. Biometrika 66, pp. 265-270
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
Extraction functions: \code{\link{coef.gets}}, \code{\link{fitted.gets}}, \code{\link{paths}}, \code{\link{print.gets}}, \code{\link{residuals.gets}},\cr
\code{\link{summary.gets}}, \code{\link{terminals}}, \code{\link{vcov.gets}}\cr

Related functions: \code{\link{arx}}, \code{\link{eqwma}}, \code{\link{leqwma}}, \code{\link{zoo}}, \code{\link{na.trim}}
}

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 100)

##Simulate four independent Gaussian regressors:
xregs <- matrix(rnorm(4*100), 100, 4)

##estimate an AR(2) with intercept and four conditioning
##regressors in the mean, and log-ARCH(3) with log(xregs^2) as
##regressors in the log-variance:
gum01 <- arx(y, mc=TRUE, ar=1:2, mxreg=xregs, arch=1:3,
  vxreg=log(xregs^2))

##GETS model selection of the mean:
meanmod01 <- getsm(gum01)

##GETS model selection of the log-variance:
varmod01 <- getsv(gum01)

##GETS model selection of the log-variance (simplified model):
varmod02 <- getsv(meanmod01)

##GETS model selection of the mean with mean intercept excluded
##from removal:
meanmod02 <- getsm(gum01, keep=1)

##GETS model selection of the mean with non-default
#serial-correlation diagnostics settings:
meanmod03 <- getsm(gum01, ar.LjungB=list(pval=0.05))

##GETS model selection of the mean with very liberal
##(i.e. 20 percent) significance levels:
meanmod04 <- getsm(gum01, t.pval=0.2, wald.pval=0.2)

##GETS model selection of log-variance with all the
##log-ARCH terms excluded from removal:
varmod03 <- getsv(gum01, keep=2:4)

##If ret is a daily (weekends excluded) financial return series,
##then the following log-variance specification includes a
##lagged volatility proxy for the week (5-day average of squared
##return) and for the month (20-day average of squared returns),
##in addition to five log-ARCH terms:
ret <- lgarchSim(1100, arch=0.1, garch=0)
gum02 <- arx(ret, arch=1:5, log.ewma=list(length=c(5,20)) )

##GETS model selection of gum02:
varmod04 <- getsv(gum02)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}
