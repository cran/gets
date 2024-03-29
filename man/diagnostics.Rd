\name{diagnostics}
\alias{diagnostics}

\title{Diagnostics tests}

\description{Auxiliary function (i.e. not intended for the average user) called by the \code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}, \code{\link{getsFun}} and \code{\link{blocksFun}} functions. The \code{diagnostics} function undertakes tests for autocorrelation, ARCH and non-normality in a residual series, and user-defined diagnostics provided via the \code{user.fun} argument (see details). The autocorrelation and ARCH tests are conducted as Ljung and Box (1979) tests for autocorrelation in the residuals and squared residuals, respectively, whereas the test for non-normality is that of Jarque and Bera (1980).}

\usage{
diagnostics(x, ar.LjungB=c(1, 0.025), arch.LjungB=c(1, 0.025),
  normality.JarqueB=NULL, verbose=TRUE, user.fun=NULL, ...)
}

\arguments{
  \item{x}{a \code{\link{list}}, for example the estimation result of \code{\link{ols}}. The tests for serial correlation, ARCH and normality look for an entry in the list named \code{std.residuals} or \code{residuals}}
  \item{ar.LjungB}{a two element vector or \code{NULL}. In the former case, the first element contains the AR-order, the second element the significance level. If \code{NULL}, then a test for autocorrelation is not conducted}
  \item{arch.LjungB}{a two element vector or \code{NULL}. In the former case, the first element contains the ARCH-order, the second element the significance level. If \code{NULL}, then a test for ARCH is not conducted}
  \item{normality.JarqueB}{\code{NULL} (the default) or a value between 0 and 1. In the latter case, a test for non-normality is conducted using a significance level equal to \code{normality.JarqueB}. If \code{NULL}, then no test for non-normality is conducted}
  \item{verbose}{logical. If \code{TRUE}, then a \code{\link{data.frame}} with the results of the diagnostics is returned. If \code{FALSE}, then the return-value is a \code{logical} that indicates whether the model passes the diagnostics (\code{TRUE} if it does, otherwise \code{FALSE}) }
  \item{user.fun}{\code{NULL} or a \code{\link{list}} with at least one entry, \code{name} (must be of class \code{character}), which should contain the name of the user-defined function. See details }
  \item{...}{further arguments (ignored) to accommodate deleted arguments from past versions of the functions}
}

\details{
  The argument \code{user.fun} enables the user to specify additional diagnostics. To do this, the argument should be a \code{\link{list}} with at least one entry, \code{name} (of class \code{character}), that contains the name of the user-defined function. The call to this function is executed with \code{\link{do.call}}, whose default value on \code{envir} is \code{parent.frame()}. Usually, this will be the global environment (\code{.GlobalEnv}), but it can be changed by adding an entry named \code{envir} to the list that indicates where the user-defined function resides. If the \code{verbose} argument is set to \code{FALSE}, then an entry named \code{pval} must be provided. This entry should contain the chosen significance level or levels, i.e. either a scalar or a vector of length equal to the number of \emph{p}-values returned by the user-defined diagnostics function (see examples). The user can also specify whether a rejection of the tests should cause the diagnostics to fail (default behaviour) or whether a rejection is desirable. For that purpose, a named entry \code{is.reject.bad} can be added that stores a logical vector of length equal to the number of tests conducted in the user diagnostics function. The first entry of the vector governs the diagnostic decision for the first row that the user diagnostics function returns, the second entry the decision for the second row etc. Additional entries in the \code{list} are passed on as arguments to the user-defined function.
 
  The user-defined function should refer to the named items of the estimation result \code{x} (see examples), and the value returned by the user-defined function should be a matrix of dimension \emph{m} x 3. Here, \emph{m} is the number of diagnostic tests performed by the user-defined function. For example, if only a single test is performed, then \emph{m} = 1 and so the returned value should be a 1 x 3 matrix (or a vector of length 3). The three columns of the \emph{m} x 3 matrix should contain, in the following order, 1) the value(s) of the test-statistic(s) (or \code{NA}), 2) the degree(s) of freedom(s) (or \code{NA}) of the tests, and 3) the \emph{p}-value(s) of the test(s). When checking whether the model passes the diagnostics or not, the \emph{p}-value(s) is(are) checked against the value(s) in the entry named \code{pval} in the \code{list} provided to \code{user.fun}. By default, a calculated \emph{p}-value below the corresponding element in \code{pval} causes the diagnostics to fail. If a named entry \code{is.reject.bad} exists, this decision rule is only applied to tests whose corresponding entry is \code{TRUE} while the decision rule is reversed for those with entry \code{FALSE}. For these tests, the diagnostics fail if the hypothesis cannot be rejected.
}

\value{
  \item{verbose=TRUE}{a \code{\link{data.frame}} that contains the diagnostics results}
  \item{verbose=FALSE}{a \code{logical} (of length one) indicating whether the residuals and/or model passes ALL the diagnostics (\code{TRUE} if it does, \code{FALSE} otherwise)}
}

\references{
C. Jarque and A. Bera (1980): 'Efficient Tests for Normality,
  Homoscedasticity and Serial Independence'. Economics Letters 6,
  pp. 255-259

G. Ljung and G. Box (1979): 'On a Measure of Lack of Fit in Time
  Series Models'. Biometrika 66, pp. 265-270
}

\author{
Genaro Sucarrat, \url{http://www.sucarrat.net/}\cr
Jonas Kurle, \url{https://www.jonaskurle.com/}
}

\seealso{\code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}, \code{\link{isat}}, \code{\link{getsFun}}, \code{\link{blocksFun}}}

\examples{
##generate some data:
set.seed(123)
vY <- rnorm(20) #the regressand
mX <- matrix(rnorm(3*20), 20, 3) #the regressors
est <- ols(vY,mX)

##return a data-frame with autocorrelation and ARCH diagnostics (default),
##and check whether they pass (the default p-value is 0.025):
diagnostics(est)
diagnostics(est, verbose=FALSE)

##add the Jarque-Bera normality test to the diagnostics (w/p-value=0.05):
diagnostics(est, normality.JarqueB=0.05)
diagnostics(est, normality.JarqueB=0.05, verbose=FALSE)

##user-defined Shapiro-Wilks test for non-normality of the residuals:
SWtest <- function(x, ...){
  tmp <- shapiro.test(x$residuals)
  return( c(tmp$statistic, NA, tmp$p.value) )
}
diagnostics(est, user.fun=list(name="SWtest", pval=0.05))
diagnostics(est, user.fun=list(name="SWtest", pval=0.05), verbose=FALSE)

##user-defined test but diagnostics fail if do not reject (illustration only)
diagnostics(est, user.fun=list(name="SWtest", pval=0.05, is.reject.bad = FALSE))
diagnostics(est, user.fun=list(name="SWtest", pval=0.05, is.reject.bad = FALSE),
  verbose=FALSE)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Climate Econometrics}
\keyword{Financial Econometrics}
