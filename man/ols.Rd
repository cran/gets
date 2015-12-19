\name{ols}
\alias{ols}

\title{OLS estimation}

\description{OLS estimation with the QR decomposition and, for some options, computation of variance-covariance matrices}

\usage{
ols(y, x, tol=1e-07, LAPACK=FALSE, method=1, user.fun=NULL, user.options=NULL)
}

\arguments{
  \item{y}{numeric vector, the regressand}
  \item{x}{numeric matrix, the regressors}
  \item{tol}{numeric value. The tolerance for detecting linear dependencies in the columns of the regressors, see \code{\link{qr}} function. Only used if LAPACK is FALSE}
  \item{LAPACK}{logical, TRUE or FALSE (default). If true use LAPACK otherwise use LINPACK, see \code{\link{qr}} function}
  \item{method}{1 (default) or 2. Method 2 returns slightly more information, which means it is slightly slower. However, the information returned can be used to speed up the computation of variance-covariance matrices}
  \item{user.fun}{the name of the user-function (a character)}
  \item{user.options}{a list with arguments (entries) that are passed on to the user-function}
}

\details{
\code{method=1} or 2 only returns the OLS coefficient estimates together with the QR-information. \code{method=1} is slightly faster than \code{method=2}. \code{method=3} returns in addition the ordinary variance-covariance matrix of the OLS estimator. \code{method=4} returns the White (1980) heteroscedasticity robust variance-covariance matrix in addition to the information returned by \code{method=3}, whereas \code{method=5} does the same except that the variance-covariance matrix now is that of Newey and West (1987).
}

\value{A list with some or all of the following elements:
  \item{qr}{}
  \item{rank}{}
  \item{qraux}{}
  \item{pivot}{}
  \item{xtxinv}{}
  \item{xtx}{}
  \item{xty}{}
  \item{coefficients}{}
  \item{fitted}{}
  \item{residuals}{}
  \item{resids2}{}
  \item{rss}{}
  \item{df}{}
  \item{sigma2}{}
  \item{omegahat}{}
  \item{vcov}{}
}

\references{
H. White (1980): 'A Heteroskedasticity-Consistent Covariance Matrix and a Direct Test for Heteroskedasticity', Econometrica 48, pp. 817-838.

W. Newey and K. West (1987): 'A Simple Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix', Econometrica 55, pp. 703-708.
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{qr}}, \code{\link{solve.qr}}}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
