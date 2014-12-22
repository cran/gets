\name{ols}
\alias{ols}

\title{OLS estimation with the QR decomposition}

\description{Accurate, robust and fast OLS estimation with the QR decomposition}

\usage{
ols(y, x, tol=1e-07, LAPACK=FALSE, method=1)
}

\arguments{
  \item{y}{numeric vector, the regressand}
  \item{x}{numeric matrix, the regressors}
  \item{tol}{numeric value. The tolerance for detecting linear dependencies in the columns of the regressors, see \code{\link{qr}} function. Only used if LAPACK is FALSE}
  \item{LAPACK}{logical, TRUE or FALSE (default). If true use LAPACK otherwise use LINPACK, see \code{\link{qr}} function}
  \item{method}{1 (default) or 2. Method 2 returns slightly more information, which means it is slightly slower. However, the information returned can be used to speed up the computation of variance-covariance matrices}
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
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{\code{\link{qr}}, \code{\link{solve.qr}}}

\keyword{Statistical Models}