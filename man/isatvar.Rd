\name{isatvar}
\alias{isatvar}

\title{Variance of the coefficient path}

\description{
Takes a 'gets' object returned by the \code{isat} function as input and returns the variance and standard errors of the constant on the coefficient path.
}

\usage{
isatvar(x, ...)
}

\arguments{
  \item{x}{a 'gets' obtained with the \code{\link{isat}} function}
  \item{\dots}{additional arguments}
}

\details{The function computes the variance and standard errors of the intercept term with structural breaks determined by \code{\link{isat}}. This permits hypothesis testing and plotting of confidence intervals for the intercept in the presence of structural breaks.}

\value{A Tx2 matrix (with T = number of observations) where the first column denotes the variance, and the second column the standard errors of the intercept term.}

\author{
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade, \url{https://sites.google.com/site/jjamesreade/}
}

\seealso{
\code{\link{isat}}, \code{\link{coef.gets}} \code{\link{plot.gets}}
}

\examples{
##Variance in presence of a break
nile <- as.zoo(Nile)
isat.nile <- isat(nile, sis=TRUE, iis=FALSE, plot=FALSE, t.pval=0.005)
var <- isatvar(isat.nile)

plot(nile)
lines(isat.nile$mean.fit, col="red")
lines(isat.nile$mean.fit + 2*var[,2], col="blue", lty=3)
lines(isat.nile$mean.fit - 2*var[,2], col="blue", lty=3)


##Variance when there is no break
set.seed(1)
x <- as.zoo(rnorm(100, 0, 1))
isat.x <- isat(x, sis=TRUE, iis=FALSE, plot=TRUE, t.pval=0.005)
var.x <- isatvar(isat.x)

plot(x)
lines(isat.x$mean.fit, col="red")
lines(isat.x$mean.fit + 2*var.x[,2], col="blue", lty=3)
lines(isat.x$mean.fit - 2*var.x[,2], col="blue", lty=3)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Financial Econometrics}

