\name{isatvar}
\alias{isatvar}

\title{Variance of the coefficient path}

\description{
Takes a 'gets' object returned by the \code{isat} function as input and returns the coefficient path of the constant (and long-run equilibrium if 'lr' is specified) together with its approximate variance and standard errors.
}

\usage{
isatvar(x, lr=FALSE)
}

\arguments{
  \item{x}{a 'gets' object obtained with the \code{\link{isat}} function}
  \item{lr}{logical. If TRUE and 'x' contains autoregressive elements, then \code{isatvar} also returns the long-run equilibrium coefficient path with its variance and standard deviation. See Pretis (2015).}
}

\details{The function computes the approximate variance and standard errors of the intercept term with structural breaks determined by \code{\link{isat}}. This permits hypothesis testing and plotting of approximate confidence intervals for the intercept in the presence of structural breaks. For dynamic autoregressive models in \code{\link{isat}} the \code{lr} argument returns the time-varying long-run equilibrium together with its approximate variance and standard errors.}

\value{If \code{lr=FALSE}: A Tx4 matrix (with T = number of observations) where the first column denotes the coefficient path relative to the full sample coefficient, the second column the coefficient path of the intercept, the third the approximate variance of the coefficient path, and the fourth column the approximate standard errors of the coefficient path. If \code{lr=TRUE}: A Tx7 matrix where the first four columns are identical to the \code{lr=FALSE} case, and the additional columns denote the long-run equilibrium coefficient path, together with the approximate variance and standard errors of the long-run equilibrium coefficient path.}

\references{

Pretis, F. (2015): 'Testing for time-varying predictive accuracy using bias-corrected indicator saturation'. Oxford Department of Economics Discussion Paper.

}
                
\author{
Felix Pretis, \url{http://www.felixpretis.org/}\cr
James Reade, \url{https://sites.google.com/site/jjamesreade/}
}

\seealso{
\code{\link{isat}}, \code{\link{coef.gets}}, \code{\link{plot.gets}}, \code{\link{biascorr}}, \code{\link{isattest}}
}

\examples{
##Variance in presence of a break
#nile <- as.zoo(Nile)
#isat.nile <- isat(nile, sis=TRUE, iis=FALSE, plot=FALSE, t.pval=0.005)
#var <- isatvar(isat.nile)

#plot(nile)
#lines(isat.nile$mean.fit, col="red")
#lines(isat.nile$mean.fit + 2*var$const.se, col="blue", lty=3)
#lines(isat.nile$mean.fit - 2*var$const.se, col="blue", lty=3)

##Variance when there is no break
#set.seed(1)
#x <- as.zoo(rnorm(100, 0, 1))
#isat.x <- isat(x, sis=TRUE, iis=FALSE, plot=TRUE, t.pval=0.005)
#var.x <- isatvar(isat.x)

#plot(x)
#lines(isat.x$mean.fit, col="red")
#lines(isat.x$mean.fit + 2*var.x[,2], col="blue", lty=3)
#lines(isat.x$mean.fit - 2*var.x[,2], col="blue", lty=3)

##Variance of the long-run equilibrium coefficient path

#nile <- as.zoo(Nile)
#isat.nile <- isat(nile, sis=TRUE, iis=FALSE, plot=TRUE, t.pval=0.005, ar=1:2)
#var <- isatvar(isat.nile, lr=TRUE)
}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
