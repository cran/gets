\name{gets-package}
\alias{gets-package}
\alias{gets}
\docType{package}
\title{General-to-Specific (GETS) Model Selection}

\description{
Automated multi-path General-to-Specific (GETS) model selection of an AR-X model with log-ARCH-X errors.

The three main functions of the package are \code{\link{arx}}, \code{\link{getsm}} and \code{\link{getsv}}. The \code{arx} function returns an object of the 'arx' class, i.e. a model with an AR-X specification in the mean and a log-ARCH-X specification of the log-variance. The second function \code{getsm} undertakes GETS model selection of the mean specification of an \code{arx} object, and returns an object of the 'gets' class. The third function \code{getsv} undertakes GETS model selection of the log-variance specification of an \code{arx} object (or alternatively of a \code{gets} object obtained with the \code{getsm} function), and returns an object of the 'gets' class.

The package also provides auxiliary functions used by \code{arx}, \code{getsm} and \code{getsv}, in addition to extraction functions (mainly S3 methods) that can be applied to objects of the 'arx' and 'gets' classes.
}
\details{
\tabular{ll}{
Package: \tab gets\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-10-14\cr
License: \tab GPL-2\cr
}
The code was originally developed in relation with G. Sucarrat and A. Escribano (2012): 'Automated Financial Model Selection: General-to-Specific Modelling of the Mean and Volatility
Specifications', Oxford Bulletin of Economics and Statistics 74, Issue 5 (October), pp. 716-735.
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\references{
G. Sucarrat and A. Escribano (2012): 'Automated Financial Model Selection: General-to-Specific Modelling of the Mean and Volatility
Specifications', Oxford Bulletin of Economics and Statistics 74, Issue 5 (October), pp. 716-735
}

\seealso{
\code{\link{arx}}, \code{\link{getsm}}, \code{\link{getsv}}
}

\examples{
##Simulate from an AR(1):
set.seed(123)
y <- arima.sim(list(ar=0.4), 100)

##Estimate an AR(2) with intercept as mean specification
##and a log-ARCH(4) as log-volatility specification:
myModel <- arx(y, mc=TRUE, ar=1:2, arch=1:4)

##GETS model selection of the mean of myModel:
simpleMean <- getsm(myModel)

##GETS model selection of the log-variance of simpleMean:
simpleVar <- getsv(simpleMean)

##results:
print(simpleMean)
print(simpleVar)

}
