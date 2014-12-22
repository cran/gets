\name{jb.test}
\alias{jb.test}

\title{Jarque-Bera test for normality}

\description{Jarque-Bera test for normality}

\usage{jb.test(x)}

\arguments{
  \item{x}{a numeric vector}
}

\value{A list with elements:
  \item{skewness}{a numeric, an estimate of the excess skewness relative to the normal}
  \item{kurtosis}{a numeric, an estimate of the excess kurtosis relative to the normal}
  \item{statistic}{a numeric, the Jarque-Bera test statistic}
  \item{df}{2, the degrees of freedom}
  \item{p.value}{numeric between 0 and 1, the p-value of the test under the null of normality}
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\keyword{Statistical Models}