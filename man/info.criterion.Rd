\name{info.criterion}
\alias{info.criterion}

\title{Computes the Value of an Information Criterion}

\description{
Given a log-likelihood, the number of observations and the number of estimated parameters, the value of a chosen information criterion is computed
}

\usage{
info.criterion(logl, n=NULL, k=NULL, method=c("sc", "aic", "aicc", "hq"))
}

\arguments{
  \item{logl}{numeric, the value of the log-likelihood}
  \item{n}{integer, number of observations}
  \item{k}{integer, number of parameters}
  \item{method}{character, either "sc" (default), "aic", "aicc" or "hq"}
}

\value{a list with elements:
  \item{method}{type of information criterion}
  \item{n}{number of observations}
  \item{k}{number of parameters}
  \item{value}{the value on the information criterion}
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\keyword{Statistical Models}