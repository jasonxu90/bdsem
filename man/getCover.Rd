\name{getCover}
\alias{getCover}
\title{Calculates coverage indicators for estimated parameters}
\usage{
  getCover(opt, trueParam)
}
\arguments{
  \item{opt}{An optim object}

  \item{trueParam}{A vector containing the true parameters}
}
\value{
  A vector of indicators, 1 indicating that the true
  parameter was included in the corresponding confidence
  interval
}
\description{
  This function takes in an object returned by \code{optim}
  and returns logicals on whether the corresponding 95%
  confidence interval of each estimate contains the true
  value.
}

