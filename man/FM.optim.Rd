\name{FM.optim}
\alias{FM.optim}
\title{Optimizes the frequent monitoring log-likelihood}
\usage{
  FM.optim(simDataList, u, initGuess)
}
\arguments{
  \item{simDataList}{A list of synthetic observed datasets,
  returned by \code{\link{makedata.simple}}}

  \item{u}{The index of the desired entry of simDataList}

  \item{initGuess}{A vector, initial guess for beta}
}
\value{
  An optim object corresponding to maximized frequent
  monitoring log-likelihood
}
\description{
  This function converts a dataset generated in the format
  returned by \code{\link{makedata.simple}} and finds the
  intervals compatible under FM using
  \code{\link{FM.data}}, then infers the MLE birth, shift,
  and death rates.
}

