\name{FM.run}
\alias{FM.run}
\title{Main function for generating observations from birth-shift-death process and inferring parameters under frequent monitoring}
\usage{
  FM.run(N, tList, lam, v, mu, initList, initGuess)
}
\arguments{
  \item{initGuess}{Vector of numbers, initial guess for
  \code{optim}}

  \item{N}{An integer specifying the number of observation
  intervals/realization from the simple BSD process.}

  \item{tList}{A list of observation interval lengths. The
  number of datasets returned is equal to the length of
  tList}

  \item{lam}{Per-particle birth rate}

  \item{v}{Per-particle shift rate}

  \item{mu}{Per-particle death rate}

  \item{initList}{A vector containing possible initial
  population sizes}
}
\value{
  A list of optim objects
}
\description{
  This function generates observation intervals from the
  simple birth-shift-death process without covariates, and
  generates a dataset for each observation time length in
  tList. Next, the frequent monitoring log-likelihood is
  maximized using \code{\link{FM.optim}} for each dataset.
}
\details{
  See accompanying vignette for example usage.
}

