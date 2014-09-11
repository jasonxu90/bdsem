\name{FM.replicate}
\alias{FM.replicate}
\title{Replicate \code{\link{FM.run}}}
\usage{
  FM.replicate(numReps, N, tList, lam, v, mu, initList,
    initGuess)
}
\arguments{
  \item{numReps}{The number of desired replications}

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

  \item{initGuess}{Vector of numbers, initial guess for
  \code{optim}}
}
\value{
  An array with entries of type returned by
  \code{\link{FM.run}}. Rows correspond to a dt value in
  tList, and columns correspond to replications
}
\description{
  Replicate \code{\link{FM.run}}
}
\examples{
tList <- c(.2,.4,.6); initList <- c(1:15)
lam = .06; v = .02; mu = .11
trueParam <- c(lam,v,mu)
N = 100; numReps = 2
example <- FM.replicate(numReps,N,tList,lam,v,mu,initList, trueParam)
}

