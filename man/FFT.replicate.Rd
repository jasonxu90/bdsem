\name{FFT.replicate}
\alias{FFT.replicate}
\title{Replicate the function \code{\link{FFT.run}}}
\usage{
  FFT.replicate(numReps, N, tList, lam, v, mu, initList,
    initGuess, s1.seq, s2.seq)
}
\arguments{
  \item{numReps}{The number of replications}

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

  \item{initGuess}{A vector containing the initial guess of
  birth, shift, and death rates}

  \item{s1.seq}{A vector of complex arguments evenly spaced
  along the unit circle}

  \item{s2.seq}{A vector of complex arguments evenly spaced
  along the unit circle}
}
\value{
  An array of optim objects in the same layout as
  \code{\link{FM.replicate}}
}
\description{
  Replicate the function \code{\link{FFT.run}}
}

