\name{FFT.run}
\alias{FFT.run}
\title{Generate synthetic data from simple BSD process and infer rates using generating function method}
\usage{
  FFT.run(N, tList, lam, v, mu, initList, initGuess,
    s1.seq, s2.seq)
}
\arguments{
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
  A list of optim objects
}
\description{
  The main function for simulation studies assessing our
  generating function approach in the the discretely
  observed simple birth-shift-death-process without
  covariates. Generates synthetic datasets using
  \code{\link{makedata.simple}}, and infers the MLE rates
  for each using \code{\link{FFT.optim}}.
}

