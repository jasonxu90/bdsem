\name{logFFT}
\alias{logFFT}
\title{Calculate log-likelihood using generating function technique}
\usage{
  logFFT(param, u, simData, initList, s1.seq, s2.seq)
}
\arguments{
  \item{param}{A vector of three numbers containing the
  birth, shift, and death rate respectively}

  \item{u}{A number, the observation interval length}

  \item{simData}{An matrix in the format of a list entry
  returned by \code{\link{makedata.simple}}}

  \item{initList}{A vector of possible initial populations}

  \item{s1.seq}{A vector of complex arguments evenly spaced
  along the unit circle}

  \item{s2.seq}{A vector of complex arguments evenly spaced
  along the unit circle}
}
\value{
  The negative log-likelihood value
}
\description{
  \code{logFFT} calculates the log-likelihood where
  transition probabilities are computed using our FFT
  generating function method
}

