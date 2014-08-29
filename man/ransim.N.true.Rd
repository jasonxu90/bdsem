\name{ransim.N.true}
\alias{ransim.N.true}
\title{Simulate N realizations from the BSD process, each with random initial number each time}
\usage{
  ransim.N.true(N, t.end, lam, v, mu, range)
}
\arguments{
  \item{N}{An integer, the number of realizations to
  simulate}

  \item{range}{A vector containing possible initial
  populations, typically a sequence of integers}

  \item{t.end}{A number giving the length of time for the
  simulation}

  \item{lam}{Per-particle birth rate}

  \item{v}{Per-particle shift rate}

  \item{mu}{Per-particle death rate}
}
\value{
  A 3 by N matrix; each \emph{i}th column corresponds to
  the \emph{i}th realization. The first row contains
  initial numbers, the second row contains the number of
  original indices still present in the population by
  t.end, and the third row contains the number of new
  indices present.
}
\description{
  Wrapper for \code{\link{simulate.one.ran}} that simply
  replicates the function N times. Each replication begins
  with a random number of initial particles uniformly
  generated from a specified range.
}
\examples{
ransim.N.true(10,2,.2,.12,.15, seq(7,15))
}

