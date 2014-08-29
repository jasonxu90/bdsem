\name{sim.N.true}
\alias{sim.N.true}
\title{Simulate N realizations from the BSD process}
\usage{
  sim.N.true(N, t.end, lam, v, mu, initNum)
}
\arguments{
  \item{N}{An integer, the number of realizations to
  simulate}

  \item{t.end}{A number giving the length of time for the
  simulation}

  \item{lam}{Per-particle birth rate}

  \item{v}{Per-particle shift rate}

  \item{mu}{Per-particle death rate}

  \item{initNum}{An integer, initial total number of
  particles}
}
\value{
  A 2 by N matrix; each \emph{i}th column is a pair of
  integers giving the number of initial indices still
  present followed by the number of new indices present in
  the population in the \emph{i}th realization
}
\description{
  Wrapper for \code{\link{simulate.one.true}} that simply
  replicates the function N times
}
\examples{
sim.N.true(10,2,.2,.12,.15,10)
}

