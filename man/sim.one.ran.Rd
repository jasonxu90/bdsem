\name{sim.one.ran}
\alias{sim.one.ran}
\title{Simulate a BSD process, returning initial number as well as end state}
\usage{
  sim.one.ran(t.end, lam, v, mu, initNum)
}
\arguments{
  \item{t.end}{A number giving the length of time for the
  simulation}

  \item{lam}{Per-particle birth rate}

  \item{v}{Per-particle shift rate}

  \item{mu}{Per-particle death rate}

  \item{initNum}{An integer, initial total number of
  particles}
}
\value{
  A vector containing the initial number of particles,
  followed by a pair of integers giving the number of
  initial indices still present followed by the number of
  new indices present in the population
}
\description{
  Wrapper for \code{\link{simulate.true}} that catches any
  possible errors and returns the initial number.
}
\examples{
sim.one.ran(2,.2,.12,.15,10)
}

