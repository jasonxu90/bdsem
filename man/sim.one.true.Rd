\name{sim.one.true}
\alias{sim.one.true}
\title{Simulate a BSD process with error catch}
\usage{
  sim.one.true(t.end, lam, v, mu, initNum)
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
  A pair of integers giving the number of initial indices
  still present followed by the number of new indices
  present in the population
}
\description{
  Wrapper for \code{\link{simulate.true}} that catches any
  possible errors
}
\examples{
sim.one.true(2,.2,.12,.15,10)
}

