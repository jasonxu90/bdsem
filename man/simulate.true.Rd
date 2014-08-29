\name{simulate.true}
\alias{simulate.true}
\title{Simulate a birth-shift-death-process}
\usage{
  simulate.true(t.end, lam, v, mu, initNum)
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
  present in the population. Otherwise, returns '999' as an
  error code
}
\description{
  This function simulates one realization of a
  birth-shift-death CTMC with per-particle birth, death,
  and shift rates lambda, nu, and mu.
}
\details{
  The process begins with initNum particles, each assigned
  a random location indexed by a number between 1 and
  100,000. Simulation occurs until time t.end. The process
  returns an ordered pair giving the number of initial
  locations (indices) still occupied, followed by the
  number of new indices present, by t.end.
}
\examples{
simulate.true(2,.2,.12,.15,10)
}

