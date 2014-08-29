\name{sim.N.eventcount}
\alias{sim.N.eventcount}
\title{Simulate N BSD processes and return sufficient statistics}
\usage{
  sim.N.eventcount(N, t.end, lam, v, mu, initNum)
}
\arguments{
  \item{N}{The number of replications, an integer}

  \item{t.end}{A number giving the length of time for the
  simulation}

  \item{lam}{Per-particle birth rate}

  \item{v}{Per-particle shift rate}

  \item{mu}{Per-particle death rate}

  \item{initNum}{An integer, initial total number of
  particles}
}
\value{
  A 4 by N matrix, where rows correspond to total copies,
  shifts, deaths, and particle time per realization,
  respectively. Each column corresponds to one realization
  of the process.
}
\description{
  Wrapper that replicates
  \code{\link{simulate.one.eventcount}} N times. Each
  replication begins with initNum particles.
}
\examples{
sim.N.eventcount(10,2,.2,.12,.15,10)
}

