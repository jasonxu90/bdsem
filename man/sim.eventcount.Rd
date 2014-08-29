\name{sim.eventcount}
\alias{sim.eventcount}
\title{Simulate a BSD process and return sufficient statistics}
\usage{
  sim.eventcount(t.end, lam, v, mu, initNum)
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
  A vector of four numbers giving total copies, shifts,
  deaths, and particle time, respectively. Returns the
  integer '999' as an error code if error occurs.
}
\description{
  Simulates a birth-shift-death process and records the
  number of births, deaths, and shifts, as well as total
  particle time, over the observation interval until time
  t.end. Used toward simulation functions that check
  accuracy of expected restricted moment calculations
}
\examples{
sim.eventcount(2,.2,.12,.15,10)
}

