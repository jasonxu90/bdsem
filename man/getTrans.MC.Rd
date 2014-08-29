\name{getTrans.MC}
\alias{getTrans.MC}
\title{Calculate empirical transition probabilities}
\usage{
  getTrans.MC(N, tList, lam, v, mu, initNum)
}
\arguments{
  \item{N}{An integer, the number of MC simulations per
  element of tList}

  \item{tList}{A list of observation interval lengths to
  simulate}

  \item{lam}{Per-particle birth rate}

  \item{v}{Per-particle shift rate}

  \item{mu}{Per-particle death rate}

  \item{initNum}{Integer giving the number of initial
  particles}
}
\value{
  A list of matrices, where each entry of the list
  corresponds to an element of tList. The i,j entry of each
  matrix in the list give the probability of the process
  ending with i type 1 particles and j type 2 particles,
  beginning with initNum type 1 particles, by the end of
  the corresponding observation length.
}
\description{
  Calculates transition probabilities empirically using
  Monte Carlo simulation
}
\details{
  Note: function can be modified to initialize matrices to
  be larger sized in the case where rates are large
}
\examples{
N = 500; tList = c(.5,1); lam = .2; v = .1; mu = .15; initNum = 10
getTrans.MC(N,tList,lam,v,mu,initNum)
}

