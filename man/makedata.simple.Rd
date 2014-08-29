\name{makedata.simple}
\alias{makedata.simple}
\title{Generate synthetic dataset for BSD process with simple rates}
\usage{
  makedata.simple(N, tList, lam, v, mu, initList)
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
}
\value{
  A list of matrices. Each entry in the list is in the
  format returned by \code{\link{ransim.N.true}}, and
  corresponds to equal observation times specified in the
  corresponding entry in tList
}
\description{
  \code{makedata.simple} creates a list of synthetic
  observed datasets, each with N observation intervals from
  the BSD process. Observation interval lengths are entries
  in tList. This is used in simulation studies checking
  inferential procedures. Simulates observations using
  \code{\link{ransim.N.true}}.
}
\examples{
makedata.simple(30,c(.2,.4,.6),.2,.12,.15,seq(8,12))
}

