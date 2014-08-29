\name{logFM}
\alias{logFM}
\title{Calculate approximate log-likelihood based on Frequent Monitoring computations}
\usage{
  logFM(param, u, FM)
}
\arguments{
  \item{param}{A vector of three numbers containing the
  birth, shift, and death rate respectively}

  \item{u}{A number, the observation interval length}

  \item{FM}{An object returned by \code{\link{FM.data}}
  containing the number of each event under FM}
}
\value{
  The negative frequent monitoring log-likelihood
}
\description{
  \code{logFM} calculates the log-likelihood where
  transition probabilities are computed under the frequent
  monitoring assumption. Used for simulation studies
  comparing results using frequent monitoring with our
  methods.
}

