\name{getTrans.FreqMon}
\alias{getTrans.FreqMon}
\title{Calculates the four transition probabilities defined under Frequent Monitoring}
\usage{
  getTrans.FreqMon(tList, lam, v, mu, initNum)
}
\arguments{
  \item{tList}{A vector of observation interval lengths}

  \item{lam}{Per-particle birth rate}

  \item{v}{Per-particle shift rate}

  \item{mu}{Per-particle death rate}

  \item{initNum}{Integer, the number of initial particles}
}
\value{
  A list of \eqn{2 \times 2} matrices, each containing the
  probabilitiy of one birth, one shift, one death, or no
  event occurring over a corresponding observation interval
  length from tList
}
\description{
  \code{getTrans.FreqMon} uses simple closed form
  expressions to compute transition probabilities available
  under the Frequent Monitoring assumption, which allows at
  most one event to occur per observation interval.
  Computes the four possible transition possibilities for
  intervals with lengths corresponding to entries in tList
}
\examples{
getTrans.FreqMon(c(.5,1,10),.3,.1,.2,10)
}

