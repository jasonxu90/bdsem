\name{de.birth}
\alias{de.birth}
\title{Expected births ODE}
\usage{
  de.birth(t, state, param)
}
\arguments{
  \item{t}{A number or vector of numbers: evaluation times
  of ODE}

  \item{state}{A named vector containing initial value of
  the state variable G, complex valued}

  \item{param}{A named vector of numbers containing the
  other arguments lam, v, mu, r, and complex number s2}
}
\value{
  The rate of change of G, the generating function for
  expected births, in list form
}
\description{
  Evaluates the ODE for generating function corresponding
  to transition probabilities. This is in a format to be
  solved using zvode from package deSolve; see deSolve
  documentation/vignette for further details
}
\examples{
t = .5; state = c(G=exp(2*1i)); param = c(lam = .2, v = .05, mu = .1, r = 3, s2 = exp(3*1i))
de.birth(t,state,param)
}

