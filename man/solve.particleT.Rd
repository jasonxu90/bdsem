\name{solve.particleT}
\alias{solve.particleT}
\title{Expected births ODE solver}
\usage{
  solve.particleT(time, dt, s1, s2, r, lam, v, mu)
}
\arguments{
  \item{time}{A number corresponding to the desired
  evaluation time of ODE}

  \item{dt}{A number giving the increment length used in
  solving the ODE}

  \item{s1}{A complex number giving the initial value of
  the ODE G}

  \item{s2}{A complex number}

  \item{lam}{Birth rate}

  \item{v}{Shift rate}

  \item{mu}{Death rate}

  \item{r}{a real number}
}
\value{
  The function value of the particle time generating
  function
}
\description{
  Wrapper that numerically solves the ODE defined in
  \code{\link{de.particleT}} using \code{zvode} from
  package \code{deSolve} given evaluation time, initial
  state values s1 and s2, and birth/shift/death rates
}
\examples{
time = 5;  dt = 1; s1 = exp(2*1i); s2 = exp(3*1i); lam = .5; v = .2; mu = .4; r = 3
solve.particleT(time,dt,s1,s2,r,lam,v,mu)
}

