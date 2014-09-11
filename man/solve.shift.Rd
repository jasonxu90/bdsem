\name{solve.shift}
\alias{solve.shift}
\title{Expected shifts ODE solver}
\usage{
  solve.shift(time, dt, s1, s2, r, lam, v, mu)
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

  \item{r}{A real number}
}
\value{
  The function value of the shifts generating function
}
\description{
  Numerically solves the ODE defined in
  \code{\link{de.shift}} using \code{zvode} from package
  \code{deSolve} given evaluation time, initial state
  values s1 and s2, and birth/shift/death rates
}
\examples{
time = 5;  dt = 1; s1 = exp(2*1i); s2 = exp(3*1i); lam = .5; v = .2; mu = .4; r = 2
solve.shift(time,dt,s1,s2, r, lam,v,mu)
}

