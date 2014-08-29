\name{makeGrid.death.partial}
\alias{makeGrid.death.partial}
\title{Evaluates partial derivative of expected deaths ODE over 2D grid}
\usage{
  makeGrid.death.partial(time, dt, s1.seq, s2.seq, lam, v,
    mu)
}
\arguments{
  \item{time}{A number corresponding to the desired
  evaluation time of ODEs}

  \item{dt}{A number giving the increment length used in
  solving the ODE}

  \item{s1.seq}{A vector of complex numbers; initial values
  of the ODE G}

  \item{s2.seq}{A vector of complex numbers as inputs of
  s2.seq}

  \item{lam}{Birth rate}

  \item{v}{Shift rate}

  \item{mu}{Death rate}
}
\value{
  A matrix of dimension length(s1.seq) by length(s2.seq) of
  the function values
}
\description{
  Wrapper that numerically partially differentiates the
  solution given in \code{\link{solve.death}} over a grid
  of input values s1, s2 at one fixed time, and r=1
}
\examples{
time = 5;  dt = 5; lam = .5; v = .2; mu = .4
gridLength = 32
s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
makeGrid.death.partial(time,dt,s1.seq,s2.seq,lam,v,mu)
}

