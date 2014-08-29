\name{getShiftMeans.initList}
\alias{getShiftMeans.initList}
\title{Compute expected shifts over a grid at a list of initial particle counts}
\usage{
  getShiftMeans.initList(u, initList, lam, v, mu, s1.seq,
    s2.seq, dt)
}
\arguments{
  \item{initList}{A vector of integers corresponding to the
  desired initial particle counts}

  \item{dt}{A number giving the increment length used in
  solving the ODE}

  \item{s1.seq}{A vector of complex numbers; initial values
  of the ODE G}

  \item{s2.seq}{A vector of complex numbers as inputs of
  s2.seq}

  \item{u}{A number giving the observation interval length,
  equivalently the time to evaluate the ODEs}

  \item{lam}{Birth rate, must be positive}

  \item{v}{Shift rate, must be positive}

  \item{mu}{Death rate}
}
\value{
  A list of matrices of dimension length(s1.seq) by
  length(s2.seq); each list entry corresponds to an initial
  number of particles from initList
}
\description{
  Wrapper that transforms the output matrices from
  \code{\link{makeGrid.shift.r1}} and
  \code{\link{makeGrid.shift.partial}} to a list of
  expected sufficient statistics matrices using \code{fft}.
  Does this at several initial particle counts from
  initList Returns a list of matrices, where each list
  entry corresponds to a number of initial particles The
  i,j entry of each matrix corresponds to the expected
  shifts when the process has i type 1 particles and j type
  2 particles, beginning with initNum type 1 particles, by
  the end of corresponding time interval.
}
\examples{
initList = c(10,11)
#gives matrices of expected sufficient statistics corresponding to 10 initial particles and 11 initial particles
u = 1; dt = 1; lam = .5; v = .2; mu = .4
gridLength = 32
s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
getShiftMeans.initList(u, initList, lam, v, mu, s1.seq, s2.seq, dt)
}

