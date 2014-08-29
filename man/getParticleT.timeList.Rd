\name{getParticleT.timeList}
\alias{getParticleT.timeList}
\title{Compute expected particle time over a grid at a list of evaluation times}
\usage{
  getParticleT.timeList(tList, lam, v, mu, initNum, s1.seq,
    s2.seq, dt)
}
\arguments{
  \item{tList}{A vector of numbers corresponding to the
  desired evaluation times}

  \item{dt}{A number giving the increment length used in
  solving the ODE}

  \item{s1.seq}{A vector of complex numbers; initial values
  of the ODE G}

  \item{s2.seq}{A vector of complex numbers as inputs of
  s2.seq}

  \item{initNum}{An integer giving the number of initial
  particles}

  \item{lam}{Birth rate, must be positive}

  \item{v}{Shift rate, must be positive}

  \item{mu}{Death rate}
}
\value{
  A list of matrices of dimension length(s1.seq) by
  length(s2.seq); each list entry corresponds to an
  evaluation time from tList
}
\description{
  Wrapper that transforms the output matrices from
  \code{\link{makeGrid.particleT.r0}} and
  \code{\link{makeGrid.particleT.partial}} to a list of
  expected sufficient statistics matrices using \code{fft}.
  Does this at several input times given in tList. Returns
  a list of matrices, where each list entry corresponds to
  a time in tList. The i,j entry of each matrix corresponds
  to the expected particle time spent with i type 1
  particles and j type 2 particles, beginning with initNum
  type 1 particles, over the corresponding time length.
}
\examples{
tList = c(1,2);  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
gridLength = 32
s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
getParticleT.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)
}

