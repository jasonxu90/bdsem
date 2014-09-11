\name{getShiftMeans.timeList}
\alias{getShiftMeans.timeList}
\title{Compute expected shifts over a grid at a list of evaluation times}
\usage{
  getShiftMeans.timeList(tList, lam, v, mu, initNum,
    s1.seq, s2.seq, dt)
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

  \item{lam}{Per-particle birth rate}

  \item{v}{Per-particle shift rate}

  \item{mu}{Per-particle death rate}
}
\value{
  A list of matrices of dimension length(s1.seq) by
  length(s2.seq); each list entry corresponds to an
  evaluation time from tList
}
\description{
  Transforms the output matrices from
  \code{\link{makeGrid.shift.r1}} and
  \code{\link{makeGrid.shift.partial}} to a list of
  expected sufficient statistics matrices using \code{fft}.
  Does this at several input times given in tList. Returns
  a list of matrices, where each list entry corresponds to
  a time in tList. The i,j entry of each matrix corresponds
  to the number of expected shifts when process has i type
  1 particles and j type 2 particles, beginning with
  initNum type 1 particles, at the end of corresponding
  time interval.
}
\examples{
tList = c(1,2);  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
gridLength = 32
s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
getShiftMeans.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)
}

