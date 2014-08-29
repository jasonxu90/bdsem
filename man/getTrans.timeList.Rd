\name{getTrans.timeList}
\alias{getTrans.timeList}
\title{Compute transition probabilities over a grid at a list of evaluation times}
\usage{
  getTrans.timeList(tList, lam, v, mu, initNum, s1.seq,
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
  length(s2.seq), each list entry corresponds to an
  evaluation time from tList, each matrix is a matrix of
  transition probabilities
}
\description{
  Wrapper that transforms the output matrix from
  \code{\link{makeGrid.trans}} to interpretable transition
  probabilities using \code{fft}. Does this at several
  input times given in tList. Returns a list of matrices,
  where each list entry corresponds to a time in tList. The
  i,j entry of each matrix corresponds to the probability
  of transitioning from initNum type 1 particles and 0 type
  2 particles to i type 1 particles, j type 2 particles,
  over the corresponding time length.
}
\examples{
tList = c(1,2);  dt = 1; lam = .5; v = .2; mu = .4; initNum = 10
gridLength = 32
s1.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
s2.seq <- exp(2*pi*1i*seq(from = 0, to = (gridLength-1))/gridLength)
getTrans.timeList(tList, lam, v, mu, initNum, s1.seq, s2.seq, dt)
}

