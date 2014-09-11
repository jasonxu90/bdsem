\name{FFT.pat.optim}
\alias{FFT.pat.optim}
\title{Optimizes the function \code{\link{logFFT.patients}} using \code{optim} package}
\usage{
  FFT.pat.optim(betaInit, t.pat, num.patients, PATIENTDATA,
    patients.design, s1.seq, s2.seq, tol, max)
}
\arguments{
  \item{betaInit}{A vector, the initial guess for the
  algorithm}

  \item{t.pat}{A number, the observation interval length}

  \item{num.patients}{An integer, number of unique
  patients}

  \item{PATIENTDATA}{A matrix in the form returned by
  \code{\link{MakePatientData}} containing the set of
  observation intervals}

  \item{patients.design}{A design matrix in the same form
  as returned by \code{\link{PatientDesignExample}}}

  \item{s1.seq}{A vector of complex arguments evenly spaced
  along the unit circle}

  \item{s2.seq}{A vector of complex arguments evenly spaced
  along the unit circle}

  \item{tol}{A number for setting the relative tolerance
  for the algorithm (the reltol argument in \code{optim})}

  \item{max}{An integer, the max number of iterations
  before termination}
}
\value{
  An \code{optim} type object
}
\description{
  Function uses Nelder-Mead optimization as implemented in
  \code{optim} to maximize the log likelihood function
  \code{\link{logFFT.patients}}.
}
\details{
  Example code is not included here because of runtime; see
  vignette for tutorial on using this function.
}

