\name{ESTEP}
\alias{ESTEP}
\title{Perform one E-step of the EM algorithm}
\usage{
  ESTEP(betaVec, t.pat, num.patients, PATIENTDATA,
    patients.design, s1.seq, s2.seq)
}
\arguments{
  \item{betaVec}{A vector, the setting of beta
  coefficients}

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
}
\value{
  A list containing a matrix of the expected sufficient
  statistics as well as the observed log likelihood value
}
\description{
  \code{ESTEP} performs one E-step of the EM algorithm,
  computing expected sufficient statistics given current
  settings of the parameters. This function is the
  "accelerated" version, meaning that intervals with no
  observed changes are computed more efficiently using
  closed form expressions, bypassing generating function
  and FFT calculations on these intervals.
}

