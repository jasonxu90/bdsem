\name{logFFT.patients}
\alias{logFFT.patients}
\title{Log likelihood for BSD process with covariates}
\usage{
  logFFT.patients(betas, t.pat, num.patients, PATIENTDATA,
    patients.design, s1.seq, s2.seq)
}
\arguments{
  \item{betas}{A vector of numbers \eqn{\mathbf\beta =
  (\mathbf\beta^\lambda, \mathbf\beta^\nu,
  \mathbf\beta^\mu)}}

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
  The negative log likelihood of the observations in
  PATIENTDATA, given rates determined by coefficient vector
  betas and covariate values in patients.design
}
\description{
  \code{logFFT.patients} evaluates the log likelihood of a
  dataset with observations corresponding to "patients" in
  the setting where rates of the process depend on
  patient-specific covariates. The transition probabilities
  given the states of the process at endpoints of each
  observation interval are computed using the
  FFT/generating function method, relying on
  \code{\link{getTrans.initList}}. The log likelihood is
  then the sum of these log transition probabilities.
}
\details{
  Note: this function is used so that MLE estimation of the
  coefficient vector can be accomplished, i.e. using
  \code{optim}, and is also used in numerically computing
  standard errors at the MLE.

  Vectors s1.seq and s2.seq should be of length greater
  than the total number of particles of either type at any
  observation interval
}

