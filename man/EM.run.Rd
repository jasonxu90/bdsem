\name{EM.run}
\alias{EM.run}
\title{Inference via EM algorithm}
\usage{
  EM.run(betaInit, t.pat, num.patients, PATIENTDATA,
    patients.design, s1.seq, s2.seq, relTol)
}
\arguments{
  \item{betaInit}{A vector, the initial guess for
  coefficients beta}

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

  \item{relTol}{A number, the relative convergence
  criterion}
}
\value{
  A list containing the log-likelihood value at
  convergence, the final beta estimate, and the number of
  iterations
}
\description{
  This function runs the EM algorithm from an initial
  guess. Infers the coefficient vector in setting where
  rates depend on patient-specific covariates. The EM
  algorithm alternates between calling \code{\link{ESTEP}}
  and \code{\link{MSTEP}} until the change in observed
  log-likelihood changes less than a specified relative
  tolerance between iterations
}

