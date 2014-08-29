\name{getStandardErrors}
\alias{getStandardErrors}
\title{Calculate numerical standard errors of estimated coefficients}
\usage{
  getStandardErrors(betas, t.pat, num.patients,
    PATIENTDATA, patients.design, s1.seq, s2.seq)
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
  A vector the length of the coefficient vector, giving
  standard errors for each estimated coefficient
}
\description{
  \code{getStandardErrors} uses \code{optimHess} to
  numerically compute the hessian at the value of the
  estimated MLE, which can be obtained using
  \code{\link{EM.run}}. The function inverts the hessian
  and takes the square root of diagonal entries to yield
  standard errors for each coefficient estimate.
}

