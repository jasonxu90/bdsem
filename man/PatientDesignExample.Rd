\name{PatientDesignExample}
\alias{PatientDesignExample}
\title{Create example design matrix in the format required for \code{\link{MakePatientRates}}}
\usage{
  PatientDesignExample(num.patients)
}
\arguments{
  \item{num.patients}{An integer giving the number of
  patients in the design matrix}
}
\value{
  A 2 by num.patients matrix: the first row is constant 1's
  corresponding to intercepts, and second row is a
  patient-specific covariate sampled uniformly between 6
  and 10
}
\description{
  This function creates one possible design matrix as an
  input for the function \code{\link{MakePatientRates}}
  This particular design matrix features one covariate
  uniformly sampled between 6 and 10, as well as a row of
  1's corresponding to an intercept term. The design matrix
  is thus 2 by number of patients
}

