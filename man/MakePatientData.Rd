\name{MakePatientData}
\alias{MakePatientData}
\title{Generates synthetic patient data for inference in discretely observed BSD process with covariates}
\usage{
  MakePatientData(t, num.patients, patients.rates)
}
\arguments{
  \item{t}{A number giving the observation interval length}

  \item{num.patients}{An integer, the number of synthetic
  patients}

  \item{patients.rates}{A matrix in the format returned by
  \code{\link{MakePatientRates}}}
}
\value{
  A \eqn{5 \times m} matrix where m is the number of total
  observation intervals generated. The first row
  corresponds to the patient ID, the second row gives the
  initial number of particles for that observation
  interval, the third through fifth column are in the
  format returned by \code{\link{ransim.N.true}}
}
\description{
  \code{MakePatientData} is the main function for
  generating a synthetic dataset with covariates. It
  simulates observations from a birth-shift-death process
  for a number of synthetic "patients", using
  \code{\link{ransim.N.true}}. This provides data for
  simulation studies toward assessing inference in a panel
  data setting, with rates depending on multiple
  covariates.
}
\details{
  Each observation interval has a common fixed length given
  by the argument t. For each patient, between 2 and 6
  observation intervals are generated, and each begins with
  between 2 and 12 initial particles. These numbers are
  initialized uniformly at random, and passed in as
  arguments to \code{\link{ransim.N.true}}, with the true
  rates set to the corresponding column of patients.rates.
}
\examples{
num.patients = 10; t = .4
patients.design <- PatientDesignExample(num.patients)
beta.lam <- c(log(8), log(.6)); beta.v <- c( log(.5), log(.7)); beta.mu <- c(log(.8), log(.8))
betas <- c(beta.lam,beta.v,beta.mu)
patients.rates <- MakePatientRates(patients.design, betas)
MakePatientData(t, num.patients, patients.rates)
}

