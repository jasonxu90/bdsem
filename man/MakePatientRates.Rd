\name{MakePatientRates}
\alias{MakePatientRates}
\title{Makes matrix of patient-specific birth, shift, and death rates}
\usage{
  MakePatientRates(patients.design, betaVec)
}
\arguments{
  \item{patients.design}{A \eqn{n \times m} matrix, where n
  is number of covariates (including intercept) in the
  model and m is number of patients}

  \item{betaVec}{The vector \eqn{\mathbf\beta} of
  regression coefficients. Note this function assumes that
  \eqn{ \mathbf\beta^\lambda, \mathbf\beta^\nu,
  \mathbf\beta^\mu} are equal in length. That is, each rate
  depends on the same number of covariates}
}
\value{
  A \eqn{3 \times m} matrix where m is the number of
  patients, and rows correspond to birth, shift, and death
  rates respectively.
}
\description{
  This function returns a matrix containing the birth,
  shift, and death rates of each patient. The rates are
  calculated based on the log-linear relationship between
  regression coefficients \eqn{\mathbf{\beta} =
  (\mathbf{\beta^\lambda}, \mathbf{\beta}^\nu,
  \mathbf{\beta}^\mu)} and covariates \eqn{\mathbf{z_p}} in
  the \emph{p}th column of the design matrix given by
  \eqn{\log(\lambda_p) = \mathbf\beta^\lambda \cdot
  \mathbf{z_p}, \log(\nu_p) = \mathbf\beta^\nu \cdot
  \mathbf{z_p}, \log(\mu_p) = \mathbf\beta^\mu \cdot
  \mathbf{z_p} }
}
\examples{
num.patients = 100
patients.design <- PatientDesignExample(num.patients)
beta.lam <- c(log(8), log(.6)); beta.v <- c( log(.5), log(.7)); beta.mu <- c(log(.8), log(.8))
betas <- c(beta.lam,beta.v,beta.mu)
MakePatientRates(patients.design, betas)
}

