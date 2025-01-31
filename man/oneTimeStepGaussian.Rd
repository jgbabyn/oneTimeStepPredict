% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneTSP.R
\name{oneTimeStepGaussian}
\alias{oneTimeStepGaussian}
\title{Find the residuals for one time step assuming the conditional distribution follows a normal distribution}
\usage{
oneTimeStepGaussian(
  ts,
  obj,
  t_index,
  data.vector.indicators,
  observation.names,
  cholesky = FALSE
)
}
\arguments{
\item{ts}{the current time step}

\item{obj}{the newly created obj with each observation as a parameter}

\item{data.vector.indicators}{the names of the DATA INDICATORS in the TMB template}

\item{observation.names}{the names corresponding to the observations}

\item{cholesky}{whether or not to use cholesky factorization or the inverse square root to do the MV normal standardization}
}
\description{
Find the residuals for one time step assuming the conditional distribution follows a normal distribution
}
