% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oneTSP.R
\name{oneTimeStepPredict}
\alias{oneTimeStepPredict}
\title{Find the standardized residuals predicting one time step at a time}
\usage{
oneTimeStepPredict(
  obj,
  observation.names = NULL,
  data.vector.indicators = NULL,
  time.step.indicators = NULL,
  method = c("oneTimeStepGaussian"),
  cholesky = FALSE
)
}
\arguments{
\item{obj}{the optimized TMB model object to find the residuals of}

\item{observation.names}{the names corresponding to the observations}

\item{data.vector.indicators}{the names of the DATA INDICATORS in the TMB template}

\item{time.step.indicators}{vector or list of vectors giving the time steps corresponding to the observations/data.vector.indicators}

\item{method}{Which method to use, currently only oneTimeStepGaussian is supported}

\item{cholesky}{whether or not to use cholesky factorization (order dependent) when standardizing the residuals.}
}
\description{
oneTimeStepPredict is similar to TMB's built-in oneStepPredict except oneTimeStepPredict
predicts all observations on the same time step simulataneously rather than one observation at a time.
With the right choice of standardization of the residuals, order doesn't matter unlike OSA residuals.
}
