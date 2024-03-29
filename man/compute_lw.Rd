% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampling.R
\name{compute_lw}
\alias{compute_lw}
\title{Calculate importance samples}
\usage{
compute_lw(
  prop_theta,
  n_particles,
  subj_est,
  dist_funcs,
  mix,
  samples,
  show = FALSE
)
}
\arguments{
\item{prop_theta}{- one proposals for theta_mu, the group level model
parameter estimates}

\item{n_particles}{- the number of particles to draw for each importance
sample.}

\item{subj_est}{- A named list containing the mean `mu` vector and covariance
matrices `sigma` of each subjects random effects as well as the number of
parameters used for the random effect estimates.}

\item{dist_funcs}{- A named list with the specific calculation for the log
density of the `prior` and `group` distribution}

\item{mix}{- The calculated values for the importance mixing vars}

\item{samples}{- The object containing original data, model design (number of
subjects, number of parameters and parameter names), the specification of
the prior and thelog likelihood function used to generate the model
estimates.}

\item{show}{- Set to TRUE to show feedback after each call}
}
\value{
The logweight of the importance samples
}
\description{
Calculate importance samples
}
