% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvrnorm_sim.R
\name{mvrnorm_corr_gen}
\alias{mvrnorm_corr_gen}
\title{Generate Multivariate Random Normal Longitudinal Data}
\usage{
mvrnorm_corr_gen(
  n,
  obs,
  t,
  mu,
  sigma,
  rho,
  corr_str = c("ar1", "compound", "ind"),
  zero_trunc = TRUE
)
}
\arguments{
\item{n}{integer scalar representing the total number of individuals}

\item{obs}{vector of length \code{n} specifying the number of observations
per indivdiual.}

\item{t}{vector corresponding to the timepoints for each individual.}

\item{mu}{vector specifying the mean value for individuals.}

\item{sigma}{scalar specifying the standard deviation for
all observations.}

\item{rho}{numeric scalar value between [0, 1] specifying the amount of
correlation between. assumes that the correlation is consistent for all
subjects.}

\item{corr_str}{character value specifying the correlation structure.
Currently available methods are \'ar1\', \'compound\', and \'ind\' which
correspond to first-order autoregressive, compound or equicorrelation,
and independence respecitvely.}

\item{zero_trunc}{logical value to specifying whether the generating
distribution should come from a multivariate zero truncated normal or an
untruncated multivariate normal. by default we assume that zero truncation
occurs since this is assummed in our microbiome setting.}
}
\value{
This function returns a list with the following objects:

\code{df} - data.frame object with complete outcome \code{Y}, subject ID,
time, group, and outcome with missing data

\code{Y} - vector of complete outcome

\code{Mu} - vector of complete mean specifications used during simulation

\code{Sigma} - block diagonal symmetric matrix of complete data used during
simulation

\code{N} - total number of observations
}
\description{
For this methodology we assume that we draw a set of \code{n} independent each
with \eqn{q_{i}} observations.
}
\examples{
size <- 15
reps <- 4
N <- size*reps
mvrnorm_corr_gen(n=size, obs=rep(reps, size), t=rep(seq_len(4), size),
mu=rep(1, N), sigma=2, rho=0.9, corr_str="ar1")

}
