% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conversion_functions.R
\name{simulate2phyloseq}
\alias{simulate2phyloseq}
\title{Convert simulated output to phyloseq object}
\usage{
simulate2phyloseq(obj, missing = FALSE)
}
\arguments{
\item{obj}{output from either \code{\link{gen_norm_microbiome}} or
\code{\link{mvrnorm_sim}}}

\item{missing}{logical indicator for objects from \code{\link{mvrnorm_sim}}.
If missing = TRUE then create MRexperiment object with \code{Y_obs} else
use \code{Y}.}
}
\value{
A phyloseq object
}
\description{
This function will convert simulated data into a \code{\link{phyloseq}}
object.
}
\examples{
bug_gen <- gen_norm_microbiome(features=6, diff_abun_features=3,
                               n_control=30, n_treat=20, control_mean=2,
                               sigma=2, num_timepoints=4, t_interval=c(0, 3),
                               rho=0.9, corr_str="compound", func_form="M",
                               beta=c(4, 3), IP=c(2, 3.3, 6),
                               missing_pct=0.2, missing_per_subject=2,
                               miss_val=0, asynch_time=TRUE)
bug_gen_phyloseq <- simulate2MRexperiment(bug_gen)
class(bug_gen_phyloseq)

}
