% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot.spaghetti.R
\name{ggplot_spaghetti}
\alias{ggplot_spaghetti}
\title{Spaghetti Plots using \code{ggplot2}}
\usage{
ggplot_spaghetti(
  y,
  id,
  time,
  alpha = 0.2,
  method = "loess",
  jit = 0,
  group = NULL
)
}
\arguments{
\item{y}{This is the y-axis parameter to specify.
Generally it is a continuous variable.}

\item{id}{This is the id parameter that identifies the unique individuals or
units.}

\item{time}{This is the time vector and must be numeric.}

\item{alpha}{Scalar value between [0,1] that specifies the transparencey
of the lineplots.}

\item{method}{Character value that specifies which type of method to use for
fitting. Optional methods come from \code{\link[ggplot2]{geom_smooth}}
function.}

\item{jit}{Scalar value that specifies how much you want to jitter each
individual observation. Useful if many of the values share the same y values
at a time point.}

\item{group}{Specifies a grouping variable to be used, and will plot it by
color on one single plot.}
}
\value{
Plots a time series data by each individual/unit with group trends
overlayed.
}
\description{
This function allows the user to create spaghetti plots for individuals
with time varying covariates. You can also break this down into subgroups
to analyze different trentds.
}
\details{
Note that the data must be in long format.
}
\examples{
library(ggplot2)
num_subjects_per_group <- 15
sim_obj <- mvrnorm_sim(n_control=num_subjects_per_group,
                       n_treat=num_subjects_per_group,
                       control_mean=5, sigma=1, num_timepoints=5,
                       t_interval = c(0, 4),
                       rho=0.95, corr_str='ar1', func_form='linear',
                       beta=c(0, 0.25),
                       missing_pct=0.6, missing_per_subject=2)

with(sim_obj$df, suppressWarnings(ggplot_spaghetti(y=Y_obs, id=ID, time=time,
                                                  jit=0.1, group=group)))+
  labs(title="Simulated Microbiome Data from Multivariate Normal",
       y="Normalized Reads", x="Time") +
  scale_linetype_manual(values=c("solid","dashed"), name="Group") +
  scale_color_manual(values=c("#F8766D", "#00BFC4"), name="Group")

}
