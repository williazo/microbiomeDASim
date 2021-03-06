% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_norm_microbiome.R
\name{final_output_gen}
\alias{final_output_gen}
\title{Generating the final combined bug output}
\usage{
final_output_gen(
  no_diff_feat,
  diff_abun_features,
  diff_Y,
  null_Y,
  diff_bugs,
  nodiff_bugs,
  final_output = NULL
)
}
\arguments{
\item{no_diff_feat}{number of non differentially abundant features}

\item{diff_abun_features}{number of differentially abundant features}

\item{diff_Y}{simulated outcome for differentially abundant features}

\item{null_Y}{simulated outcome for non differentially abundant features}

\item{diff_bugs}{sample information for differentially abundant features}

\item{nodiff_bugs}{sample information for non differentially abundant
features}

\item{final_output}{final object that will store the simulated data}
}
\value{
final output list with the OTU table and corresponding bug feature data.frame
}
\description{
Generating the final combined bug output
}
\keyword{internal}
