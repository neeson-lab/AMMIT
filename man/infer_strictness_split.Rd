% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer_strictness.R
\name{infer_strictness_split}
\alias{infer_strictness_split}
\title{Infer the AMMIT strictness parameter required to derive manual threshold values from unmixed marker intensity data}
\usage{
infer_strictness_split(
  object,
  markers = "all",
  return = "strictness",
  subset = NULL
)
}
\arguments{
\item{object}{A list of SpatialExperiment objects, each with assay 'data' present.}

\item{markers}{Character; one or multiple markers found in the SpatialExperiment object for which to infer strictness, or "all" (default) to do so for all markers.}

\item{return}{Character; either 'strictness' (default) to return a list of named vectors of the inferred strictness for each marker, or 'object' to return the list of SpatialExperiment objects with this information embedded in the rowData.}

\item{subset}{Character; Which of the split subsets for which to infer strictness. If left NULL (default), all subsets will have strictness inferred.}
}
\value{
A list of: Numeric vectors of inferred strictness values (if 'strictness' specified for return), or SpatialExperiment objects with this information embedded in rowData (if 'object' specified).
}
\description{
Infer the AMMIT strictness parameter required to derive manual threshold values from unmixed marker intensity data
}
\examples{
data("ammit_spe_dummy")
spe_list <- split_spe(ammit_spe_dummy)
spe_list <- infer_manual_thresholds(spe_list, markers = "M7")
spe_list <- unmix_intensities(spe_list, markers="M7", k=3)
infer_strictness(spe_list, markers = "M7", return = "strictness")
}
