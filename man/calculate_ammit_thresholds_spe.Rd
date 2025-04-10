% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_ammit_thresholds.R
\name{calculate_ammit_thresholds_spe}
\alias{calculate_ammit_thresholds_spe}
\title{Derive the AMMIT threshold from unmixed marker intensity data and a given strictness parameter.}
\usage{
calculate_ammit_thresholds_spe(object, markers, strictness, return = "object")
}
\arguments{
\item{object}{A SpatialExperiment object, with assay 'data' present.}

\item{markers}{Character; one or multiple markers found in the SpatialExperiment object for which to derive thresholds, or "all" (default) to do so for all markers.}

\item{strictness}{Numeric; one or multiple strictness values matching to each of the marker(s) specified.}

\item{return}{Character; either 'thresholds' to return a named vector of the derived thresholds for each marker, or 'object' (default) to return the SpatialExperiment object with this information embedded in the rowData.}
}
\value{
A numeric vector of derived thresholds (if 'thresholds' specified for return), or a SpatialExperiment object with this information embedded in rowData (if 'object' specified).
}
\description{
Derive the AMMIT threshold from unmixed marker intensity data and a given strictness parameter.
}
\examples{
data("ammit_spe_dummy")
ammit_spe_dummy <- infer_manual_thresholds(ammit_spe_dummy, markers = "M7")
ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=3)
strictness <- infer_strictness(ammit_spe_dummy, markers = "M7", return = "strictness")
calculate_ammit_thresholds(ammit_spe_dummy, markers = "M7", strictness = strictness, return = "thresholds")
}
