#' Title
#'
#' @param spe A SpatialExperiment object with assays 'data' and 'reference' present.
#' @param markers A character vector of markers matching, at least partially, to the rownames of your SpatialExperiment object.
#'
#' @return The input SpatialExperiment object, with manual thresholds for each of the elements in markers added into the rowData under manual_thresholds
#' @export
#'
#' @examples
infer_manual_thresholds <- function(spe, markers) {

  if (!methods::is(spe, "SpatialExperiment")) {
    stop("Your `spe` object is not a SpatialExperiment object!")
    }
  if (!"reference" %in% names(SummarizedExperiment::assays(spe))) {
    stop("No reference data found in your SpatialExperiment object. Have you added reference data?")
    }
  if (!"data" %in% names(SummarizedExperiment::assays(spe))) {
    stop("No intensity data found in your SpatialExperiment object. Have you added intensity data?")
    }
  if (!any(markers %in% rownames(spe))) {
    stop("None of the specified markers were found in your SpatialExperiment object!")
    }
  if (!all(rownames(spe) %in% markers)) {
    warning("The following markers were not specified, and manual thresholds for them will not be calculated: ", rownames(spe)[!rownames(spe) %in% markers])
    }

  manual_thresholds <- sapply(markers, function (marker) {
    positives <- SummarizedExperiment::assay(spe, "reference")[marker,]==1
    threshold <- min(SummarizedExperiment::assay(spe, "data")[marker,positives])
    return(threshold)
  })

  output <- data.frame(row.names = rownames(spe))
  output[match(markers, rownames(spe)),"manual_threshold"] <- manual_thresholds

  if ("manual_threshold" %in% colnames(SummarizedExperiment::rowData(spe))) {
    SummarizedExperiment::rowData(spe)$manual_threshold <- dplyr::coalesce(output$manual_threshold, SummarizedExperiment::rowData(spe)$manual_threshold)
  } else {
    SummarizedExperiment::rowData(spe)$manual_threshold <- output$manual_threshold
  }

}
