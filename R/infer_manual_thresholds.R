#' Infer manual thresholds for marker intensities from pre-existing classifications
#'
#' @param object A SpatialExperiment object with assays 'data' and 'reference' present, or a list of SpatialExperiment objects.
#' @param markers A character vector of markers matching, at least partially, to the rownames of your SpatialExperiment object(s), or "all" (default) for all markers.
#' @param ... Arguments to pass to infer_manual_thresholds_split (if a list is input)
#'
#' @return The input SpatialExperiment object, with manual thresholds for all markers specified inserted in the rowData under manual_thresholds.
#' @export
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment rowData<-
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assays
#'
#' @examples
#' data("ammit_spe_dummy")
#' markers <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7")
#' ammit_spe_dummy <- infer_manual_thresholds(ammit_spe_dummy, markers=markers)
#' SummarizedExperiment::rowData(ammit_spe_dummy)
infer_manual_thresholds <- function(object, markers, ...) {

  if (!methods::is(object, "SpatialExperiment") && !methods::is(object, "list")) {
    stop("Your object is not a SpatialExperiment object, or list of SPE objects!")
  }

  if (methods::is(object, "list")) {

    object <- infer_manual_thresholds_split(object=object, markers=markers, subset=NULL)

    } else {
    if (!"reference" %in% names(SummarizedExperiment::assays(object))) {
      stop("No reference data found in your SpatialExperiment object. Have you added reference data?")
      }
    if (!"data" %in% names(SummarizedExperiment::assays(object))) {
      stop("No intensity data found in your SpatialExperiment object. Have you added intensity data?")
      }
    if (length(markers)==1 && markers=="all") {
        markers <- rownames(object)
      } else if (!all(markers %in% rownames(object))) {
        stop("Not all markers found in your SPE object.")
      }
      if (!all(rownames(object) %in% markers)) {
      warning("The following markers were not specified, and manual thresholds for them will not be calculated: ", rownames(object)[!rownames(object) %in% markers])
      }

    manual_thresholds <- sapply(markers, function (marker) {
      positives <- SummarizedExperiment::assay(object, "reference")[marker,]==1
      threshold <- min(SummarizedExperiment::assay(object, "data")[marker,positives])
      return(threshold)
    })

    output <- data.frame(row.names = rownames(object))
    output[match(markers, rownames(object)),"manual_threshold"] <- manual_thresholds

    if ("manual_threshold" %in% colnames(SummarizedExperiment::rowData(object))) {
      SummarizedExperiment::rowData(object)$manual_threshold <- dplyr::coalesce(output$manual_threshold, SummarizedExperiment::rowData(object)$manual_threshold)
    } else {
      SummarizedExperiment::rowData(object)$manual_threshold <- output$manual_threshold
    }
  }

  return(object)
}

#' Infer manual thresholds for marker intensities from pre-existing classifications on a list of SpatialExperiment objects
#'
#' @param object A list of SpatialExperiment objects, each with assays 'data' and 'reference' present.
#' @param markers A character vector of markers matching, at least partially, to the rownames of your SpatialExperiment objects, or "all" (default) for all markers.
#' @param subset A subset of the names of the list of objects for which you wish to infer manual thresholds.
#'
#' @return The input list of SpatialExperiment objects, each with manual thresholds for all markers specified inserted in the rowData under manual_thresholds.
#' @export
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' data("ammit_spe_dummy")
#' spe_list <- split_spe(ammit_spe_dummy, split.by="Analysis.Region")
#' markers <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7")
#' spe_list <- infer_manual_thresholds(spe_list, markers=markers)
#' lapply(spe_list, SummarizedExperiment::rowData)
infer_manual_thresholds_split <- function(object, markers, subset=NULL) {

  if (is.null(subset)) {
    subset <- names(object)
  } else if (!all(subset %in% names(object))) {
    stop("Not all of the specified subsets are found in your list of objects.")
  }

  for (region in subset) {
    if (!"reference" %in% names(assays(object[[region]]))) {
      stop("No reference data found in your SpatialExperiment object, in region: ", region, ". Have you added reference data?")
    }
    if (!"data" %in% names(SummarizedExperiment::assays(object[[region]]))) {
      stop("No intensity data found in your SpatialExperiment object, in region: ", region, ". Have you added intensity data?")
    }
    if (length(markers)==1 && markers=="all") {
      markers <- rownames(object[[region]])
    } else if (!all(markers %in% rownames(object[[region]]))) {
      stop("None of the specified markers were found in your SpatialExperiment object, in region: ", region, "!")
    }
    if (!all(rownames(object[[region]]) %in% markers)) {
      warning("The following markers were not specified, in region: ", region, ", and manual thresholds for them will not be calculated: ",
              paste(rownames(object[[region]])[!rownames(object[[region]]) %in% markers], collapse=" "))
    }
  }

  temp <- lapply(subset, \(region){

    object <- object[[region]]

    manual_thresholds <- sapply(markers, function (marker) {
      positives <- assay(object, "reference")[marker,]==1
      threshold <- min(assay(object, "data")[marker,positives])
      return(threshold)
    })

    output <- data.frame(row.names = rownames(object))
    output[match(markers, rownames(object)),"manual_threshold"] <- manual_thresholds

    if ("manual_threshold" %in% colnames(rowData(object))) {
      rowData(object)$manual_threshold <- dplyr::coalesce(output$manual_threshold, rowData(object)$manual_threshold)
    } else {
      rowData(object)$manual_threshold <- output$manual_threshold
    }

    return(object)

  })

  names(temp) <- subset

  for (region in subset) { object[[region]] <-  temp[[region]] }

  rm(temp)
  return(object)

}
