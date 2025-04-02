#' Infer the AMMIT strictness parameter required to derive manual threshold values from unmixed marker intensity data
#'
#' @param object A SpatialExperiment object with assay 'data' present, or a list of SpatialExperiment objects.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object for which to infer strictness, or "all" (default) to do so for all markers.
#' @param return Character; either 'strictness' (default) to return a named vector of the inferred strictness for each marker, or 'object' to return the SpatialExperiment object with this information embedded in the rowData.
#' @param ... Arguments to pass to infer_strictness_split (subset) for a list input
#'
#' @return A single, or list of: Numeric vector(s) of inferred strictness values (if 'strictness' specified for return), or SpatialExperiment object(s) with this information embedded in rowData (if 'object' specified).
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' ammit_spe_dummy <- infer_manual_thresholds(ammit_spe_dummy, markers = "M7")
#' ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=3)
#' infer_strictness(ammit_spe_dummy, markers = "M7", return = "strictness")
infer_strictness <- function(object,
                             markers="all",
                             return="strictness",
                             ...) {
  if (!methods::is(object, "SpatialExperiment") && !methods::is(object, "list")) {
    stop("Your object is not a SpatialExperiment object, or list of SPE objects!")
  }
  if (!methods::is(markers, "character")) {
    stop("`markers` should be a character!")
  }

  if (methods::is(object, "list")) {
    object <- infer_strictness_split(object = object,
                                     markers = markers,
                                     return = return,
                                     ...)
  }

  if (methods::is(object, "SpatialExperiment")) {
    object <- infer_strictness_spe(object = object,
                                   markers = markers,
                                   return = return,
                                   ...)
  }
  return(object)
}

#' Infer the AMMIT strictness parameter required to derive manual threshold values from unmixed marker intensity data
#'
#' @param object A list of SpatialExperiment objects, each with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object for which to infer strictness, or "all" (default) to do so for all markers.
#' @param return Character; either 'strictness' (default) to return a list of named vectors of the inferred strictness for each marker, or 'object' to return the list of SpatialExperiment objects with this information embedded in the rowData.
#' @param subset Character; Which of the split subsets for which to infer strictness. If left NULL (default), all subsets will have strictness inferred.
#'
#' @return A list of: Numeric vectors of inferred strictness values (if 'strictness' specified for return), or SpatialExperiment objects with this information embedded in rowData (if 'object' specified).
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' spe_list <- split_spe(ammit_spe_dummy)
#' spe_list <- infer_manual_thresholds(spe_list, markers = "M7")
#' spe_list <- unmix_intensities(spe_list, markers="M7", k=3)
#' infer_strictness(spe_list, markers = "M7", return = "strictness")
infer_strictness_split <- function(object,
                                   markers="all",
                                   return="strictness",
                                   subset=NULL) {
  if (!methods::is(object, "list")) {
    stop("Your object is not a list of SpatialExperiment objects!")
  }

  if (is.null(subset)) {
    subset <- names(object)
  } else if (!all(subset %in% names(object))) {
    stop("Not all of the specified subsets are found in your list of objects.")
  }

  if (return=="strictness") {
    result <- lapply(subset, function(region) {
      infer_strictness_spe(object = object[[region]],
                           markers = markers,
                           return = "strictness")
    })
    names(result) <- subset
    return(result)
  } else if (return=="object") {
    object[subset] <- lapply(subset, function(region) {
      infer_strictness_spe(object = object[[region]],
                            markers = markers,
                            return = "object")
    })
    names(object) <- subset
    return(object)
  } else { stop("Unrecognised option supplied for `return`") }

}

#' Infer the AMMIT strictness parameter required to derive manual threshold values from unmixed marker intensity data
#'
#' @param object A SpatialExperiment object, with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object for which to infer strictness, or "all" (default) to do so for all markers.
#' @param return Character; either 'strictness' (default) to return a named vector of the inferred strictness for each marker, or 'object' to return the SpatialExperiment object with this information embedded in the rowData.
#'
#' @return A numeric vector of inferred strictness values (if 'strictness' specified for return), or a SpatialExperiment object with this information embedded in rowData (if 'object' specified).
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' ammit_spe_dummy <- infer_manual_thresholds(ammit_spe_dummy, markers = "M7")
#' ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=3)
#' infer_strictness(ammit_spe_dummy, markers = "M7", return = "strictness")
infer_strictness_spe <- function(object,
                                 markers="all",
                                 return="strictness") {
  if (length(markers)==1 && markers=="all") {
    markers <- rownames(object)
  } else if (!all(markers %in% rownames(object))) {
    stop("Not all markers found in your SPE object.")
  }

  if (!"manual_threshold" %in% colnames(SummarizedExperiment::rowData(object))) {
    stop("Manual thresholds not found! Please enter them into rowData as 'manual_threshold', or run infer_manual_thresholds on this object!")
  }
  if (!"unmixed" %in% colnames(SummarizedExperiment::rowData(object))) {
    stop("Unmixing has not been completed for this object! Please run unmix_intensities before inferring strictness.")
  } else if (!all(SummarizedExperiment::rowData(object)[markers,"unmixed"])) {
    stop("Unmixing has not completed for all the markers you specified! Please run unmix_intensities for the relevant markers.")
  }

  strictness <- sapply(markers, function(marker){

    rowdata <- SummarizedExperiment::rowData(object)[marker,]
    manual_threshold <- rowdata$manual_threshold
    manual_threshold <- transform_intensities(manual_threshold, method=unique(rowdata$unmix_transformation))
    rowdata <- data.frame(unlist(rowdata))
    colnames(rowdata) <- gsub(colnames(rowdata), pattern=paste0("\\.", marker, "$"), replacement="")
    k <- unique(rowdata$k)
    rowdata$idx <- 1:k
    rowdata <- rowdata |>
      dplyr::mutate(delta = shape / (sqrt(1 + (shape^2)))) |>
      dplyr::mutate(coeff = sqrt(2 / pi)) |>
      dplyr::mutate(scale = sqrt(sigma2)) |>
      dplyr::mutate(mean=mu+scale*delta*coeff) |>
      dplyr::arrange(mean)

    negative_distributions <- 1:(k-1)

    sums <- do.call(sum, lapply(negative_distributions, FUN=\(id){
      rowdata$pii[id]*sn::psn(x=manual_threshold,
                              xi=rowdata$mu[id],
                              omega=rowdata$scale[id],
                              alpha=rowdata$shape[id])
    }))
    strict <- -log10(1-sums/sum(rowdata$pii[negative_distributions]))
    return(strict)
  })

  if (return=="strictness") {
    return(strictness)
  } else if (return == "object") {
    add <- data.frame(row.names=rownames(object))
    add$inferred_strictness <- NA
    add$inferred_strictness[match(names(strictness), rownames(object))] <- strictness
    new_rowdata <- cbind(SummarizedExperiment::rowData(object), add)
    SummarizedExperiment::rowData(object) <- new_rowdata
    return(object)
  } else { stop("Unrecognised option supplied for `return`") }
}
