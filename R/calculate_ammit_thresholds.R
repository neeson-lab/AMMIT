#' Derive the AMMIT threshold from unmixed marker intensity data and a given strictness parameter.
#'
#' @param object A SpatialExperiment object with assay 'data' present, or a list of SpatialExperiment objects.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object for which to derive thresholds, or "all" (default) to do so for all markers.
#' @param strictness Numeric; one or multiple strictness values matching to each of the marker(s) specified.
#' @param return Character; either 'thresholds' to return a named vector of the derived thresholds for each marker, or 'object' (default) to return the SpatialExperiment object with this information embedded in the rowData.
#' @param ... Arguments to pass to infer_strictness_split (subset) for a list input
#'
#' @return A single, or list of: Numeric vector(s) of derived thresholds (if 'thresholds' specified for return), or SpatialExperiment object(s) with this information embedded in rowData (if 'object' specified).
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' ammit_spe_dummy <- infer_manual_thresholds(ammit_spe_dummy, markers = "M7")
#' ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=3)
#' strictness <- infer_strictness(ammit_spe_dummy, markers = "M7", return = "strictness")
#' calculate_ammit_thresholds(ammit_spe_dummy, markers = "M7", strictness = strictness, return = "thresholds")
calculate_ammit_thresholds <- function (object,
                                  markers,
                                  strictness,
                                  return = "object",
                                  ...) {
  if (!methods::is(object, "SpatialExperiment") && !methods::is(object, "list")) {
    stop("Your object is not a SpatialExperiment object, or list of SPE objects!")
  }
  if (!methods::is(markers, "character")) {
    stop("`markers` should be a character!")
  }

  if (methods::is(object, "list")) {
    object <- calculate_ammit_thresholds_split(object = object,
                                     markers = markers,
                                     strictness = strictness,
                                     return = return,
                                     ...)
  }

  if (methods::is(object, "SpatialExperiment")) {
    object <- calculate_ammit_thresholds_spe(object = object,
                                   markers = markers,
                                   strictness = strictness,
                                   return = return,
                                   ...)
  }
  return(object)
}

#' Derive the AMMIT threshold from unmixed marker intensity data and a given strictness parameter.
#'
#' @param object A list of SpatialExperiment objects, each with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object for which to derive thresholds, or "all" (default) to do so for all markers.
#' @param strictness Numeric; one or multiple strictness values matching to each of the marker(s) specified.
#' @param return Character; either 'thresholds' to return a list of named vectors of the derived thresholds for each marker, or 'object' (default) to return the list of SpatialExperiment objects with this information embedded in the rowData.
#' @param subset Character; Which of the split subsets for which to derive thresholds. If left NULL (default), all subsets will have strictness inferred.
#'
#' @return A list of: Numeric vectors of derived thresholds (if 'thresholds' specified for return), or SpatialExperiment objects with this information embedded in rowData (if 'object' specified).
#'
#' @examples
#' data("ammit_spe_dummy")
#' spe_list <- split_spe(ammit_spe_dummy)
#' spe_list <- infer_manual_thresholds(spe_list, markers = "M7")
#' spe_list <- unmix_intensities(spe_list, markers="M7", k=3)
#' strictness <- infer_strictness(spe_list, markers = "M7", return = "strictness")
#' calculate_ammit_thresholds(spe_list, markers = "M7", strictness = strictness, return = "thresholds")
calculate_ammit_thresholds_split <- function (object,
                                        markers,
                                        strictness,
                                        return = "object",
                                        subset=NULL) {
  if (!methods::is(object, "list")) {
    stop("Your object is not a list of SpatialExperiment objects!")
  }

  if (is.null(subset)) {
    subset <- names(object)
  } else if (!all(subset %in% names(object))) {
    stop("Not all of the specified subsets are found in your list of objects.")
  }

  if (return=="thresholds") {
    result <- lapply(subset, function(region) {
      if (methods::is(strictness, "list")) {
        strictness <- strictness[[region]]
      }
      calculate_ammit_thresholds_spe(object = object[[region]],
                           markers = markers,
                           strictness = strictness,
                           return = "thresholds")
    })
    names(result) <- subset
    return(result)
  } else if (return=="object") {
    object[subset] <- lapply(subset, function(region) {
      if (methods::is(strictness, "list")) {
        strictness <- strictness[[region]]
      }
      calculate_ammit_thresholds_spe(object = object[[region]],
                           markers = markers,
                           strictness = strictness,
                           return = "object")
    })
    names(object) <- subset
    return(object)
  } else { stop("Unrecognised option supplied for `return`") }

}

#' Derive the AMMIT threshold from unmixed marker intensity data and a given strictness parameter.
#'
#' @param object A SpatialExperiment object, with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object for which to derive thresholds, or "all" (default) to do so for all markers.
#' @param strictness Numeric; one or multiple strictness values matching to each of the marker(s) specified.
#' @param return Character; either 'thresholds' to return a named vector of the derived thresholds for each marker, or 'object' (default) to return the SpatialExperiment object with this information embedded in the rowData.
#'
#' @return A numeric vector of derived thresholds (if 'thresholds' specified for return), or a SpatialExperiment object with this information embedded in rowData (if 'object' specified).
#'
#' @examples
#' data("ammit_spe_dummy")
#' ammit_spe_dummy <- infer_manual_thresholds(ammit_spe_dummy, markers = "M7")
#' ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=3)
#' strictness <- infer_strictness(ammit_spe_dummy, markers = "M7", return = "strictness")
#' calculate_ammit_thresholds(ammit_spe_dummy, markers = "M7", strictness = strictness, return = "thresholds")
calculate_ammit_thresholds_spe <- function (object,
                                      markers,
                                      strictness,
                                      return = "object") {
  if (length(markers)==1 && markers=="all") {
    markers <- rownames(object)
  } else if (!all(markers %in% rownames(object))) {
    stop("Not all markers found in your SPE object.")
  }

  if (!"unmixed" %in% colnames(SummarizedExperiment::rowData(object))) {
    stop("Unmixing has not been completed for this object! Please run unmix_intensities before deriving the AMMIT threshold.")
  } else if (!all(SummarizedExperiment::rowData(object)[markers,"unmixed"])) {
    stop("Unmixing has not completed for all the markers you specified! Please run unmix_intensities for the relevant markers.")
  }

  if (length(markers)!=length(strictness)) {
    stop("Length of markers and strictness values does not match!")
  }


  ammit_thresholds <- sapply(markers, \(marker) {
    rowdata <- SummarizedExperiment::rowData(object)[marker, ]

    strict <- strictness[match(marker, markers)]

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

    p_thresh <- 1-10^-strict

    thresh_get <- function(x) {
      distribution_sums <- do.call(sum, lapply(negative_distributions, FUN=\(id){
        rowdata$pii[id] * sn::psn(x=x,
                                xi=rowdata$mu[id],
                                omega=rowdata$scale[id],
                                alpha=rowdata$shape[id])
      }))
      difference <- (distribution_sums/sum(rowdata$pii[negative_distributions])) - p_thresh
      return(difference)
    }

    thresh <- uniroot(thresh_get, c(0,20))$root

    if (unique(rowdata$unmix_transformation) == "asinh") {
      thresh <- sinh(thresh)
    } else if (unique(rowdata$unmix_transformation) == "sqrt") {
      thresh <- thresh^2
    }

    return(thresh)
  })

  if (return=="thresholds") {
    return(ammit_thresholds)
  } else if (return == "object") {
    add <- data.frame(row.names=rownames(object))
    add$ammit_threshold <- NA
    add$ammit_threshold[match(names(ammit_thresholds), rownames(object))] <- ammit_thresholds
    if ("ammit_threshold" %in% colnames(SummarizedExperiment::rowData(object))) {
      merge <- cbind(add, "old"=SummarizedExperiment::rowData(object)$ammit_threshold)
      add <- merge |>
        dplyr::mutate(ammit_threshold=dplyr::coalesce(ammit_threshold, old)) |>
        dplyr::select(ammit_threshold)
    }
    SummarizedExperiment::rowData(object)[,"ammit_threshold"] <- add
    return(object)
  } else { stop("Unrecognised option supplied for `return`") }

}
