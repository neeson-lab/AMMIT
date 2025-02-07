#' Unmix marker intensity data from an AMMIT SpatialExperiment object into multiple skew-normal distributions
#'
#' @param object A SpatialExperiment object with assay 'data' present, or a list of SpatialExperiment objects.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to unmix, or "all" (default) to unmix all markers.
#' @param k Numeric; How many mixtures to derive from the intensity data. Will be recycled to match the number of markers. Markers can have different k values!
#' @param downsample Numeric; If provided, how many cells to downsample intensity data to before unmixing. Default value is 15000. If NULL or non-numeric, no downsampling will be done.
#' @param transformation Character; Passed to transform_intensities. Either "none" to use raw values, or one of "asinh" or "sqrt" to apply the respective transformation to the data. Defaults to "asinh".
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before unmixing. Recommended; defaults to TRUE.
#' @param ... Arguments to pass to unmix_intensites_split (subset) for a list input; or to smsn.mix, such as nu and iter.max.
#'
#' @return The input object(s), with model parameters for the unmixed distributions for each marker located in their respective rowData.
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=2)
#' SummarizedExperiment::rowData(ammit_spe_dummy)["M7", ]
unmix_intensities <- function(object,
                              markers,
                              k,
                              downsample = 15000,
                              transformation = "asinh",
                              censor_zeroes = TRUE,
                              ...) {

  if (!methods::is(object, "SpatialExperiment") && !methods::is(object, "list")) {
    stop("Your object is not a SpatialExperiment object, or list of SPE objects!")
  }
  if (!methods::is(markers, "character")) {
    stop("`markers` should be a character!")
  }
  if (!transformation %in% c("none", "asinh", "sqrt")) {
    stop("`transformation` should be one of: 'none', 'asinh', or 'sqrt'!")
  }


  if (methods::is(object, "list")) {
    object <- unmix_intensities_split(object = object,
                            markers = markers,
                            k = k,
                            downsample = downsample,
                            transformation = transformation,
                            censor_zeroes = censor_zeroes,
                           ...)
  }

  if (methods::is(object, "SpatialExperiment")) {
    object <- unmix_intensities_spe(object = object,
                          markers = markers,
                          k = k,
                          downsample = downsample,
                          transformation = transformation,
                          censor_zeroes = censor_zeroes,
                          ...)
  }

  return(object)

}

#' Unmix marker intensity data from a list of SpatialExperiment objects into multiple skew-normal distributions
#'
#' @param object A list of SpatialExperiment objects, each with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to unmix, or "all" (default) to unmix all markers.
#' @param k Numeric; How many mixtures to derive from the intensity data. Will be recycled to match the number of markers. Markers can have different k values!
#' @param downsample Numeric; If provided, how many cells to downsample intensity data to before unmixing. Default value is 15000. If NULL or non-numeric, no downsampling will be done.
#' @param transformation Character; Passed to transform_intensities. Either "none" to use raw values, or one of "asinh" or "sqrt" to apply the respective transformation to the data. Defaults to "asinh".
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before unmixing. Recommended; defaults to TRUE.
#' @param subset Character; Which of the split subsets to unmix. If left NULL (default), all subsets will be plotted.
#' @param ... Arguments to pass to smnsn.mix; such as nu and iter.max.
#'
#' @return The input objects, with model parameters for the unmixed distributions for each marker located in their respective rowData.
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' spe_list <- split_spe(ammit_spe_dummy, split.by="Analysis.Region")
#' spe_list <- unmix_intensities(spe_list, markers="M7", k=2, subset="Tumour")
#' SummarizedExperiment::rowData(spe_list[["Tumour"]])["M7", ]
unmix_intensities_split <- function(object,
                                    markers,
                                    k,
                                    downsample = 15000,
                                    transformation = "asinh",
                                    censor_zeroes = TRUE,
                                    subset=NULL,
                                    ...) {
  if (!methods::is(object, "list")) {
    stop("Your object is not a list of SpatialExperiment objects!")
  }

  if (is.null(subset)) {
    subset <- names(object)
  } else if (!all(subset %in% names(object))) {
    stop("Not all of the specified subsets are found in your list of objects.")
  }

  object[subset] <- lapply(subset, \(region) {
    unmix_intensities_spe(object = object[[region]],
                          markers = markers,
                          k = k,
                          downsample = downsample,
                          transformation = transformation,
                          censor_zeroes = censor_zeroes,
                          ...)
  })
  return(object)
}

#' Unmix marker intensity data from a SpatialExperiment object into multiple skew-normal distributions
#'
#' @param object A SpatialExperiment object, with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to unmix, or "all" (default) to unmix all markers.
#' @param k Numeric; How many mixtures to derive from the intensity data. Will be recycled to match the number of markers. Markers can have different k values!
#' @param downsample Numeric; If provided, how many cells to downsample intensity data to before unmixing. Default value is 15000. If NULL or non-numeric, no downsampling will be done.
#' @param transformation Character; Passed to transform_intensities. Either "none" to use raw values, or one of "asinh" or "sqrt" to apply the respective transformation to the data. Defaults to "asinh".
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before unmixing. Recommended; defaults to TRUE.
#' @param nu Passed to smsns.mix. Defaults to 0.
#' @param iter.max Passed to smsns.mix. Defaults to 2000.
#' @param ... Arguments to pass to smnsn.mix.
#'
#' @return The input object, with model parameters for the unmixed distributions for each marker located in its rowData.
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=3)
#' SummarizedExperiment::rowData(ammit_spe_dummy)["M7", ]
unmix_intensities_spe <- function(object,
                                  markers,
                                  k,
                                  downsample = 15000,
                                  transformation = "asinh",
                                  censor_zeroes = TRUE,
                                  nu=0,
                                  iter.max=2000,
                                  ...) {
  if (length(markers)==1 && markers=="all") {
    markers <- rownames(object)
  } else if (!all(markers %in% rownames(object))) {
    stop("Not all markers found in your SPE object.")
  }
  if (length(k)!=length(markers)) {
    warning("The lengths of `k` and `markers` do not match. k will be recycled!")
    k <- rep_len(k, length.out=length(markers))
  }

  k.frame <- cbind(data.frame(row.names=markers), k)

  models <- sapply(markers, function(marker) {

    matrix <- SummarizedExperiment::assay(object, "data")[marker,]
    if (censor_zeroes) {
      matrix <- matrix[matrix>0]
    }
    if (!is.null(downsample) && methods::is(downsample, "numeric")) {
      if (length(matrix)>downsample) {
        matrix <- sample(matrix, size=downsample)
      }
    } else { warning("Not dowsampling your intensities before unmixing. This could increase runtime significantly...") }
    matrix <- transform_intensities(matrix, method=transformation)

    k.actual <- k.frame[marker,]
    unmix <- mixsmsn::smsn.mix(matrix, nu = nu, g = k.actual, family = "Skew.normal", calc.im = F, iter.max = iter.max, ...)
    unmix$k <- k.actual
    unmix$unmix_transformation <- transformation
    unmix$unmixed <- TRUE
    return(unmix)
  }, simplify = T)

  add <- data.frame(t(models))
  missing.rows <- setdiff(rownames(object), rownames(add))
  add$unmix_transformation <- transformation
  add$unmixed <- TRUE
  add <- dplyr::bind_rows(add, data.frame(row.names=missing.rows, unmixed=rep(FALSE, length(missing.rows))))
  add <- add[rownames(object),]

  new_rowdata <- cbind(SummarizedExperiment::rowData(object), add)

  SummarizedExperiment::rowData(object) <- new_rowdata
  return(object)
}
