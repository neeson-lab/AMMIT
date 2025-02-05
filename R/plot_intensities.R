#' Plot the marker intensities of an AMMIT SpatialExperiment Object
#'
#' @param object A SpatialExperiment object with assay 'data' present, or a list of SpatialExperiment objects.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to plot, or "all" (default) to plot all markers.
#' @param thresholds Character; "none" to plot no threshold values, "manual" to plot manual thresholds (see infer_manual_thresholds), "ammit" to plot ammit thresholds, or "both" for both. Defaults to "none".
#' @param transformation Character; Passed to transform_intensities. Either "none" to use raw values, or one of "asinh" or "sqrt" to apply the respective transformation to the data. Defaults to "asinh".
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before plotting. Recommended; defaults to TRUE.
#' @param ... Arguments to be passed to plot_intensities_split in the case of a list of SpatialExperiment objects.
#'
#' @return A ggplot density plot for each of the markers selected
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' plot_intensities(ammit_spe_dummy,
#'                  markers = "M7",
#'                  thresholds = "none",
#'                  transformation = "asinh",
#'                  censor_zeroes = TRUE)
plot_intensities <- function(object,
                             markers="all",
                             thresholds = "none",
                             transformation = "asinh",
                             censor_zeroes = TRUE,
                             ...) {
  if (!methods::is(object, "SpatialExperiment") && !methods::is(object, "list")) {
    stop("Your object is not a SpatialExperiment object, or list of SPE objects!")
  }
  if (!methods::is(markers, "character")) {
    stop("`markers` should be a character!")
  }
  if (!thresholds %in% c("none", "manual", "ammit", "both")) {
    stop("`thresholds` should be one of: 'none', 'manual', 'ammit', or 'both'!")
  }
  if (!transformation %in% c("none", "asinh", "sqrt")) {
    stop("`transformation` should be one of: 'none', 'asinh', or 'sqrt'!")
  }

  if (methods::is(object, "list")) {
    plot_intensities_split(object = object,
                           markers = markers,
                           thresholds = thresholds,
                           transformation = transformation,
                           censor_zeroes = censor_zeroes,
                           ...)
  }

  if (methods::is(object, "SpatialExperiment")) {
    plot_intensities_spe(object = object,
                         markers = markers,
                         thresholds = thresholds,
                         transformation = transformation,
                         censor_zeroes = censor_zeroes)
  }
}

#' Plot the marker intensities for a list of SpatialExperiment Objects
#'
#' @param object A list of SpatialExperiment objects, each with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to plot, or "all" (default) to plot all markers.
#' @param thresholds Character; "none" to plot no threshold values, "manual" to plot manual thresholds (see infer_manual_thresholds), "ammit" to plot ammit thresholds, or "both" for both. Defaults to "none".
#' @param transformation Character; Passed to transform_intensities. Either "none" to use raw values, or one of "asinh" or "sqrt" to apply the respective transformation to the data. Defaults to "asinh".
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before plotting. Recommended; defaults to TRUE.
#' @param subset Character; Which of the split subsets to plot. If left NULL (default), all subsets will be plotted.
#'
#' @return A list of patchworked ggplot density plots with each of the markers selected, for each subset of the list.
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' spe_list <- split_spe(ammit_spe_dummy, split.by="Analysis.Region")
#' plot_intensities_split(spe_list,
#'                  markers = "M7",
#'                  thresholds = "none",
#'                  transformation = "asinh",
#'                  censor_zeroes = TRUE,
#'                  subset = "Tumour")
plot_intensities_split <- function(object,
                                   markers="all",
                                   thresholds = "none",
                                   transformation="asinh",
                                   censor_zeroes = TRUE,
                                   subset=NULL) {
  if (!methods::is(object, "list")) {
    stop("Your object is not a list of SpatialExperiment objects!")
  }

  if (is.null(subset)) {
    subset <- names(object)
  } else if (!all(subset %in% names(object))) {
    stop("Not all of the specified subsets are found in your list of objects.")
  }

  lapply(subset, \(region) {
    object <- object[[region]]
    plot_intensities_spe(object = object,
                         markers = markers,
                         thresholds = thresholds,
                         transformation = transformation,
                         censor_zeroes = censor_zeroes)
  })
}

#' Plot the marker intensities of an AMMIT SpatialExperiment Object
#'
#' @param object A SpatialExperiment object with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to plot, or "all" (default) to plot all markers.
#' @param thresholds Character; "none" to plot no threshold values, "manual" to plot manual thresholds (see infer_manual_thresholds), "ammit" to plot ammit thresholds, or "both" for both. Defaults to "none".
#' @param transformation Character; Passed to transform_intensities. Either "none" to use raw values, or one of "asinh" or "sqrt" to apply the respective transformation to the data. Defaults to "asinh".
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before plotting. Recommended; defaults to TRUE.
#'
#' @return A ggplot density plot for each of the markers selected.
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' plot_intensities_spe(ammit_spe_dummy,
#'                  markers = "M7",
#'                  thresholds = "none",
#'                  transformation = "asinh",
#'                  censor_zeroes = TRUE)
plot_intensities_spe <- function(object,
                                 markers="all",
                                 thresholds = "none",
                                 transformation="asinh",
                                 censor_zeroes = TRUE) {
  if (!methods::is(object, "SpatialExperiment")) {
    stop("Your object is not a SpatialExperiment object!")
  }

  if (markers=="all") {
    markers <- rownames(object)
  } else if (!all(markers %in% rownames(object))) {
    stop("Not all markers found in your SPE object.")
  }

  lapply(markers, \(marker){
    plot_intensities_single(object = object,
                            marker = marker,
                            thresholds = thresholds,
                            transformation = transformation,
                            censor_zeroes = censor_zeroes)
  }) |>
    patchwork::wrap_plots(guides = "collect")

}

#' Plot a single marker intensity from an AMMIT SpatialExperiment Object
#'
#' @param object A SpatialExperiment object with assay 'data' present.
#' @param marker Character; A single marker present in the rownames of the SpatialExperiment Object, to be plotted.
#' @param thresholds Character; "none" to plot no threshold values, "manual" to plot manual thresholds (see infer_manual_thresholds), "ammit" to plot ammit thresholds, or "both" for both. Defaults to "none".
#' @param transformation Character; Passed to transform_intensities. Either "none" to use raw values, or one of "asinh" or "sqrt" to apply the respective transformation to the data. Defaults to "asinh".
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before plotting. Recommended; defaults to TRUE.
#'
#' @return A single ggplot density plot of the intensities for the marker selected.
#' @export
#'
#' @examples
#' #' data("ammit_spe_dummy")
#' plot_intensities_spe(ammit_spe_dummy,
#'                  marker = "M7",
#'                  thresholds = "none",
#'                  transformation = "asinh",
#'                  censor_zeroes = TRUE)
plot_intensities_single <- function(object,
                                    marker,
                                    thresholds = "none",
                                    transformation = "asinh",
                                    censor_zeroes = TRUE) {

  plot.data <- data.frame("Intensity"=SummarizedExperiment::assay(object, "data")[marker,])

  if (censor_zeroes) {
    plot.data <- data.frame("Intensity"=plot.data[plot.data$Intensity>0,])
  }

  plot.data$Intensity <- transform_intensities(plot.data$Intensity, method = transformation)
  trans.label <- c("asinh-transformed", "sqrt-transformed", "raw")[match(transformation, c("asinh", "sqrt", "none"))]

  p <- ggplot2::ggplot(data = plot.data, mapping = ggplot2::aes(x = Intensity)) +
    ggplot2::geom_density() +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste(marker, trans.label, "intensity")) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10))

  if (thresholds=="manual") {
    thresh = SummarizedExperiment::rowData(object)[marker, "manual_threshold"] |>
      transform_intensities(method = transformation)
    manual_threshold <- data.frame(thresh = thresh, color = "Manual")
    p <- p +
      ggplot2::geom_vline(data = manual_threshold, mapping = ggplot2::aes(color = "Manual", xintercept = thresh), lty = 2) +
      ggplot2::scale_color_manual(values = c("Manual" = "orange"), name = "Threshold")
  }
  if (thresholds=="ammit") {
    athresh = SummarizedExperiment::rowData(object)[marker, "ammit_threshold"] |>
      transform_intensities(method = transformation)
    ammit_threshold <- data.frame(thresh = athresh, color = "AMMIT")
    p <- p +
      ggplot2::geom_vline(data = ammit_threshold, mapping = ggplot2::aes(color = "AMMIT", xintercept = thresh), lty = 1) +
      ggplot2::scale_color_manual(values = c("AMMIT" = "darkgreen"), name = "Threshold")
  }
  if (thresholds=="both") {
    thresh = SummarizedExperiment::rowData(object)[marker, "manual_threshold"] |>
      transform_intensities(method = transformation)
    athresh = SummarizedExperiment::rowData(object)[marker, "ammit_threshold"] |>
      transform_intensities(method = transformation)
    manual_threshold <- data.frame(thresh = thresh, color = "Manual")
    ammit_threshold <- data.frame(thresh = thresh, color = "AMMIT")
    p <- p +
      ggplot2::geom_vline(data = manual_threshold, mapping = ggplot2::aes(color = "Manual", xintercept = thresh), lty = 2) +
      ggplot2::geom_vline(data = ammit_threshold, mapping = ggplot2::aes(color = "AMMIT", xintercept = thresh), lty = 1) +
      ggplot2::scale_color_manual(values = c("Manual" = "orange", "AMMIT" = "darkgreen"), name = "Threshold")
  }

  return(p)

}
