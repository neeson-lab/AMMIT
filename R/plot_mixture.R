#' Plot the unmixed distributions of marker intensities of an AMMIT SpatialExperiment Object
#'
#' @param object A SpatialExperiment object with assay 'data' present, or a list of SpatialExperiment objects.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to plot, or "all" (default) to plot all markers.
#' @param thresholds Character; "none" to plot no threshold values, "manual" to plot manual thresholds (see infer_manual_thresholds), "ammit" to plot ammit thresholds, or "both" for both. Defaults to "none".
#' @param combine_negatives Logical; Whether to combine the negative distributions, if there are multiple, to plot them as one curve (on which the ammit threshold would have been calculated)
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before plotting. Recommended; defaults to TRUE.
#' @param ... Arguments to be passed to plot_mixture_split in the case of a list of SpatialExperiment objects.
#'
#' @return A ggplot density plot for each of the markers selected
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' ammit_spe_dummy <- infer_manual_thresholds(ammit_spe_dummy, markers = "M7")
#' ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=2)
#' ammit_spe_dummy <- get_ammit_thresholds(ammit_spe_dummy,
#'                                         markers = "M7",
#'                                         strictness = 0.1,
#'                                         return = "object")
#' plot_mixture(ammit_spe_dummy,
#'              markers = "M7",
#'              thresholds = "both",
#'              combine_negative = FALSE,
#'              censor_zeroes = TRUE)
plot_mixture <- function(object,
                         markers="all",
                         thresholds = "none",
                         combine_negatives = TRUE,
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

  if (methods::is(object, "list")) {
    plot_mixture_split(object = object,
                           markers = markers,
                           thresholds = thresholds,
                           combine_negatives = combine_negatives,
                           censor_zeroes = censor_zeroes,
                           ...)
  }

  if (methods::is(object, "SpatialExperiment")) {
    plot_mixture_spe(object = object,
                         markers = markers,
                         thresholds = thresholds,
                     combine_negatives = combine_negatives,
                     censor_zeroes = censor_zeroes)
  }
}

#' Plot the unmixed distributions of marker intensities for a list of SpatialExperiment Objects
#'
#' @param object A list of SpatialExperiment objects, each with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to plot, or "all" (default) to plot all markers.
#' @param thresholds Character; "none" to plot no threshold values, "manual" to plot manual thresholds (see infer_manual_thresholds), "ammit" to plot ammit thresholds, or "both" for both. Defaults to "none".
#' @param combine_negatives Logical; Whether to combine the negative distributions, if there are multiple, to plot them as one curve (on which the ammit threshold would have been calculated)
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before plotting. Recommended; defaults to TRUE.
#' @param subset Character; Which of the split subsets to plot. If left NULL (default), all subsets will be plotted.
#'
#' @return A list of patchworked ggplot density plots with each of the markers selected, for each subset of the list.
plot_mixture_split <- function(object,
                               markers="all",
                               thresholds = "none",
                               combine_negatives = TRUE,
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
    plot_mixture_spe(object = object,
                     markers = markers,
                     thresholds = thresholds,
                     combine_negatives = combine_negatives,
                     censor_zeroes = censor_zeroes)
  })
}

#' Plot the unmixed distributions of marker intensities of an AMMIT SpatialExperiment Object
#'
#' @param object A SpatialExperiment object with assay 'data' present.
#' @param markers Character; one or multiple markers found in the SpatialExperiment object to plot, or "all" (default) to plot all markers.
#' @param thresholds Character; "none" to plot no threshold values, "manual" to plot manual thresholds (see infer_manual_thresholds), "ammit" to plot ammit thresholds, or "both" for both. Defaults to "none".
#' @param combine_negatives Logical; Whether to combine the negative distributions, if there are multiple, to plot them as one curve (on which the ammit threshold would have been calculated)
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before plotting. Recommended; defaults to TRUE.
#'
#' @return A ggplot density plot for each of the markers selected.
plot_mixture_spe <- function(object,
                             markers="all",
                             thresholds = "none",
                             combine_negatives = TRUE,
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
    plot_mixture_single(object = object,
                        marker = marker,
                        thresholds = thresholds,
                        combine_negatives = combine_negatives,
                        censor_zeroes = censor_zeroes)
  }) |>
    patchwork::wrap_plots(guides = "collect")

}


#' Plot the AMMIT-unmixed distributions of a marker's intensity values
#'
#' @param object A SpatialExperiment object with assay 'data' present.
#' @param marker Character; A marker, found in the SpatialExperiment object, to plot.
#' @param thresholds Character; "none" to plot no threshold values, "manual" to plot manual thresholds (see infer_manual_thresholds), "ammit" to plot ammit thresholds, or "both" for both. Defaults to "none".
#' @param combine_negatives Logical; Whether to combine the negative distributions, if there are multiple, to plot them as one curve (on which the ammit threshold would have been calculated)
#' @param censor_zeroes Logical; Whether to remove zeroes from intensity data before plotting. Recommended; defaults to TRUE.
#'
#' @return A single ggplot histogram plot of the intensities for the marker selected, with the unmixed distributions plotted on top.
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' ammit_spe_dummy <- infer_manual_thresholds(ammit_spe_dummy, markers = "M7")
#' ammit_spe_dummy <- unmix_intensities(ammit_spe_dummy, markers="M7", k=2)
#' ammit_spe_dummy <- get_ammit_thresholds(ammit_spe_dummy,
#'                                         markers = "M7",
#'                                         strictness = 0.1,
#'                                         return = "object")
#' plot_mixture_single(ammit_spe_dummy,
#'                     marker = "M7",
#'                     thresholds = "both",
#'                     combine_negative = FALSE,
#'                     censor_zeroes = TRUE)
plot_mixture_single <- function(object,
                                marker,
                                thresholds = "none",
                                combine_negatives = TRUE,
                                censor_zeroes = TRUE) {

  plot.data <- data.frame("intensity"=SummarizedExperiment::assay(object, "data")[marker,])
  if (censor_zeroes) {
    plot.data <- data.frame("intensity"=plot.data[plot.data$intensity>0,])
  }

  rowdata <- SummarizedExperiment::rowData(object)[marker, ]


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

  negative_distributions <- rowdata$idx[1:(k-1)]
  negative_distributions <- 1:(k-1)

  transformation <- unique(rowdata$unmix_transformation)
  plot.data$intensity <- transform_intensities(plot.data$intensity, method = transformation)

  ## base plot (histogram)

  p <- ggplot2::ggplot(data = plot.data, mapping = ggplot2::aes(x = intensity)) +
    ggplot2::geom_histogram(bins=100, color="dimgrey", fill="grey", ggplot2::aes(y = ggplot2::after_stat(density) )) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste(marker, "Unmixing results")) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 10))

  ## add function lines for unmixed density curves

  if (combine_negatives) {
    neg_line <- ggplot2::geom_function(fun = \(x) {
      density_combine <- \(x) { do.call(sum,
                                        lapply(negative_distributions,
                                               \(id) {
                                                 rowdata$pii[id] * sn::dsn(x,
                                                                           xi = rowdata$mu[id],
                                                                           omega = rowdata$scale[id],
                                                                           alpha = rowdata$shape[id])
                                               }))
      }
      Vectorize(density_combine)(x)
    },
    color="red")
    pos_line <- ggplot2::geom_function(fun = \(x) rowdata$pii[k] * sn::dsn(x,
                                                                           xi = rowdata$mu[k],
                                                                           omega = rowdata$scale[k],
                                                                           alpha = rowdata$shape[k]),
                                       color="green")
    p <- p + neg_line + pos_line
  } else {
    distr.colours <- c(c("red", "blue", "purple", "orange", "violet", "brown")[1:(k-1)], "green")
    names(distr.colours) <- rowdata$idx
    lines <- lapply(1:k, \(id) {
      ggplot2::geom_function(fun = \(x) rowdata$pii[id] * sn::dsn(x,
                                                                  xi = rowdata$mu[id],
                                                                  omega = rowdata$scale[id],
                                                                  alpha = rowdata$shape[id]),
                             color = distr.colours[[id]])
    })
    p <- p + lines
  }

  ## add thresholds

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
    ammit_threshold <- data.frame(athresh = athresh, color = "AMMIT")
    p <- p +
      ggplot2::geom_vline(data = manual_threshold, mapping = ggplot2::aes(color = "Manual", xintercept = thresh), lty = 2) +
      ggplot2::geom_vline(data = ammit_threshold, mapping = ggplot2::aes(color = "AMMIT", xintercept = athresh), lty = 1) +
      ggplot2::scale_color_manual(values = c("Manual" = "orange", "AMMIT" = "darkgreen"), name = "Threshold")
  }

  return(p)

}
