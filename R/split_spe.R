#' Split a SpatialExperiment object based on a variable
#'
#' @param spe A SpatialExperiment object to be split.
#' @param split.by Character; A column name present in the colData of spe, to split upon. Defaults to "Analysis.Region".
#'
#' @importFrom SummarizedExperiment colData
#'
#' @return A named list of SpatialExperiment objects, where each object is the subset of the original for one value of the splitting variable.
#' @export
#'
#' @examples
#' data("ammit_spe_dummy")
#' split_spe(ammit_spe_dummy, split.by="Analysis.Region")
split_spe <- function(spe,
                            split.by="Analysis.Region") {
  split.by <- as.character(split.by)
  regions <- unique(colData(spe)[, split.by])
  split <- lapply(regions, \(region) {
    subset <- colData(spe)[, split.by] == region
    spe[,subset]
  } )
  names(split) <- regions
  return(split)
}
