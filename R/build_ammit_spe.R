#' Title
#'
#' @param method
#' @param markers
#' @param reference
#' @param filter.dapi
#' @param HALO_filepath
#' @param locations
#' @param intensity_matrix
#' @param reference_matrix
#' @param spatial_coords
#' @param manual_thresholds
#'
#' @return
#' @export
#'
#' @examples
build_ammit_spe = function(method=c("HALO", "manual"), markers, reference, filter.dapi,
                           HALO_filepath, locations=NULL,
                           intensity_matrix=NULL, reference_matrix=NULL, spatial_coords=NULL,
                           manual_thresholds=NULL) {

  if (method=="manual") {
    if (is.null(intensity_matrix)) { stop("Please specify an intensity_matrix (matrix of markers x cells, where data is intensity of a marker for a specific cell") }
    if (reference & is.null(reference_matrix)) { stop("Please specify a reference matrix (binary matrix of markers x cells, determining whether a cell is considered positive for a marker).") }
    if (is.null(spatial_coords)) { stop("Please specify spatial coordinates for your cells.") }
    }

  if (method=="HALO") {
    if (!file.exists(HALO_filepath)) { stop("Cannot find file specified in HALO_filepath.") }
    if (!is.null(locations) & length(locations)!=length(markers)) { stop("`locations` and `markers` should be character vectors of the same length. ") }
  }

  if (!is(reference, "logical")) { stop("reference should be set to TRUE or FALSE depending on whether your data includes pre-calculated positivity classification for markers") }

  if (method=="HALO") { spe = build_ammit_spe_HALO(markers=markers, reference=reference, HALO_filepath=HALO_filepath, locations=locations, filter_dapi=filter.dapi, manual_thresholds=NULL) }

  if (method=="manual") { spe = SpatialExperiment(assays=list("data"=intensity_matrix, "reference"=reference_matrix), colData = col_data, spatialCoords = spatial_coords) }

  return(spe)
}
