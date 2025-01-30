#' Title
#'
#' @param method Character; "HALO" if reading from a HALO output file; otherwise, "manual" if you are providing intensity matrix, reference matrix, spatial coords, etc. directly.
#' @param ... Arguments to be passed to build_ammit_spe_HALO, if method "HALO" is selected.
#' @param intensity_matrix Numeric; A matrix of markers (rows) x cells (columns) where values are the intensity of each marker.
#' @param reference_matrix Numeric; A matrix of markers (rows) x cells (columns) where values are the binary classification for whether a cell is positive for a marker.
#' @param spatial_coords Numeric; A matrix of cells (rows) x coordinates (columns), columns should be named Cell.X.Position and Cell.Y.Position.
#' @param col_data Data frame; colData to be attached to the SpatialExperiment object; i.e. metadata for each cell. This may include Analysis Region, for example.
#'
#' @return A SpatialExperiment object with assay 'data' containing the intensities per cell for each marker, spatial coordinates, metadata (in colData), and positive classification (if provided) in the assay 'reference'.
#' @export
#'
#' @examples
build_ammit_spe = function(method=c("HALO", "manual"),
                           ...,
                           intensity_matrix=NULL,
                           reference_matrix=NULL,
                           spatial_coords=NULL,
                           col_data=NULL) {

  if (method=="manual") {
    if (is.null(intensity_matrix)) {
      stop("Please specify an intensity_matrix (matrix of markers x cells, where data is intensity of a marker for a specific cell")
      }
    if (reference & is.null(reference_matrix)) {
      stop("Please specify a reference matrix (binary matrix of markers x cells, determining whether a cell is considered positive for a marker).")
      }
    if (is.null(spatial_coords)) {
      stop("Please specify spatial coordinates for your cells.")
      }
    }


  if (method=="HALO") { spe = build_ammit_spe_HALO(...) }

  if (method=="manual") { spe = SpatialExperiment::SpatialExperiment(assays=list("data"=intensity_matrix, "reference"=reference_matrix), colData = col_data, spatialCoords = spatial_coords) }

  return(spe)
}
