#' Title
#'
#' @param markers
#' @param reference
#' @param HALO_filepath
#' @param locations
#' @param filter_dapi
#' @param manual_thresholds
#'
#' @return
#' @export
#'
#' @examples
build_ammit_spe_HALO = function(markers, reference, HALO_filepath, locations, filter_dapi,
                                manual_thresholds=NULL) {


  if (!is.null(HALO_filepath)) {
    message("Reading from file...")
    if (grepl(pattern="\\.csv$", HALO_filepath, ignore.case=T)) {
      data <- read.csv(HALO_filepath)
    } else if (grepl(pattern="\\.tsv$", HALO_filepath, ignore.case=T)) {
      data <- read.tsv(HALO_filepath)
    } else { stop("If specifying `HALO_filepath`, the file in question should have a .csv or .tsv extension") }
  }

  # find intensity columns
  marker.matches <- lapply(markers, function(x) grep(colnames(data), pattern=paste0(x, ".*Intensity")))
  if (any(sapply(marker.matches, length)>1)) {
    if (is.null(locations)) {
      warning("The following markers matched to more than one intensity column: ", markers[sapply(marker.matches, length)>1])
      warning("Taking the first option in these cases. To avoid this behaviour, specify location.")
      marker.matches <- sapply(marker.matches, min)
    } else {
      marker.matches <- lapply(markers, function(x) grep(colnames(data), pattern=paste0(x, ".*", locations, ".*Intensity")))
      marker.matches <- unlist(marker.matches)
    }
  } else {
    marker.matches <- unlist(marker.matches)
  }
  intensity_matrix <- data[,marker.matches]
  names(intensity_matrix) <- markers
  intensity_matrix <- t(intensity_matrix)

  if(reference) {
    # find classification columns
    marker.matches <- lapply(markers, function(x) grep(colnames(data), pattern=paste0(x, ".*Classification")))
    if (any(sapply(marker.matches, length)>1)) {
      if (is.null(locations)) {
        warning("The following markers matched to more than one classification column: ", markers[sapply(marker.matches, length)>1])
        warning("Taking the first option in these cases. To avoid this behaviour, specify location.")
        marker.matches <- sapply(marker.matches, min)
      } else {
        marker.matches <- lapply(markers, function(x) grep(colnames(data), pattern=paste0(x, ".*", locations, ".*Classification")))
        marker.matches <- unlist(marker.matches)
      }
    } else {
      marker.matches <- unlist(marker.matches)
    }
    reference_matrix <- data[,marker.matches]
    names(reference_matrix) <- markers
    reference_matrix <- t(reference_matrix)
  }


  # find meta columns
  meta_selection <- grep(colnames(data), pattern="^Object|DAPI.*Classification", value = T)
  if (filter_dapi) { dapi.column <- grep(colnames(data), pattern="DAPI.*Classification", value = T) }
  col_data <- data[,meta_selection]

  # find spatial coords
  spatial_coords <- data[,c("XMin","XMax", "YMin", "YMax")]
  spatial_coords <- spatial_coords %>% mutate(Cell.X.Position=(XMin+XMax)/2, Cell.Y.Position=(YMin+YMax)/2) %>% select(Cell.X.Position, Cell.Y.Position) %>% as.matrix


  spe <- SpatialExperiment(assays=list("data"=intensity_matrix, "reference"=reference_matrix), colData = col_data, spatialCoords = spatial_coords)

  if (filter_dapi) { spe <- spe[,data[,dapi.column] == 1] }

  return(spe)

}
