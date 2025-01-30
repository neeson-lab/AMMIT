#' Title
#'
#' @param HALO_filepath Character; path to your HALO output file (.csv or .tsv), generally named 'Total_Object_Results.csv'
#' @param markers Character; a vector of the markers (i.e., proteins being measured) in the HALO output
#' @param reference Logical; whether the input file has manual thresholds computed, i.e. whether there are columns named "marker.Positive.Classification" or similar.
#' @param locations Character; a vector of the locations associated with each marker (e.g. Cytoplasm, Nucleus) for selecting the desired intensity data, if there are multiple intensity columns per marker.
#' @param filter_dapi Logical (Default TRUE); whether to filter the input to only DAPI-positive cells.
#'
#' @return A SpatialExperiment object with assay 'data' containing the intensities per cell for each marker, spatial coordinates, metadata (in colData), and positive classification (if provided) in the assay 'reference'.
#' @export
#'
#' @examples
build_ammit_spe_HALO = function (HALO_filepath,
                                 markers,
                                 reference,
                                 locations=NULL,
                                 filter_dapi=TRUE) {
  if (!file.exists(HALO_filepath)) {
    stop("Cannot find the file specified by `HALO_filepath`.")
    }
  if (!is.null(locations) && length(locations)!=length(markers)) {
    stop("`locations` and `markers` should be character vectors of the same length.")
    }
  if (!methods::is(reference, "logical")) {
    stop("`reference` should be set to TRUE or FALSE depending on whether your data includes pre-calculated positivity classification for markers")
    }
  if (grepl(pattern="\\.csv$", HALO_filepath, ignore.case=T)) {
    data <- utils::read.csv(HALO_filepath)
  } else if (grepl(pattern="\\.tsv$", HALO_filepath, ignore.case=T)) {
    data <- utils::read.delim(HALO_filepath)
  } else {
    stop("If specifying `HALO_filepath`, the file in question should have a .csv or .tsv extension")
  }

  # find intensity columns
  marker.matches <- lapply(markers, \ (x) grep(colnames(data), pattern=paste0(x, ".*Intensity")))
  if (any(sapply(marker.matches, length)>1)) {
    if (is.null(locations)) {
      warning("The following markers matched to more than one intensity column: ", markers[sapply(marker.matches, length) > 1])
      warning("Taking the first option in these cases. To avoid this behaviour, specify location.")
      marker.matches <- sapply(marker.matches, min)
    } else {
      marker.matches <- lapply(markers, \ (x) grep(colnames(data), pattern=paste0(x, ".*", locations, ".*Intensity")))
      marker.matches <- unlist(marker.matches)
    }
  } else {
    marker.matches <- unlist(marker.matches)
  }
  intensity_matrix <- data[, marker.matches]
  names(intensity_matrix) <- markers
  intensity_matrix <- t(intensity_matrix)

  if (reference) {
    marker.matches <- lapply(markers, \ (x) grep(colnames(data), pattern=paste0(x, ".*Classification")))
    if (any(sapply(marker.matches, length) > 1)) {
      if (is.null(locations)) {
        warning("The following markers matched to more than one classification column: ", markers[sapply(marker.matches, length) > 1])
        warning("Taking the first option in these cases. To avoid this behaviour, specify location.")
        marker.matches <- sapply(marker.matches, min)
      } else {
        marker.matches <- lapply(markers, \ (x) grep(colnames(data), pattern=paste0(x, ".*", locations, ".*Classification")))
        marker.matches <- unlist(marker.matches)
      }
    } else {
      marker.matches <- unlist(marker.matches)
    }
    reference_matrix <- data[, marker.matches]
    names(reference_matrix) <- markers
    reference_matrix <- t(reference_matrix)
  }

  meta_selection <- grep(colnames(data), pattern="^Object|DAPI.*Classification", value = TRUE)
  if (filter_dapi) {
    dapi.column <- grep(colnames(data), pattern="DAPI.*Classification", value = TRUE)
    }
  col_data <- data[, meta_selection]

  spatial_coords <- data[, c("XMin","XMax", "YMin", "YMax")]
  spatial_coords <- spatial_coords |>
    dplyr::mutate(Cell.X.Position=(XMin+XMax) / 2, Cell.Y.Position=(YMin+YMax) / 2) |>
    dplyr::select(Cell.X.Position, Cell.Y.Position) |>
    as.matrix()

  spe <- SpatialExperiment::SpatialExperiment(assays=list("data"=intensity_matrix, "reference"=reference_matrix), colData = col_data, spatialCoords = spatial_coords)

  if (filter_dapi) {
    spe <- spe[,data[,dapi.column] == 1]
    }

  return(spe)
}
