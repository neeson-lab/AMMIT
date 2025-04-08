#' Build a SpatialExperiment object suitable for AMMIT
#'
#' @param method Character; "HALO" if reading from a HALO output file; otherwise, "manual" if you are providing intensity matrix, reference matrix, spatial coords, etc. directly.
#' @param reference Logical; whether the input file has manual thresholds computed, i.e. whether the data for reference_matrix is present and should be included in the output.
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
#' HALO_filepath <- system.file("extdata", "halo_output.csv", package="AMMIT")
#' markers <- c("M1", "M2", "M3", "M4",
#'              "M5", "M6", "M7")
#' locations <-  c("Nucleus", "Cytoplasm", "Nucleus", "Cytoplasm",
#'                 "Cytoplasm", "Cytoplasm", "Cytoplasm")
#' spe <- build_ammit_spe(method="HALO",
#'                        HALO_filepath = HALO_filepath,
#'                        markers = markers,
#'                        reference = TRUE,
#'                        locations = locations,
#'                        filter_dapi = TRUE)
build_ammit_spe = function(method=c("HALO", "manual"),
                           reference,
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


  if (method=="HALO") { spe = build_ammit_spe_HALO(reference=reference, ...) }

  if (method=="manual") {
      if (reference) {
        spe = SpatialExperiment::SpatialExperiment(assays=list("data"=intensity_matrix, "reference"=reference_matrix), colData = col_data, spatialCoords = spatial_coords)
      } else {
        spe = SpatialExperiment::SpatialExperiment(assays=list("data"=intensity_matrix), colData = col_data, spatialCoords = spatial_coords)
      }
    }

  return(spe)
}

#' Build a SpatialExperiment object from a HALO output suitable for AMMIT
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
#' HALO_filepath <- system.file("extdata", "halo_output.csv", package="AMMIT")
#' markers <- c("M1", "M2", "M3", "M4",
#'              "M5", "M6", "M7")
#' locations <-  c("Nucleus", "Cytoplasm", "Nucleus", "Cytoplasm",
#'                 "Cytoplasm", "Cytoplasm", "Cytoplasm")
#' spe <- build_ammit_spe_HALO(HALO_filepath = HALO_filepath,
#'                             markers = markers,
#'                             reference = TRUE,
#'                             locations = locations,
#'                             filter_dapi = TRUE)
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
  if (any(sapply(marker.matches, length)==0)) { stop("One or more markers not found!") }
  if (any(sapply(marker.matches, length)>1)) {
    if (is.null(locations)) {
      warning("The following markers matched to more than one intensity column: ", paste(markers[sapply(marker.matches, length) > 1], collapse=" "))
      warning("Taking the first option in these cases. To avoid this behaviour, specify `location`.")
      marker.matches <- sapply(marker.matches, min)
    } else {
      marker.matches <- lapply(markers, \ (x) grep(x = colnames(data), pattern=paste0(x, ".*", locations[match(x, markers)], ".*Intensity")))
      marker.matches <- unlist(marker.matches)
    }
  } else {
    marker.matches <- unlist(marker.matches)
  }
  intensity_matrix <- data[, marker.matches]
  names(intensity_matrix) <- markers
  intensity_matrix <- t(intensity_matrix)

  # check for NAs in intensity matrix
  if (any(apply(intensity_matrix, MARGIN=2, \(x) any(is.na(x))))) {
    warning("Some of your intensity values are NA. Removing these 'cells'... ")
    rows_to_remove <- which(apply(intensity_matrix, MARGIN=2, \(x) any(is.na(x))))
    warning("'cells' removed: ", length(rows_to_remove))
    intensity_matrix <- intensity_matrix[,-rows_to_remove]
    data <- data[-rows_to_remove,]
  }


  meta_selection <- grep(colnames(data), pattern=paste(markers, collapse="|"), value = TRUE, invert = TRUE)
  if (filter_dapi) {
    dapi.column <- grep(colnames(data), pattern="DAPI.*Classification", value = TRUE)[1]
  }
  col_data <- data[, meta_selection]


  spatial_coords <- data.frame("Cell.X.Position"=(data$XMin+data$XMax)/2,
                               "Cell.Y.Position"=(data$YMin+data$YMax)/2) |>
    as.matrix()

  if (reference) {
    marker.matches <- lapply(markers, \ (x) grep(colnames(data), pattern=paste0(x, ".*Classification")))
    if (any(sapply(marker.matches, length) > 1)) {
      if (is.null(locations)) {
        warning("The following markers matched to more than one classification column: ", markers[sapply(marker.matches, length) > 1])
        warning("Taking the first option in these cases. To avoid this behaviour, specify location.")
        marker.matches <- sapply(marker.matches, min)
      } else {
        marker.matches <- lapply(markers, \ (x) grep(x = colnames(data), pattern=paste0(x, ".*", locations[match(x, markers)], ".*Classification")))
        marker.matches <- unlist(marker.matches)
      }
    } else {
      marker.matches <- unlist(marker.matches)
    }
    reference_matrix <- data[, marker.matches]
    names(reference_matrix) <- markers
    reference_matrix <- t(reference_matrix)
    spe <- SpatialExperiment::SpatialExperiment(assays=list("data"=intensity_matrix, "reference"=reference_matrix), colData = col_data, spatialCoords = spatial_coords)
  } else {
    spe <- SpatialExperiment::SpatialExperiment(assays=list("data"=intensity_matrix), colData = col_data, spatialCoords = spatial_coords)
  }

  if (filter_dapi) {
    spe <- spe[,data[,dapi.column] == 1]
  }

  return(spe)
}

