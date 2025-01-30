infer_manual_thresholds <- function(spe, markers) {

  if (class(spe)!="SpatialExperiment") { stop("Your `spe` object is not a SpatialExperiment object!") }
  if (!"reference" %in% names(assays(spe))) { stop("No reference data found in your SpatialExperiment object. Have you added reference data?") }
  if (!"data" %in% names(assays(spe))) { stop("No intensity data found in your SpatialExperiment object. Have you added intensity data?") }
  if (!any(markers %in% rownames(spe))) { stop("None of the specified markers were found in your SpatialExperiment object!")}
  if (!all(rownames(spe) %in% markers)) { warning("The following markers were not specified, and manual thresholds for them will not be calculated: ", rownames(se)[!rownames(se) %in% markers])}

  manual_thresholds <- sapply(markers, function(marker) {
    positives <- assay(spe, "reference")[marker,]==1
    threshold <- min(assay(spe, "data")[marker,positives])
    return(threshold)
  })

  output <- data.frame(row.names = rownames(spe))

  output[match(markers, rownames(spe)),"manual_threshold"] <- manual_thresholds

  if ("manual_threshold" %in% colnames(rowData(spe))) {
    rowData(spe)$manual_threshold <- coalesce(output$manual_threshold, rowData(spe)$manual_threshold)
  } else {
    rowData(spe)$manual_threshold <- output$manual_threshold
  }

}
