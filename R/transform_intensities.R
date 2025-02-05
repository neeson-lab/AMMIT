#' Transform marker intensities for unmixing or plotting
#'
#' @param intensities Numeric; a vector of marker intensities (or an intensity threshold)
#' @param method Character; one of: 'none' for no transformation, 'asinh' for arcsinh transformation, or 'sqrt' for square-root transformation.
#'
#' @return A transformed numeric vector.
#' @export
#'
#' @examples
#' intensities <- c(1,3,3,4,5,3,3)
#' transform_intensities(intensities, method = "asinh")
transform_intensities <- function(intensities, method) {
  if (!method %in% c("none", "asinh", "sqrt")) {
    stop("`method` should be one of: 'none', 'asinh', or 'sqrt'!")
  }
  if (method == "asinh") {
    intensities <- asinh(intensities)
  }
  if (method == "sqrt") {
    intensities <- sqrt(intensities)
  }
  return(intensities)
}
