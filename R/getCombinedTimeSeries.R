#' @title Get combined time series
#' @description This is a utility function to sum the raw data time series that correspond to the specified compartments
#' @param model The model from which to get the data values
#' @param compartments The vector of compartment names that contain the values to combine
#' @return The combined time series
#' @keywords internal
getCombinedTimeSeries <- function(model, compartments) {
  compartmentLength = length(model$parameters$populationFractions)
  if (length(compartments) > 1) {
    if (compartmentLength > 1) {
      return(sapply(1:compartmentLength, function(x){rowSums(model$rawOutput[, paste0(compartments, x), drop = FALSE])}))
    } else {
      return(as.matrix(rowSums(model$rawOutput[, compartments, drop = FALSE])))
    }
  } else {
    if (compartmentLength > 1) {
      return(model$rawOutput[, paste0(compartments, 1:compartmentLength), drop = FALSE])
    } else {
      return(model$rawOutput[, compartments, drop = FALSE])
    }
  }
}