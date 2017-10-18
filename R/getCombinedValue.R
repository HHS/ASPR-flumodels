#' @title Get combined values
#' @description This is a utility function to sum the raw data values that correspond to the specified compartments at the given index
#' @param model The model from which to get the data values
#' @param compartments The vector of compartment names that contain the values to combine
#' @param index The index in the array from which data is desired
#' @return The combined value
#' @keywords internal
getCombinedValue <- function(model, compartments, index) {
  compartmentLength = length(model$parameters$populationFractions)
  if (length(compartments) > 1) {
    if (compartmentLength > 1) {
      return(sapply(1:compartmentLength, function(x){sum(model$rawOutput[index, paste0(compartments, x)])}))
    } else {
      return(sum(model$rawOutput[index, compartments]))
    }
  } else {
    if (compartmentLength > 1) {
      return(model$rawOutput[index, paste0(compartments, 1:compartmentLength)])
    } else {
      return(model$rawOutput[index, compartments])
    }
  }
}