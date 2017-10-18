#' @title Get labels
#' @description This is a utility function to obtain the numbered compartment names
#' associated with the specified compartment string
#' @param compartmentName The base string name of the compartment
#' @param compartmentLength The length of the vector compartment
#' @keywords internal
getLabels <- function(compartmentName, compartmentLength) {
  if(compartmentLength > 1) {
    return(paste0(compartmentName, 1:compartmentLength))
  } else {
    return(compartmentName)
  }
}
