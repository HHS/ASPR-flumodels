#' @title Get day of peak
#' @description Gets the day of peak infections from the model
#' @details This is the generic function, for class-specific versions see:
#' \itemize{
#' 	\item \code{\link{getDayOfPeak.SEIRModel}} for SEIR-type models
#' }
#' @param model A flumodels model object
#' @param ... Other parameters for the class-specific functions
#' @export
getDayOfPeak <- function(model, ...) {
  UseMethod("getDayOfPeak", model)
}

#SEIR
#' @title Get day of peak
#' @description Gets the day of peak infections from the model
#' @param model The model from which to get the data
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @return The day of peak infections
#' @method getDayOfPeak SEIRModel
#' @keywords internal
#' @export
getDayOfPeak.SEIRModel <- function(model, incidence = FALSE, symptomatic = FALSE) {
  dailyTotalInfections <- getInfectionTimeSeries(model, byGroup = FALSE, incidence = incidence, symptomatic = symptomatic,
                                                 fractionSymptomatic = 1) #Fraction symptomatic is irrelevant
  return(which.max(dailyTotalInfections) - 1) #Find the index of the maximum and subtract 1, since time starts at 0
}

#SEAIR
#' @title Get day of peak
#' @description Gets the day of peak infections from the model
#' @param model The model from which to get the data
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @return The day of peak infections
#' @method getDayOfPeak SEIRModel
#' @keywords internal
#' @export
getDayOfPeak.SEAIRTVModel <- function(model, incidence = FALSE, symptomatic = FALSE) {
  dailyTotalInfections <- getInfectionTimeSeries(model, byGroup = FALSE, incidence = incidence, symptomatic = symptomatic) #Fraction symptomatic is irrelevant
  return(which.max(dailyTotalInfections) - 1) #Find the index of the maximum and subtract 1, since time starts at 0
}
