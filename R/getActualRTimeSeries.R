#' @title Get actual reproductive number time series
#' @description Gets the time series of the actual observed effective reproductive 
#'   numbers for the model.
#' @param model A flumodels model object
#' @keywords internal
getActualRTimeSeries <- function(model) {
  UseMethod("getActualRTimeSeries", model)
}

#SEIR
#' @rdname getActualRTimeSeries
#' @method getActualRTimeSeries SEIRModel
#' @keywords internal
#' @export
getActualRTimeSeries.SEIRModel <- function(model) {
  return(computeActualRTimeSeries(model, reconstructState.SEIR, getDerivative.SEIR))
}

#SEIRT
#' @rdname getActualRTimeSeries
#' @method getActualRTimeSeries SEIRTModel
#' @keywords internal
#' @export
getActualRTimeSeries.SEIRTModel <- function(model) {
  return(computeActualRTimeSeries(model, reconstructState.SEIR, getDerivative.SEIRT))
}

#SEIRV
#' @rdname getActualRTimeSeries
#' @method getActualRTimeSeries SEIRVModel
#' @keywords internal
#' @export
getActualRTimeSeries.SEIRVModel <- function(model) {
  return(computeActualRTimeSeries(model, reconstructState.SEIRV, getDerivative.SEIRV))
}

#SEIRTV
#' @rdname getActualRTimeSeries
#' @method getActualRTimeSeries SEIRTVModel
#' @keywords internal
#' @export
getActualRTimeSeries.SEIRTVModel <- function(model) {
  return(computeActualRTimeSeries(model, reconstructState.SEIRV, getDerivative.SEIRTV))
}

#SEIRV2Dose
#' @rdname getActualRTimeSeries
#' @method getActualRTimeSeries SEIRV2DoseModel
#' @keywords internal
#' @export
getActualRTimeSeries.SEIRV2DoseModel <- function(model) {
  return(computeActualRTimeSeries(model, reconstructState.SEIRV2Dose, getDerivative.SEIRV2Dose))
}

#SEIRTV2Dose
#' @rdname getActualRTimeSeries
#' @method getActualRTimeSeries SEIRTV2DoseModel
#' @keywords internal
#' @export
getActualRTimeSeries.SEIRTV2DoseModel <- function(model) {
  return(computeActualRTimeSeries(model, reconstructState.SEIRV2Dose, getDerivative.SEIRTV2Dose))
}

#' @title Compute actual reproductive number time series
#' @description This is a helper function that implements the calculations in
#'   the various 'getActualRTimeSeries' model class-specific functions
#' @param model A flumodels model object
#' @param reconstructState A function that reconstructs the state of the model
#' @param getDerivative A function that computes the derivative of the model
#' @param model A flumodels model object
#' @keywords internal
#' @export
computeActualRTimeSeries <- function(model, reconstructState, getDerivative) {
  susceptibleCompartments <- getCompartments(model, "S")
  infectiousCompartments <- getCompartments(model, "I")
  return(apply(model$rawOutput, 1,
               function(row) {
                 state <- reconstructState(row[-1])
                 derivatives <- reconstructState(unlist(getDerivative(t = row[1],
                                                                      state = row[-1],
                                                                      parameters = model$parameters)))
                 infections <- do.call(sum, state[infectiousCompartments])
                 return(ifelse(infections > 0,
                               -do.call(sum, derivatives[susceptibleCompartments]) / infections 
                               / model$parameters$gamma,
                               NA))
               }))
}