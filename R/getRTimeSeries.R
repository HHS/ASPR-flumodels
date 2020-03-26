#' @title Get effective reproductive number time series
#' @description Gets the time series of effective reproductive numbers for the model.
#' @param model A flumodels model object
#' @param actual A boolean value that determines whether to calculate the actual observed
#'   effective R value (actual = TRUE), or the theoretical effective R value.
#' @return The effective R time series
#' @export
getRTimeSeries <- function(model, actual = FALSE) {
  if (actual) { #Calculate from derivatives
    return(getActualRTimeSeries(model))
  } else {
    susceptibilityVectorTimeSeries <- getSusceptibilityVectorTimeSeries(model)
    infectiousnessVectorTimeSeries <- getInfectiousnessVectorTimeSeries(model)
    infectiousPeriod <- if ("SEIRModel" %in% class(model)) {
      1 / model$parameters$gamma
    } else {
      1 / model$parameters$lambda2 + 1 / model$parameters$gamma
    }

    compartmentLength <- length(model$parameters$populationFractions)
    contactMatrixArray <- array(model$parameters$contactMatrix,
                                dim = c(compartmentLength,
                                        compartmentLength,
                                        nrow(susceptibilityVectorTimeSeries)))
    if (model$parameters$useCommunityMitigation) {
      contactMatrixArray[, , (model$parameters$communityMitigationStartDay + 1):
                             model$parameters$communityMitigationEndDay] <-
        model$parameters$communityMitigationMultiplier * model$parameters$contactMatrix
    }

    effectiveRTimeSeries <-  unlist(lapply(1:nrow(susceptibilityVectorTimeSeries), function(x) {
        effectiveTransmissionMatrix <- t(t(susceptibilityVectorTimeSeries[x, ] * contactMatrixArray[, , x]) 
                                         * infectiousnessVectorTimeSeries[x, ])
        return(Mod(eigen(effectiveTransmissionMatrix, symmetric = FALSE, only.values = TRUE)$values)[1])
      })) * model$parameters$beta * infectiousPeriod
    return(effectiveRTimeSeries)
  }
}


#' @title Get susceptibility vector time series
#' @description Gets the time series of effective fractions of the population
#'   by age group that are susceptible to infection.
#' @param model A flumodels model object
#' @keywords internal
getSusceptibilityVectorTimeSeries <- function(model) {
  UseMethod("getSusceptibilityVectorTimeSeries", model)
}

#' @rdname getSusceptibilityVectorTimeSeries
#' @method getSusceptibilityVectorTimeSeries SEIRModel
#' @keywords internal
#' @export
getSusceptibilityVectorTimeSeries.SEIRModel <- function(model) {
  return(sweep(getCombinedTimeSeries(model, "S"),
               2, model$parameters$populationFractions, "/"))
}

#' @rdname getSusceptibilityVectorTimeSeries
#' @method getSusceptibilityVectorTimeSeries SEIRVModel
#' @keywords internal
#' @export
getSusceptibilityVectorTimeSeries.SEIRVModel <- function(model) {
  STimeSeries <- getCombinedTimeSeries(model, "S")
  SvTimeSeries <- getCombinedTimeSeries(model, "Sv")
  return(sweep(STimeSeries + (1 - model$parameters$VEs) * SvTimeSeries,
               2, model$parameters$populationFractions, "/"))
}

#' @rdname getSusceptibilityVectorTimeSeries
#' @method getSusceptibilityVectorTimeSeries SEIRVModel
#' @keywords internal
#' @export
getSusceptibilityVectorTimeSeries.SEAIRTVModel <- function(model) {
  return(getSusceptibilityVectorTimeSeries.SEIRVModel(model))
}

#' @rdname getSusceptibilityVectorTimeSeries
#' @method getSusceptibilityVectorTimeSeries SEIRV2DoseModel
#' @keywords internal
#' @export
getSusceptibilityVectorTimeSeries.SEIRV2DoseModel <- function(model) {
  STimeSeries <- getCombinedTimeSeries(model, "S")
  SvTimeSeries <- getCombinedTimeSeries(model, "Sv")
  SvbTimeSeries <- getCombinedTimeSeries(model, "Svb")
  return(sweep(STimeSeries + 
               (1 - model$parameters$VEs1) * SvTimeSeries +
               (1 - model$parameters$VEs2) * SvbTimeSeries,
               2, model$parameters$populationFractions, "/"))
}

#' @title Get infectiousness vector time series
#' @description Gets the time series of effective fractions of the population
#'   by age group that are infectious.
#' @param model A flumodels model object
#' @keywords internal
getInfectiousnessVectorTimeSeries <- function(model) {
  UseMethod("getInfectiousnessVectorTimeSeries", model)
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEIRModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEIRModel <- function(model) {
  return(matrix(1, nrow = nrow(model$rawOutput),
                ncol = length(model$parameters$populationFractions)))
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEIRTModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEIRTModel <- function(model) {
  return(matrix(1 - model$parameters$AVEi.eff, nrow = nrow(model$rawOutput),
                ncol = length(model$parameters$populationFractions), byrow = TRUE))
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEIRVModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEIRVModel <- function(model) {
  ITimeSeries <- getCombinedTimeSeries(model, "I")
  IvTimeSeries <- getCombinedTimeSeries(model, "Iv")
  return(ifelse(IvTimeSeries > 0,
                1 - (model$parameters$VEi * IvTimeSeries / (ITimeSeries + IvTimeSeries)),
                1))
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEIRTVModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEIRTVModel <- function(model) {
  ITimeSeries <- getCombinedTimeSeries(model, "I")
  IvTimeSeries <- getCombinedTimeSeries(model, "Iv")
  return(sweep(ifelse(IvTimeSeries > 0,
                      1 - (model$parameters$VEi * IvTimeSeries / (ITimeSeries + IvTimeSeries)),
                      1),
               2, 1 - model$parameters$AVEi.eff, "*"))
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEAIRTVModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEAIRTVModel <- function(model) {
  ITimeSeries <- getCombinedTimeSeries(model, "A") + getCombinedTimeSeries(model, "I")
  IvTimeSeries <- getCombinedTimeSeries(model, "Av") + getCombinedTimeSeries(model, "Iv")
  return(sweep(ifelse(IvTimeSeries > 0,
                      1 - (model$parameters$VEi * IvTimeSeries / (ITimeSeries + IvTimeSeries)),
                      1),
               2, 1 - model$parameters$AVEi.eff, "*"))
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEIRV2DoseModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEIRV2DoseModel <- function(model) {
  ITimeSeries <- getCombinedTimeSeries(model, "I")
  IvTimeSeries <- getCombinedTimeSeries(model, "Iv")
  IvbTimeSeries <- getCombinedTimeSeries(model, "Ivb")
  ITotalTimeSeries <- ITimeSeries + IvTimeSeries + IvbTimeSeries
  return(ifelse(ITotalTimeSeries > 0,
                1 - (model$parameters$VEi1 * IvTimeSeries / ITotalTimeSeries)
                  - (model$parameters$VEi2 * IvbTimeSeries / ITotalTimeSeries),
                1))
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEIRTV2DoseModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEIRTV2DoseModel <- function(model) {
  ITimeSeries <- getCombinedTimeSeries(model, "I")
  IvTimeSeries <- getCombinedTimeSeries(model, "Iv")
  IvbTimeSeries <- getCombinedTimeSeries(model, "Ivb")
  ITotalTimeSeries <- ITimeSeries + IvTimeSeries + IvbTimeSeries
  return(sweep(ifelse(ITotalTimeSeries > 0,
                      1 - (model$parameters$VEi1 * IvTimeSeries / ITotalTimeSeries)
                        - (model$parameters$VEi2 * IvbTimeSeries / ITotalTimeSeries),
                      1),
               2, 1 - model$parameters$AVEi.eff, "*"))
}

#' @rdname getSusceptibilityVectorTimeSeries
#' @method getSusceptibilityVectorTimeSeries SEIRVMonoModel
#' @keywords internal
#' @export
getSusceptibilityVectorTimeSeries.SEIRVMonoModel <- function(model) {
  STimeSeries <- getCombinedTimeSeries(model, "S")
  SvTimeSeries <- getCombinedTimeSeries(model, "Sv")
  SvMTimeSeries <- getCombinedTimeSeries(model, "SvM")
  return(sweep(STimeSeries + 
                 (1 - model$parameters$VEs1) * SvTimeSeries +
                 (1 - model$parameters$VEs2) * SvMTimeSeries,
               2, model$parameters$populationFractions, "/"))
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEIRVMonoModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEIRVMonoModel <- function(model) {
  ITimeSeries <- getCombinedTimeSeries(model, "I")
  IvTimeSeries <- getCombinedTimeSeries(model, "Iv")
  IvMTimeSeries <- getCombinedTimeSeries(model, "IvM")
  ITotalTimeSeries <- ITimeSeries + IvTimeSeries + IvMTimeSeries
  return(ifelse(ITotalTimeSeries > 0,
                1 - (model$parameters$VEi1 * IvTimeSeries / ITotalTimeSeries)
                - (model$parameters$VEiM * IvMTimeSeries / ITotalTimeSeries),
                1))
}

#' @rdname getSusceptibilityVectorTimeSeries
#' @method getSusceptibilityVectorTimeSeries SEIRVPrimeBoostModel
#' @keywords internal
#' @export
getSusceptibilityVectorTimeSeries.SEIRVPrimeBoostModel <- function(model) {
    getSusceptibilityVectorTimeSeries.SEIRV2DoseMosel(model)
}

#' @rdname getInfectiousnessVectorTimeSeries
#' @method getInfectiousnessVectorTimeSeries SEIRVPrimeBoostModel
#' @keywords internal
#' @export
getInfectiousnessVectorTimeSeries.SEIRVPrimeBoostModel <- function(model) {
    getInfectiousnessVectorTimeSeries.SEIRV2DoseModel(model)
}
