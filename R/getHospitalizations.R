#' @title Get hospitalizations
#' @description Gets the total hospitalizations from the given model
#' @details This is the generic function, for class-specific versions see:
#' \itemize{
#'  \item \code{\link{getHospitalizations.SEIRModel}} for SEIR-type models
#'  \item \code{\link{getHospitalizations.SEIRTModel}} for SEIRT-type models
#' }
#' @param model A flumodels model object
#' @param ... Other parameters for the class-specific functions
#' @export
getHospitalizations <- function(model, ...) {
  UseMethod("getHospitalizations", model)
}

#SEIR
#' @title Get hospitalizations
#' @description Gets the total hospitalizations from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; single fraction or vector of fractions by population group; must be specified
#' @param caseHospitalizationRatio The fraction of cases that are hospitalized; single fraction or vector of fractions by population group; must be specified
#' @param timeSeries Whether to return a time series or a final cumulative value
#' @return A vector of hospitalizations or hospitalization rates
#' @method getHospitalizations SEIRModel
#' @keywords internal
#' @export
getHospitalizations.SEIRModel <- function(model, byGroup = TRUE, asRate = FALSE, fractionSymptomatic, caseHospitalizationRatio, timeSeries = FALSE) {
  if (missing(fractionSymptomatic)) {
    stop("fractionSymptomatic must be specified.")
  }
  checkBetween0and1(fractionSymptomatic)
  checkDimensionsMatch(fractionSymptomatic, model$parameters$populationFractions)
  if (missing(caseHospitalizationRatio)) {
    stop("caseHospitalizationRatio must be specified.")
  }
  checkBetween0and1(caseHospitalizationRatio)
  checkDimensionsMatch(caseHospitalizationRatio, model$parameters$populationFractions)

  if (timeSeries) {
    hospitalizations <- t(t(getInfectionTimeSeries(model, byGroup = TRUE, asRate = asRate, symptomatic = TRUE, fractionSymptomatic = fractionSymptomatic, incidence = TRUE)) * 
      caseHospitalizationRatio)
  }
  else {
    hospitalizations <- caseHospitalizationRatio * getInfections(model, byGroup = TRUE, asRate = asRate, symptomatic = TRUE, fractionSymptomatic = fractionSymptomatic)
  }                  
      
  if (asRate) {
    names(hospitalizations) <- getLabels("hospitalizationRate", length(model$parameters$populationFractions)) 
  } else { #Number
    names(hospitalizations) <- getLabels("hospitalizations", length(model$parameters$populationFractions))
  }
  if (byGroup) {
    return(hospitalizations)
  } else {
    if (timeSeries) {
      return(rowSums(hospitalizations))
    }
    else {
      return(sum(hospitalizations))
    }
  }
}

#SEIRT
#' @title Get hospitalizations
#' @description Gets the total hospitalizations from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param caseHospitalizationRatio The fraction of cases that are hospitalized; single fraction or vector of fractions by population group; must be specified
#' @param timeSeries Whether to return a time series or a final cumulative value
#' @return A vector of hospitalizations or hospitalization rates
#' @method getHospitalizations SEIRTModel
#' @keywords internal
#' @export
getHospitalizations.SEIRTModel <- function(model, byGroup = TRUE, asRate = FALSE, caseHospitalizationRatio, timeSeries = FALSE) {
  if (missing(caseHospitalizationRatio)) {
    stop("caseHospitalizationRatio must be specified.")
  }
  AVEp.outpatient.eff <- model$parameters$AVEp * 
                           model$parameters$fractionAdhere * 
                           model$parameters$fractionDiagnosedAndPrescribedOutpatient * 
                           model$parameters$fractionSeekCare
  
  if (timeSeries) {
    hospitalizations <- t(t(getInfectionTimeSeries(model, byGroup = TRUE, asRate = asRate, symptomatic = TRUE, incidence = TRUE)) * 
      (1 - AVEp.outpatient.eff) * caseHospitalizationRatio)
  }
  else {
    hospitalizations <- getInfections(model, byGroup = TRUE, asRate = asRate, symptomatic = TRUE) * (1 - AVEp.outpatient.eff) *
                          caseHospitalizationRatio
  }
  
  if (asRate) {
    names(hospitalizations) <- getLabels("hospitalizationRate", length(model$parameters$populationFractions)) 
  } else { #Number
    names(hospitalizations) <- getLabels("hospitalizations", length(model$parameters$populationFractions))
  }
  if (byGroup) {
    return(hospitalizations)
  } else {
    if (timeSeries) {
      return(rowSums(hospitalizations))
    }
    else {
      return(sum(hospitalizations))
    }
  }
}

#SEAIRTV
#' @title Get hospitalizations
#' @description Gets the total hospitalizations from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param caseHospitalizationRatio The fraction of cases that are hospitalized; single fraction or vector of fractions by population group; must be specified
#' @param timeSeries Whether to return a time series or a final cumulative value
#' @return A vector of hospitalizations or hospitalization rates
#' @method getHospitalizations SEAIRTVModel
#' @keywords internal
#' @export
getHospitalizations.SEAIRTVModel <- function(model, byGroup = TRUE, asRate = FALSE, caseHospitalizationRatio, timeSeries = FALSE) {
  getHospitalizations.SEIRTModel(model, byGroup, asRate, caseHospitalizationRatio, timeSeries)
}