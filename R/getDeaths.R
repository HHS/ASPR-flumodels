#' @title Get deaths
#' @description Gets the total deaths from the given model
#' @details This is the generic function, for class-specific versions see:
#' \itemize{
#'  \item \code{\link{getDeaths.SEIRModel}} for SEIR-type models
#'  \item \code{\link{getDeaths.SEIRTModel}} for SEIRT-type models
#' }
#' @param model A flumodels model object
#' @param ... Other parameters for the class-specific functions
#' @export
getDeaths <- function(model, ...) {
  UseMethod("getDeaths", model)
}

#SEIR
#' @title Get deaths
#' @description Gets the total deaths from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; single fraction or vector of fractions by population group; must be specified
#' @param caseFatalityRatio The fraction of cases that result in fatalities; single fraction or vector of fractions by population group; must be specified
#' @param timeSeries Whether to return a time series or a final cumulative value
#' @return A vector of deaths or death rates
#' @method getDeaths SEIRModel
#' @keywords internal
#' @export
getDeaths.SEIRModel <- function(model, byGroup = TRUE, asRate = FALSE, fractionSymptomatic, caseFatalityRatio, timeSeries = FALSE) {
  if (missing(fractionSymptomatic)) {
    stop("fractionSymptomatic must be specified.")
  }
  checkBetween0and1(fractionSymptomatic)
  checkDimensionsMatch(fractionSymptomatic, model$parameters$populationFractions)
  if (missing(caseFatalityRatio)) {
    stop("caseFatalityRatio must be specified.")
  }
  checkBetween0and1(caseFatalityRatio)
  checkDimensionsMatch(caseFatalityRatio, model$parameters$populationFractions)
  
  if (timeSeries) {
    deaths <- t(t(getInfectionTimeSeries(model, byGroup = TRUE, asRate = asRate, symptomatic = TRUE, fractionSymptomatic = fractionSymptomatic, incidence = TRUE)) * 
      caseFatalityRatio)
  }
  else {
    deaths <- caseFatalityRatio * getInfections(model, byGroup = TRUE, asRate = asRate, symptomatic = TRUE, fractionSymptomatic = fractionSymptomatic)
  }
  
  if (asRate) {
    names(deaths) <- getLabels("deathRate", length(model$parameters$populationFractions)) 
  } else { #Number
    names(deaths) <- getLabels("deaths", length(model$parameters$populationFractions))
  }
  if (byGroup) {
    return(deaths)
  } else {
    if (timeSeries) {
      return(rowSums(deaths))
    }
    else {
    return(sum(deaths))
    }
  }
}

#SEIRT
#' @title Get deaths
#' @description Gets the total deaths from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param caseFatalityRatio The fraction of cases that result in fatalities; single fraction or vector of fractions by population group; must be specified
#' @param timeSeries Whether to return a time series or a final cumulative value
#' @return A vector of deaths or death rates
#' @method getDeaths SEIRTModel
#' @keywords internal
#' @export
getDeaths.SEIRTModel <- function(model, byGroup = TRUE, asRate = FALSE, caseFatalityRatio, timeSeries = FALSE) {
  if (missing(caseFatalityRatio)) {
    stop("caseFatalityRatio must be specified.")
  }
  checkBetween0and1(caseFatalityRatio)
  checkDimensionsMatch(caseFatalityRatio, model$parameters$populationFractions)

  AVEp.outpatient.eff <- model$parameters$AVEp * 
                           model$parameters$fractionAdhere * 
                           model$parameters$fractionDiagnosedAndPrescribedOutpatient * 
                           model$parameters$fractionSeekCare
  AVEp.inpatient.eff <-  model$parameters$AVEp * model$parameters$fractionDiagnosedAndPrescribedInpatient *
                           model$parameters$fractionAdmitted * (1 - AVEp.outpatient.eff)
  
  if (timeSeries) {
    deaths <- t(t(getInfectionTimeSeries(model, byGroup = TRUE, asRate = asRate, symptomatic = TRUE, incidence = TRUE)) * 
      (1 - AVEp.outpatient.eff - AVEp.inpatient.eff) * caseFatalityRatio)
  }
  else {
    deaths <- caseFatalityRatio * (1 - AVEp.outpatient.eff - AVEp.inpatient.eff) * 
              getInfections(model, byGroup = TRUE, asRate = asRate, symptomatic = TRUE)
  }
  if (asRate) {
    names(deaths) <- getLabels("deathRate", length(model$parameters$populationFractions)) 
  } else { #Number
    names(deaths) <- getLabels("deaths", length(model$parameters$populationFractions))
  }
  if (byGroup) {
    return(deaths)
  } else {
    if (timeSeries) {
      return(rowSums(deaths))
    }
    else {
    return(sum(deaths))
    }
  }
}

#SEAIRTV
#' @title Get deaths
#' @description Gets the total deaths from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param caseFatalityRatio The fraction of cases that result in fatalities; single fraction or vector of fractions by population group; must be specified
#' @param timeSeries Whether to return a time series or a final cumulative value
#' @return A vector of deaths or death rates
#' @method getDeaths SEAIRTVModel
#' @keywords internal
#' @export
getDeaths.SEAIRTVModel <- function(model, byGroup = TRUE, asRate = FALSE, caseFatalityRatio, timeSeries = FALSE) {
  getDeaths.SEIRTModel(model, byGroup, asRate, caseFatalityRatio, timeSeries)
}