#' @title Get infections
#' @description Gets the total infections from the given model
#' @details This is the generic function, for class-specific versions see:
#' \itemize{
#'  \item \code{\link{getInfections.SEIRModel}} for SEIR-type models
#'  \item \code{\link{getInfections.SEIRVModel}} for SEIRV-type models
#'  \item \code{\link{getInfections.SEIRV2DoseModel}} for SEIRV2Dose-type models
#'  \item \code{\link{getInfections.SEIRTModel}} for SEIRT-type models
#'  \item \code{\link{getInfections.SEIRTVModel}} for SEIRTV-type models
#'  \item \code{\link{getInfections.SEAIRTVModel}} for SEAIRTV-type models
#'  \item \code{\link{getInfections.SEIRTV2DoseModel}} for SEIRTV2Dose-type models
#' }
#' @param model A flumodels model object
#' @param ... Other parameters for the class-specific functions
#' @export
getInfections <- function(model, ...) {
  UseMethod("getInfections", model)
}

#SEIR
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @return A vector of infections 
#' @method getInfections SEIRModel
#' @keywords internal
#' @export
getInfections.SEIRModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE, fractionSymptomatic) {
  compartmentLength <- length(model$parameters$populationFractions)
  compartments <- getCompartments(model, type = "R")
  if (symptomatic) {
    if (missing(fractionSymptomatic)) {
      stop("fractionSymptomatic must be specified.")
    }
    attackRateByGroup <- (getCombinedValue(model, compartments, index = nrow(model$rawOutput)) -
                          getCombinedValue(model, compartments, index = 1)) * fractionSymptomatic
    if (asRate) {
      infectionsByGroup <- attackRateByGroup
      names(infectionsByGroup) <- getLabels("symptomaticAttackRate", compartmentLength)
    } else { #Number
      infectionsByGroup <- attackRateByGroup * model$parameters$population
      names(infectionsByGroup) <- getLabels("symptomaticCases", compartmentLength)
    } 
  } else {
    attackRateByGroup <- getCombinedValue(model, compartments, index = nrow(model$rawOutput)) -
                         getCombinedValue(model, compartments, index = 1)
    if (asRate) {
      infectionsByGroup <- attackRateByGroup
      names(infectionsByGroup) <- getLabels("serologicAttackRate", compartmentLength)
    } else { #Number
      infectionsByGroup <- attackRateByGroup * model$parameters$population
      names(infectionsByGroup) <- getLabels("serologicCases", compartmentLength)
    } 
  }
  if (byGroup) {
    return(infectionsByGroup)
  } else {
    return(sum(infectionsByGroup))
  }
}

#SEIRV
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @return A vector of infections 
#' @method getInfections SEIRVModel
#' @keywords internal
#' @export
getInfections.SEIRVModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE, fractionSymptomatic) {
  compartmentLength <- length(model$parameters$populationFractions)
  compartments <- getCompartments(model, type = "R")
  if (symptomatic) {
    if (missing(fractionSymptomatic)) {
      stop("fractionSymptomatic must be specified.")
    }
    
    # Those who get VEp
    attackRateByGroup <- (getCombinedValue(model, "Rv", index = nrow(model$rawOutput)) -
                            getCombinedValue(model, "Rv", index = 1)) * fractionSymptomatic * (1 - model$parameters$VEp)
    # And those who don't
    attackRateByGroup <- attackRateByGroup + (getCombinedValue(model, "R", index = nrow(model$rawOutput)) -
                            getCombinedValue(model, "R", index = 1)) * fractionSymptomatic
    if (asRate) {
      infectionsByGroup <- attackRateByGroup
      names(infectionsByGroup) <- getLabels("symptomaticAttackRate", compartmentLength)
    } else { #Number
      infectionsByGroup <- attackRateByGroup * model$parameters$population
      names(infectionsByGroup) <- getLabels("symptomaticCases", compartmentLength)
    } 
  } else {
    attackRateByGroup <- getCombinedValue(model, compartments, index = nrow(model$rawOutput)) -
      getCombinedValue(model, compartments, index = 1)
    if (asRate) {
      infectionsByGroup <- attackRateByGroup
      names(infectionsByGroup) <- getLabels("serologicAttackRate", compartmentLength)
    } else { #Number
      infectionsByGroup <- attackRateByGroup * model$parameters$population
      names(infectionsByGroup) <- getLabels("serologicCases", compartmentLength)
    } 
  }
  if (byGroup) {
    return(infectionsByGroup)
  } else {
    return(sum(infectionsByGroup))
  }
}

#SEIRV2Dose
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @return A vector of infections 
#' @method getInfections SEIRV2DoseModel
#' @keywords internal
#' @export
getInfections.SEIRV2DoseModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE, fractionSymptomatic) {
  compartmentLength <- length(model$parameters$populationFractions)
  compartments <- getCompartments(model, type = "R")
  if (symptomatic) {
    if (missing(fractionSymptomatic)) {
      stop("fractionSymptomatic must be specified.")
    }
    
    # Those who get VEp
    attackRateByGroup <- (getCombinedValue(model, c("Rv"), index = nrow(model$rawOutput)) -
                            getCombinedValue(model, c("Rv"), index = 1)) * fractionSymptomatic * (1 - model$parameters$VEp1) +
      (getCombinedValue(model, c("Rvb"), index = nrow(model$rawOutput)) -
         getCombinedValue(model, c("Rvb"), index = 1)) * fractionSymptomatic * (1 - model$parameters$VEp2)
    
    # And those who don't
    attackRateByGroup <- attackRateByGroup + (getCombinedValue(model, "R", index = nrow(model$rawOutput)) -
                                                getCombinedValue(model, "R", index = 1)) * fractionSymptomatic
    if (asRate) {
      infectionsByGroup <- attackRateByGroup
      names(infectionsByGroup) <- getLabels("symptomaticAttackRate", compartmentLength)
    } else { #Number
      infectionsByGroup <- attackRateByGroup * model$parameters$population
      names(infectionsByGroup) <- getLabels("symptomaticCases", compartmentLength)
    } 
  } else {
    attackRateByGroup <- getCombinedValue(model, compartments, index = nrow(model$rawOutput)) -
      getCombinedValue(model, compartments, index = 1)
    if (asRate) {
      infectionsByGroup <- attackRateByGroup
      names(infectionsByGroup) <- getLabels("serologicAttackRate", compartmentLength)
    } else { #Number
      infectionsByGroup <- attackRateByGroup * model$parameters$population
      names(infectionsByGroup) <- getLabels("serologicCases", compartmentLength)
    } 
  }
  if (byGroup) {
    return(infectionsByGroup)
  } else {
    return(sum(infectionsByGroup))
  }
}

#SEIRT
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @return A vector of infections
#' @method getInfections SEIRTModel
#' @keywords internal
#' @export
getInfections.SEIRTModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE) {
  getInfections.SEIRModel(model, byGroup, asRate, symptomatic, fractionSymptomatic = model$parameters$fractionSymptomatic)
}

#SEIRTV
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @return A vector of infections
#' @method getInfections SEIRTVModel
#' @keywords internal
#' @export
getInfections.SEIRTVModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE) {
  getInfections.SEIRVModel(model, byGroup, asRate, symptomatic, fractionSymptomatic = model$parameters$fractionSymptomatic)
}

#SEAIRTV
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @return A vector of infections
#' @method getInfections SEAIRTVModel
#' @keywords internal
#' @export
getInfections.SEAIRTVModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE) {
  getInfections.SEIRVModel(model, byGroup, asRate, symptomatic, fractionSymptomatic = model$parameters$fractionSymptomatic)
}

#SEIRTV2Dose
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @return A vector of infections
#' @method getInfections SEIRTV2DoseModel
#' @keywords internal
#' @export
getInfections.SEIRTV2DoseModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE) {
  return(getInfections.SEIRV2DoseModel(model, byGroup, asRate, symptomatic, fractionSymptomatic = model$parameters$fractionSymptomatic))
}

#SEIRVMono
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @return A vector of infections
#' @method getInfections SEIRVMonoModel
#' @keywords internal
#' @export
getInfections.SEIRVMonoModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE, fractionSymptomatic) {
  compartmentLength <- length(model$parameters$populationFractions)
  compartments <- getCompartments(model, type = "R")
  if (symptomatic) {
    if (missing(fractionSymptomatic)) {
      stop("fractionSymptomatic must be specified.")
    }
    
    # Those who get VEp
    attackRateByGroup <- (getCombinedValue(model, c("Rv"), index = nrow(model$rawOutput)) -
                            getCombinedValue(model, c("Rv"), index = 1)) * fractionSymptomatic * (1 - model$parameters$VEp1) +
      (getCombinedValue(model, c("RvM"), index = nrow(model$rawOutput)) -
         getCombinedValue(model, c("RvM"), index = 1)) * fractionSymptomatic * (1 - model$parameters$VEpM)
    
    # And those who don't
    attackRateByGroup <- attackRateByGroup + (getCombinedValue(model, "R", index = nrow(model$rawOutput)) -
                                                getCombinedValue(model, "R", index = 1)) * fractionSymptomatic
    if (asRate) {
      infectionsByGroup <- attackRateByGroup
      names(infectionsByGroup) <- getLabels("symptomaticAttackRate", compartmentLength)
    } else { #Number
      infectionsByGroup <- attackRateByGroup * model$parameters$population
      names(infectionsByGroup) <- getLabels("symptomaticCases", compartmentLength)
    } 
  } else {
    attackRateByGroup <- getCombinedValue(model, compartments, index = nrow(model$rawOutput)) -
      getCombinedValue(model, compartments, index = 1)
    if (asRate) {
      infectionsByGroup <- attackRateByGroup
      names(infectionsByGroup) <- getLabels("serologicAttackRate", compartmentLength)
    } else { #Number
      infectionsByGroup <- attackRateByGroup * model$parameters$population
      names(infectionsByGroup) <- getLabels("serologicCases", compartmentLength)
    } 
  }
  if (byGroup) {
    return(infectionsByGroup)
  } else {
    return(sum(infectionsByGroup))
  }
}

#SEIRVPrimeBoost
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @return A vector of infections
#' @method getInfections SEIRVPrimeBoostModel
#' @keywords internal
#' @export
getInfections.SEIRVPrimeBoostModel <- getInfections.SEIRV2DoseModel

#SEIRTVPrimeBoost
#' @title Get infections
#' @description Gets the total infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @return A vector of infections
#' @method getInfections SEIRTVPrimeBoostModel
#' @keywords internal
#' @export
getInfections.SEIRTVPrimeBoostModel <- function(model, byGroup = TRUE, asRate = FALSE, symptomatic = FALSE) {
  return(getInfections.SEIRV2DoseModel(model, byGroup, asRate, symptomatic, fractionSymptomatic = model$parameters$fractionSymptomatic))
}
