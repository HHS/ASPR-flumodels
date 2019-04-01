#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @details This is the generic function, for class-specific versions see:
#' \itemize{
#'  \item \code{\link{getInfectionTimeSeries.SEIRModel}} for SEIR-type models
#'  \item \code{\link{getInfectionTimeSeries.SEIRVModel}} for SEIRV-type models
#'  \item \code{\link{getInfectionTimeSeries.SEIRV2DoseModel}} for SEIRV2Dose-type models
#'  \item \code{\link{getInfectionTimeSeries.SEIRTModel}} for SEIRT-type models
#'  \item \code{\link{getInfectionTimeSeries.SEIRTVModel}} for SEIRTV-type models
#'  \item \code{\link{getInfectionTimeSeries.SEIRTV2DoseModel}} for SEIRTV2Dose-type models
#' }
#' @param model A flumodels model object
#' @param ... Other parameters for the class-specific functions
#' @export
getInfectionTimeSeries <- function(model, ...) {
  UseMethod("getInfectionTimeSeries", model)
}

#SEIR
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @method getInfectionTimeSeries SEIRModel
#' @keywords internal
#' @export
getInfectionTimeSeries.SEIRModel <- function(model, byGroup = TRUE, asRate = FALSE, incidence = FALSE,
                                             symptomatic = FALSE, fractionSymptomatic, byWeek = FALSE) {
  compartmentLength = length(model$parameters$populationFractions)
  if (incidence) { #Incidence
    if (symptomatic) {
      if (missing(fractionSymptomatic)) {
        stop("fractionSymptomatic must be specified.")
      }
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append(getCompartments(model, type = "I"),
                                                                        getCompartments(model, type = "R")))
      dailyInfectionsByGroup <- rbind(rep(0, compartmentLength),
                                      cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                        cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic
    } else {
      dailySusceptiblesByGroup <- getCombinedTimeSeries(model, getCompartments(model, type = "S"))
      dailyInfectionsByGroup <- rbind(rep(0, compartmentLength),
                                      dailySusceptiblesByGroup[1:(nrow(dailySusceptiblesByGroup) - 1), , drop = FALSE] -
                                        dailySusceptiblesByGroup[2:nrow(dailySusceptiblesByGroup), , drop = FALSE])
    }
    dailyInfectionsByGroup[model$parameters$seedStartDay + 1, ] <- NA #Ignore the seeded infections
    colnames(dailyInfectionsByGroup) <- getLabels("incidence", compartmentLength)
  } else { #Prevalence
    if (symptomatic) {
      if (missing(fractionSymptomatic)) {
        stop("fractionSymptomatic must be specified.")
      }
      dailyInfectionsByGroup <- getCombinedTimeSeries(model, getCompartments(model, type = "I")) * fractionSymptomatic
    } else {
      dailyInfectionsByGroup <- getCombinedTimeSeries(model, append(getCompartments(model, type = "E"),
                                                                    getCompartments(model, type = "I")))
    }
    colnames(dailyInfectionsByGroup) <- getLabels("prevalence", compartmentLength)
  }
  makeOutputFromDailyInfectionsByGroup(model, dailyInfectionsByGroup, byGroup, asRate, byWeek)
}

#SEIRV
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @method getInfectionTimeSeries SEIRVModel
#' @keywords internal
#' @export
getInfectionTimeSeries.SEIRVModel <- function(model, byGroup = TRUE, asRate = FALSE, incidence = FALSE,
                                              symptomatic = FALSE, fractionSymptomatic, byWeek = FALSE) {
  compartmentLength = length(model$parameters$populationFractions)
  if (incidence) { #Incidence
    if (symptomatic) {
      if (missing(fractionSymptomatic)) {
        stop("fractionSymptomatic must be specified.")
      }
      
      # Those who get VEp
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append("Iv", "Rv"))
      dailyInfectionsByGroup <- rbind(rep(0, compartmentLength),
                                      cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                        cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic * (1 - model$parameters$VEp)
      # Those who do not get VEp
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append("I", "R"))
      dailyInfectionsByGroup <- dailyInfectionsByGroup + rbind(rep(0, compartmentLength),
                                                               cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                                                 cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic
    } else {
      dailySusceptiblesByGroup <- getCombinedTimeSeries(model, getCompartments(model, type = "S"))
      dailyInfectionsByGroup <- rbind(rep(0, compartmentLength),
                                      dailySusceptiblesByGroup[1:(nrow(dailySusceptiblesByGroup) - 1), , drop = FALSE] -
                                        dailySusceptiblesByGroup[2:nrow(dailySusceptiblesByGroup), , drop = FALSE])
    }
    dailyInfectionsByGroup[model$parameters$seedStartDay + 1, ] <- NA #Ignore the seeded infections
    colnames(dailyInfectionsByGroup) <- getLabels("incidence", compartmentLength)
  } else { #Prevalence
    if (symptomatic) {
      if (missing(fractionSymptomatic)) {
        stop("fractionSymptomatic must be specified.")
      }
      # Those who get VEp
      dailyInfectionsByGroup <- getCombinedTimeSeries(model, "Iv") * fractionSymptomatic * (1 - model$parameters$VEp)
      # Those who don't
      dailyInfectionsByGroup <- dailyInfectionsByGroup + getCombinedTimeSeries(model, "I") * fractionSymptomatic
    } else {
      dailyInfectionsByGroup <- getCombinedTimeSeries(model, append(getCompartments(model, type = "E"),
                                                                    getCompartments(model, type = "I")))
    }
    colnames(dailyInfectionsByGroup) <- getLabels("prevalence", compartmentLength)
  }
  
  makeOutputFromDailyInfectionsByGroup(model, dailyInfectionsByGroup, byGroup, asRate, byWeek)
}

#SEIRV2Dose
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @method getInfectionTimeSeries SEIRV2DoseModel
#' @keywords internal
#' @export
getInfectionTimeSeries.SEIRV2DoseModel <- function(model, byGroup = TRUE, asRate = FALSE, incidence = FALSE,
                                                   symptomatic = FALSE, fractionSymptomatic, byWeek = FALSE) {
  compartmentLength = length(model$parameters$populationFractions)
  if (incidence) { #Incidence
    if (symptomatic) {
      if (missing(fractionSymptomatic)) {
        stop("fractionSymptomatic must be specified.")
      }
      
      # Those who get VEp1
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append("Iv", "Rv"))
      dailyInfectionsByGroup <- rbind(rep(0, compartmentLength),
                                      cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                        cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic * (1 - model$parameters$VEp1)
      # Those who get VEp2
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append("Ivb", "Rvb"))
      dailyInfectionsByGroup <- dailyInfectionsByGroup + rbind(rep(0, compartmentLength),
                                                               cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                                                 cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic * (1 - model$parameters$VEp2)
      # Those who do not get VEp
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append("I", "R"))
      dailyInfectionsByGroup <- dailyInfectionsByGroup + rbind(rep(0, compartmentLength),
                                                               cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                                                 cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic
    } else {
      dailySusceptiblesByGroup <- getCombinedTimeSeries(model, getCompartments(model, type = "S"))
      dailyInfectionsByGroup <- rbind(rep(0, compartmentLength),
                                      dailySusceptiblesByGroup[1:(nrow(dailySusceptiblesByGroup) - 1), , drop = FALSE] -
                                        dailySusceptiblesByGroup[2:nrow(dailySusceptiblesByGroup), , drop = FALSE])
    }
    dailyInfectionsByGroup[model$parameters$seedStartDay + 1, ] <- NA #Ignore the seeded infections
    colnames(dailyInfectionsByGroup) <- getLabels("incidence", compartmentLength)
  } else { #Prevalence
    if (symptomatic) {
      if (missing(fractionSymptomatic)) {
        stop("fractionSymptomatic must be specified.")
      }
      # Those who get VEp1
      dailyInfectionsByGroup <- getCombinedTimeSeries(model, "Iv") * fractionSymptomatic * (1 - model$parameters$VEp1)
      # Those who get VEp2
      dailyInfectionsByGroup <- dailyInfectionsByGroup + getCombinedTimeSeries(model, "Ivb") *
        fractionSymptomatic * (1 - model$parameters$VEp2)
      # Those who don't
      dailyInfectionsByGroup <- dailyInfectionsByGroup + getCombinedTimeSeries(model, "I") * fractionSymptomatic
    } else {
      dailyInfectionsByGroup <- getCombinedTimeSeries(model, append(getCompartments(model, type = "E"),
                                                                    getCompartments(model, type = "I")))
    }
    colnames(dailyInfectionsByGroup) <- getLabels("prevalence", compartmentLength)
  }
  
  makeOutputFromDailyInfectionsByGroup(model, dailyInfectionsByGroup, byGroup, asRate, byWeek)
}

#SEIRT
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @method getInfectionTimeSeries SEIRTModel
#' @keywords internal
#' @export
getInfectionTimeSeries.SEIRTModel <- function(model, byGroup = TRUE, asRate = FALSE, incidence = FALSE,
                                              symptomatic = FALSE, byWeek = FALSE) {
  return(getInfectionTimeSeries.SEIRModel(model, byGroup, asRate, incidence, symptomatic,
                                          fractionSymptomatic = model$parameters$fractionSymptomatic,
                                          byWeek = FALSE))
}

#SEIRTV
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @method getInfectionTimeSeries SEIRTVModel
#' @keywords internal
#' @export
getInfectionTimeSeries.SEIRTVModel <- function(model, byGroup = TRUE, asRate = FALSE, incidence = FALSE,
                                               symptomatic = FALSE, byWeek = FALSE) {
  return(getInfectionTimeSeries.SEIRVModel(model, byGroup, asRate, incidence, symptomatic,
                                           fractionSymptomatic = model$parameters$fractionSymptomatic,
                                           byWeek = FALSE))
}

#SEIRTV2Dose
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @method getInfectionTimeSeries SEIRTV2DoseModel
#' @keywords internal
#' @export
getInfectionTimeSeries.SEIRTV2DoseModel <- function(model, byGroup = TRUE, asRate = FALSE, incidence = FALSE,
                                                    symptomatic = FALSE, byWeek = FALSE) {
  return(getInfectionTimeSeries.SEIRV2DoseModel(model, byGroup, asRate, incidence, symptomatic,
                                                fractionSymptomatic = model$parameters$fractionSymptomatic,
                                                byWeek = FALSE))
}


#SEIRVMono
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @method getInfectionTimeSeries SEIRVMonoModel
#' @keywords internal
#' @export
getInfectionTimeSeries.SEIRVMonoModel <- function(model, byGroup = TRUE, asRate = FALSE, incidence = FALSE,
                                                  symptomatic = FALSE, fractionSymptomatic, byWeek = FALSE) {
  compartmentLength = length(model$parameters$populationFractions)
  if (incidence) { #Incidence
    if (symptomatic) {
      if (missing(fractionSymptomatic)) {
        stop("fractionSymptomatic must be specified.")
      }
      
      # Those who get VEp1
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append("Iv", "Rv"))
      dailyInfectionsByGroup <- rbind(rep(0, compartmentLength),
                                      cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                        cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic * (1 - model$parameters$VEp1)
      # Those who get VEpM
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append("IvM", "RvM"))
      dailyInfectionsByGroup <- dailyInfectionsByGroup + rbind(rep(0, compartmentLength),
                                                               cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                                                 cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic * (1 - model$parameters$VEpM)
      # Those who do not get VEp
      cumulativeIncidenceByGroup <- getCombinedTimeSeries(model, append("I", "R"))
      dailyInfectionsByGroup <- dailyInfectionsByGroup + rbind(rep(0, compartmentLength),
                                                               cumulativeIncidenceByGroup[2:nrow(cumulativeIncidenceByGroup), , drop = FALSE] -
                                                                 cumulativeIncidenceByGroup[1:(nrow(cumulativeIncidenceByGroup) - 1), , drop = FALSE]) *
        fractionSymptomatic
    } else {
      dailySusceptiblesByGroup <- getCombinedTimeSeries(model, getCompartments(model, type = "S"))
      dailyInfectionsByGroup <- rbind(rep(0, compartmentLength),
                                      dailySusceptiblesByGroup[1:(nrow(dailySusceptiblesByGroup) - 1), , drop = FALSE] -
                                        dailySusceptiblesByGroup[2:nrow(dailySusceptiblesByGroup), , drop = FALSE])
    }
    dailyInfectionsByGroup[model$parameters$seedStartDay + 1, ] <- NA #Ignore the seeded infections
    colnames(dailyInfectionsByGroup) <- getLabels("incidence", compartmentLength)
  } else { #Prevalence
    if (symptomatic) {
      if (missing(fractionSymptomatic)) {
        stop("fractionSymptomatic must be specified.")
      }
      # Those who get VEp1
      dailyInfectionsByGroup <- getCombinedTimeSeries(model, "Iv") * fractionSymptomatic * (1 - model$parameters$VEp1)
      # Those who get VEpM
      dailyInfectionsByGroup <- dailyInfectionsByGroup + getCombinedTimeSeries(model, "IvM") *
        fractionSymptomatic * (1 - model$parameters$VEpM)
      # Those who don't
      dailyInfectionsByGroup <- dailyInfectionsByGroup + getCombinedTimeSeries(model, "I") * fractionSymptomatic
    } else {
      dailyInfectionsByGroup <- getCombinedTimeSeries(model, append(getCompartments(model, type = "E"),
                                                                    getCompartments(model, type = "I")))
    }
    colnames(dailyInfectionsByGroup) <- getLabels("prevalence", compartmentLength)
  }
  
  makeOutputFromDailyInfectionsByGroup(model, dailyInfectionsByGroup, byGroup, asRate, byWeek)
}

#SEIRVPrimeBoost
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param incidence If true, returns infection incidence, otherwise returns infection prevalence; defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @method getInfectionTimeSeries SEIRVPrimeBoostModel
#' @keywords internal
#' @export
getInfectionTimeSeries.SEIRVPrimeBoostModel <- getInfectionTimeSeries.SEIRV2DoseModel

#Convenience function to produce the final step of output
#' @title Get infection time series
#' @description Gets the time series of infections from the given model
#' @param model The model from which to get the data
#' @param dailyInfectionsByGroup Calculated daily infection prevalence (either as incidence or prevalence)
#' @param byGroup Whether or not to return data by population group; defaults to TRUE
#' @param asRate Whether to return results as a rate (fraction of population) or else a number; defaults to FALSE
#' @param byWeek If true, returns the output on a weekly, not daily basis; defaults to FALSE
#' @return A matrix that contains the infections by simulation day (or week, if byWeek is selected)
#' @keywords internal
#' @author Matt Clay clay.matt@gmail.com
makeOutputFromDailyInfectionsByGroup <- function(model, dailyInfectionsByGroup, byGroup, asRate, byWeek) {
  if (!asRate) {
    dailyInfectionsByGroup <- dailyInfectionsByGroup * model$parameters$population
  }
  
  if (byWeek) {
    # Shrink dailyInfectionsByGroup by the fake first 0
    # Sum by week using rowsum. While not as elegant as some other options, it is robust for single-column matrices & for weeks
    # that only have a single day in them.
    dailyInfectionsByGroup <- dailyInfectionsByGroup[2:nrow(dailyInfectionsByGroup),]
    
    if (class(dailyInfectionsByGroup) == "matrix") {
      daily.length <- nrow(dailyInfectionsByGroup)
    } else {
      daily.length <- length(dailyInfectionsByGroup)
    }

    group <- rep(1:ceiling(daily.length / 7), each = 7)
    group <- group[1:daily.length]
    
    weeklyInfectionsByGroup <- rowsum(dailyInfectionsByGroup, group, na.rm = TRUE)
    
    if (byGroup) {
      return(weeklyInfectionsByGroup)
    } else {
      return(rowSums(weeklyInfectionsByGroup))
    }
  } else {
    if (byGroup) {
      return(dailyInfectionsByGroup)
    } else {
      return(rowSums(dailyInfectionsByGroup))
    }
  }
}
