#' @title SEIR+V Model
#' @description A model for influenza that uses a SEIR framework with age 
#'   structure and that allows for vaccination
#' @param population Size of population; defaults to 1
#' @param populationFractions Vector of population fractions (all non-negative, 
#'   sum to 1); defaults to 1, representing a single population group
#' @param contactMatrix A matrix whose (row i, column j) entry denotes the 
#'   number of potentially infectious contacts a single individual from group j 
#'   has with individuals from group i each day; defaults to proportional mixing
#' @param R0 Average number of secondary cases from a single infected individual
#'   in a completely susceptible population; must be specified
#' @param latentPeriod Latent period in days; must be specified
#' @param infectiousPeriod Infectious period in days; must be specified
#' @param seedInfections Fraction of the population to seed with infections; 
#'   single fraction or vector of fractions by population group; defaults to 0
#' @param priorImmunity Fraction of the population with prior immunity; single 
#'   fraction, or vector of fractions by population group; defaults to 0
#' @param useCommunityMitigation Whether or not to use community mitigation
#'   implemented by modulation of the contact matrix; defaults to FALSE
#' @param communityMitigationStartDay If using community mitigation, day of the
#'   simulation on which to start mitigation; must be specified if applicable
#' @param communityMitigationDuration If using community mitigation, duration of
#'   time during which mitigation is in effect; must be specified if applicable
#' @param communityMitigationMultiplier If using community mitigation, the
#'   non-negative matrix of multipliers that will be used to modulate the contact
#'   matrix by elementwise multiplication; must be specified if applicable
#' @param vaccineAdministrationRatePerDay Vaccine administration rate each day; 
#'   defaults to 0
#' @param vaccineAvailabilityByDay Vector that contains the amount of vaccine 
#'   available each day; defaults to 0
#' @param vaccineUptakeMultiplier Vector of multipliers that determines the
#'   relative rate at which vaccine is given to each age group; defaults to
#'   vaccine being allotted proportionally by population
#' @param VEs Vaccine efficacy: protection for vaccinated susceptible 
#'   individuals; single fraction or vector of fractions by population group;
#'   defaults to 0
#' @param VEi Vaccine efficacy: prevention of transmission from vaccinated 
#'   infected individuals; single fraction or vector of fractions by population
#'   group; defaults to 0
#' @param VEp Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals; single fraction or vector of fractions by population
#'   group; defaults to 0
#' @param vaccineEfficacyDelay Delay in days between administration of dose and 
#'   onset of protection; defaults to 7
#' @param simulationLength Number of days to simulate after seeding infections; 
#'   defaults to 240
#' @param seedStartDay Day on which to seed initial infections; defaults to 0
#' @param tolerance Absolute tolerance for numerical integration; defaults to 
#'   1e-8
#' @param method Which integration method to use. Defaults to lsoda
#' @return a SEIRVModel object
#' @export
SEIRVModel2 <- function(population, populationFractions, contactMatrix, R0,
                       latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                       useCommunityMitigation, communityMitigationStartDay,
                       communityMitigationDuration, communityMitigationMultiplier,
                       vaccineAdministrationRatePerDay, vaccineAvailabilityByDay,
                       vaccineUptakeMultiplier, VEs, VEi, VEp, vaccineEfficacyDelay,
                       simulationLength, seedStartDay, tolerance, method) {
  #Check inputs #TODO: Add checks for all inputs
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIRV", argumentList) #Get parameters from checked inputs

  initialState <- with(parameters, {
    c(S  = (1 - priorImmunity) * populationFractions,
      E  = 0 * populationFractions,
      I  = 0 * populationFractions,
      R  = priorImmunity * populationFractions,
      Sv = 0 * populationFractions,
      Ev = 0 * populationFractions,
      Iv = 0 * populationFractions,
      Rv = 0 * populationFractions,
      V  = 0 * populationFractions,
      vaccinating = rep(1, length(populationFractions)))
  })

  rootFunction <- function(t, state, parameters) {
    stateList <- reconstructState.SEIRV2(state)
    with(append(stateList, parameters), {

      return(ifelse(vaccinating > 0, populationFractions - V - tolerance, 1))
      })
  }
  eventFunction <- function(t, state, parameters) {
    stateList <- reconstructState.SEIRV2(state)
    with(append(stateList, parameters), {
      state[getLabels("vaccinating", length(populationFractions))] <- 
          ifelse(populationFractions - V > tolerance, 1, 0)
      return(state)
      })
  }
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIRV2,
                              seedFunction = doSeed.SEIRV2,
                              rootFunction = rootFunction,
                              eventFunction = eventFunction,
                              method = method)
  
  #Build the SEIRVModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- c("SEIRVModel", "SEIRModel")
  return(model)
}

#' @title Check SEIR+V inputs
#' @description Checks the input parameters for the SEIR+V model
#' @return List of parameters for the SEIR+V model
#' @keywords internal
checkInputs.SEIRV <- function(population, populationFractions, contactMatrix, R0,
                              latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                              useCommunityMitigation, communityMitigationStartDay,
                              communityMitigationDuration, communityMitigationMultiplier, 
                              vaccineAdministrationRatePerDay, vaccineAvailabilityByDay,
                              vaccineUptakeMultiplier, VEs, VEi, VEp, vaccineEfficacyDelay,
                              simulationLength, seedStartDay, tolerance, method) {
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  SEIRParameters <- do.call("checkInputs.SEIR", argumentList)
  #Update arguments passed to checkInputs.Vaccine using SEIRParameters
  argumentList$population <- SEIRParameters$population
  argumentList$populationFractions <- SEIRParameters$populationFractions
  argumentList$seedStartDay <- SEIRParameters$seedStartDay
  argumentList$simulationLength <- SEIRParameters$simulationLength
  vaccineParameters <- do.call("checkInputs.Vaccine", argumentList)
  #Return the parameters
  return(c(SEIRParameters, vaccineParameters))
}

#' @title Check vaccine inputs
#' @description Checks vaccine inputs and computes the vaccine parameters
#' @return List of vaccine parameters
#' @keywords internal
checkInputs.Vaccine <- function(population, populationFractions, seedStartDay, simulationLength,
                                vaccineAdministrationRatePerDay = 0, vaccineAvailabilityByDay = 0,
                                vaccineUptakeMultiplier = 1, VEs = 0, VEi = 0, VEp = 0, 
                                vaccineEfficacyDelay = 7, ...) {
  #Validate vaccine parameters
  #vaccineAdministrationRatePerDay
  checkNonNegativeNumber(vaccineAdministrationRatePerDay)
  #vaccineAvailabilityByDay
  checkNonNegative(vaccineAvailabilityByDay)
  #vaccineUptakeMultiplier
  checkNonNegative(vaccineUptakeMultiplier)
  checkDimensionsMatch(vaccineUptakeMultiplier, populationFractions)
  #VEs
  checkBetween0and1(VEs)
  checkDimensionsMatch(VEs, populationFractions)
  #VEi
  checkBetween0and1(VEi)
  checkDimensionsMatch(VEi, populationFractions)
  #VEp
  checkBetween0and1(VEp)
  checkDimensionsMatch(VEp, populationFractions)
  #vaccineEfficacyDelay
  checkNonNegativeNumber(vaccineEfficacyDelay)


  #Compute the daily vaccination rate
  totalSimulationLength <- seedStartDay + simulationLength
  vaccinationRateByDay <- rep(0, totalSimulationLength)
  currentVaccineAvailability <- 0
  for (i in 1:totalSimulationLength) {
    if (i <= length(vaccineAvailabilityByDay)){
      currentVaccineAvailability <- currentVaccineAvailability + vaccineAvailabilityByDay[i]
    }
    vaccinationRateByDay[i] <- min(vaccineAdministrationRatePerDay, currentVaccineAvailability)
    currentVaccineAvailability <- currentVaccineAvailability - vaccinationRateByDay[i]
  }
  vaccinationRateByDay <- vaccinationRateByDay / population #Normalize
  
  #Define vaccination rate function
  vaccinationRate <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(vaccinationRateByDay[floor(t - vaccineEfficacyDelay + 1)])
    }
  }

  #Compute the vaccination rate age multiplier
  vaccinationRateAgeMultiplier <- vaccineUptakeMultiplier * populationFractions
  totalMultiplier <- sum(vaccinationRateAgeMultiplier)
  if (totalMultiplier > 0) {
    vaccinationRateAgeMultiplier <- vaccinationRateAgeMultiplier / totalMultiplier
  } else {
    warning("vaccineUptakeMultiplier prevents vaccination from occurring.", call. = FALSE)
  }

  #Return the parameters
  return(list(vaccinationRate = vaccinationRate, vaccinationRateAgeMultiplier = vaccinationRateAgeMultiplier,
              VEs = VEs, VEi = VEi, VEp = VEp, vaccineEfficacyDelay = vaccineEfficacyDelay))
}

#This is a utility function that reconstructs the model state as a list so that equations can refer to compartments by name
reconstructState.SEIRV2 <- function(state) {
  numberOfClasses <- length(state) / 10 #Each of the 10 classes are vectors of the same length
  S  <- state[                       1 :     numberOfClasses ]
  E  <- state[    (numberOfClasses + 1):(2 * numberOfClasses)]
  I  <- state[(2 * numberOfClasses + 1):(3 * numberOfClasses)]
  R  <- state[(3 * numberOfClasses + 1):(4 * numberOfClasses)]
  Sv <- state[(4 * numberOfClasses + 1):(5 * numberOfClasses)]
  Ev <- state[(5 * numberOfClasses + 1):(6 * numberOfClasses)]
  Iv <- state[(6 * numberOfClasses + 1):(7 * numberOfClasses)]
  Rv <- state[(7 * numberOfClasses + 1):(8 * numberOfClasses)]
  V  <- state[(8 * numberOfClasses + 1):(9 * numberOfClasses)]
  vaccinating <- state[(9 * numberOfClasses + 1):(10 * numberOfClasses)]

  return(as.list(environment()))
}

#This function implements the multivariate derivative of the SEIR+V model
#parameters should define populationFractions, contactMatrix, beta, lambda, gamma,
#VEs, VEi, and the function vaccinationRate(t)
#Note that the total population is normalized to be 1
getDerivative.SEIRV2 <- function(t, state, parameters) {
  stateList <- reconstructState.SEIRV2(state)
  with(append(stateList, parameters), {
    if (useCommunityMitigation) {
      if ((t >= communityMitigationStartDay) && (t < communityMitigationEndDay)) {
        contactMatrix <- communityMitigationMultiplier * contactMatrix
      } 
    }

    isVaccinatingByAge <- (vaccinating > 0) & (populationFractions - V > 0)
    effectiveVaccinationMultiplier <- sum(ifelse(isVaccinatingByAge, 1, 0) * vaccinationRateAgeMultiplier)
    if (effectiveVaccinationMultiplier > 0) {
      vaccinationRateByAge <- vaccinationRate(t) * vaccinationRateAgeMultiplier /
        effectiveVaccinationMultiplier
    } else {
      vaccinationRateByAge <- 0
    }
   
    #Flows
    forceOfInfection <- beta / populationFractions * (contactMatrix %*% (I + ((1 - VEi) * Iv)))

    S_to_E <- S * forceOfInfection
    E_to_I <- lambda * E
    I_to_R <- gamma * I
    
    Sv_to_Ev <- Sv * (1 - VEs) * forceOfInfection
    Ev_to_Iv <- lambda * Ev
    Iv_to_Rv <- gamma * Iv
    
    S_to_Sv <-  ifelse(isVaccinatingByAge, vaccinationRateByAge * S / (populationFractions - V), 0)
    
    #Derivatives
    #Non-vaccinated compartments
    dS <- -S_to_E - S_to_Sv
    dE <- S_to_E - E_to_I
    dI <- E_to_I - I_to_R
    dR <- I_to_R
    #Vaccinated compartments
    dSv <- -Sv_to_Ev + S_to_Sv
    dEv <- Sv_to_Ev - Ev_to_Iv
    dIv <- Ev_to_Iv - Iv_to_Rv
    dRv <- Iv_to_Rv
    #Auxiliary vaccinated compartment
    dV <- ifelse(isVaccinatingByAge, vaccinationRateByAge, 0)

    #Return derivative
    return(list(c(dS, dE, dI, dR, dSv, dEv, dIv, dRv, dV, 0 * populationFractions)))
  })
}

#This function implements seeding infections in the SEIR+V model
#parameters should define seedInfections, lambda, and gamma
doSeed.SEIRV2 <- function(state, parameters) {
  stateList <- reconstructState.SEIRV2(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions
    E <- E + seedInfectionsFractions / (1 + lambda / gamma)
    I <- I + seedInfectionsFractions / (1 + gamma / lambda)
    #Return derivative
    return(c(S, E, I, R, Sv, Ev, Iv, Rv, V, vaccinating))
  })
}