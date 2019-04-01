#' @title SEIR+V Prime-Boost Model
#' @description A model for influenza that uses a SEIR framework with age 
#'   structure and that allows for vaccination with a prime-boost vaccine
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
#'   single fraction, or vector of fractions by population group; defaults to 0
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
#' @param vaccineAvailabilityByDayPrime Vector that contains the amount of Prime vaccine 
#'   available each day; defaults to 0
#' @param vaccineAvailabilityByDayBoost Vector that contains the amount of boost vaccine 
#'   available each day; defaults to 0
#' @param vaccineUptakeMultiplier Vector of multipliers that determines the
#'   relative rate at which vaccine is given to each age group; defaults to
#'   vaccine being allotted proportionally by population
#' @param boostDelay Minimum days between prime and boost doses; defaults to 14
#' @param VEs1 Vaccine efficacy: protection for vaccinated susceptible 
#'   individuals after priming dose; single fraction or vector of fractions by
#'   population group; defaults to 0
#' @param VEs2 Vaccine efficacy: protection for vaccinated susceptible 
#'   individuals after priming and boost doses; single fraction or vector of fractions by
#'   population group; defaults to 0
#' @param VEi1 Vaccine efficacy: prevention of transmission from vaccinated 
#'   infected individuals after priming dose; single fraction or vector of fractions
#'   by population group; defaults to 0
#' @param VEi2 Vaccine efficacy: prevention of transmission from vaccinated 
#'   infected individuals after priming and boost doses; single fraction or vector of fractions
#'   by population group; defaults to 0
#' @param VEp1 Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals after priming dose; single fraction or vector of fractions 
#'   by population group; defaults to 0
#' @param VEp2 Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals after priming and boost doses; single fraction or vector of fractions 
#'   by population group; defaults to 0
#' @param vaccineEfficacyDelay Delay in days between administration of dose and 
#'   onset of protection; defaults to 7
#' @param simulationLength Number of days to simulate after seeding infections; 
#'   defaults to 240
#' @param seedStartDay Day on which to seed initial infections; defaults to 0
#' @param tolerance Absolute tolerance for numerical integration; defaults to 
#'   1e-8
#' @param method Which integration method to use. Defaults to lsoda
#' @return a SEIRVPrimeBoostModel object
#' @export
SEIRVPrimeBoostModel <- function(population, populationFractions, contactMatrix, R0,
                            latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                            useCommunityMitigation, communityMitigationStartDay,
                            communityMitigationDuration, communityMitigationMultiplier,
                            vaccineAdministrationRatePerDay, vaccineAvailabilityByDayPrime,
                            vaccineAvailabilityByDayBoost, vaccineUptakeMultiplier, boostDelay,
                            VEs1, VEs2, VEi1, VEi2, VEp1, VEp2, vaccineEfficacyDelay,
                            simulationLength, seedStartDay, tolerance, method) {
  #Check inputs #TODO: Add checks for all inputs
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIRVPrimeBoost", argumentList) #Get parameters from checked inputs
  
  initialState <- with(parameters, {
    c(S   = (1 - priorImmunity) * populationFractions,
      E   = 0 * populationFractions,
      I   = 0 * populationFractions,
      R   = priorImmunity * populationFractions,
      Sv  = 0 * populationFractions,
      Ev  = 0 * populationFractions,
      Iv  = 0 * populationFractions,
      Rv  = 0 * populationFractions,
      Svb = 0 * populationFractions,
      Evb = 0 * populationFractions,
      Ivb = 0 * populationFractions,
      Rvb = 0 * populationFractions,
      V   = 0 * populationFractions,
      Vb  = 0 * populationFractions,
      vaccinatingPrime = rep(1, length(populationFractions)),
      vaccinatingBoost = rep(0, length(populationFractions)))
  })
  
  rootFunction <- function(t, state, parameters) {
    stateList <- reconstructState.SEIRVPrimeBoost(state)
    with(append(stateList, parameters), {
      return(c(ifelse(vaccinatingPrime > 0, populationFractions - V - tolerance, 1),
               ifelse(vaccinatingBoost > 0, V - Vb - tolerance, 1),
               ifelse(V - Vb - tolerance > tolerance & vaccinatingBoost <= 0, 0, 1)))
      })
  }
  eventFunction <- function(t, state, parameters) {
    stateList <- reconstructState.SEIRVPrimeBoost(state)
    with(append(stateList, parameters), {
      state[getLabels("vaccinatingPrime", length(populationFractions))] <-
          ifelse(populationFractions - V > tolerance, 1, 0)
      state[getLabels("vaccinatingBoost", length(populationFractions))] <-
          ifelse(V - Vb > tolerance, 1, 0)
      return(state)
      })
  }

  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIRVPrimeBoost,
                              seedFunction = doSeed.SEIRVPrimeBoost,
                              rootFunction = rootFunction,
                              eventFunction = eventFunction,
                              method = method)
  
  #Build the SEIRVPrimeBoostModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- c("SEIRVPrimeBoostModel", "SEIRVModel", "SEIRModel")
  return(model)
}

#' @title Check SEIR+V prime-boost inputs
#' @description Checks the input parameters for the SEIR+V prime-boost model
#' @return List of parameters for the SEIR+V prime-boost model
#' @keywords internal
checkInputs.SEIRVPrimeBoost <- function(population, populationFractions, contactMatrix, R0,
                                        latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                                        useCommunityMitigation, communityMitigationStartDay,
                                        communityMitigationDuration, communityMitigationMultiplier,
                                        vaccineAdministrationRatePerDay, vaccineAvailabilityByDayPrime,
                                        vaccineAvailabilityByDayBoost, vaccineUptakeMultiplier, boostDelay,
                                        VEs1, VEs2, VEi1, VEi2, VEp1, VEp2, vaccineEfficacyDelay,
                                        simulationLength, seedStartDay, tolerance, method) {
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  SEIRParameters <- do.call("checkInputs.SEIR", argumentList)
  #Update call to checkInputs.Vaccine2Dose using SEIRParameters
  argumentList$population <- SEIRParameters$population
  argumentList$populationFractions <- SEIRParameters$populationFractions
  argumentList$seedStartDay <- SEIRParameters$seedStartDay
  argumentList$simulationLength <- SEIRParameters$simulationLength
  vaccinePrimeBoostParameters <- do.call("checkInputs.VaccinePrimeBoost", argumentList)
  #Return the parameters
  return(c(SEIRParameters, vaccinePrimeBoostParameters))
}

#' @title Check prime-boost vaccine inputs
#' @description Checks prime-boost vaccine inputs and computes the vaccine parameters
#' @return List of prime-boost vaccine parameters
#' @keywords internal
checkInputs.VaccinePrimeBoost <- function(population, populationFractions, seedStartDay, simulationLength,
                                          vaccineAdministrationRatePerDay = 0, vaccineAvailabilityByDayPrime = 0,
                                          vaccineAvailabilityByDayBoost = 0, vaccineUptakeMultiplier = 1,
                                          boostDelay = 14, VEs1 = 0, VEs2 = 0, VEi1 = 0, VEi2 = 0,
                                          VEp1 = 0, VEp2 = 0, vaccineEfficacyDelay = 7, ...) {
  #Validate vaccine parameters
  #vaccineAdministrationRatePerDay
  checkNonNegativeNumber(vaccineAdministrationRatePerDay)
  #vaccineAvailabilityByDayPrime
  checkNonNegative(vaccineAvailabilityByDayPrime)
  #vaccineAvailabilityByDayBoost
  checkNonNegative(vaccineAvailabilityByDayBoost)
  #vaccineUptakeMultiplier
  checkNonNegative(vaccineUptakeMultiplier)
  checkDimensionsMatch(vaccineUptakeMultiplier, populationFractions)
  #boostDelay
  checkNonNegativeNumber(boostDelay)
  #VEs1
  checkBetween0and1(VEs1)
  checkDimensionsMatch(VEs1, populationFractions)
  #VEi1
  checkBetween0and1(VEi1)
  checkDimensionsMatch(VEi1, populationFractions)
  #VEp1
  checkBetween0and1(VEp1)
  checkDimensionsMatch(VEp1, populationFractions)
   #VEs2
  checkBetween0and1(VEs2)
  checkDimensionsMatch(VEs2, populationFractions)
  #VEi2
  checkBetween0and1(VEi2)
  checkDimensionsMatch(VEi2, populationFractions)
  #VEp2
  checkBetween0and1(VEp2)
  checkDimensionsMatch(VEp2, populationFractions)
  #vaccineEfficacyDelay
  checkNonNegativeNumber(vaccineEfficacyDelay)

  #Compute the daily vaccination rates for priming and boost doses
  totalSimulationLength <- seedStartDay + simulationLength
  vaccinationRatePrimeByDay <- rep(0, totalSimulationLength)
  vaccinationRateBoostByDay <- rep(0, totalSimulationLength)
  currentVaccineAvailabilityPrime <- 0
  currentVaccineAvailabilityBoost <- 0
  newPrimedToBoost <- rep(0, totalSimulationLength)
  currentPrimedToBoost <- 0
  for (i in 1:totalSimulationLength) {
    if (i <= length(vaccineAvailabilityByDayPrime)){
      currentVaccineAvailabilityPrime <- currentVaccineAvailabilityPrime + vaccineAvailabilityByDayPrime[i] 
    }
    if (i <= length(vaccineAvailabilityByDayBoost)){
      currentVaccineAvailabilityBoost <- currentVaccineAvailabilityBoost + vaccineAvailabilityByDayBoost[i] 
    }
    currentPrimedToBoost <- currentPrimedToBoost + newPrimedToBoost[i]
    #Prioritize boosting over priming
    vaccinationRateBoostByDay[i] <- min(vaccineAdministrationRatePerDay, currentVaccineAvailabilityBoost,
                                        currentPrimedToBoost)
    vaccinationRatePrimeByDay[i] <- max(min(vaccineAdministrationRatePerDay, currentVaccineAvailabilityPrime) -
                                              vaccinationRateBoostByDay[i], 0)
    newPrimedToBoost[i + boostDelay] <- vaccinationRatePrimeByDay[i] #Works if boostDelay >=1
    currentPrimedToBoost <- currentPrimedToBoost - vaccinationRateBoostByDay[i]
    currentVaccineAvailabilityPrime <- currentVaccineAvailabilityPrime - vaccinationRatePrimeByDay[i]
    currentVaccineAvailabilityBoost <- currentVaccineAvailabilityBoost - vaccinationRateBoostByDay[i]
  }
  vaccinationRatePrimeByDay <- vaccinationRatePrimeByDay / population #Normalize
  vaccinationRateBoostByDay <- vaccinationRateBoostByDay / population
  
  #Define vaccination rate functions
  vaccinationRatePrime <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(vaccinationRatePrimeByDay[floor(t - vaccineEfficacyDelay + 1)])
    }
  }
  
  vaccinationRateBoost <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(vaccinationRateBoostByDay[floor(t - vaccineEfficacyDelay + 1)])
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
  return(list(vaccinationRatePrime = vaccinationRatePrime,
              vaccinationRateBoost = vaccinationRateBoost,
              vaccinationRateAgeMultiplier = vaccinationRateAgeMultiplier,
              VEs1 = VEs1, VEs2 = VEs2, VEi1 = VEi1, VEi2 = VEi2,
              VEp1 = VEp1, VEp2 = VEp2, vaccineEfficacyDelay = vaccineEfficacyDelay))
}

#This is a utility function that reconstructs the model state as a list so that equations can refer to compartments by name
reconstructState.SEIRVPrimeBoost <- function(state) {
  numberOfClasses <- length(state) / 16 #Each of the 16 classes are vectors of the same length
  S  <-  state[                        1 :     numberOfClasses ]
  E  <-  state[     (numberOfClasses + 1):(2 * numberOfClasses)]
  I  <-  state[ (2 * numberOfClasses + 1):(3 * numberOfClasses)]
  R  <-  state[ (3 * numberOfClasses + 1):(4 * numberOfClasses)]
  Sv <-  state[ (4 * numberOfClasses + 1):(5 * numberOfClasses)]
  Ev <-  state[ (5 * numberOfClasses + 1):(6 * numberOfClasses)]
  Iv <-  state[ (6 * numberOfClasses + 1):(7 * numberOfClasses)]
  Rv <-  state[ (7 * numberOfClasses + 1):(8 * numberOfClasses)]
  Svb <- state[ (8 * numberOfClasses + 1):(9 * numberOfClasses)]
  Evb <- state[ (9 * numberOfClasses + 1):(10 * numberOfClasses)]
  Ivb <- state[(10 * numberOfClasses + 1):(11 * numberOfClasses)]
  Rvb <- state[(11 * numberOfClasses + 1):(12 * numberOfClasses)]
  V  <-  state[(12 * numberOfClasses + 1):(13 * numberOfClasses)]
  Vb <-  state[(13 * numberOfClasses + 1):(14 * numberOfClasses)]
  vaccinatingPrime <-  state[(14 * numberOfClasses + 1):(15 * numberOfClasses)]
  vaccinatingBoost <-  state[(15 * numberOfClasses + 1):(16 * numberOfClasses)]
  return(as.list(environment()))
}

#This function Implements the multivariate derivative of the SEIR+V model with prime-boost dosing
#parameters should define populationFractions, contactMatrix, beta, lambda, gamma,
#VEs1, VEi1, VEs2, VEi2, VEp1, VEp2, and the functions vaccinationRatePrime(t) and vaccinationRateBoost(t)
#Note that the total population is normalized to be 1
getDerivative.SEIRVPrimeBoost <- function(t, state, parameters) {
  stateList <- reconstructState.SEIRVPrimeBoost(state)
  with(append(stateList, parameters), {
    if (useCommunityMitigation) {
      if ((t >= communityMitigationStartDay) && (t < communityMitigationEndDay)) {
        contactMatrix <- communityMitigationMultiplier * contactMatrix
      } 
    }

    isVaccinatingPrimeByAge <- (vaccinatingPrime > 0) & (populationFractions - V > 0)
    effectiveVaccinationPrimeMultiplier <- sum(ifelse(isVaccinatingPrimeByAge, 1, 0) * vaccinationRateAgeMultiplier)
    if (effectiveVaccinationPrimeMultiplier > 0) {
      vaccinationRatePrimeByAge <- vaccinationRatePrime(t) * vaccinationRateAgeMultiplier /
                                     effectiveVaccinationPrimeMultiplier
    } else {
      vaccinationRatePrimeByAge <- 0
    }

    isVaccinatingBoostByAge <- (vaccinatingBoost > 0) & (V - Vb > 0)
    effectiveVaccinationBoostMultiplier <- sum(ifelse(isVaccinatingBoostByAge, 1, 0) * vaccinationRateAgeMultiplier)
    if (effectiveVaccinationBoostMultiplier > 0) {
      vaccinationRateBoostByAge <- vaccinationRateBoost(t) * vaccinationRateAgeMultiplier /
                                     effectiveVaccinationBoostMultiplier
    } else {
      vaccinationRateBoostByAge <- 0
    }

    #Flows
    forceOfInfection <- beta / populationFractions * (contactMatrix %*% (I + ((1 - VEi1) * Iv) + ((1 - VEi2) * Ivb)))
    
    S_to_E <- S * forceOfInfection
    E_to_I <- lambda * E
    I_to_R <- gamma * I
    
    Sv_to_Ev <- Sv * (1 - VEs1) * forceOfInfection
    Ev_to_Iv <- lambda * Ev
    Iv_to_Rv <- gamma * Iv
    
    Svb_to_Evb <- Svb * (1 - VEs2) * forceOfInfection
    Evb_to_Ivb <- lambda * Evb
    Ivb_to_Rvb <- gamma * Ivb

    S_to_Sv <-  ifelse(isVaccinatingPrimeByAge, vaccinationRatePrimeByAge * S / (populationFractions - V), 0)
    Sv_to_Svb <-  ifelse(isVaccinatingBoostByAge, vaccinationRateBoostByAge * Sv / (V - Vb), 0)
    
    #Derivatives
    #Non-vaccinated compartments
    dS <- -S_to_E - S_to_Sv
    dE <- S_to_E - E_to_I
    dI <- E_to_I - I_to_R
    dR <- I_to_R
    #Vaccinated compartments
    dSv <- -Sv_to_Ev + S_to_Sv - Sv_to_Svb
    dEv <- Sv_to_Ev - Ev_to_Iv
    dIv <- Ev_to_Iv - Iv_to_Rv
    dRv <- Iv_to_Rv

    dSvb <- -Svb_to_Evb + Sv_to_Svb
    dEvb <- Svb_to_Evb - Evb_to_Ivb
    dIvb <- Evb_to_Ivb - Ivb_to_Rvb
    dRvb <- Ivb_to_Rvb
    #Auxiliary vaccinated compartment - people with a priming dose
    dV <- ifelse(isVaccinatingPrimeByAge, vaccinationRatePrimeByAge, 0)
    #Auxiliary vaccinated compartment - people with a boosting dose
    dVb <- ifelse(isVaccinatingBoostByAge, vaccinationRateBoostByAge, 0)

    zeroVec <- 0 * populationFractions
    
    #Return derivative
    return(list(c(dS, dE, dI, dR, dSv, dEv, dIv, dRv, dSvb, dEvb, dIvb, dRvb, dV, dVb, zeroVec, zeroVec)))
  })
}

#This function implements seeding infections in the SEIR+V 2 dose model
#parameters should define seedInfections, lambda, and gamma
doSeed.SEIRVPrimeBoost <- function(state, parameters) {
  stateList <- reconstructState.SEIRVPrimeBoost(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions
    E <- E + seedInfectionsFractions / (1 + lambda / gamma)
    I <- I + seedInfectionsFractions / (1 + gamma / lambda)
    #Return derivative
    return(c(S, E, I, R, Sv, Ev, Iv, Rv, Svb, Evb, Ivb, Rvb, V, Vb, vaccinatingPrime, vaccinatingBoost))
  })
}