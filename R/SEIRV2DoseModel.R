#' @title SEIR+V 2-dose Model
#' @description A model for influenza that uses a SEIR framework with age 
#'   structure and that allows for vaccination with a 2-dose vaccine
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
#' @param vaccineAvailabilityByDay Vector that contains the amount of vaccine 
#'   available each day; defaults to 0
#' @param vaccineUptakeMultiplier Vector of multipliers that determines the
#'   relative rate at which vaccine is given to each age group; defaults to
#'   vaccine being allotted proportionally by population
#' @param dose2Delay Days between doses 1 and 2; defaults to 0
#' @param VEs1 Vaccine efficacy: protection for vaccinated susceptible 
#'   individuals after 1 dose; single fraction or vector of fractions by
#'   population group; defaults to 0
#' @param VEs2 Vaccine efficacy: protection for vaccinated susceptible 
#'   individuals after 2 doses; single fraction or vector of fractions by
#'   population group; defaults to 0
#' @param VEi1 Vaccine efficacy: prevention of transmission from vaccinated 
#'   infected individuals after 1 dose; single fraction or vector of fractions
#'   by population group; defaults to 0
#' @param VEi2 Vaccine efficacy: prevention of transmission from vaccinated 
#'   infected individuals after 2 doses; single fraction or vector of fractions
#'   by population group; defaults to 0
#' @param VEp1 Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals after 1 dose; single fraction or vector of fractions 
#'   by population group; defaults to 0
#' @param VEp2 Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals after 2 doses; single fraction or vector of fractions 
#'   by population group; defaults to 0
#' @param vaccineEfficacyDelay Delay in days between administration of dose and 
#'   onset of protection; defaults to 7
#' @param simulationLength Number of days to simulate after seeding infections; 
#'   defaults to 240
#' @param seedStartDay Day on which to seed initial infections; defaults to 0
#' @param tolerance Absolute tolerance for numerical integration; defaults to 
#'   1e-8
#' @param method Which integration method to use. Defaults to lsoda
#' @return a SEIRV2DoseModel object
#' @export
SEIRV2DoseModel <- function(population, populationFractions, contactMatrix, R0,
                            latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                            useCommunityMitigation, communityMitigationStartDay,
                            communityMitigationDuration, communityMitigationMultiplier,
                            vaccineAdministrationRatePerDay, vaccineAvailabilityByDay,
                            vaccineUptakeMultiplier, dose2Delay, VEs1, VEs2, VEi1, VEi2, VEp1, VEp2, vaccineEfficacyDelay,
                            simulationLength, seedStartDay, tolerance, method) {
  #Check inputs #TODO: Add checks for all inputs
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIRV2Dose", argumentList) #Get parameters from checked inputs
  
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
      Vb  = 0 * populationFractions)
  })
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIRV2Dose,
                              seedFunction = doSeed.SEIRV2Dose, method = method)
  
  #Build the SEIRV2DoseModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- c("SEIRV2DoseModel", "SEIRVModel", "SEIRModel")
  return(model)
}

#' @title Check SEIR+V 2-dose inputs
#' @description Checks the input parameters for the SEIR+V 2-dose model
#' @return List of parameters for the SEIR+V 2-dose model
#' @keywords internal
checkInputs.SEIRV2Dose <- function(population, populationFractions, contactMatrix, R0,
                                   latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                                   useCommunityMitigation, communityMitigationStartDay,
                                   communityMitigationDuration, communityMitigationMultiplier,
                                   vaccineAdministrationRatePerDay, vaccineAvailabilityByDay,
                                   vaccineUptakeMultiplier, dose2Delay, VEs1, VEs2, VEi1, VEi2,
                                   VEp1, VEp2, vaccineEfficacyDelay, simulationLength, seedStartDay, tolerance, method) {
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  SEIRParameters <- do.call("checkInputs.SEIR", argumentList)
  #Update call to checkInputs.Vaccine2Dose using SEIRParameters
  argumentList$population <- SEIRParameters$population
  argumentList$populationFractions <- SEIRParameters$populationFractions
  argumentList$seedStartDay <- SEIRParameters$seedStartDay
  argumentList$simulationLength <- SEIRParameters$simulationLength
  vaccine2DoseParameters <- do.call("checkInputs.Vaccine2Dose", argumentList)
  #Return the parameters
  return(c(SEIRParameters, vaccine2DoseParameters))
}

#' @title Check 2-Dose vaccine inputs
#' @description Checks 2-Dose vaccine inputs and computes the vaccine parameters
#' @return List of 2-Dose vaccine parameters
#' @keywords internal
checkInputs.Vaccine2Dose <- function(population, populationFractions, seedStartDay, simulationLength,
                                     vaccineAdministrationRatePerDay = 0, vaccineAvailabilityByDay = 0,
                                     vaccineUptakeMultiplier = 1, dose2Delay = 0, VEs1 = 0, VEs2 = 0,
                                     VEi1 = 0, VEi2 = 0, VEp1 = 0, VEp2 = 0, vaccineEfficacyDelay = 7, ...) {
  #Validate vaccine parameters
  #vaccineAdministrationRatePerDay
  checkNonNegativeNumber(vaccineAdministrationRatePerDay)
  #vaccineAvailabilityByDay
  checkNonNegative(vaccineAvailabilityByDay)
  #vaccineUptakeMultiplier
  checkNonNegative(vaccineUptakeMultiplier)
  checkDimensionsMatch(vaccineUptakeMultiplier, populationFractions)
  #dose2Delay
  checkNonNegativeNumber(dose2Delay)
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

  #Compute the daily vaccination rates for both doses
  totalSimulationLength <- seedStartDay + simulationLength
  vaccinationRateDose1ByDay <- rep(0, totalSimulationLength)
  vaccinationRateDose2ByDay <- rep(0, totalSimulationLength)
  currentVaccineAvailability <- 0
  for (i in 1:totalSimulationLength) {
    if (i <= length(vaccineAvailabilityByDay)){
      currentVaccineAvailability <- currentVaccineAvailability + vaccineAvailabilityByDay[i] 
    }
    if (dose2Delay > 0) {
      vaccinationRateDose1ByDay[i] <- min(vaccineAdministrationRatePerDay - 
        vaccinationRateDose2ByDay[i], currentVaccineAvailability / 2, na.rm = TRUE)
    } else {
      vaccinationRateDose1ByDay[i] <- min(vaccineAdministrationRatePerDay/2,
                                          currentVaccineAvailability / 2, na.rm = TRUE)
    }
    if (i + dose2Delay <= totalSimulationLength) { #Schedule 2nd dose
      vaccinationRateDose2ByDay[i + dose2Delay] <- vaccinationRateDose1ByDay[i]
    }
    currentVaccineAvailability <- currentVaccineAvailability - 2 * vaccinationRateDose1ByDay[i]
  }
  vaccinationRateDose1ByDay <- vaccinationRateDose1ByDay / population #Normalize
  vaccinationRateDose2ByDay <- vaccinationRateDose2ByDay / population
  
  #Define vaccination rate functions
  vaccinationRateDose1 <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(vaccinationRateDose1ByDay[floor(t - vaccineEfficacyDelay + 1)])
    }
  }
  
  vaccinationRateDose2 <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(vaccinationRateDose2ByDay[floor(t - vaccineEfficacyDelay + 1)])
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
  return(list(vaccinationRateDose1 = vaccinationRateDose1,
              vaccinationRateDose2 = vaccinationRateDose2,
              vaccinationRateAgeMultiplier = vaccinationRateAgeMultiplier,
              dose2Delay = dose2Delay, VEs1 = VEs1, VEs2 = VEs2,
              VEi1 = VEi1, VEi2 = VEi2, VEp1 = VEp1, VEp2 = VEp2, vaccineEfficacyDelay = vaccineEfficacyDelay))
}

#This is a utility function that reconstructs the model state as a list so that equations can refer to compartments by name
reconstructState.SEIRV2Dose <- function(state) {
  numberOfClasses <- length(state) / 14 #Each of the 14 classes are vectors of the same length
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
  return(as.list(environment()))
}

#This function Implements the multivariate derivative of the SEIR+V model with two doses
#parameters should define populationFractions, contactMatrix, beta, lambda, gamma,
#VEs1, VEi1, VEs2, VEi2, VEp1, VEp2, dose2Delay, and the functions vaccinationRateDose1(t) and vaccinationRateDose2(t)
#Note that the total population is normalized to be 1
getDerivative.SEIRV2Dose <- function(t, state, parameters) {
  stateList <- reconstructState.SEIRV2Dose(state)
  with(append(stateList, parameters), {
    if (useCommunityMitigation) {
      if ((t >= communityMitigationStartDay) && (t < communityMitigationEndDay)) {
        contactMatrix <- communityMitigationMultiplier * contactMatrix
      } 
    }

    effectiveVaccinationDose1Multiplier <- sum(ifelse(V < populationFractions, 1, 0) * vaccinationRateAgeMultiplier)
    if (effectiveVaccinationDose1Multiplier > 0) {
      vaccinationRateDose1ByAge <- vaccinationRateDose1(t) * vaccinationRateAgeMultiplier /
                                     effectiveVaccinationDose1Multiplier
    } else {
      vaccinationRateDose1ByAge <- 0
    }

    effectiveVaccinationDose2Multiplier <- sum(ifelse(Vb < V, 1, 0) * vaccinationRateAgeMultiplier)
    if (effectiveVaccinationDose2Multiplier > 0) {
      vaccinationRateDose2ByAge <- vaccinationRateDose2(t) * vaccinationRateAgeMultiplier /
                                     effectiveVaccinationDose2Multiplier
    } else {
      vaccinationRateDose2ByAge <- 0
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

    S_to_Sv <-  ifelse(V < populationFractions, vaccinationRateDose1ByAge * S / (populationFractions - V), 0)
    Sv_to_Svb <-  ifelse(Vb < V, vaccinationRateDose2ByAge * Sv / (V - Vb), 0)
    
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
    #Auxiliary vaccinated compartment - people with at least the 1st dose
    dV <- ifelse(V < populationFractions, vaccinationRateDose1ByAge, 0)
    #Auxiliary vaccinated compartment - people with 2 doses
    dVb <- ifelse(Vb < V, vaccinationRateDose2ByAge, 0)
    
    #Return derivative
    return(list(c(dS, dE, dI, dR, dSv, dEv, dIv, dRv, dSvb, dEvb, dIvb, dRvb, dV, dVb)))
  })
}

#This function implements seeding infections in the SEIR+V 2 dose model
#parameters should define seedInfections, lambda, and gamma
doSeed.SEIRV2Dose <- function(state, parameters) {
  stateList <- reconstructState.SEIRV2Dose(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions
    E <- E + seedInfectionsFractions / (1 + lambda / gamma)
    I <- I + seedInfectionsFractions / (1 + gamma / lambda)
    #Return derivative
    return(c(S, E, I, R, Sv, Ev, Iv, Rv, Svb, Evb, Ivb, Rvb, V, Vb))
  })
}