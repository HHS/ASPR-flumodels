#' @title SEIR+V+Mono Model
#' @description A model for influenza that uses a SEIR framework with age 
#'   structure and that allows for vaccination and subsequent monovalent vaccination
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
#' @param vaccineAvailabilityDose1ByDay Vector that contains the amount of first dose of vaccine 
#'   available each day; defaults to 0
#' @param vaccineAvailabilityDoseMByDay Vector that contains the amount of monovalent vaccine 
#'   available each day; defaults to 0
#' @param vaccineUptakeMultiplier Vector of multipliers that determines the
#'   relative rate at which vaccine is given to each age group; defaults to
#'   vaccine being allotted proportionally by population
#' @param timeOfExclusiveMonovalent Time when a monovalent vaccine is exclusively administered. Defaults to never.
#' @param VEs1 Vaccine efficacy: protection for vaccinated susceptible 
#'   individuals after initial dose; single fraction or vector of fractions by
#'   population group; defaults to 0
#' @param VEsM Vaccine efficacy: protection for vaccinated susceptible 
#'   individuals after monovalent dose; single fraction or vector of fractions by
#'   population group; defaults to 0
#' @param VEi1 Vaccine efficacy: prevention of transmission from vaccinated 
#'   infected individuals after initial dose; single fraction or vector of fractions
#'   by population group; defaults to 0
#' @param VEiM Vaccine efficacy: prevention of transmission from vaccinated 
#'   infected individuals after monovalent dose; single fraction or vector of fractions
#'   by population group; defaults to 0
#' @param VEp1 Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals after initial dose; single fraction or vector of fractions 
#'   by population group; defaults to 0
#' @param VEpM Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals after monovalent dose; single fraction or vector of fractions 
#'   by population group; defaults to 0
#' @param vaccineEfficacyDelay Delay in days between administration of dose and 
#'   onset of protection; defaults to 7
#' @param simulationLength Number of days to simulate after seeding infections; 
#'   defaults to 240
#' @param seedStartDay Day on which to seed initial infections; defaults to 0
#' @param tolerance Absolute tolerance for numerical integration; defaults to 
#'   1e-8
#' @param method Which integration method to use. Defaults to lsoda
#' @return a SEIRVMonoModel object
#' @export
SEIRVMonoModel <- function(population, populationFractions, contactMatrix, R0,
                            latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                            useCommunityMitigation, communityMitigationStartDay,
                            communityMitigationDuration, communityMitigationMultiplier,
                            vaccineAdministrationRatePerDay, vaccineAvailabilityDose1ByDay, vaccineAvailabilityDoseMByDay,
                            vaccineUptakeMultiplier, timeOfExclusiveMonovalent, VEs1, VEsM, VEi1, VEiM, VEp1, VEpM, vaccineEfficacyDelay,
                            simulationLength, seedStartDay, tolerance, method) {
  #Check inputs 
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIRVMono", argumentList) #Get parameters from checked inputs
  
  initialState <- with(parameters, {
    c(S   = (1 - priorImmunity) * populationFractions,
      E   = 0 * populationFractions,
      I   = 0 * populationFractions,
      R   = priorImmunity * populationFractions,
      Sv  = 0 * populationFractions,
      Ev  = 0 * populationFractions,
      Iv  = 0 * populationFractions,
      Rv  = 0 * populationFractions,
      SvM = 0 * populationFractions,
      EvM = 0 * populationFractions,
      IvM = 0 * populationFractions,
      RvM = 0 * populationFractions,
      V   = 0 * populationFractions,
      VM  = 0 * populationFractions)
  })
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIRVMono,
                              seedFunction = doSeed.SEIRVMono, method = method)
  
  #Build the SEIRVMonoModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- c("SEIRVMonoModel", "SEIRVModel", "SEIRModel")
  return(model)
}


#' @title Check SEIR+VM inputs
#' @description Checks the input parameters for the SEIR+V model
#' @return List of parameters for the SEIR+V model
#' @keywords internal
checkInputs.SEIRVMono <- function(population, populationFractions, contactMatrix, R0,
                                   latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                                   useCommunityMitigation, communityMitigationStartDay,
                                   communityMitigationDuration, communityMitigationMultiplier,
                                   vaccineAdministrationRatePerDay, vaccineAvailabilityDose1ByDay, vaccineAvailabilityDoseMByDay,
                                   vaccineUptakeMultiplier, timeOfExclusiveMonovalent, VEs1, VEsM, VEi1, VEiM,
                                   VEp1, VEpM, vaccineEfficacyDelay, simulationLength, seedStartDay, tolerance, method) {
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  SEIRParameters <- do.call("checkInputs.SEIR", argumentList)
  #Update call to checkInputs.VaccineMDose using SEIRParameters
  argumentList$population <- SEIRParameters$population
  argumentList$populationFractions <- SEIRParameters$populationFractions
  argumentList$seedStartDay <- SEIRParameters$seedStartDay
  argumentList$simulationLength <- SEIRParameters$simulationLength
  vaccineMDoseParameters <- do.call("checkInputs.VaccineMonovalent", argumentList)
  #Return the parameters
  return(c(SEIRParameters, vaccineMDoseParameters))
}

#' @title Check vaccine inputs
#' @description Checks vaccine inputs and computes the vaccine parameters
#' @return List of vaccine parameters
#' @keywords internal
checkInputs.VaccineMonovalent <- function(population, populationFractions, seedStartDay, simulationLength,
                                          vaccineAdministrationRatePerDay = 0, vaccineAvailabilityDose1ByDay = 0,
                                          vaccineAvailabilityDoseMByDay = 0,
                                          vaccineUptakeMultiplier = 1, timeOfExclusiveMonovalent = Inf, VEs1 = 0, VEsM = 0,
                                          VEi1 = 0, VEiM = 0, VEp1 = 0, VEpM = 0, vaccineEfficacyDelay = 7, ...) {
  #Validate vaccine parameters
  #vaccineAdministrationRatePerDay
  checkNonNegativeNumber(vaccineAdministrationRatePerDay)
  #vaccineAvailabilityDose1ByDay
  checkNonNegative(vaccineAvailabilityDose1ByDay)
  #vaccineAvailabilityDoseMByDay
  checkNonNegative(vaccineAvailabilityDoseMByDay)
  #vaccineUptakeMultiplier
  checkNonNegative(vaccineUptakeMultiplier)
  checkDimensionsMatch(vaccineUptakeMultiplier, populationFractions)
  #timeOfExclusiveMonovalent
  checkNonNegativeNumber(timeOfExclusiveMonovalent)
  #VEs1
  checkBetween0and1(VEs1)
  checkDimensionsMatch(VEs1, populationFractions)
  #VEi1
  checkBetween0and1(VEi1)
  checkDimensionsMatch(VEi1, populationFractions)
  #VEp1
  checkBetween0and1(VEp1)
  checkDimensionsMatch(VEp1, populationFractions)
  #VEsM
  checkBetween0and1(VEsM)
  checkDimensionsMatch(VEsM, populationFractions)
  #VEiM
  checkBetween0and1(VEiM)
  checkDimensionsMatch(VEiM, populationFractions)
  #VEpM
  checkBetween0and1(VEpM)
  checkDimensionsMatch(VEpM, populationFractions)
  #vaccineEfficacyDelay
  checkNonNegativeNumber(vaccineEfficacyDelay)
  
  #Compute the daily vaccination rates for both doses
  totalSimulationLength <- seedStartDay + simulationLength
  vaccinationRateDose1ByDay <- rep(0, totalSimulationLength)
  vaccinationRateDoseMByDay <- rep(0, totalSimulationLength)
  currentVaccineAvailabilityDose1 <- 0
  currentVaccineAvailabilityDoseM <- 0
  for (i in 1:totalSimulationLength) {
    if (i <= length(vaccineAvailabilityDose1ByDay)){
      currentVaccineAvailabilityDose1 <- currentVaccineAvailabilityDose1 + vaccineAvailabilityDose1ByDay[i]
    }
    if (i <= length(vaccineAvailabilityDoseMByDay)){
      currentVaccineAvailabilityDoseM <- currentVaccineAvailabilityDoseM + vaccineAvailabilityDoseMByDay[i]
    }

    vaccinationRateDoseMByDay[i] <- min(vaccineAdministrationRatePerDay, currentVaccineAvailabilityDoseM)
    # Presume that we give dose 1 along with the monovalent
    currentVaccineAvailabilityDose1 <- currentVaccineAvailabilityDose1 - vaccinationRateDoseMByDay[i]
    # Give dose 1 if we don't have enough monovalent to give    
    vaccinationRateDose1ByDay[i] <- ifelse(i < timeOfExclusiveMonovalent, 
                                           max(min(vaccineAdministrationRatePerDay - vaccinationRateDoseMByDay[i], 
                                               currentVaccineAvailabilityDose1), 0), 0)
    currentVaccineAvailabilityDose1 <- currentVaccineAvailabilityDose1 - vaccinationRateDose1ByDay[i]
    currentVaccineAvailabilityDoseM <- currentVaccineAvailabilityDoseM - vaccinationRateDoseMByDay[i]
  }
  vaccinationRateDose1ByDay <- vaccinationRateDose1ByDay / population #Normalize
  vaccinationRateDoseMByDay <- vaccinationRateDoseMByDay / population #Normalize
  
  #Define vaccination rate functions
  vaccinationRateDose1 <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(vaccinationRateDose1ByDay[floor(t - vaccineEfficacyDelay + 1)])
    }
  }
  
  vaccinationRateDoseM <- function(t) {
    if ((t < vaccineEfficacyDelay) || (t >= totalSimulationLength + vaccineEfficacyDelay)) {
      return(0)
    } else {
      return(vaccinationRateDoseMByDay[floor(t - vaccineEfficacyDelay + 1)])
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
              vaccinationRateDoseM = vaccinationRateDoseM,
              vaccinationRateAgeMultiplier = vaccinationRateAgeMultiplier,
              timeOfExclusiveMonovalent = timeOfExclusiveMonovalent, VEs1 = VEs1, VEsM = VEsM,
              VEi1 = VEi1, VEiM = VEiM, VEp1 = VEp1, VEpM = VEpM, vaccineEfficacyDelay = vaccineEfficacyDelay))
}

#This is a utility function that reconstructs the model state as a list so that equations can refer to compartments by name
reconstructState.SEIRVMono <- function(state) {
  numberOfClasses <- length(state) / 14 #Each of the 9 classes are vectors of the same length
  S  <-  state[                        1 :     numberOfClasses ]
  E  <-  state[     (numberOfClasses + 1):(2 * numberOfClasses)]
  I  <-  state[ (2 * numberOfClasses + 1):(3 * numberOfClasses)]
  R  <-  state[ (3 * numberOfClasses + 1):(4 * numberOfClasses)]
  Sv <-  state[ (4 * numberOfClasses + 1):(5 * numberOfClasses)]
  Ev <-  state[ (5 * numberOfClasses + 1):(6 * numberOfClasses)]
  Iv <-  state[ (6 * numberOfClasses + 1):(7 * numberOfClasses)]
  Rv <-  state[ (7 * numberOfClasses + 1):(8 * numberOfClasses)]
  SvM <- state[ (8 * numberOfClasses + 1):(9 * numberOfClasses)]
  EvM <- state[ (9 * numberOfClasses + 1):(10 * numberOfClasses)]
  IvM <- state[(10 * numberOfClasses + 1):(11 * numberOfClasses)]
  RvM <- state[(11 * numberOfClasses + 1):(12 * numberOfClasses)]
  V  <-  state[(12 * numberOfClasses + 1):(13 * numberOfClasses)]
  VM  <- state[(13 * numberOfClasses + 1):(14 * numberOfClasses)]
  return(as.list(environment()))
}

#This function implements the multivariate derivative of the SEIR+VM model
#parameters should define populationFractions, contactMatrix, beta, lambda, gamma,
#VEs, VEi, and the function vaccinationRate(t)
#Note that the total population is normalized to be 1
getDerivative.SEIRVMono <- function(t, state, parameters) {
  stateList <- reconstructState.SEIRVMono(state)
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
    
    effectiveVaccinationDoseMMultiplier <- sum(ifelse(VM < populationFractions, 1, 0) * vaccinationRateAgeMultiplier)
    if (effectiveVaccinationDoseMMultiplier > 0) {
      vaccinationRateDoseMByAge <- vaccinationRateDoseM(t) * vaccinationRateAgeMultiplier /
        effectiveVaccinationDoseMMultiplier
    } else {
      vaccinationRateDoseMByAge <- 0
    }
   
    #Flows
    forceOfInfection <- beta / populationFractions * (contactMatrix %*% (I + ((1 - VEi1) * Iv) + ((1 - VEiM) * IvM)))

    S_to_E <- S * forceOfInfection
    E_to_I <- lambda * E
    I_to_R <- gamma * I
    
    Sv_to_Ev <- Sv * (1 - VEs1) * forceOfInfection
    Ev_to_Iv <- lambda * Ev
    Iv_to_Rv <- gamma * Iv
    
    SvM_to_EvM <- SvM * (1 - VEsM) * forceOfInfection
    EvM_to_IvM <- lambda * EvM
    IvM_to_RvM <- gamma * IvM
    
    # Presume that only new people get monovalent vaccine
    S_to_Sv <-  ifelse(V < populationFractions, vaccinationRateDose1ByAge * S / (populationFractions - V), 0)
#     S_to_SvM <- ifelse(VM < populationFractions, vaccinationRateDoseMByAge * S^2 / ((S + Sv) * (populationFractions - VM - V)), 0)
#     Sv_to_SvM <-  ifelse(VM < populationFractions, vaccinationRateDoseMByAge * Sv^2 / ((S + Sv) * (V - VM)), 0)
    S_to_SvM <- ifelse(VM + V < populationFractions, vaccinationRateDoseMByAge * S / ( (populationFractions - VM - V)), 0)
    Sv_to_SvM <- 0
    
    #Derivatives
    #Non-vaccinated compartments
    dS <- -S_to_E - S_to_Sv - S_to_SvM
    dE <- S_to_E - E_to_I
    dI <- E_to_I - I_to_R
    dR <- I_to_R
    #Vaccinated compartments
    dSv <- -Sv_to_Ev + S_to_Sv - Sv_to_SvM
    dEv <- Sv_to_Ev - Ev_to_Iv
    dIv <- Ev_to_Iv - Iv_to_Rv
    dRv <- Iv_to_Rv
    
    dSvM <- -SvM_to_EvM + Sv_to_SvM + S_to_SvM
    dEvM <- SvM_to_EvM - EvM_to_IvM
    dIvM <- EvM_to_IvM - IvM_to_RvM
    dRvM <- IvM_to_RvM
    #Auxiliary vaccinated compartment
    # Would theoretically give V to everyone who gets VM, but it's difficult to track that, so only track V by itself
    dV <- ifelse(V < populationFractions, vaccinationRateDose1ByAge, 0) #+ ifelse(VM < populationFractions, vaccinationRateDoseMByAge, 0)
    #Auxiliary vaccinated compartment - people with monovalent dose
    dVM <- ifelse(VM + V < populationFractions, vaccinationRateDoseMByAge, 0)
    
    #Return derivative
    return(list(c(dS, dE, dI, dR, dSv, dEv, dIv, dRv, dSvM, dEvM, dIvM, dRvM, dV, dVM)))
  })
}

#This function implements seeding infections in the SEIR+V model
#parameters should define seedInfections, lambda, and gamma
doSeed.SEIRVMono <- function(state, parameters) {
  stateList <- reconstructState.SEIRVMono(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions
    E <- E + seedInfectionsFractions / (1 + lambda / gamma)
    I <- I + seedInfectionsFractions / (1 + gamma / lambda)
    #Return derivative
    return(c(S, E, I, R, Sv, Ev, Iv, Rv, SvM, EvM, IvM, RvM, V, VM))
  })
}