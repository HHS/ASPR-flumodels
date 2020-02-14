#' @title SEIR+TV Model
#' @description A model for influenza that uses a SEIR framework with age 
#'   structure and that allows for antiviral treatment and vaccination
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
#' @param fractionSymptomatic Fraction of the infections that are symptomatic 
#'   cases; single fraction or vector of fractions by population group;  defaults
#'   to 0.5
#' @param fractionSeekCare Fraction of the symptomatic cases that seek medical
#'   care; single fraction or vector of fractions by population group; defaults
#'   to 0.6
#' @param fractionDiagnosedAndPrescribedOutpatient Fraction of the outpatient
#'   medical care seeking cases that are diagnosed and prescribed antiviral
#'   drugs; single fraction or vector of fractions by population group; defaults
#'   to 0.7
#' @param fractionAdhere Fraction of the cases that are prescribed antiviral
#'   drugs that adhere to the regimen; single fraction or vector of fractions by
#'   population group; defaults to 0.8
#' @param fractionAdmitted Fraction of the cases that require hospitalization
#'   that are admitted; single fraction or vector of fractions by population
#'   group; defaults to 1
#' @param fractionDiagnosedAndPrescribedInpatient Fraction of the hospitalized
#'   cases that are diagnosed and prescribed antiviral drugs; single fraction or
#'   vector of fractions by population group; defaults to 1
#' @param AVEi Antiviral efficacy: prevention of transmission from infected 
#'   individuals taking antiviral drugs; defaults to 0
#' @param AVEp Antiviral efficacy: probability that antiviral treatment averts 
#'   hospitalization and/or death; defaults to 0
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
#' @return a SEIRTVModel object
#' @export
SEIRTVModel <- function(population, populationFractions, contactMatrix, R0,
                        latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                        useCommunityMitigation, communityMitigationStartDay,
                        communityMitigationDuration, communityMitigationMultiplier,
                        fractionSymptomatic,  fractionSeekCare, fractionDiagnosedAndPrescribedOutpatient,
                        fractionAdhere, fractionAdmitted, fractionDiagnosedAndPrescribedInpatient, AVEi, AVEp,
                        vaccineAdministrationRatePerDay, vaccineAvailabilityByDay,
                        vaccineUptakeMultiplier, VEs, VEi, VEp, vaccineEfficacyDelay,
                        simulationLength, seedStartDay, tolerance, method) {
  #Check inputs #TODO: Add checks for all inputs
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIRTV", argumentList) #Get parameters from checked inputs

  initialState <- with(parameters, {
    c(S  = (1 - priorImmunity) * populationFractions,
      E  = 0 * populationFractions,
      I  = 0 * populationFractions,
      R  = priorImmunity * populationFractions,
      Sv = 0 * populationFractions,
      Ev = 0 * populationFractions,
      Iv = 0 * populationFractions,
      Rv = 0 * populationFractions,
      V  = 0 * populationFractions)
  })
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIRTV,
                              seedFunction = doSeed.SEIRV)
  
  #Build the SEIRVModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- c("SEIRTVModel", "SEIRTModel", "SEIRVModel", "SEIRModel")
  return(model)
}

#' @title Check SEIR+TV inputs
#' @description Checks the input parameters for the SEIR+V model
#' @return List of parameters for the SEIR+V model
#' @keywords internal
checkInputs.SEIRTV <- function(population, populationFractions, contactMatrix, R0,
                               latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                               useCommunityMitigation, communityMitigationStartDay,
                               communityMitigationDuration, communityMitigationMultiplier,
                               fractionSymptomatic,  fractionSeekCare, fractionDiagnosedAndPrescribedOutpatient,
                               fractionAdhere, fractionAdmitted, fractionDiagnosedAndPrescribedInpatient, AVEi, AVEp,
                               vaccineAdministrationRatePerDay, vaccineAvailabilityByDay,
                               vaccineUptakeMultiplier, VEs, VEi, VEp, vaccineEfficacyDelay,
                               simulationLength, seedStartDay, tolerance, method) {
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  SEIRParameters <- do.call("checkInputs.SEIR", argumentList)
  antiviralParameters <- do.call("checkInputs.Antiviral", argumentList)
  #Update arguments passed to checkInputs.Vaccine using SEIRParameters
  argumentList$population <- SEIRParameters$population
  argumentList$populationFractions <- SEIRParameters$populationFractions
  argumentList$seedStartDay <- SEIRParameters$seedStartDay
  argumentList$simulationLength <- SEIRParameters$simulationLength
  vaccineParameters <- do.call("checkInputs.Vaccine", argumentList)
  #Return the parameters
  return(c(SEIRParameters, antiviralParameters, vaccineParameters))
}

#This function implements the multivariate derivative of the SEIR+TV model
#parameters should define populationFractions, contactMatrix, beta, lambda, gamma,
#AVEi.eff, VEs, VEi, and the function vaccinationRate(t)
#Note that the total population is normalized to be 1
getDerivative.SEIRTV <- function(t, state, parameters) {
  stateList <- reconstructState.SEIRV(state)
  with(append(stateList, parameters), {
    if (useCommunityMitigation) {
      if ((t >= communityMitigationStartDay) && (t < communityMitigationEndDay)) {
        contactMatrix <- communityMitigationMultiplier * contactMatrix
      } 
    }

    effectiveVaccinationMultiplier <- sum(ifelse(V < populationFractions, 1, 0) * vaccinationRateAgeMultiplier)
    if (effectiveVaccinationMultiplier > 0) {
      vaccinationRateByAge <- vaccinationRate(t) * vaccinationRateAgeMultiplier /
        effectiveVaccinationMultiplier
    } else {
      vaccinationRateByAge <- 0
    }
    
    #Flows
    # forceOfInfection <- beta / populationFractions * (contactMatrix %*% ((1 - AVEi.eff) * (I + ((1 - VEi) * Iv))))
    # Adjusted to account for VEp, which reduces the impact of AVEi since it, in essence, reduces fractionSymptomatic
    forceOfInfection <- beta / populationFractions * (contactMatrix %*% ((1 - AVEi.eff) * (I + (1 - VEi) * Iv) + 
                                                                           VEp *AVEi.eff * (1 - VEi) * Iv ))
    
    S_to_E <- S * forceOfInfection
    E_to_I <- lambda * E 
    I_to_R <- gamma * I
    
    Sv_to_Ev <- Sv * (1 - VEs) * forceOfInfection
    Ev_to_Iv <- lambda * Ev
    Iv_to_Rv <- gamma * Iv
    
    S_to_Sv <-  ifelse(V < populationFractions, vaccinationRateByAge * S / (populationFractions - V), 0)
    
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
    dV <- ifelse(V < populationFractions, vaccinationRateByAge, 0)
    
    #Return derivative
    return(list(c(dS, dE, dI, dR, dSv, dEv, dIv, dRv, dV)))
  })
}