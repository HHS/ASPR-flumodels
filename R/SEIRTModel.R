#' @title SEIR+T Model
#' @description A model for influenza that uses a SEIR framework with age 
#'   structure and that allows for antiviral treatment
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
#'   individuals taking antiviral drugs; single fraction or vector of fractions
#'   by population group; defaults to 0
#' @param AVEp Antiviral efficacy: probability that antiviral treatment averts 
#'   hospitalization and/or death; single fraction or vector of fractions by
#'   population group; defaults to 0
#' @param simulationLength Number of days to simulate after seeding infections; 
#'   defaults to 240
#' @param seedStartDay Day on which to seed initial infections; defaults to 0
#' @param tolerance Absolute tolerance for numerical integration; defaults to 
#'   1e-8
#' @param method Which integration method to use. Defaults to lsoda
#' @return a SEIRTModel object
#' @export
SEIRTModel <- function(population, populationFractions, contactMatrix, R0,
                       latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                       useCommunityMitigation, communityMitigationStartDay,
                       communityMitigationDuration, communityMitigationMultiplier,
                       fractionSymptomatic,  fractionSeekCare, fractionDiagnosedAndPrescribedOutpatient,
                       fractionAdhere, fractionAdmitted, fractionDiagnosedAndPrescribedInpatient,
                       AVEi, AVEp, simulationLength, seedStartDay, tolerance, method) {
  #Check inputs
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIRT", argumentList) #Get parameters from checked inputs

  initialState <- with(parameters, {
    c(S  = (1 - priorImmunity) * populationFractions,
      E  = 0 * populationFractions,
      I  = 0 * populationFractions,
      R  = priorImmunity * populationFractions)
  })
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIRT,
                              seedFunction = doSeed.SEIR)
  
  #Build the SEIRModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- c("SEIRTModel", "SEIRModel")
  return(model)
}

#' @title Check SEIR+T inputs
#' @description Checks the input parameters for the SEIR+T model
#' @return List of parameters for the SEIR+T model
#' @keywords internal
checkInputs.SEIRT <- function(population, populationFractions, contactMatrix, R0,
                              latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                              useCommunityMitigation, communityMitigationStartDay,
                              communityMitigationDuration, communityMitigationMultiplier,
                              fractionSymptomatic,  fractionSeekCare, fractionDiagnosedAndPrescribedOutpatient,
                              fractionAdhere, fractionAdmitted, fractionDiagnosedAndPrescribedInpatient,
                              AVEi, AVEp, simulationLength, seedStartDay, tolerance, method) {
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  SEIRParameters <- do.call("checkInputs.SEIR", argumentList)
  antiviralParameters <- do.call("checkInputs.Antiviral", argumentList)
  #Return the parameters
  return(c(SEIRParameters, antiviralParameters))
}

#' @title Check antiviral treatment inputs
#' @description Checks antiviral treatment inputs and computes the antiviral treatment parameters
#' @return List of antiviral treatment parameters
#' @keywords internal
checkInputs.Antiviral <- function(populationFractions, fractionSymptomatic = 0.5,  fractionSeekCare = 0.6,
                                  fractionDiagnosedAndPrescribedOutpatient = 0.7, fractionAdhere = 0.8,
                                  fractionAdmitted = 1, fractionDiagnosedAndPrescribedInpatient = 1, AVEi = 0, AVEp = 0, ...) {
  #Validate parameters
  checkBetween0and1(fractionSymptomatic)
  checkDimensionsMatch(fractionSymptomatic, populationFractions)
  checkBetween0and1(fractionSeekCare)
  checkDimensionsMatch(fractionSeekCare, populationFractions)
  checkBetween0and1(fractionDiagnosedAndPrescribedOutpatient)
  checkDimensionsMatch(fractionDiagnosedAndPrescribedOutpatient, populationFractions)
  checkBetween0and1(fractionAdhere)
  checkDimensionsMatch(fractionAdhere, populationFractions)
  checkBetween0and1(fractionAdmitted)
  checkDimensionsMatch(fractionAdmitted, populationFractions)
  checkBetween0and1(fractionDiagnosedAndPrescribedInpatient)
  checkDimensionsMatch(fractionDiagnosedAndPrescribedInpatient, populationFractions)
  checkBetween0and1(AVEi)
  checkDimensionsMatch(AVEi, populationFractions)
  checkBetween0and1(AVEp)
  checkDimensionsMatch(AVEp, populationFractions)

  AVEi.eff <- AVEi * fractionAdhere * fractionDiagnosedAndPrescribedOutpatient * fractionSeekCare * fractionSymptomatic
  return(list(fractionSymptomatic = fractionSymptomatic, fractionSeekCare = fractionSeekCare, 
              fractionDiagnosedAndPrescribedOutpatient = fractionDiagnosedAndPrescribedOutpatient, fractionAdhere = fractionAdhere,
              fractionAdmitted = fractionAdmitted, fractionDiagnosedAndPrescribedInpatient = fractionDiagnosedAndPrescribedInpatient,
              AVEi.eff = AVEi.eff, AVEp = AVEp))
}

#This function implements the multivariate derivative of the SEIR+T model
#parameters should define populationFractions, normalizedContactMatrix,
#beta, lambda, gamma, and AVEi.eff
#Note that the total population is normalized to be 1
getDerivative.SEIRT <- function(t, state, parameters) {
  stateList <- reconstructState.SEIR(state)
  with(append(stateList, parameters), {
    if (useCommunityMitigation) {
      if ((t >= communityMitigationStartDay) && (t < communityMitigationEndDay)) {
        contactMatrix <- communityMitigationMultiplier * contactMatrix
      } 
    }
    #Flows
    forceOfInfection <- beta / populationFractions * (contactMatrix %*% ((1 - AVEi.eff) * I))
    
    S_to_E <- S * forceOfInfection
    E_to_I <- lambda * E
    I_to_R <- gamma * I
    
    #Derivatives
    dS <- -S_to_E
    dE <- S_to_E - E_to_I
    dI <- E_to_I - I_to_R
    dR <- I_to_R
    
    #Return derivative
    return(list(c(dS, dE, dI, dR)))
  })
}