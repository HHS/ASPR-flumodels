#' @title SEIRPrep Model
#' @description A model for influenza that uses a SEIR framework with age 
#'   structure
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
#' @param seedInfections Number of infections to seed; single number or vector
#'   of numbers by population group; defaults to 0
#' @param priorImmunity Fraction of the population with prior immunity; single 
#'   fraction, or vector of fractions by population group; defaults to 0
#' @param useCommunityMitigation Whether or not to use community mitigation
#'   implemented by modulation of the contact matrix; defaults to FALSE
#' @param communityMitigationStartDay If using community mitigation, days of the
#'   simulation on which each cycle of mitigation starts; must be specified 
#'   sequentially in time if applicable
#' @param communityMitigationDuration If using community mitigation, duration of
#'   time during which each cycle of mitigation is in effect (different cycles of
#'   mitigation cannot overlap); must be specified if applicable
#' @param communityMitigationMultiplier The communityMitigationMultiplier is the set
#'  of multipliers that will be used to modulate the contact
#'   matrix by elementwise multiplication and can be:
#' (i) a single non-negative scalar, if mitigation is the same for all subgroups 
#' and all mitigation cycles; or, (ii) an atomic vector of non-negative entries 
#' with length equal to the number of mitigation rounds, if mitigation is the 
#' same for all subgroups but varies by  mitigation cycles; or
#' (iii) a single 2D array of non-negative entries whose dimensions are identical
#'  to those of the contact matrix, if mitigation varies by subgroup but is the 
#'  same across mitigation cycles; or (iv) a list of 2D arrays, each of which has
#'  dimensions that are identical to those of the contact matrix;
#' the length of the list must match the number of mitigation cycles. This is 
#' the most general case, in which mitigation varies by subgroup and by mitigation cycle.
#' Must be specified if applicable.
#' @param prepDiagnosisRate Rate at which PrEP is diagnosed; defaults to 0
#' @param fractionOnPrEP Fraction of the population that takes PrEP; single fraction 
#' or vector of fractions by population group;  defaults to 0.0
#' @param timeBetweenPrEPDoses Time (in days) elapsed between PrEP doses; single number or 
#' or vector by population group; defaults to 1.0
#' @param AVEs Antiviral efficacy: prevention of transmission to susceptible
#'   individuals taking PrEP; single fraction or vector of fractions
#'   by population group; defaults to 0
#' @param AVEi Antiviral efficacy: prevention of transmission from infected 
#'   individuals taking PrEP; single fraction or vector of fractions
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
#' @return a SEIRModel object
#' @export
SEIRPrEPModel <- function(population, populationFractions, contactMatrix, R0,
                      latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                      useCommunityMitigation, communityMitigationStartDay,
                      communityMitigationDuration, communityMitigationMultiplier,
                      prepDiagnosisRate, fractionOnPrEP, timeBetweenPrEPDoses, AVEs, AVEi, AVEp,
                      simulationLength, seedStartDay, tolerance, method) {
  #Check inputs
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIRPrEP", argumentList) #Get parameters from checked inputs
  
  initialState <- with(parameters, {
    c(S  = (1 - priorImmunity) * populationFractions,
      E  = 0 * populationFractions,
      I  = 0 * populationFractions,
      R  = priorImmunity * populationFractions,
      Sprep = 0 * populationFractions,
      Eprep = 0 * populationFractions,
      Iprep = 0 * populationFractions,
      Rprep = 0 * populationFractions,
      PrEP = 0 * populationFractions)
  })
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIRPrEP,
                              seedFunction = doSeed.SEIRPrEP)
  
  #Build the SEIRModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- "SEIRPrEPModel"
  return(model)
}


#' @title Check SEIR+V inputs
#' @description Checks the input parameters for the SEIR+V model
#' @return List of parameters for the SEIR+V model
#' @keywords internal
checkInputs.SEIRPrEP <- function(population, populationFractions, contactMatrix, R0,
                                 latentPeriod, infectiousPeriod, seedInfections, priorImmunity,
                                 useCommunityMitigation, communityMitigationStartDay,
                                 communityMitigationDuration, communityMitigationMultiplier, 
                                 prepDiagnosisRate, fractionOnPrEP, timeBetweenPrEPDoses, AVEs, AVEi, AVEp,
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
  prepParameters <- do.call("checkInputs.PrEP", argumentList)
  #Return the parameters
  return(c(SEIRParameters, prepParameters))
}



#' @title Check PrEP inputs
#' @description Checks PrEP inputs and computes the parameters
#' @return List of PrEP parameters
#' @keywords internal
checkInputs.PrEP <- function(populationFractions, prepDiagnosisRate = 0, fractionOnPrEP = 0, timeBetweenPrEPDoses = 1, AVEs = 0, AVEi = 0, AVEp = 0, ...) {
  #Validate PrEP parameters
  #prepDiagnosisRate
  checkPositiveNumber(prepDiagnosisRate)
  #fractionOnPrEP
  checkBetween0and1(fractionOnPrEP)
  checkDimensionsMatch(fractionOnPrEP, populationFractions)
  #timeBetweenPrEPDoses
  checkPositiveNumber(timeBetweenPrEPDoses)
  checkDimensionsMatch(timeBetweenPrEPDoses, populationFractions)
  #AVEs
  checkBetween0and1(AVEs)
  checkDimensionsMatch(AVEs, populationFractions)
  #AVEi
  checkBetween0and1(AVEi)
  checkDimensionsMatch(AVEi, populationFractions)
  #AVEp
  checkBetween0and1(AVEp)
  checkDimensionsMatch(AVEp, populationFractions)
  
  return(list(prepDiagnosisRate = prepDiagnosisRate,
              fractionOnPrEP = fractionOnPrEP, timeBetweenPrEPDoses = timeBetweenPrEPDoses,
              AVEs = AVEs, AVEi = AVEi, AVEp = AVEp))
}




#This is a utility function that reconstructs the model state as a list so that equations can refer to compartments by name
reconstructState.SEIRPrEP <- function(state) {
  numberOfClasses <- length(state) / 9 #Each of the 9 classes are vectors of the same length
  S  <- state[                       1 :     numberOfClasses ]
  E  <- state[    (numberOfClasses + 1):(2 * numberOfClasses)]
  I  <- state[(2 * numberOfClasses + 1):(3 * numberOfClasses)]
  R  <- state[(3 * numberOfClasses + 1):(4 * numberOfClasses)]
  Sprep <- state[(4 * numberOfClasses + 1):(5 * numberOfClasses)]
  Eprep <- state[(5 * numberOfClasses + 1):(6 * numberOfClasses)]
  Iprep <- state[(6 * numberOfClasses + 1):(7 * numberOfClasses)]
  Rprep <- state[(7 * numberOfClasses + 1):(8 * numberOfClasses)]
  PrEP  <- state[(8 * numberOfClasses + 1):(9 * numberOfClasses)]
  return(as.list(environment()))
}

#This function implements the multivariate derivative of the SEIRPrEP model
#parameters should define populationFractions, normalizedContactMatrix, beta, lambda, and gamma
#Note that the total population is normalized to be 1
getDerivative.SEIRPrEP <- function(t, state, parameters) {
  stateList <- reconstructState.SEIRPrEP(state)
  with(append(stateList, parameters), {
    if (useCommunityMitigation) {
      # Is the simulation in the midst of a mitigation cycle?
      if ( any((t >= communityMitigationStartDay) & (t < communityMitigationEndDay)) ) {
        # If so, which one?
        communityMitigationCycle <- 
          which(communityMitigationStartDay <= t & t < communityMitigationEndDay)
        
        # communityMitigationMultiplier can be: (i) a single number for each mitigation cycle and for
        # all subgroups; (ii) a 2D matrix specifying a mitigation across subgroups that does not vary between
        # mitigation cycles; or (iii) or a list of 2D arrays specifying mitigation intensity across subgroups 
        # at each mitigation cycle
        if ( is.vector(communityMitigationMultiplier, mode = "numeric") ){
          contactMatrix <- communityMitigationMultiplier[communityMitigationCycle] * contactMatrix
          
        } else if ( is.matrix(communityMitigationMultiplier) ) {
          contactMatrix <- communityMitigationMultiplier * contactMatrix
          
        } else if ( is.list(communityMitigationMultiplier) ) {
          contactMatrix <- communityMitigationMultiplier[[communityMitigationCycle]] * contactMatrix
          
        }
      }
    }
    
    #Flows
    forceOfInfection <- beta / populationFractions * (contactMatrix %*% (I + (1 - AVEi) * Iprep))
    
    S_to_E <- S * forceOfInfection
    E_to_I <- lambda * E
    I_to_R <- gamma * I
    
    # S_to_Sv <-  ifelse(V < populationFractions, vaccinationRateByAge * S / (populationFractions - V), 0)
    S_to_Sprep <-  ifelse(Sprep + Eprep + Iprep + Rprep < populationFractions * fractionOnPrEP, S * prepDiagnosisRate, 0)
    
    Sprep_to_Eprep <- Sprep * (1 - AVEs) * forceOfInfection
    Eprep_to_Iprep <- lambda * Eprep
    Iprep_to_Rprep <- gamma * Iprep
    
    
    #Derivatives
    dS <- - S_to_E - S_to_Sprep
    dE <- S_to_E - E_to_I
    dI <- E_to_I - I_to_R
    dR <- I_to_R
    
    #PrEP compartments
    dSprep <- - Sprep_to_Eprep + S_to_Sprep
    dEprep <- Sprep_to_Eprep - Eprep_to_Iprep
    dIprep <- Eprep_to_Iprep - Iprep_to_Rprep
    dRprep <- Iprep_to_Rprep
    
    #Compartment for tracking doses
    dPrEP <- (Sprep + Eprep + Iprep + Rprep) / timeBetweenPrEPDoses
    
    #Return derivative
    return(list(c(dS, dE, dI, dR, dSprep, dEprep, dIprep, dRprep, dPrEP)))
  })
}

#This function implements seeding infections in the SEIR model
#parameters should define seedInfections, lambda, and gamma
#Note that the total population is normalized to be 1
doSeed.SEIRPrEP <- function(state, parameters) {
  stateList <- reconstructState.SEIRPrEP(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions * S / (S + Sprep)
    E <- E + seedInfectionsFractions * S / (S + Sprep) / (1 + lambda / gamma)
    I <- I + seedInfectionsFractions * S / (S + Sprep) / (1 + gamma / lambda)
    
    Sprep <- Sprep - seedInfectionsFractions * Sprep / (S + Sprep)
    Eprep <- Eprep + seedInfectionsFractions * Sprep / (S + Sprep) / (1 + lambda / gamma)
    Iprep <- Iprep + seedInfectionsFractions * Sprep / (S + Sprep) / (1 + gamma / lambda)
    #Return derivative
    return(c(S, E, I, R, Sprep, Eprep, Iprep, Rprep, PrEP))
  })
}
