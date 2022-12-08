#' @title SEIR Model
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
#' @param simulationLength Number of days to simulate after seeding infections; 
#'   defaults to 240
#' @param seedStartDay Day on which to seed initial infections; defaults to 0
#' @param tolerance Absolute tolerance for numerical integration; defaults to 
#'   1e-8
#' @param method Which integration method to use. Defaults to lsoda
#' @return a SEIRModel object
#' @export
SEIRModel <- function(population,
                      populationFractions,
                      contactMatrix,
                      R0,
                      latentPeriod,
                      infectiousPeriod,
                      seedInfections,
                      priorImmunity,
                      useCommunityMitigation,
                      communityMitigationStartDay,
                      communityMitigationDuration,
                      communityMitigationMultiplier,
                      simulationLength,
                      seedStartDay,
                      tolerance,
                      method)
{
  #Check inputs
  specifiedArguments <- names(match.call())[-1]
  argumentList <- lapply(specifiedArguments, as.name)
  names(argumentList) <- specifiedArguments
  parameters <- do.call("checkInputs.SEIR", argumentList) #Get parameters from checked inputs
  
  initialState <- with(parameters,
                       {
                         c(S  = (1 - priorImmunity) * populationFractions,
                           E  = 0 * populationFractions,
                           I  = 0 * populationFractions,
                           R  = priorImmunity * populationFractions)
                       })
  
  rawOutput <- integrateModel(initialState = initialState,
                              parameters = parameters,
                              derivativeFunction = getDerivative.SEIR,
                              seedFunction = doSeed.SEIR)
  
  #Build the SEIRModel object to return
  model <- list(parameters = parameters,
                rawOutput = rawOutput)
  class(model) <- "SEIRModel"
  return(model)
}

#' @title Check SEIR inputs
#' @description Checks the input parameters for the SEIR model
#' @return List of parameters for the SEIR model
#' @keywords internal
checkInputs.SEIR <- function(population = 1,
                             populationFractions = 1,
                             contactMatrix,
                             R0,
                             latentPeriod,
                             infectiousPeriod,
                             seedInfections = 0,
                             priorImmunity = 0,
                             useCommunityMitigation = FALSE,
                             communityMitigationStartDay, 
                             communityMitigationDuration,
                             communityMitigationMultiplier, 
                             simulationLength = 240,
                             seedStartDay = 0,
                             tolerance = 1e-8,
                             method = "lsoda", ...)
{
  #population
  checkPositiveNumber(population)
  #populationFractions
  if (!((abs(sum(populationFractions) - 1) < tolerance) && all(populationFractions >= 0)))
  {
    stop("populationFractions must be positive and sum to 1.", call. = FALSE)
  }
  populationFractions <- as.vector(populationFractions) #Drop names to prevent errors in internal data storage
  #contactMatrix
  if (missing(contactMatrix))
  { #If contact matrix is not supplied, default to proportional mixing
    contactMatrix <- as.matrix(populationFractions)[ , rep(1, length(populationFractions)), drop = FALSE]
  }
  else
  {
    if (!(all(contactMatrix >=0) && (ncol(contactMatrix) == nrow(contactMatrix))))
    {
      stop("contactMatrix must be square and non-negative.", call. = FALSE)
    }
    if (length(populationFractions) != ncol(contactMatrix))
    {
      stop("contactMatrix dimensions do not match populationFractions.", call. = FALSE)
    }
  } 
  #R0
  if (missing(R0))
  {
    stop("R0 must be specified.", call. = FALSE)
  }
  checkPositive(R0)
  #latentPeriod
  if (missing(latentPeriod))
  {
    stop("latentPeriod must be specified.", call. = FALSE)
  }
  checkPositive(latentPeriod)
  #infectiousPeriod
  if (missing(infectiousPeriod))
  {
    stop("infectiousPeriod must be specified.", call. = FALSE)
  }
  checkPositive(infectiousPeriod)
  #seedInfections
  checkNonNegative(seedInfections)
  checkDimensionsMatch(seedInfections, populationFractions)
  if (length(seedInfections) == 1)
  {
    if (seedInfections > population)
    {
      stop("seedInfections can not exceed the population.", call. = FALSE)
    }
    seedInfections <- seedInfections * populationFractions #Distribute seed infections among population groups proportionately
  } else if (!all(seedInfections <= (population * populationFractions)))
  {
    stop("seedInfections can not exceed the population by group.", call. = FALSE)
  }
  #priorImmunity
  checkBetween0and1(priorImmunity)
  checkDimensionsMatch(priorImmunity, populationFractions)
  #Community Mitigation
  if (useCommunityMitigation) 
  {
    if (missing(communityMitigationStartDay))
    {
      stop("communityMitigationStartDay must be specified when using community mitigation.", 
           call. = FALSE)
    }
    checkNonNegative(communityMitigationStartDay)
    if (any(sort(communityMitigationStartDay) != communityMitigationStartDay)) {
      stop("communityMitigationStartDays must be ordered sequentially in time.", 
           call. = FALSE)
    }
    if (missing(communityMitigationDuration)) {
      stop("communityMitigationDuration must be specified when using community mitigation.", 
           call. = FALSE)
    }
    checkNonNegative(communityMitigationDuration)
    if ((length(communityMitigationStartDay) != length(communityMitigationDuration)) &
        (length(communityMitigationDuration) != 1)) {
      stop("communityMitigationStartDay and communityMitigationDuration must have the same length or
           length(communityMitigationDuration) = 1 for it to be recycled.", 
           call. = FALSE)
    }
    if ( any(tail(communityMitigationStartDay, -1) < head(communityMitigationStartDay + communityMitigationDuration, -1))){
      stop("There cannot be overlap between mitigation time windows.", 
           call. = FALSE)
    }
    if (missing(communityMitigationMultiplier)) {
      stop("communityMitigationMultiplier must be specified when using community mitigation.", 
           call. = FALSE)
    }
    if ( is.vector(communityMitigationMultiplier, mode = "numeric") ){
      checkNonNegative(communityMitigationMultiplier)
      
      if ( !(length(communityMitigationMultiplier) %in% c(1, length(communityMitigationStartDay))) ){
        stop("If communityMitigationMultiplier is the same for all subgroups it must either be:
              (i) a scalar; or,
              (ii) an atomic vector with length equal to the number of mitigation rounds.", 
             call. = FALSE)   
      }
    }
    if ( is.matrix(communityMitigationMultiplier) ){
      checkNonNegative(communityMitigationMultiplier)
      
      if ( !all(dim(communityMitigationMultiplier) == dim(contactMatrix)) ) {
        stop("If communityMitigationMultiplier is a matrix, then its dimensions
             must match those of the contact matrix.", 
             call. = FALSE)
      }  
    }
    if ( is.list(communityMitigationMultiplier) ){
      if ( !all(lapply(communityMitigationMultiplier, class) == "matrix") ){
        stop("If communityMitigationMultiplier is a list, then its individual elements must be matrices
              with dimensions matching those of the contact matrix.", 
             call. = FALSE)
      } else {
        sapply(communityMitigationMultiplier, checkNonNegative)
        
        if ( !all(sapply(communityMitigationMultiplier, function(l) { all(dim(l) == dim(contactMatrix)) })) ){
          stop("If communityMitigationMultiplier is a list, then its individual elements must be matrices
              with dimensions matching those of the contact matrix.", 
               call. = FALSE)
        } 
      }
    }

  }
  
  #Collect and return the parameters
  parameters <- list(population = population,
                     populationFractions = populationFractions,
                     contactMatrix = contactMatrix,
                     beta = R0 / infectiousPeriod / max(Mod(eigen(contactMatrix, symmetric = FALSE, only.values = TRUE)$values)),
                     gamma = 1 / infectiousPeriod,
                     lambda = 1 / latentPeriod,
                     seedInfections = seedInfections,
                     priorImmunity = priorImmunity,
                     useCommunityMitigation = useCommunityMitigation,
                     simulationLength = simulationLength,
                     seedStartDay = seedStartDay,
                     tolerance = tolerance,
                     method = method)
  
  if (useCommunityMitigation)
  {
    parameters <- append(parameters,
                         list(communityMitigationStartDay = communityMitigationStartDay, 
                              communityMitigationEndDay = communityMitigationStartDay + communityMitigationDuration,
                              communityMitigationMultiplier = communityMitigationMultiplier))
  }
  return(parameters)
}

#This is a utility function that reconstructs the model state as a list so that equations can refer to compartments by name
reconstructState.SEIR <- function(state)
{
  numberOfClasses <- length(state)/4 #Each of the 4 classes are vectors of the same length
  S  <- state[                       1 :     numberOfClasses ]
  E  <- state[    (numberOfClasses + 1):(2 * numberOfClasses)]
  I  <- state[(2 * numberOfClasses + 1):(3 * numberOfClasses)]
  R  <- state[(3 * numberOfClasses + 1):(4 * numberOfClasses)]
  return(as.list(environment()))
}

#This function implements the multivariate derivative of the SEIR model
#parameters should define populationFractions, normalizedContactMatrix, beta, lambda, and gamma
#Note that the total population is normalized to be 1
getDerivative.SEIR <- function(t,
                               state,
                               parameters)
{
  stateList <- reconstructState.SEIR(state)
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
    forceOfInfection <- beta / populationFractions * (contactMatrix %*% I)
    
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

#This function implements seeding infections in the SEIR model
#parameters should define seedInfections, lambda, and gamma
#Note that the total population is normalized to be 1
doSeed.SEIR <- function(state, parameters)
{
  stateList <- reconstructState.SEIR(state)
  with(append(stateList, parameters), {
    seedInfectionsFractions <- seedInfections / population
    S <- S - seedInfectionsFractions
    E <- E + seedInfectionsFractions / (1 + lambda / gamma)
    I <- I + seedInfectionsFractions / (1 + gamma / lambda)
    #Return derivative
    return(c(S, E, I, R))
  })
}
