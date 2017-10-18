#' @title Integrate model
#' @description Calculates and returns the raw output from the deSolve ODE integration
#' @param initialState The initial state of the model's ODE system prior to seeding
#' @param parameters The model's parameters
#' @param derivativeFunction A function the computes the derivatives for the model's ODE system
#' @param seedFunction A function that modifies the state to seed infections
#' @param method Which integration method to use. Defaults to lsoda
#' @return A matrix of deSolve ode output for the model
#' @importFrom deSolve ode
#' @keywords internal
integrateModel <- function(initialState, parameters, derivativeFunction, seedFunction, method = "lsoda") {
  with(parameters, {
    #Integrate the model, being mindful of the seedStartDay
    if (parameters$seedStartDay > 0) { #If seeding is not immediate
      #Compute first seedStartDay's worth of data
      times <- seq(0, parameters$seedStartDay, by = 1)
      preSeedData <- as.matrix(deSolve::ode(y = initialState, times = times, func = derivativeFunction, parms = parameters,
                                   atol = tolerance, method = method))
      
      #Next, integrate model starting from seed date
      initialState <- seedFunction(preSeedData[nrow(preSeedData), -1], parameters) #Copy state from pre-seed data
      times <- seq(parameters$seedStartDay, parameters$seedStartDay + parameters$simulationLength, by = 1)
      postSeedData <- as.matrix(deSolve::ode(y = initialState, times = times, func = derivativeFunction, parms = parameters,
                                    atol = parameters$tolerance, method = method))
      
      #Bind data together and return
      return(rbind(preSeedData[-nrow(preSeedData), ], postSeedData)) #postSeedData contains the data for day seedStartDay
      
    } else { #Seeding is immediate
      #Set up the initial state
      initialState <- seedFunction(initialState, parameters)
      times <- seq(0, parameters$simulationLength, by = 1)
      #Integrate and return the model output
      return(as.matrix(deSolve::ode(y = initialState, times = times, func = derivativeFunction, parms = parameters,
                           atol = parameters$tolerance, method = method)))
    }
  })
}