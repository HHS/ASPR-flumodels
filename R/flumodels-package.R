#' @title Deterministic compartmental models for influenza with mitigations.
#' @description Provides several compartmental models for influenza that enable
#' the consideration of various mitigation strategies including vaccines and
#' antiviral drugs.
#' @details \tabular{ll}{
#' Package: \tab flumodels \cr
#' Type: \tab Package \cr
#' Version: \tab 1.0.9 \cr
#' Date: \tab 2019-06-12 \cr
#' License: \tab BSD-2-Clause \cr
#' }
#' @name flumodels
#' @aliases flumodels
#' @docType package
#' @author Jason Asher \email{jason.m.asher@@gmail.com}
#' @import deSolve
#' @keywords package
#' @examples
#' model <- SEIRModel(R0 = 1.3,
#'                    latentPeriod = 1.5,
#'                    infectiousPeriod = 2.5,
#'                    seedInfections = 0.0001 * 310e6,
#'                    population = 310e6)
#' plot(model)
#'
#' model2 <- SEIRModel(R0 = 1.3,
#'                     latentPeriod = 1.5,
#'                     infectiousPeriod = 2.5,
#'                     seedInfections = 0.0001 * 310e6,
#'                     population = 310e6,
#'                     populationFractions = c(0.25, 0.75),
#'                     contactMatrix = matrix(c(18,  3,
#'                                               9, 12), ncol=2, byrow = TRUE))
#' plot(model2)
NULL
