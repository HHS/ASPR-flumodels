#' @title Print
#' @description Prints a brief summary of a SEIRModel object.
#' @rdname print
#' @method print SEIRModel
#' @keywords internal
#' @export
print.SEIRModel <- function(x, ...) {
  R0 <- x$parameters$beta / x$parameters$gamma * max(Mod(eigen(x$parameters$contactMatrix)$values))
  cat(class(x)[1], "\n")
  latentPeriod <- 1 / x$parameters$lambda
  recoveryPeriod <- 1 / x$parameters$gamma
  generationTime <- latentPeriod + recoveryPeriod
  if (length(x$parameters$populationFractions) > 1) {
    cat("The model has ", length(x$parameters$populationFractions), "population groups\n" )
  }
  cat("R0:", R0, "\n")
  cat("Generation time:", generationTime, "days (latent period:",
      latentPeriod, "days, recovery period:", recoveryPeriod, "days)\n")
  cat("Serologic attack rate:", paste0(round(100*getInfections(x, byGroup = FALSE, asRate = TRUE), 2), "%"), "\n")
}