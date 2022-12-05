#' @name plot
#' @title Plot
#' @description Plot the model
#' @details This is the generic function, for class-specific versions see:
#' \itemize{
#'  \item \code{\link{plot.SEIRModel}} for SEIR-type models
#'  \item \code{\link{plot.SEIRTModel}} for SEIRT-type models
#'  \item \code{\link{plot.SEAIRTVModel}} for SEAIRTV-type models
#' }
#' @usage plot(x, y, ...)
#' @param x The model to plot
#' @param y Unused
#' @param ... Other parameters for the class-specific functions
#' @rdname plot
NULL

#' @title Plot
#' @description Plot the model
#' @param x The model to plot
#' @param y Unused
#' @param asRate Whether to plot rates (fractions of population) or else
#'   numbers; defaults to FALSE unless model population is 1
#' @param incidence If true, plots incidence, otherwise plots prevalence;
#'   defaults to FALSE
#' @param symptomatic Whether or not to only report symptomatic infections; defaults to FALSE
#' @param fractionSymptomatic The fraction of cases that are symptomatic; must be specified if symptomatic = TRUE
#' @param normalizeRate Whether or not to normalize rates by population;
#'   defaults to FALSE
#' @param populationLabels Labels for the population groups; optional
#' @method plot SEIRModel
#' @importFrom graphics axis legend lines plot title points mtext
#' @keywords internal
#' @export
plot.SEIRModel <- function(x, y, ..., asRate = FALSE, incidence = FALSE, symptomatic = FALSE, fractionSymptomatic,
                           normalizeRate = FALSE, populationLabels) {
  if (x$parameters$population == 1) { #If population is 1, assume that all numbers are rates
    asRate <- TRUE
  }
  if (symptomatic) {
    if (missing(fractionSymptomatic)) {
      stop("fractionSymptomatic must be specified.")
    }
    infectionsByGroup <- getInfectionTimeSeries(x, byGroup = TRUE, asRate = asRate, incidence = incidence,
                                                symptomatic = symptomatic, fractionSymptomatic = fractionSymptomatic)
  } else {
    infectionsByGroup <- getInfectionTimeSeries(x, byGroup = TRUE, asRate = asRate, incidence = incidence)
  }
  overallInfections <- rowSums(infectionsByGroup)
  if (asRate && normalizeRate) {
    infectionsByGroup <- sweep(infectionsByGroup, 2, x$parameters$populationFractions, "/")
  }
  times <- x$rawOutput[, "time"]
  scalePoints <- pretty(c(0, max(overallInfections, infectionsByGroup, na.rm = TRUE)))
  ylim <- c(0, max(scalePoints))
  graphics::plot(times, overallInfections, 
       type = 'l', xlab = "Time (days)", ylab = "", ylim = ylim, yaxt = "n", lwd = 3)
  makeLegend <- FALSE
  legendLabels <- NULL
  legendCol <- NULL
  legendLty <- NULL
  legendPch <- NULL
  if (incidence) {
    mainTitle <- "Infection Incidence"
    if (asRate) {
      seedInfections <- sum(x$parameters$seedInfections)
    } else {
      seedInfections <- sum(x$parameters$seedInfections) * x$parameters$population
    }
    graphics::points(x$parameters$seedStartDay, seedInfections,
           type = "p", pch = 20, cex = 1)
    makeLegend <- TRUE
    legendLabels <- "Seed Cases"
    legendCol <- 1
    legendLty <- NA
    legendPch <- 20
  } else {
    mainTitle <- "Infection Prevalence"
  }
  if (symptomatic) {
    mainTitle <- paste("Symptomatic", mainTitle);
  }
  graphics::title(main = mainTitle)
  #Plot population group data
  if (ncol(infectionsByGroup) > 1) {
    for (i in 1:ncol(infectionsByGroup)) {
      graphics::lines(times, infectionsByGroup[, i], col = i + 1, lwd = 3)
    }
    if (missing(populationLabels)) {
      legendLabels <- append(append("Overall", paste0("Group ", 1:ncol(infectionsByGroup))),
                             legendLabels)
    } else {
      legendLabels <- append("Overall", populationLabels)
    }
    makeLegend <- TRUE
    legendCol <- append(1:(ncol(infectionsByGroup) + 1), legendCol)
    legendLty <- append(rep(1, ncol(infectionsByGroup) + 1), legendLty)
    legendPch <- append(rep(NA, ncol(infectionsByGroup) + 1), legendPch)
  }
  #Scale and label y-axis
  if (asRate) {
    graphics::axis(2, at = scalePoints, las = 1,
         labels = paste0(scalePoints * 100, "%"))
  } else {
    exponent <- floor(log10(ylim[2]) / 3)
    if (exponent < 2) { #Hundreds of thousands or fewer
      graphics::axis(2, at = scalePoints, las = 1, labels = format(scalePoints, big.mark=",", scientific=FALSE))
    } else if (exponent < 3) { #Millions
      graphics::axis(2, at = scalePoints, las = 1, labels = scalePoints / 1e6)
      graphics::mtext("Millions", line = 3, side = 2, adj = 1)
    } else if (exponent < 4) { #Billions
      graphics::axis(2, at = scalePoints, las = 1, labels = scalePoints / 1e9)
      graphics::mtext("Billions", line = 3, side = 2, adj = 1)
    } else {
      graphics::axis(2, at = scalePoints, las = 1, labels = format(scalePoints, scientific=TRUE))
    }
  }
  #Make legend
  if (makeLegend) {
    graphics::legend("topright", legend = legendLabels,
           col = legendCol, lty = legendLty, pch = legendPch,
           bty = "n", lwd = 1)
  }
}

#' @title Plot
#' @description Plot the model
#' @param x The model to plot
#' @param y Unused
#' @param asRate Whether to plot rates (fractions of population) or else
#'   numbers; defaults to FALSE unless model population is 1
#' @param incidence If true, plots incidence, otherwise plots prevalence;
#'   defaults to FALSE
#' @param symptomatic Whether or not to only plot symptomatic infections; defaults to FALSE
#' @param normalizeRate Whether or not to normalize rates by population;
#'   defaults to FALSE
#' @param populationLabels Labels for the population groups; optional
#' @method plot SEIRTModel
#' @keywords internal
#' @export
plot.SEIRTModel <- function(x, ..., asRate = FALSE, incidence = FALSE, symptomatic = FALSE,
                            normalizeRate = FALSE, populationLabels) {
  plot.SEIRModel(x, asRate, incidence, symptomatic, fractionSymptomatic = x$parameters$fractionSymptomatic,
                 normalizeRate, populationLabels)
}

#' @title Plot
#' @description Plot the model
#' @param x The model to plot
#' @param y Unused
#' @param asRate Whether to plot rates (fractions of population) or else
#'   numbers; defaults to FALSE unless model population is 1
#' @param incidence If true, plots incidence, otherwise plots prevalence;
#'   defaults to FALSE
#' @param symptomatic Whether or not to only plot symptomatic infections; defaults to FALSE
#' @param normalizeRate Whether or not to normalize rates by population;
#'   defaults to FALSE
#' @param populationLabels Labels for the population groups; optional
#' @method plot SEAIRTVModel
#' @keywords internal
#' @export
plot.SEAIRTVModel <- function(x, ..., asRate = FALSE, incidence = FALSE, symptomatic = FALSE,
                            normalizeRate = FALSE, populationLabels) {
  plot.SEIRTModel(x, ..., asRate = FALSE, incidence = FALSE, symptomatic = FALSE,
                                normalizeRate = FALSE, populationLabels)
}