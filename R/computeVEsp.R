#' @title computeVEp
#' @description Compute VEp (Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals) for a given VEsp (computed VEp + VEs). The presumed relationship is:
#'   (1 - VEsp) = (1 - VEp) * (1 - VEs).
#' @param VEsp Vaccine efficacy for combined protection of and prevention of symptomatic illness in vaccinated individuals.
#'    A single fraction or vector of fractions by population group.
#' @param VEsRatio Ratio of VEs / VEp. VEsRatio = -1 is equivilant to VEsRatio = infinity (VEp = 0). One of VEsRatio, VEpRatio,
#'    and VEs must be specified.
#' @param VEpRatio of VEp / VEs. VEpRatio = -1 is equivilant to VEpRatio = infinity (VEs = 0). One of VEsRatio, VEpRatio,
#'    and VEs must be specified.
#' @param VEs Vaccine efficacy: protection for vaccinated susceptible
#'   individuals; single fraction or vector of fractions by population group;
#'   One of VEsRatio, VEpRatio, and VEs must be specified.
#' @return A numeric or vector of numerics of same length as VEsp.
#' @author Matthew Clay <clay.matt@gmail.com>
#' @export
computeVEp <- function(VEsp, VEsRatio, VEpRatio, VEs) {

  checkBetween0and1(VEsp)

  if (missing(VEsRatio) & missing(VEpRatio) & missing(VEs))
    stop("One of VEsRatio, VEpRatio, and VEs must be specified")

  if (length(match.call()[-1]) > 2)
    stop("Only one of VEsRatio, VEpRatio, and VEs should be specified")

  if (!missing(VEsRatio)) {
    if ( VEsRatio < 0 & VEsRatio != -1 )
      stop("VEsRatio must be > 0 (or equal to -1 to represent infinity)")
    if (VEsRatio == -1)
      return( rep(0, length(VEsp) ))
    else if (VEsRatio == 0)
      return (VEsp)
    return( ( ( VEsRatio + 1 ) - sqrt( 1 +2 * VEsRatio + VEsRatio^2 - 4 * VEsp * VEsRatio ) ) /
              (2 * VEsRatio) )
  }

  if (!missing(VEpRatio)) {
    if ( VEpRatio < 0 & VEpRatio != -1 )
      stop("VEpRatio must be > 0 (or equal to -1 to represent infinity)")
    if (VEpRatio == -1)
      return( VEsp )
    else if (VEpRatio == 0)
      return ( rep(0, length(VEsp) ) )
    return( 0.5 * ( VEpRatio + 1  - sqrt( 1 + 2 * VEpRatio + VEpRatio^2 - 4 * VEsp * VEpRatio ) ) )
  }

  # We must be dealing with VEs specified
  checkBetween0and1(VEs)
  if (length(VEsp) != 1 & length(VEs) != 1)
    checkDimensionsMatch(VEsp, VEs)
  if ( sum(VEsp < VEs) > 0 )
    stop("VEs can not be larger than VEsp")
  return( 1 - (1 - VEsp) / (1- VEs) )
}


#' @title computeVEs
#' @description Compute VEs (Vaccine efficacy: protection for vaccinated susceptible
#'   individuals) for a given VEsp (computed VEp + VEs). The presumed relationship is:
#'   (1 - VEsp) = (1 - VEp) * (1 - VEs).
#' @param VEsp Vaccine efficacy for combined protection of and prevention of symptomatic illness in vaccinated individuals.
#'    A single fraction or vector of fractions by population group.
#' @param VEsRatio Ratio of VEs / VEp. VEsRatio = -1 is equivilant to VEsRatio = infinity (VEp = 0). One of VEsRatio, VEpRatio,
#'    and VEs must be specified.
#' @param VEpRatio Ratio of VEp / VEs. VEpRatio = -1 is equivilant to VEpRatio = infinity (VEs = 0). One of VEsRatio, VEpRatio,
#'    and VEs must be specified.
#' @param VEp Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals; single fraction or vector of fractions by population group.
#'   One of VEsRatio, VEpRatio, and VEs must be specified.
#' @return A numeric or vector of numerics of same length as VEsp.
#' @author Matthew Clay <clay.matt@gmail.com>
#' @export
computeVEs <- function(VEsp, VEsRatio, VEpRatio, VEp) {

  checkBetween0and1(VEsp)

  if (missing(VEsRatio) & missing(VEpRatio) & missing(VEp))
    stop("One of VEsRatio, VEpRatio, and VEs must be specified")

  if (length(match.call()[-1]) > 2)
    stop("Only one of VEsRatio, VEpRatio, and VEs should be specified")

  if (!missing(VEsRatio)) {
    if ( VEsRatio < 0 & VEsRatio != -1 )
      stop("VEsRatio must be > 0 (or equal to -1 to represent infinity)")
    if (VEsRatio == -1)
      return( VEsp )
    else if (VEsRatio == 0)
      return ( rep(0, length(VEsp) ) )
    return( 0.5 * ( VEsRatio + 1 - sqrt( 1 + 2 * VEsRatio + VEsRatio^2 - 4 * VEsp * VEsRatio ) ) )
  }

  if (!missing(VEpRatio)) {
    if ( VEpRatio < 0 & VEpRatio != -1 )
      stop("VEpRatio must be > 0 (or equal to -1 to represent infinity)")
    if (VEpRatio == -1)
      return( rep(0, length(VEsp) ) )
    else if (VEpRatio == 0)
      return ( VEsp )
    return( ( ( VEpRatio + 1 ) - sqrt( 1 +2 * VEpRatio + VEpRatio^2 - 4 * VEsp * VEpRatio ) ) /
              (2 * VEpRatio) )
  }

  # We must be dealing with VEp specified
  checkBetween0and1(VEp)
  if (length(VEsp) != 1 & length(VEp) != 1)
    checkDimensionsMatch(VEsp, VEp)
  if ( sum(VEsp < VEp) > 0 )
    stop("VEp can not be larger than VEsp")
  return( ifelse(VEp == VEsp, 0, 1 - (1 - VEsp) / (1- VEp) ))
}


#' @title computeVEsp
#' @description Compute VEsp (Vaccine efficacy for combined protection of and prevention of symptomatic illness in
#' vaccinated individuals) for a given VEp and VEs.
#' @param VEp Vaccine efficacy: prevention of symptomatic illness in
#'   infected indivduals; single fraction or vector of fractions by population group;
#'   defaults to 0. One of VEsRatio, VEpRatio, and VEs must be specified.
#' @param VEs Vaccine efficacy: protection for vaccinated susceptible
#'   individuals; single fraction or vector of fractions by population group;
#'   defaults to 0. One of VEsRatio, VEpRatio, and VEs must be specified.
#' @return A numeric or vector of numerics of same length as VEp & VEs.
#' @author Matthew Clay <clay.matt@gmail.com>
#' @export
computeVEsp <- function(VEp = 0, VEs = 0) {
  checkBetween0and1(VEp)
  checkBetween0and1(VEs)
  if (length(VEp) != 1 & length(VEs) != 1)
    checkDimensionsMatch(VEp, VEs)

  return ( 1 - (1 - VEp) * (1 - VEs))
}
