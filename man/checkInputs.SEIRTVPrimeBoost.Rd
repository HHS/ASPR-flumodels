% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEIRTVPrimeBoostModel.R
\name{checkInputs.SEIRTVPrimeBoost}
\alias{checkInputs.SEIRTVPrimeBoost}
\title{Check SEIR+V prime-boost inputs}
\usage{
checkInputs.SEIRTVPrimeBoost(
  population,
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
  fractionSymptomatic,
  fractionSeekCare,
  fractionDiagnosedAndPrescribedOutpatient,
  fractionAdhere,
  fractionAdmitted,
  fractionDiagnosedAndPrescribedInpatient,
  AVEi,
  AVEp,
  vaccineAdministrationRatePerDay,
  vaccineAvailabilityByDayPrime,
  vaccineAvailabilityByDayBoost,
  vaccineUptakeMultiplier,
  boostDelay,
  VEs1,
  VEs2,
  VEi1,
  VEi2,
  VEp1,
  VEp2,
  vaccineEfficacyDelay,
  simulationLength,
  seedStartDay,
  tolerance,
  method
)
}
\value{
List of parameters for the SEIR+V prime-boost model
}
\description{
Checks the input parameters for the SEIR+V prime-boost model
}
\keyword{internal}
