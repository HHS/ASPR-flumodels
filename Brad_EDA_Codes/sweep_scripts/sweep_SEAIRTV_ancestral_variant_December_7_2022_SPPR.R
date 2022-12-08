library(tidyverse)
library(magrittr)
library(purrr)
library(flumodels)
library(flumodelsutil)
list.files("../../R/",full.names = TRUE) %>% sapply(FUN = source)
#Get the census data first.
k1 <- read_csv(file = "../../Input_Data/fipsCensusSingleYearData/county_population_by_age_StandardCovidVaccineAgeGroups_2022-06-09.csv",
               col_types = "cfi")

k1 %<>% mutate(age_group = fct_recode(age_group,
                                      `12-17` = "12-15",
                                      `12-17` = "16-17",
                                      `18-49` = "18-24",
                                      `18-49` = "25-39",
                                      `18-49` = "40-49",
                                      `50-64` = "50-64",
                                      `65+` = "65-74",
                                      `65+` = "75-99")) %>% group_by(fips_code,age_group) %>% summarize(population = sum(population))

load("../../Input_Data/stateFipsCrosswalk/stateFipsCrosswalk.RData")

g1 <- stateFipsCrosswalk
g1 %<>% mutate(fips_code = as.integer(fips_code))
g1 %<>% mutate(fips_code = as.character(fips_code))
states <- g1$state
names(states) <- (stateFipsCrosswalk$fips_code %>% as.integer %>% as.character)
k1$state <- states[k1$fips_code]

k1 %<>% ungroup()
k1 %<>% group_by(age_group,state)
k1 %<>% summarize(population = sum(population))

k1 %<>% ungroup
k1 %<>% rename(Location = state)
k1 %<>% rename(AgeRange = age_group)
k1 %<>% mutate(AgeRange = as.ordered(AgeRange))
k1 %<>% ungroup()

pop_by_age_bracket <- aggregate(k1$population,by = list(AgeRange = k1$AgeRange),FUN = sum)
names(pop_by_age_bracket) <- c("AgeRange","Population")
pop_by_age_bracket$Fraction <- pop_by_age_bracket$Population/sum(pop_by_age_bracket$Population)
pop_by_age_bracket %<>% as_tibble()

model_list <- list()

population_list <- c(3.3e8)
populationFraction_list <- pop_by_age_bracket %>% pull(Fraction) %>% list()
R0_list <- c(2.5)
latentPeriod_list <- 5.5
infectiousPeriod_list <- c(3.0)
fractionLatentThatIsInfectious_list <- c(0.4)
relativeInfectivityAsymptomatic_list <- list(rep(0.75,6L))
seedInfections_list <- c(0.0001)
priorImmunity_list <- c(0)
useCommunityMitigation_list <- TRUE
communityMitigationStartDay_list <- c(14)
communityMitigationDuration_list <- c(0.75)
communityMitigationMultiplier_list <- c(0.75)

fractionSymptomatic_list <- list(c(0.61, #0-4
                                   0.61, #5-11
                                   0.71, #12-17
                                   0.77, #18-49
                                   0.77, #50-64
                                   0.77)) #65+

fractionSeekCare_list <- c(0.6)
fractionDiagnosedAndPrescribedOutpatient <- c(0.4)
fractionAdhere_list <- c(1.0)
fractionAdmittted <- c(1.0)
fractionDiagnosedAndPrescribedInpatient <- c(0.4)
AVEi_list <- c(0.0)
AVEp_list <- c(0.7)
vaccineAdministrationRatePerDay_list <- c(0.0)
vaccineAvailabilityByDay_list <- c(3.3e8)
VEs_list <- c(0.9)
VEi_list <- c(0.0)
VEp_list <- c(0.0)
vaccineEfficacyDelay_list <- c(14)
simulationLength <- c(360)
seedStartDay <- 0
tolerance <- 1e-8
method <- "default"

parameter_frame <- expand.grid(population_list,
                               populationFraction_list,
                               R0_list,
                               latentPeriod_list,
                               infectiousPeriod_list,
                               fractionLatentThatIsInfectious_list,
                               relativeInfectivityAsymptomatic_list,
                               seedInfections_list,
                               priorImmunity_list,
                               useCommunityMitigation_list,
                               communityMitigationStartDay_list,
                               communityMitigationDuration_list,
                               communityMitigationMultiplier_list,
                               fractionAdhere_list,
                               fractionSymptomatic_list,
                               fractionSeekCare_list,
                               fractionDiagnosedAndPrescribedOutpatient,
                               AVEi_list,
                               AVEp_list,
                               vaccineAdministrationRatePerDay_list,
                               vaccineAvailabilityByDay_list,
                               VEs_list,
                               VEi_list,
                               VEp_list,
                               vaccineEfficacyDelay_list) %>% as_tibble()

names(parameter_frame) <- c("population",
                            "populationFraction",
                            "R0",
                            "latentPeriod",
                            "infectiousPeriod",
                            "fractionLatentThatIsInfectious",
                            "relativeInfectivityAsymptomatic",
                            "seedInfections",
                            "priorImmunity",
                            "useCommunityMitigation",
                            "communityMitigationStartDay",
                            "communityMitigationDuration",
                            "communityMitigationMultiplier",
                            "fractionAdhere",
                            "fractionSymptomatic",
                            "fractionSeekCare",
                            "fractionDiagnosedAndPrescribedOutpatient",
                            "AVEi",
                            "AVEp",
                            "vaccineAdministrationRatePerDay",
                            "vaccineAvailabilityByDay",
                            "VEs",
                            "VEi",
                            "VEp",
                            "vaccineEfficacyDelay")

num_rows <- nrow(parameter_frame)
model_list <- list()
for(i in 1:num_rows)
{
  print(sprintf("Running sweep %d out of %d.",i,num_rows))
  data_row <- parameter_frame[i,]
  
  model <- SEAIRTVModel(population = data_row$population,
                        populationFraction = data_row$populationFraction %>% unlist(),
                        contactMatrix = makeContactMatrix(ages = c(4,11,17,49,65)),
                        R0 = data_row$R0,
                        latentPeriod = data_row$latentPeriod,
                        infectiousPeriod = data_row$infectiousPeriod,
                        fractionLatentThatIsInfectious = data_row$fractionLatentThatIsInfectious,
                        relativeInfectivityAsymptomatic = data_row$relativeInfectivityAsymptomatic %>% unlist(),
                        seedInfections = data_row$seedInfections,
                        priorImmunity = data_row$priorImmunity,
                        fractionAdhere = data_row$fractionAdhere,
                        useCommunityMitigation = data_row$useCommunityMitigation,
                        communityMitigationStartDay = data_row$communityMitigationStartDay,
                        communityMitigationDuration = data_row$communityMitigationDuration,
                        communityMitigationMultiplier = data_row$communityMitigationMultiplier,
                        fractionSymptomatic = data_row$fractionSymptomatic %>% unlist(),
                        fractionSeekCare = data_row$fractionSeekCare,
                        fractionDiagnosedAndPrescribedOutpatient  = data_row$fractionDiagnosedAndPrescribedOutpatient,
                        AVEi = data_row$AVEi,
                        AVEp = data_row$AVEp,
                        vaccineAdministrationRatePerDay = data_row$vaccineAdministrationRatePerDay,
                        vaccineAvailabilityByDay = data_row$vaccineAvailabilityByDay,
                        VEs = data_row$VEs,
                        VEi = data_row$VEi,
                        VEp = data_row$VEp,
                        vaccineEfficacyDelay = data_row$vaccineEfficacyDelay,
                        simulationLength = simulationLength,
                        seedStartDay  = seedStartDay,
                        tolerance  = 1e-8,
                        method = "lsoda")

model_list[[i]] <- model
}

full_tibble <- bind_cols(parameter_frame,tibble(model = model_list),index = 1:num_rows)
write_rds(x = full_tibble,file = "sweeps/sweep_SEIRV_RDS_November_30_2022/model_data.rds")