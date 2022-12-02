library(tidyverse)
library(magrittr)
library(purrr)
library(flumodels)

R0_list <- c(5,10,25)

model_list <- list()
qq <- 0L
latentPeriod_list <- c(7,14)
infectiousPeriod_list <- c(7,14)
VEs_list <- c(0.3,0.9)
VEi_list <- c(0.3,0.9)
VEp_list <- c(0.3,0.9)

par_tibble_list <- list()
vaccineEfficacyDelay_list <- c(7,14,30) #1 week, 2 weeks, 1 month)
for(R0_value in R0_list)
{
  for(latentPeriod_value in latentPeriod_list)
  {
    for(infectiousPeriod_value in infectiousPeriod_list)
    {
      for(VEs_value in VEs_list)
      {
        for(VEi_value in VEi_list)
        {
          for(VEp_value in VEp_list)
          {
            for(vaccineEfficacyDelay_value in vaccineEfficacyDelay_list)
            {
              qq <- qq + 1L
              model <- SEIRVModel(R0 = R0_value,
                                  latentPeriod = latentPeriod_value,
                                  infectiousPeriod = infectiousPeriod_value,
                                  population = 330e6,
                                  seedInfections = 10000,
                                  vaccineAvailabilityByDay = c(60e6, rep(0, 59), rep(15e6/7, 7*16)),
                                  vaccineAdministrationRatePerDay = 15e6/7,
                                  VEs = VEs_value,
                                  VEi = VEi_value,
                                  VEp = VEp_value,
                                  vaccineEfficacyDelay = vaccineEfficacyDelay_value,
                                  simulationLength = 365)
              
              par_tibble_list[[qq]] <- tibble(index = qq,
                                              R0 = R0_value,
                                              latentPeriod = latentPeriod_value,
                                              infectiousPeriod = infectiousPeriod_value,
                                              VEs = VEs_value,
                                              VEi = VEi_value,
                                              VEp = VEp_value,
                                              vaccineEfficacyDelay = vaccineEfficacyDelay_value,
                                              model= list(model))
              print(qq)
            }
          }
        }
      }
    }
  }
}

par_tibble <- par_tibble_list %>% bind_rows()

write_rds(x = par_tibble,file = "sweeps/sweep_SEIRV_RDS_November_30_2022/model_data.rds")