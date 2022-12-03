library(tidyverse)
library(magrittr)
library(purrr)
library(flumodels)

model_list <- list()
R0_list <- c(2,5,10)
latentPeriod_list <- c(7,14)
infectiousPeriod_list <- c(7,14)
VEs1_list <- c(0.1,0.5,0.9)
VEi1_list <- c(0.1,0.5,0.9)
VEp1_list <- c(0.1,0.5,0.9)
VEs2_list <- c(0.1,0.5,0.9)
VEi2_list <- c(0.1,0.5,0.9)
VEp2_list <- c(0.1,0.5,0.9)
vaccineEfficacyDelay_list <- c(7,14,30) #1 week, 2 weeks, 1 month)
boostDelay_list <- c(7,14,30)

par_grid <- expand.grid(R0 = R0_list,
                        latentPeriod_list,
                        infectiousPeriod_list,
                        VEs1_list,
                        VEs2_list,
                        VEi1_list,
                        VEi2_list,
                        VEp1_list,
                        VEp2_list,
                        vaccineEfficacyDelay_list,
                        boostDelay_list)

par_grid %<>% as_tibble()
names(par_grid) <- c("R0",
                     "latentPeriod",
                     "infectiousPeriod",
                     "VEs1",
                     "VEs2",
                     "VEi1",
                     "VEi2",
                     "VEp1",
                     "VEp2",
                     "vaccineEfficacyDelay",
                     "boostDelay")

num_rows <- nrow(par_grid)
par_tibble_list <- list()


for(i in 1:num_rows)
{
  row_selected <- par_grid[i,]
  
  R0_value <- row_selected$R0
  latentPeriod_value <- row_selected$latentPeriod
  infectiousPeriod_value <- row_selected$infectiousPeriod
  VEs1_value <- row_selected$VEs1
  VEs2_value <- row_selected$VEs2
  VEi1_value <- row_selected$VEi1
  VEi2_value <- row_selected$VEi2
  VEp1_value <- row_selected$VEp1
  VEp2_value <- row_selected$VEp2
  vaccineEfficacyDelay_value <- row_selected$vaccineEfficacyDelay
  boostDelay_value <- row_selected$boostDelay
  
  
  model <- SEIRVPrimeBoostModel(R0 = R0_value,
                                latentPeriod = latentPeriod_value,
                                infectiousPeriod = infectiousPeriod_value,
                                population = 330e6,
                                seedInfections = 10000,
                                VEs1 = VEs1_value,
                                VEs2 = VEs2_value,
                                VEi1 = VEi1_value,
                                VEi2 = VEi2_value,
                                VEp1 = VEp1_value,
                                VEp2 = VEp2_value,
                                boostDelay = boostDelay_value,
                                vaccineEfficacyDelay = vaccineEfficacyDelay_value,
                                simulationLength = 365)
  
  par_tibble_list[[i]] <- tibble(index = i,
                                 R0 = R0_value,
                                 latentPeriod = latentPeriod_value,
                                 infectiousPeriod = infectiousPeriod_value,
                                 VEs1 = VEs1_value,
                                 VEs2 = VEs2_value,
                                 VEi1 = VEi1_value,
                                 VEi2 = VEi2_value,
                                 VEp1 = VEp1_value,
                                 VEp2 = VEp2_value,
                                 vaccineEfficacyDelay = vaccineEfficacyDelay_value,
                                 boostDelay = boostDelay_value,
                                 model= list(model))
  print(i)
}
par_tibble <- par_tibble_list %>% bind_rows()

write_rds(x = par_tibble,file = "sweeps/sweep_SEIRVPrimeBoost_RDS_November_30_2022/model_data.rds")