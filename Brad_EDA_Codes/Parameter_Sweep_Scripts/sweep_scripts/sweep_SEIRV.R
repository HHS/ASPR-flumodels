library(tidyverse)
library(magrittr)
library(purrr)
library(flumodels)


model_list <- list()

R0_list <- c(7,15,25)
VEs_list <- c(0.2,0.5,0.9)
VEi_list <- c(0.2,0.5,0.9)
VEp_list <- c(0.2,0.5,0.9)



parameter_frame <- expand.grid(R0_list,VEs_list,VEi_list,VEp_list) %>% as_tibble()
names(parameter_frame) <- c("R0","VEs","VEi","VEp")
num_rows <- nrow(parameter_frame)

model_list <- list()

for(i in 1:num_rows)
{
  print(sprintf("Running sweep %d out of %d.",i,num_rows))
  data_row <- parameter_frame[i,]
  R0_value <- data_row$R0
  VEs_value <- data_row$VEs
  VEi_value <- data_row$VEi
  VEp_value <- data_row$VEp

  model <- SEIRVModel(R0 = R0_value,
                      population = 330e6,
                      seedInfections = 10000,
                      latentPeriod = 7,
                      infectiousPeriod = 7,
                      vaccineAvailabilityByDay = c(60e6, rep(0, 59), rep(15e6/7, 7*16)),
                      vaccineAdministrationRatePerDay = 15e6/7,
                      VEs = VEs_value,
                      VEi = VEi_value,
                      VEp = VEp_value,
                      simulationLength = 365)
  
  model_list[[i]] <- model
}

full_tibble <- bind_cols(parameter_frame,tibble(model = model_list),index = 1:num_rows)
write_rds(x = full_tibble,file = "sweeps/sweep_SEIRV_RDS_November_30_2022/model_data.rds")