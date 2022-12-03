library(tidyverse)
library(magrittr)
library(purrr)
library(flumodels)

R0_list <- c(1.5,2.0,3,5,7,8,10,15,25)

model_list <- list()
qq <- 0L
latentPeriod_list <- c(1,2,3,5,7,10,12,14)
infectiousPeriod_list <- c(1,2,3,5,7,10,12,14)
par_tibble_list <- list()
for(R0_value in R0_list)
{
  for(latentPeriod_value in latentPeriod_list)
  {
    for(infectiousPeriod_value in infectiousPeriod_list)
    {
      qq <- qq + 1L
      model <- SEIRModel(R0 = R0_value,
                         latentPeriod = latentPeriod_value,
                         infectiousPeriod = infectiousPeriod_value,
                         seedInfections = 0.0001,
                         simulationLength = 365)
      
      #print(model$rawOutput[100:150,])


      par_tibble_list[[qq]] <- tibble(index = qq,
                                      R0 = R0_value,
                                      latentPeriod = latentPeriod_value,
                                      infectiousPeriod = infectiousPeriod_value,
                                      model= list(model))
      print(qq)
    }
  }

  
}
par_tibble_list <- par_tibble_list %>% bind_rows()

write_rds(x = par_tibble_list,file = "sweeps/sweep_SEIR_RDS_November_30_2022/model_data.rds")