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
      print(qq)
    }
  }

  
}
par_tibble <- par_tibble_list %>% bind_rows()
#write_csv(x = par_tibble,file = "CSV_par_tibble_sweep_1_R0_restricted_2 to 10_November_30_2022/par_tibble.csv")
write_rds(x = par_tibble,file = "sweep_1_RDS_mark2_November_30_2022/model_data.rds")