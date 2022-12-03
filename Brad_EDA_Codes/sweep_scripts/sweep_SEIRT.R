library(tidyverse)
library(magrittr)
library(purrr)
library(flumodels)

R0_list <- c(1,1.25,1.5,2.0,3,5,7,8,10,25,50)

model_list <- list()
qq <- 0L
latentPeriod_list <- c(0.5,1,2,3,5,7,14)
infectiousPeriod_list <- c(0.5,1,2,3,5,7,14)
AVEi_list <- c(0.01,0.1,0.2,0.5,0.75,0.9,0.99)
AVEp_list <- c(0.01,0.1,0.2,0.5,0.75,0.9,0.99)
par_tibble_list <- list()
for(R0_value in R0_list)
{
  for(latentPeriod_value in latentPeriod_list)
  {
    for(infectiousPeriod_value in infectiousPeriod_list)
    {
      for(AVEi_value in AVEi_list)
      {
        for(AVEp_value in AVEp_list)
        {
          qq <- qq + 1L
          model <- SEIRTModel(R0 = R0_value,
                              latentPeriod = latentPeriod_value,
                              infectiousPeriod = infectiousPeriod_value,
                              seedInfections = 0.0001,
                              AVEi = AVEi_value,
                              AVEp = AVEp_value,
                              simulationLength = 365)
          # tibble_transformation <- model[["rawOutput"]] %>% as_tibble() %>% #Original model output is of class 'deSolve'; not useful for ggplot2 plotting.
          #   mutate(across(time:R,.fns = as.double)) %>% #Change all to doubles.
          #   pivot_longer(cols = S:R,
          #                values_to = "value",
          #                names_to = "state_variable") %>% #Now groupable for easy plotting
          #   arrange(state_variable) %>%
          #   mutate(state_variable = as.factor(state_variable)) #Change from character to factor() to work with ggplot2.
          
          par_tibble_list[[qq]] <- tibble(index = qq,
                                          R0 = R0_value,
                                          latentPeriod = latentPeriod_value,
                                          infectiousPeriod = infectiousPeriod_value,
                                          AVEi = AVEi_value,
                                          AVEp = AVEp_value,
                                          model= list(model))
          print(qq)
        }
      }

    }
  }
}
par_tibble <- par_tibble_list %>% bind_rows()

write_rds(x = par_tibble,file = "sweeps/sweep_SEIRT_RDS_November_30_2022/model_data.rds")