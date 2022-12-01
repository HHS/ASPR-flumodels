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
      
      tibble_transformation <- model[["rawOutput"]] %>% as_tibble() %>% #Original model output is of class 'deSolve'; not useful for ggplot2 plotting.
        mutate(across(time:R,.fns = as.double)) %>% #Change all to doubles.
        pivot_longer(cols = S:R,
                     values_to = "value",
                     names_to = "state_variable") %>% #Now groupable for easy plotting
        arrange(state_variable) %>%
        mutate(state_variable = as.factor(state_variable)) #Change from character to factor() to work with ggplot2.
      
      # p1 <- ggplot(data = tibble_transformation,mapping = aes(x = time,y = value,color = state_variable,linetype=state_variable))
      # p1 <- p1 + theme_bw()
      # p1 <- p1 + scale_x_continuous(limits = c(0,240),breaks = seq(from = 0,to = 240,by = 20),minor_breaks = seq(from = 0, to = 240,by = 4))
      # p1 <- p1 + scale_y_continuous(limits = c(0,1),breaks = seq(from = 0,to = 1,by = 0.2),minor_breaks = seq(from = 0, to = 1,0.04))
      # p1 <- p1 + geom_line(alpha=0.6,size=2)
      # LINES <- c("S" = "solid","E" = "longdash","I" = "twodash","R" = "dotted")
      # 
      # p1 <- p1 + scale_linetype_manual(values = LINES)
      # cols <- c("S" = "red", "E" = "blue", "I" = "green", "R" = "magenta")
      # p1 <- p1 + scale_colour_manual(values = cols)
      # 
      # 
      # p1 <- p1 + theme(legend.position = "top")
      # 
      # title_string <- sprintf("SEIRModel, R0 = %0.5g, latentPeriod = %0.5g, infectiousPeriod = %0.5g",R0_value,latentPeriod_value,infectiousPeriod_value)
      # p1 <- p1 + ggtitle(title_string)
      # p1 <- p1 + xlab("Time [days]")
      # p1 <- p1 + ylab("Fraction in State S, E, I, or R")
      # print(p1)
      # output_filename <- sprintf("sweep_1_R0_restricted_2 to 10_November_30_2022/SEIRModel_%d.png",qq)
      # ggsave(filename = output_filename,device = "png",dpi = 300,width = 11,height = 8.5)
      
      #write_csv(x = tibble_transformation,file = sprintf("CSV_sweep_1_R0_restricted_2 to 10_November_30_2022/SEIRModel_%d.csv",qq))
      par_tibble_list[[qq]] <- tibble(index = qq,R0 = R0_value,latentPeriod = latentPeriod_value,infectiousPeriod = infectiousPeriod_value,model= list(model))
      print(qq)
    }
  }

  
}
par_tibble <- par_tibble_list %>% bind_rows()
#write_csv(x = par_tibble,file = "CSV_par_tibble_sweep_1_R0_restricted_2 to 10_November_30_2022/par_tibble.csv")
write_rds(x = par_tibble,file = "sweep_1_RDS_mark2_November_30_2022/model_data.rds")