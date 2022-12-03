library(flumodels)
library(magrittr)
library(tidyverse)

source("convert_model_raw_output_to_tibble.R")
source("find_peak_of_state.R")
source("time_to_reach_50_percent_for_given_state.R")

x <- readRDS("sweeps/sweep_SEIRV_RDS_November_30_2022/model_data.rds")
x %<>% as_tibble()
x %<>% group_by(index)
x <- x %>% group_split()

qq <- 0L
my_list <- list()
for(xrow in x)
{
  qq <-qq + 1L
  xrow %<>% convert_model_raw_output_to_tibble()
  state_names <- xrow$tidy_raw_output[[1]] %>% pull(state_variable) %>% as.character() %>% unique()
  
  for(state_name in state_names)
  {
    xrow %<>% find_peak_of_state(state_variable_name = state_name)
    xrow %<>% time_to_reach_50_percent_for_given_state(state_variable_name = state_name)
  }
my_list[[qq]] <- xrow
}

my_list %<>% bind_rows()