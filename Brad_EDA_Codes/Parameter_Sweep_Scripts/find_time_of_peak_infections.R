find_time_of_peak_for_given_state <- function(x,state_variable_name)
{
  require(tidyverse)
  require(magrittr)
  
  df <- x %>% as_tibble %>% pull(model) %>% unlist(recursive = FALSE)
  df <- df$rawOutput
  
  df %<>% as_tibble()
  state_vector <- df %>% pull(state_variable_name)
  max_index <- which.max(state_vector)
  time_vector <- df %>% pull(time)
  time_to_max_name <- sprintf("t_%s_max",state_variable_name)
  x[[time_to_max_name]]<- time_vector[max_index]
  return(x)
}
