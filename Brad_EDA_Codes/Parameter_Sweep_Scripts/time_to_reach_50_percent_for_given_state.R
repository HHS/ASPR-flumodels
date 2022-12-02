time_to_reach_50_percent_for_given_state <- function(x,state_variable_name)
{
  require(tidyverse)
  require(magrittr)
  require(modelbased)
  
  df <- x %>% as_tibble %>% pull(model) %>% unlist(recursive = FALSE)
  df <- df$rawOutput
  
  df %<>% as_tibble()
  state_vector <- df %>% pull(state_variable_name)
  half_index <- zero_crossings(state_vector - 0.5) %>% floor() %>% as.integer()
  time_to_half_name <- sprintf("t_%s_half",state_variable_name)
  print(half_index)
  if(any(is.na(half_index)))
  {
    x[[time_to_half_name]]<- NA
    return(x)
  }
  time_vector <- df %>% pull(time)

  x[[time_to_half_name]]<- min(time_vector[half_index])
  return(x)
}
