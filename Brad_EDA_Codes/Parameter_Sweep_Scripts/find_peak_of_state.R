find_peak_of_state <- function(x,state_variable_name)
{
  require(tidyverse)
  require(magrittr)
  
  df <- x %>% as_tibble %>% pull(model) %>% unlist(recursive = FALSE)
  df <- df$rawOutput
  
  df %<>% as_tibble()
  state_vector <- df %>% pull(state_variable_name)
  max_state_name <- sprintf("%s_max",state_variable_name)
  x[[max_state_name]] <- max(state_vector)
  return(x)
}
