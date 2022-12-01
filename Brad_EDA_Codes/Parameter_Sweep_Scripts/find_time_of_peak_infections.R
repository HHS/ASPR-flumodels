find_time_of_peak_infections <- function(x)
{
  require(tidyverse)
  require(magrittr)
  
  df <- x %>% as_tibble %>% pull(model) %>% unlist(recursive = FALSE)
  df <- df$rawOutput
  
  df %<>% as_tibble()
  I_state <- df %>% pull(I)
  max_index <- which.max(I_state)
  time_vector <- df %>% pull(time)
  
  x$time_of_peak_infections <- time_vector[max_index]
  return(x)
}
