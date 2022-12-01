convert_model_raw_output_to_tibble <- function(x,...)
{
  require(tidyverse)
  require(magrittr)
  
  df <- x %>% as_tibble %>% pull(model) %>% unlist(recursive = FALSE)
  df <- df$rawOutput
  
  df %<>% as_tibble()
  df %<>% mutate(across(time:R,.fns = as.double))
  
  df %<>% pivot_longer(cols = S:R,names_to = "state_variable",values_to = "value") %>% arrange(state_variable,time)
  df %<>% mutate(state_variable = as.factor(state_variable))
  x$tidy_raw_output <- list(df)
  return(x)
}
