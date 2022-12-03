convert_model_raw_output_to_tibble <- function(x)
{
  require(tidyverse)
  require(magrittr)
  require(assertthat)
  
  x %<>% verify_sweep_row()
  
  model <- x %>% extract2("model") %>% extract2(1)
  
  model %<>% verify_SEIRModel()
  
  df <- model %>% extract2("rawOutput")
  
  N <- ncol(df)
  df %<>% as_tibble() %>%
    mutate(across(.cols = everything(),.fns = as.double)) %>%
    pivot_longer(cols = 2:N,names_to = "state_variable",values_to = "value") %>% arrange(state_variable,time) %>%
    mutate(state_variable = as.factor(state_variable))
  x$tidy_raw_output <- list(df)
  return(x)
}
