find_peak_infections <- function(x)
{
  require(tidyverse)
  require(magrittr)
  
  df <- x %>% as_tibble %>% pull(model) %>% unlist(recursive = FALSE)
  df <- df$rawOutput
  
  df %<>% as_tibble()
  I_state <- df %>% pull(I)
  
  x$peak_infections <- max(I_state)
  return(x)
}
