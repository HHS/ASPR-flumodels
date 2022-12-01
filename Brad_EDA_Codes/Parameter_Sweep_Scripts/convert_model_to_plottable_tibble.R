
df <- x$model[[1]]$rawOutput

df %<>% as_tibble()
df %<>% mutate(across(time:R,.fns = as.double))

df %<>% pivot_longer(cols = S:R,names_to = "state_variable",values_to = "value") %>% arrange(state_variable,time)
df %<>% mutate(state_variable = as.factor(state_variable))
x$tibble <- list(df)