plot_SEIR <- function(x)
{
  library(tidyverse)
  library(magrittr)
  
  model_data <- x %>% pull(model) %>% unlist(recursive = FALSE)
  simulation_time <- model_data$parameters$simulationLength
  
  df <- x %>% pull(tidy_raw_output) %>% unlist(recursive = FALSE) %>% as_tibble()
  
  p1 <- ggplot(data = df,mapping = aes(x = time,y = value,color = state_variable,linetype=state_variable))
  p1 <- p1 + theme_bw()
  p1 <- p1 + scale_x_continuous(limits = c(0,simulation_time),
                                breaks = seq(from = 0,to = simulation_time,by = 30),
                                minor_breaks = seq(from = 0, to = simulation_time,by = 5))
  
  p1 <- p1 + scale_y_continuous(limits = c(0,1),
                                breaks = seq(from = 0,to = 1,by = 0.1),
                                minor_breaks = seq(from = 0, to = 1,0.02))
  
  p1 <- p1 + geom_line(alpha=0.6,size=1.0)
  
  LINES <- c("S" = "solid","E" = "longdash","I" = "twodash","R" = "dotted")
  p1 <- p1 + scale_linetype_manual(values = LINES)
  COLORS <- c("S" = "red", "E" = "blue", "I" = "green", "R" = "magenta")
  p1 <- p1 + scale_colour_manual(values = COLORS)

  month_marks <- seq(from = 0,to = simulation_time,by = 30)
  for(x_val in month_marks)
  {
    p1 <- p1 + geom_vline(xintercept = x_val,size=0.5,alpha=0.3,color="black")
  }
  
  p1 <- p1 + geom_hline(yintercept = 0.5,color="black",alpha=0.5,size = 0.5,linetype="dotted")
  title_string <- sprintf("SEIRModel, R0 = %0.5g, latentPeriod = %0.5g, infectiousPeriod = %0.5g",
                          x$R0,
                          x$latentPeriod,
                          x$infectiousPeriod)
  
  
  p1 <- p1 + ggtitle(title_string)
  p1 <- p1 + xlab("Time [days]")
  p1 <- p1 + ylab("Fraction in State S, E, I, or R")
  p1 <- p1 + theme(legend.position = "top")
  p1 <- p1 + theme(axis.title.x = element_text(face = "bold"))
  p1 <- p1 + theme(axis.title.y = element_text(face = "bold"))
  
  x$ggplot <- list(p1)
  return(x)
}
