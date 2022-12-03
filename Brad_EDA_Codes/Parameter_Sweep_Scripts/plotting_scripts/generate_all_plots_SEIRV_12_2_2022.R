invisible(x <- readRDS("sweeps/sweep_SEIRV_RDS_November_30_2022/model_data.rds"))
source("convert_model_raw_output_to_tibble.R")
source("plot_SEIRV.R")
library(tidyverse)
library(magrittr)
library(purrr)
library(flumodels)

x %<>% group_by(index)

x %<>% group_split()
k1 <- list()
qq <- 0L
for(smallx in x)
{
  print(qq)
  qq <- qq + 1L
  smallx %<>% convert_model_raw_output_to_tibble() %>% bind_rows()
  k1[[qq]] <- plot_SEIRV(smallx)
  ggsave(filename = sprintf("images/images_sweep_SEIRV_RDS_December_2_2022/%d.png",qq),width = 10,height = 7.5,device = "png",dpi = 300)
}
x <- k1 %>% bind_rows
