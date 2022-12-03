invisible(x <- readRDS("sweeps/sweep_SEIR_RDS_November_30_2022/model_data.rds"))

source("find_peak_infections.R")
source("convert_model_raw_output_to_tibble.R")
source("find_time_of_peak_infections.R")
source("plot_SEIR.R")

x %<>% group_by(index)

x %<>% group_split()
k1 <- list()
qq <- 0L
for(smallx in x)
{
  qq <- qq + 1L
  smallx %<>% convert_model_raw_output_to_tibble() %>% bind_rows()
  smallx %<>% find_peak_infections()
  smallx %<>% find_time_of_peak_infections()
  k1[[qq]] <- plot_SEIR(smallx)
  ggsave(filename = sprintf("images/images_sweep_SEIR_RDS_November_30_2022/%d.png",qq),width = 10,height = 7.5,device = "png",dpi = 300)
}
x <- k1 %>% bind_rows
