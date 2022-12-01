invisible(x <- readRDS("sweeps/sweep_SEIR_RDS_November_30_2022/model_data.rds"))
library(tidyverse)
library(magrittr)
library(flumodels)
library(scales)

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
  k1[[qq]] <- smallx
}
x <- k1 %>% bind_rows

x %<>% filter(R0 > 1.25)


R0_names <- c(
  `1.5` = "R0 = 1.5",
  `2` = "R0 = 2",
  `3` = "R0 = 3",
  `5` = "R0 = 5",
  `7` = "R0 = 7",
  `8` = "R0 = 8",
  `10` = "R0 = 10",
  `15` = "R0 = 15",
  `25` = "R0 = 25",
  `50` = "R0 = 50"
)


p1 <- ggplot(data = x,aes(x = latentPeriod,y = infectiousPeriod,color = peak_infections)) + facet_wrap(~R0,nrow=3,labeller = as_labeller(R0_names)) + geom_point(size=3) + theme_bw()
p1 <- p1 + theme(legend.position = "right")
p1 <- p1 + xlab("tL, Latent Period [days]")
p1 <- p1 + ylab("tI, Infectious Period [days]")
p1 <- p1 + ggtitle("Peak Infection Value as a Function of ti, tL, and R0")

p1 <-p1 + scale_color_gradientn(colours = rev(c("magenta","red","orange","yellow","green","blue")), 
                                values = rescale(c(0,0.2,0.4,0.6,0.8,1.0)),
                                guide = "colorbar", limits=c(0,1),name= "Fraction")


p1 <- p1 + scale_x_continuous(limits = c(0,14),breaks = seq(from = 0,by = 2,to = 14))
p1 <- p1 + scale_y_continuous(limits = c(0,14),breaks = seq(from = 0,by = 2,to = 14))
p1 <- p1 + theme(axis.title.x = element_text(face = "bold"))
p1 <- p1 + theme(axis.title.y = element_text(face = "bold"))
p1 <- p1 + theme(title = element_text(face = "bold"))
p1 <- p1 + geom_abline(intercept = 0, slope = 1, size = 0.5,color="black",alpha=0.4,linetype="solid")
p1 <- p1 + theme(aspect.ratio = 1.0)
p1 <- p1 + theme(legend.position = "top")
ggsave(filename = "images/SEIR_dot_tile_plot_12_1_2022_sourcing_November_30_2022_data/peak_infection_tiledot_plot.png",device = "png",width = 7.5,height = 7.5,dpi = 300)
print(p1)



p1 <- ggplot(data = x,aes(x = latentPeriod,y = infectiousPeriod,color = time_of_peak_infections)) + facet_wrap(~R0,nrow=3,labeller = as_labeller(R0_names)) + geom_point(size=3) + theme_bw()
p1 <- p1 + theme(legend.position = "right")
p1 <- p1 + xlab("tL, Latent Period [days]")
p1 <- p1 + ylab("tI, Infectious Period [days]")
#p1 <- p1 + theme(legend.position = "none")
p1 <- p1 + ggtitle("Time-to-Peak Infection as a Function of ti, tL, and R0")
# p1 <- p1 + scale_color_gradient2(midpoint=0.5, low="red", mid="blue",
#                                  high="green", space ="Lab" )

p1 <-p1 + scale_color_gradientn(colours = c("magenta","red","yellow","orange","green","cyan","blue"), 
                             values = rescale(c(0,14,30,60,90,180,365)),
                             guide = "colorbar", limits=c(0,365),name= "Days")

p1 <- p1 + scale_x_continuous(limits = c(0,14),breaks = seq(from = 0,by = 2,to = 14))
p1 <- p1 + scale_y_continuous(limits = c(0,14),breaks = seq(from = 0,by = 2,to = 14))
p1 <- p1 + theme(axis.title.x = element_text(face = "bold"))
p1 <- p1 + theme(axis.title.y = element_text(face = "bold"))
p1 <- p1 + theme(title = element_text(face = "bold"))
p1 <- p1 + theme(legend.position = "top")
p1 <- p1 + theme(aspect.ratio = 1.0)
p1 <- p1 + geom_abline(intercept = 0, slope = 1, size = 0.5,color="black",alpha=0.4,linetype="solid")
ggsave(filename = "images/SEIR_dot_tile_plot_12_1_2022_sourcing_November_30_2022_data/day_of_occurrence_of_peak_infection_tiledot_plot.png",device = "png",width = 7.5,height = 7.5,dpi = 300)
print(p1)