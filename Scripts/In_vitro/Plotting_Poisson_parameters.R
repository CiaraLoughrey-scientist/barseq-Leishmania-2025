#plot poisson parameters#
#this script must be run after the Infection_Poisson_fitting_analysis script, as this generates the data frame used here#

#set wd#
setwd("") #add working directory details here#

#pakcages#
library(tidyverse)

#load data#
poisson_params <- read.csv("poisson_params.csv")


#with error bars of 2x SD#
poisson_params2 <- poisson_params
poisson_params2 <- poisson_params2 %>%
  group_by(Time) %>%
  mutate(mean_lambda = mean(Lambda))
poisson_params2 <- poisson_params2 %>%
  group_by(Time) %>%
  mutate(sd_lambda = sd(Lambda))

#plots#

ggplot(data = poisson_params2) +
  geom_jitter(aes(x = Time, y = Lambda, colour = as.factor(Barcode)), size = 5,
              position = position_jitterdodge(dodge.width = .15, jitter.width = 0.1)) +
  geom_errorbar(aes(x=Time, ymin=mean_lambda-(2*sd_lambda), ymax=mean_lambda+(2*sd_lambda)), width = 0.2) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4))+
  labs(title = "Distribution lambda changes\n", x="\nTime-point", y="Poisson\nlambda parameter\n", colour = "Barcode:")+
  theme_bw() +  theme(legend.title = element_text(size = 20),
                      legend.position = "bottom",
                      plot.title = element_text(hjust = 0.7, size = 28),
                      panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.title = element_text(size=26), 
                      axis.line = 
                        element_line(colour = "black"), 
                      axis.text.x = element_text(size = 20, colour = "black"),
                      axis.text.y = element_text(size = 20, colour = "black"),
                      legend.text = element_text(size = 20))+
  scale_colour_viridis_d() +
  scale_fill_viridis_d()

ggsave("poisson_param_plot.png", width = 6, height = 8)
ggsave("poisson_param_plot.pdf", width = 5.5, height = 8)
