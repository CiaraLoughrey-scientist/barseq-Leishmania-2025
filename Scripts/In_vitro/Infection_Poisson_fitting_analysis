#in vitro parasite burden analysis#

#packages#
library(fitdistrplus)
library(tidyverse)

#set directory#
setwd("") #add working directory details here#


#load data#
Infection_burden_b1_to_5<- read.csv("Data/Infection_burden_b1_to_5.csv") #save dataset into Data folder#

#add timepoints#
Infection_burden_b1_to_5_1hr <- subset(Infection_burden_b1_to_5, Timepoint_hr == 1)
Infection_burden_b1_to_5_3hr <- subset(Infection_burden_b1_to_5, Timepoint_hr == 3)


#fit poisson distribution#
pois_1hr_2 <- data.frame(c(), c())
list_params <- list()

#1 hr infection#
for (barcode in unique(Infection_burden_b1_to_5_1hr$Barcode_ID)){
  pois_1hr_df <- subset(Infection_burden_b1_to_5_1hr, subset = Barcode_ID == barcode, select = Parasite_burden)
  list_params[[barcode]] <- fitdistrplus::fitdist(as.numeric(pois_1hr_df$Parasite_burden), distr = "pois")
  pois_1hr_2 <- rbind(pois_1hr_2, c(list_params[[barcode]]$estimate[["lambda"]]))
}
pois_1hr_2 <- mutate(pois_1hr_2, barcode = c(1, 2, 3, 4, 5))
colnames(pois_1hr_2) <-c("Lambda","Barcode")
pois_1hr_2$Time <- "1hr"

pois_3hr_2 <- data.frame(c(), c())
list_params <- list()

#3 hr infection#
for (barcode in unique(Infection_burden_b1_to_5_3hr$Barcode_ID)){
  pois_3hr_df <- subset(Infection_burden_b1_to_5_3hr, subset = Barcode_ID == barcode, select = Parasite_burden)
  list_params[[barcode]] <- fitdistrplus::fitdist(as.numeric(pois_3hr_df$Parasite_burden), distr = "pois")
  pois_3hr_2 <- rbind(pois_3hr_2, c(list_params[[barcode]]$estimate[["lambda"]]))
}
pois_3hr_2 <- mutate(pois_3hr_2, barcode = c(1, 2, 3, 4, 5))
colnames(pois_3hr_2) <-c("Lambda","Barcode")
pois_3hr_2$Time <- "3hr"


#save parameters#
poisson_params_df <- bind_rows(pois_3hr_2, pois_1hr_2)
write.csv(poisson_params_df, "poisson_params.csv")
