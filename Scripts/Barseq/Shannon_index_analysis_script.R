#####Shannon calculation for reinfection experiment#################

#set working directory#
setwd("") #add working directory details here#

#create subdirectories#
dir.create("Results")
dir.create("Results/Diversity_analysis")



#package loading#
library(tidyverse)
library(ggsci)
library(ggpubr)
library(vegan)

#load data, reinf3#
barcode_counts <- read.csv("Data/2025-JCM-002-barcodes.csv") #save dataset into Data folder#


#subset for mice only#
barcode_counts <- barcode_counts[-c(2:7)]
#will not remove low read samples before Shannon analysis#

#remove B59#
barcode_counts <- subset(barcode_counts, X != c("TGGATCTTAGGACGCAACAT"))
#strip out barcodes#
df_mice <- barcode_counts
df_mice <- df_mice[,-1]

#shannon index#
#samples list#
list_samples <- names(df_mice)
#transform data#
barcode_t <- t(df_mice)

shannon_data <- data.frame(Sample = list_samples, 
                               Shannon_index = vegan::diversity(x = barcode_t, index = 'shannon'))

#samples list#
samples <- names(df_mice)
tissue<- gsub('[0-9]+', '', sub(".*?(\\d+.*)", "\\1", samples))

#make dataframe#
shannon_data$mouse_ID <- parse_number(samples)

#calculate Shannon index for all samples#
shannon_data <- mutate(shannon_data, 
                           tissue = tissue, 
                       expt = "3")
SI_df <- data.frame(shannon_data)

#save data#
write.csv(SI_df, "Results/Diversity_analysis/SHANNON_DATA_reinf3.csv")
