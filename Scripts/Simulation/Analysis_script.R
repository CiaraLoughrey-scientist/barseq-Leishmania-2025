#####Simulation analysis#################
#this script uses the simulation results as the input dataframes#

#set working directory#
setwd("") #add working directory details here#

#create subdirectory#
dir.create("Results")

#package loading#
library(tidyverse)
library(ggsci)
library(ggpubr)
library(vegan)



#load datasets#
#these are in the simulation folder in Data_files#

basic_sim <- read.csv("sim_output_10_to_500_barcodes.csv")
variance_sim <- read.csv("sim_output_10_to_500_barcodes_variance.csv")
clonal_rare_sim <- read.csv("sim_output_10_to_500_barcodes_clonal_expansion_1in1000000.csv")
clonal_medium_sim <- read.csv("sim_output_10_to_500_barcodes_clonal_1in10000.csv")
clonal_common_sim <- read.csv("sim_output_10_to_500_barcodes_clonal_1in100.csv")


#add groups to each#
clonal_common_sim$type <- "clonal_c_common_sim"
clonal_medium_sim$type <- "clonal_b_medium_sim"
clonal_rare_sim$type <- "clonal_a_rare_sim"
variance_sim$type <- "bb_variance_sim"
basic_sim$type <- "basic_sim"


#merge by binding together#
combo_sims_df <- rbind(clonal_common_sim, clonal_medium_sim, clonal_rare_sim, basic_sim, variance_sim)

#correlation of FP estimate (STAMP) vs real counted FP#

#############making the plots for paper figure################################

ggplot(combo_sims_df, aes(x = log10(True_FP), y= log10(Est_FP), colour = as.factor(Barcode_number))) + 
  geom_point(size = 2) + 
  labs(title = "", colour = "Barcode number:", x= "\nlog10(True FP)", y = "log10(STAMP-estimated FP)\n") +
  geom_abline(intercept = 0, slope = 1, linewidth = 1, alpha = 0.3) +
  scale_colour_viridis_d(option = "C") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() + theme(legend.title = element_text(size = 24),
                     legend.position = "bottom",
                     plot.title = element_text(hjust = 0.5, size = 24),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size=24), 
                     axis.line = 
                       element_line(colour = "black"),
                     axis.text.x = element_text(size = 20, colour = "black"),
                     axis.text.y = element_text(size = 20, colour = "black"),
                     legend.text = element_text(size = 20, colour = "black"),
                     strip.text.x = element_text(size = 24)) +
  guides(colour = guide_legend(nrow=1)) +
  facet_wrap(~type,  nrow = 1, labeller = as_labeller(
    c(basic_sim = "Basic", clonal_c_common_sim = "Clonal Expansion\n1 in 100", clonal_a_rare_sim = "Clonal Expansion\n1 in 1000000", 
      clonal_b_medium_sim = "Clonal Expansion\n1 in 10000", bb_variance_sim = "Growth rate\nvariance")))

ggsave(paste0(file = "Results/facet_simulation_plot_all_barcodes_figure.png"), height = 6, width=20)
ggsave(paste0(file = "Results/facet_simulation_plot_all_barcodes_figure.pdf"), height = 6, width=20)



#correlation#

#make data table#
titles <- c("Simulation_type", "Barcode_number", "tau", "p_value", "Z_statistic")
data_stats <- data.frame(matrix(nrow = 0, ncol = length(titles)))


#test correlation of whole dataset#
cor.test(basic_sim$True_FP, basic_sim$Est_FP, method = "kendall")
#significant correlation (p<2.2e-16, z=174.35, tau = 0.9378)#


#subset for barcode number-specific correlation#
corr_val <- with(subset(basic_sim, Barcode_number == c("10")),
     cor.test(True_FP, Est_FP, method = "kendall"))
#with 10 barcodes, tau = 0.9176, z=51.389, p<2.2e-16#
new_data <- c("Basic", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)
colnames(data_stats) <- titles

corr_val <- with(subset(basic_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)


corr_val <- with(subset(basic_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

#variance simulation#

cor.test(variance_sim$True_FP, variance_sim$Est_FP, method = "kendall")
#significant correlation (p<2.2e-16, z=151.03, tau = 0.8118)#


#subset for barcode number-specific correlation#
corr_val <- with(subset(variance_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

#clonal expansion simulation#
#rare#

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

#medium#

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

#common#

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)



#save dataset#
#supp table 1a#
write.csv(data_stats, "Results/Correlation_values.csv")



######################assess correlation specifically for FP (true) below 5000 as this fits our actual data##################

#new dataframes subsetted#
basic_sim <- subset(basic_sim, True_FP < 5001)
variance_sim <- subset(variance_sim, True_FP < 5001)
clonal_rare_sim <- subset(clonal_rare_sim, True_FP < 5001)
clonal_medium_sim <- subset(clonal_medium_sim, True_FP < 5001)
clonal_common_sim <- subset(clonal_common_sim, True_FP < 5001)


#make data table#
titles <- c("Simulation_type", "Barcode_number", "tau", "p_value", "Z_statistic")
data_stats <- data.frame(matrix(nrow = 0, ncol = length(titles)))


#test correlation of whole dataset#
cor.test(basic_sim$True_FP, basic_sim$Est_FP, method = "kendall")
#significant correlation (p<2.2e-16, z=174.35, tau = 0.9378)#


#subset for barcode number-specific correlation#
corr_val <- with(subset(basic_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
#with 10 barcodes, tau = 0.9176, z=51.389, p<2.2e-16#
new_data <- c("Basic", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)
colnames(data_stats) <- titles

corr_val <- with(subset(basic_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)


corr_val <- with(subset(basic_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(basic_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Basic", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

#variance simulation#

cor.test(variance_sim$True_FP, variance_sim$Est_FP, method = "kendall")
#significant correlation (p<2.2e-16, z=151.03, tau = 0.8118)#


#subset for barcode number-specific correlation#
corr_val <- with(subset(variance_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(variance_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("Variance", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

#clonal expansion simulation#
#rare#

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_rare_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("rare_clonal", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

#medium#

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_medium_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("medium_clonal", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

#common#

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("10")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "10", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("20")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "20", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("30")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "30", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("40")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "40", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("50")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "50", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("70")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "70", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("100")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "100", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("140")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "140", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("210")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "210", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("350")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "350", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)

corr_val <- with(subset(clonal_common_sim, Barcode_number == c("500")),
                 cor.test(True_FP, Est_FP, method = "kendall"))
new_data <- c("common_clonal", "500", corr_val$estimate, corr_val$p.value ,corr_val$statistic)
data_stats <- rbind(data_stats, new_data)



#save dataset#
#supp table 1b#
write.csv(data_stats, "Results/Correlation_values_below_5000FP.csv")

