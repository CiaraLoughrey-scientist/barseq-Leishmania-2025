#GD plots and stats#
#you must run the GD_calculation_script before running this one, as it relies on datasets created from that script#

#set working directory#
setwd("") #add working directory details here#

#create subdirectories#
dir.create("Results")
dir.create("Results/Figures_for_paper")
dir.create("Results/Stats")

#package loading#
library(tidyverse)
library(ggsci)
library(ggpubr)
library(vegan)

#load data#
reinf3_df <- read.csv("Results/GD_data/Reinf3_all_GD_data.csv")

#groups and tissue IDs#
reinf3_df$Tissue2 <- reinf3_df$Sample
reinf3_df$Tissue2[reinf3_df$Sample %in% c("LL", "LM", "LS")] <- "L" 
reinf3_df$Tissue2[reinf3_df$Sample %in% c("LN1", "LN2")] <- "LN" 
reinf3_df$Tissue2[reinf3_df$Sample %in% c("BM1", "BM2")] <- "BM" 
reinf3_df$Tissue2[reinf3_df$Sample %in% c("SKIN01", "SKIN02", "SKIN03", "SKIN04", "SKIN05", "SKIN06", "SKIN07",
                                           "SKIN08", "SKIN09", "SKIN10", "SKIN11", "SKIN12")] <- "SKIN"

reinf3_df$comparison2 <- reinf3_df$comparison
reinf3_df$comparison2[reinf3_df$comparison %in% c("LL", "LM", "LS")] <- "L" 
reinf3_df$comparison2[reinf3_df$comparison %in% c("LN1", "LN2")] <- "LN" 
reinf3_df$comparison2[reinf3_df$comparison %in% c("BM1", "BM2")] <- "BM" 
reinf3_df$comparison2[reinf3_df$comparison %in% c("SKIN01", "SKIN02", "SKIN03", "SKIN04", "SKIN05", "SKIN06", "SKIN07",
                                          "SKIN08", "SKIN09", "SKIN10", "SKIN11", "SKIN12")] <- "SKIN"
reinf3_df$Tissue2[reinf3_df$Tissue2 == "S"] <- "SPLEEN"


#add animal groups#
#add groups#
reinf3_df$group <- NULL
reinf3_df$group[reinf3_df$mouse_ID %in% c("1", "2", "3", "4", "5")] <- "1single_04wk"
reinf3_df$group[reinf3_df$mouse_ID %in% c("12", "13", "14")] <- "1single_10wk"
reinf3_df$group[reinf3_df$mouse_ID %in% c("21", "22", "23", "24", "25")] <- "1single_12wk"
reinf3_df$group[reinf3_df$mouse_ID %in% c("15", "16", "17", "18", "19", "20")] <- "2reinf_2wk"
reinf3_df$group[reinf3_df$mouse_ID %in% c( "6", "7", "8", "9", "10", "11")] <- "2reinf_2d"


reinf3_df$group2 <- NULL
reinf3_df$group2[reinf3_df$mouse_ID %in% c("1", "2", "3", "4", "5")] <- "1single_early"
reinf3_df$group2[reinf3_df$mouse_ID %in% c("12", "13", "14")] <- "1single_late"
reinf3_df$group2[reinf3_df$mouse_ID %in% c("21", "22", "23", "24", "25")] <- "1single_late"
reinf3_df$group2[reinf3_df$mouse_ID %in% c("15", "16", "17", "18", "19", "20")] <- "2reinf_2wk"
reinf3_df$group2[reinf3_df$mouse_ID %in% c( "6", "7", "8", "9", "10", "11")] <- "2reinf_2d"



#Bone marrow plot#


ggplot(subset(reinf3_df, comparison2 == "BM" & group %in% c("1single_04wk", "1single_10wk", "1single_12wk") & GD != 0 & Tissue2 != c("SKIN")), 
       aes(x = Tissue2, y= GD)) + 
  geom_jitter(aes(colour = group2),
              size = 2, stroke = 2,
              position = position_jitterdodge(dodge.width = .75, jitter.width = 0.1)) + 
  geom_boxplot(aes(fill = group2), alpha = 0.3,
               outlier.shape = NA, 
               width = 0.75, lwd = 0.8, fatten = 1.5, show.legend = FALSE) +
  labs(title = "Bone Marrow\n", colour = "Infection phase:", x= "\nTissue comparison", y = "Genetic Distance\n") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0, 1)) +
  scale_x_discrete(labels = c(BM = "Bone\nMarrow", GUT = "Gut", LIVER = "Liver", LN =  "Lymph\nNode", LUNG = "Lung", SPLEEN = "Spleen")) +
  theme_bw() + theme(legend.title = element_text(size = 16, hjust = 0.5),
                     legend.position = "bottom",
                     legend.title.position = "top",
                     plot.title = element_text(hjust = 0.5, size = 26),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size=24), 
                     axis.line = 
                       element_line(colour = "black"), 
                     axis.text.x = element_text(size = 16, colour = "black"),
                     axis.text.y = element_text(size = 16, colour = "black"),
                     legend.text = element_text(size = 16, colour = "black"),
                     strip.text.x = element_text(size = 20)) +
  scale_fill_manual(values = c("#21918c", "#440154"), labels = c("Early", "Chronic")) +
  scale_colour_manual(values = c("#21918c", "#440154"), labels = c("Early", "Chronic"))+
  guides(fill="none") +
  stat_pwc(aes(y=GD, x=Tissue2, group = group2),
           hide.ns = "p", step.increase = 0.15, vjust = -0.25,
           bracket.nudge.y = 0.1, tip.length = 0,  size = 0.5, label.size = 4.5, label = "p", 
           method = "wilcox_test")

ggsave(file = "Results/Figures_for_paper/reinf3_GD_plot_final_BM_ONLY.png", height = 6, width=6)
ggsave(file = "Results/Figures_for_paper/reinf3_GD_plot_final_BM_ONLY.pdf", height = 6, width=6)


#skin to tissues#

ggplot(subset(reinf3_df, comparison2 == c("SKIN") & group %in% c("1single_04wk", "1single_10wk", "1single_12wk") & GD != 0 & Tissue2 != "SKIN"), 
       aes(x = group2, y= GD)) + 
  geom_jitter(aes(colour = group2, shape = as.factor(mouse_ID)),
              size = 1.5, stroke = 1,
              position = position_jitterdodge(dodge.width = .75, jitter.width = 0.05)) + 
  geom_boxplot(aes(fill = group2), alpha = 0.3,
               outlier.shape = NA, 
               width = 0.75, lwd = 0.8, fatten = 1.5, show.legend = FALSE) +
  scale_shape_manual(values = c(1, 2, 3, 4, 15, 16, 17, 18, 10,7,6, 8)) +
  labs(title = "Visceral tissues to Skin\n", colour = "Infection phase:", x= "\nInfection phase", y = "Genetic Distance\n") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0, 1)) +
  scale_x_discrete(labels = c("Early", "Chronic")) +
  theme_bw() + theme(legend.title = element_text(size = 16, hjust = 0.5),
                     legend.position = "bottom",
                     legend.title.position = "top",
                     plot.title = element_text(hjust = 0.5, size = 26),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size=24), 
                     axis.line = 
                       element_line(colour = "black"), 
                     axis.text.x = element_text(size = 16, colour = "black"),
                     axis.text.y = element_text(size = 16, colour = "black"),
                     legend.text = element_text(size = 16, colour = "black"),
                     strip.text.x = element_text(size = 20)) +
  scale_fill_manual(values = c("#21918c", "#440154"), labels = c("Early", "Chronic")) +
  scale_colour_manual(values = c("#21918c", "#440154"), labels = c("Early", "Chronic"))+
  guides(fill="none", shape = "none") +
  stat_pwc(aes(y=GD, x=group2, group = group2),
           hide.ns = "p", step.increase = 0.15, vjust = -0.25,
           bracket.nudge.y = 0.1, tip.length = 0,  size = 0.5, label.size = 4.5, label = "p", 
           method = "wilcox_test") +
  facet_wrap(~Tissue2, nrow = 1, labeller = as_labeller(
    c(BM = "Bone Marrow", GUT = "Gut", LIVER = "Liver", LN =  "Lymph Node", LUNG = "Lung", SPLEEN = "Spleen", SKIN = "Skin")))

ggsave(file = "Results/Figures_for_paper/reinf3_GD_plot_final_skin.png", height = 6, width=18)
ggsave(file = "Results/Figures_for_paper/reinf3_GD_plot_final_skin.pdf", height = 6, width=18)


#skin to skin#

ggplot(subset(reinf3_df, comparison2 == c("SKIN") & group %in% c("1single_04wk", "1single_10wk", "1single_12wk") & GD != 0 & Tissue2 == "SKIN"), 
       aes(x = group2, y= GD)) + 
  geom_jitter(aes(colour = group2, shape = as.factor(mouse_ID)),
              size = 1.5, stroke = 1,
              position = position_jitterdodge(dodge.width = .75, jitter.width = 0.1)) + 
  geom_boxplot(aes(fill = group2), alpha = 0.3,
               outlier.shape = NA, 
               width = 0.75, lwd = 0.8, fatten = 1.5, show.legend = FALSE) +
  scale_shape_manual(values = c(1, 2, 3, 4, 15, 16, 17, 18, 10,7,6, 8)) +
  labs(title = "Skin to Skin\n", colour = "Infection phase:", x= "\nInfection phase", y = "Genetic Distance\n") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(0, 1)) +
  scale_x_discrete(labels = c("Early", "Chronic")) +
  theme_bw() + theme(legend.title = element_text(size = 16, hjust = 0.5),
                     legend.position = "bottom",
                     legend.title.position = "top",
                     plot.title = element_text(hjust = 0.5, size = 26),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size=24), 
                     axis.line = 
                       element_line(colour = "black"), 
                     axis.text.x = element_text(size = 16, colour = "black"),
                     axis.text.y = element_text(size = 16, colour = "black"),
                     legend.text = element_text(size = 16, colour = "black"),
                     strip.text.x = element_text(size = 20)) +
  scale_fill_manual(values = c("#21918c", "#440154"), labels = c("Early", "Chronic")) +
  scale_colour_manual(values = c("#21918c", "#440154"), labels = c("Early", "Chronic"))+
  guides(fill="none", shape = "none") +
  stat_pwc(aes(y=GD, x=group2, group = group2),
           hide.ns = "p", step.increase = 0.15, vjust = -0.25,
           bracket.nudge.y = 0.1, tip.length = 0,  size = 0.5, label.size = 4.5, label = "p", 
           method = "wilcox_test")

ggsave(file = "Results/Figures_for_paper/reinf3_GD_plot_final_skin_to_skin.png", height = 6, width=6)
ggsave(file = "Results/Figures_for_paper/reinf3_GD_plot_final_skin_to_skin.pdf", height = 6, width=6)


#mann whitney test#
#supp table 3#
df_p_val <- subset(reinf3_df, group %in% c("1single_04wk", "1single_10wk", "1single_12wk") & GD != 0) %>%
  rstatix::group_by(Tissue2, comparison2) %>% 
  rstatix::wilcox_test(GD ~ group2, detailed = T, paired = F) 
write.csv(df_p_val, file = "Results/Stats/mann_whitney_GD_chronic_vs_early.csv")
