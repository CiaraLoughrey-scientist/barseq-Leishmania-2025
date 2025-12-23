
#GD calculation#
#must be run before any other GD scripts, as these rely on datasets created here#

#set working directory#
setwd("") #add working directory details here#
#create subdirectories#
dir.create("Results")
dir.create("Results/GD_heatmaps")
dir.create("Results/GD_networks")
dir.create("Results/GD_data")
dir.create("Results/GD_metrics")

#package loading#
library(tidyverse)
library(ggsci)
library(ggpubr)
library(vegan)
library(ggraph)
library(igraph)
library(naniar)
library(Cairo)

#Genetic Distance (GD) defined as Cavalli-Sforza chord distance#
chord_dist = function(count1, count2){
  freq1 = count1/sum(unlist(count1))
  freq2 = count2/sum(unlist(count2))
  costheta = sum(sqrt(freq1*freq2))
  chord_dist = (2*(sqrt(2))/pi)*(sqrt(1-costheta))
  chord_dist
}

#load dataset#

barcode_counts <- read.csv("Data/2025-JCM-002-barcodes.csv") #save dataset into Data folder#


#subset for mice only and label the mouse ID#
barcode_counts <- barcode_counts[-c(2:7)]
#will not remove low read samples before starting analysis as difference in total read frequency is not done here, only read freq per barcode#

######################################################################################################################

#remove B59#
barcode_counts <- subset(barcode_counts, X != c("TGGATCTTAGGACGCAACAT"))

#data formatting to include mouse/tissue attributes#
long_mice <- pivot_longer(barcode_counts,
                          names_to = "Samples",
                          values_to = "barcode_count",
                          cols = colnames(barcode_counts[-c(1,1)])) #general rule will work for all dataframes created using our code for amalgamating barcodes#

samples = long_mice$Samples
long_mice$mouse_ID <- NULL
long_mice$mouse_ID <- parse_number(samples)

#name tissues#
Tissue <- sub('^M[0-9]+', '', samples)
long_mice <- mutate(long_mice, "Tissue" = Tissue)

long_mice$Tissue2 <- long_mice$Tissue
long_mice$Tissue2[long_mice$Tissue2 %in% c("LL", "LM", "LS", "L")] <- "LIVER" 
long_mice$Tissue[long_mice$Tissue %in% c("LL", "LM", "LS", "L")] <- "LIVER" 
long_mice$Tissue2[long_mice$Tissue2 %in% c("LN1", "LN2")] <- "LN" 
long_mice$Tissue2[long_mice$Tissue2 %in% c("BM1", "BM2")] <- "BM" 
long_mice$Tissue2[long_mice$Tissue2 %in% c("SKIN1", "SKIN2", "SKIN3", "SKIN4", "SKIN5", "SKIN6", "SKIN7",
                                           "SKIN8", "SKIN9", "SKIN10", "SKIN11", "SKIN12")] <- "SKIN" 
long_mice$Tissue[long_mice$Tissue %in% c("SKIN1")] <- "SKIN01" 
long_mice$Tissue[long_mice$Tissue %in% c("SKIN2")] <- "SKIN02"
long_mice$Tissue[long_mice$Tissue %in% c("SKIN3")] <- "SKIN03"
long_mice$Tissue[long_mice$Tissue %in% c("SKIN4")] <- "SKIN04"
long_mice$Tissue[long_mice$Tissue %in% c("SKIN5")] <- "SKIN05"
long_mice$Tissue[long_mice$Tissue %in% c("SKIN6")] <- "SKIN06"
long_mice$Tissue[long_mice$Tissue %in% c("SKIN7")] <- "SKIN07"
long_mice$Tissue[long_mice$Tissue %in% c("SKIN8")] <- "SKIN08"
long_mice$Tissue[long_mice$Tissue %in% c("SKIN9")] <- "SKIN09"
#long_mice$Tissue2[long_mice$Tissue2 %in% c("S")] <- "SPLEEN" 
#long_mice$Tissue[long_mice$Tissue %in% c("S")] <- "SPLEEN" #this causes spleen to separate from other tissues in heatmaps, we will rename in the heatmap instead to avoid this#



#reinf3 dataset#
#create dataframes for bet and harm and eigen centrality metrics#
bet_data <- data.frame()
harm_data <- data.frame()
eigen_data <- data.frame()
degree_data <- data.frame()
all_GD_data <- data.frame()

#find the mice in current study
mice_in_study <- unique(((long_mice$mouse_ID)))
#remove all nas
mice_in_study <- sort(mice_in_study[!is.na(mice_in_study)])
#create wide format for each mouse here
#create list that will store all dfs
list_DFs_per_mouse <- list()
#run loop for all mice
for (mouse in mice_in_study){
  #back to wide format for each mouse#
  list_DFs_per_mouse[[mouse]] <- pivot_wider(subset(long_mice, mouse_ID == mouse),
                                             id_cols = c("X"),
                                             names_from = Tissue, 
                                             values_from = barcode_count)
  print(paste0("For mouse: ", mouse))
  print(head(list_DFs_per_mouse[[mouse]]))
  
  #now run code on subdivided dataframe for each mouse individually#
  #for GD calculation#
  
  df <- list_DFs_per_mouse[[mouse]][-c(1)]
  paired <- combn(colnames(df), 2)
  len = length(colnames(df))
  chord_df <- matrix(ncol = len, nrow = len)
  chord_df[,] <- 0
  colnames(chord_df) <- colnames(df)
  rownames(chord_df) <- colnames(df)
  
  for (p in seq(1, ncol(paired))) {
    colA <- paired[1,p]
    colB <- paired[2,p]
    cd <- chord_dist(df[colA], df[colB])
    chord_df[colA,colB] <- cd
    chord_df[colB,colA] <- cd
  }
  
  chord_df_data <- data.frame(chord_df)
  chord_df_data <- mutate(chord_df_data, "comparison" = rownames(chord_df_data))
  
  
  df_long <- pivot_longer(
    chord_df_data,
    cols = colnames(df),
    names_to = "Sample",
    values_to = "GD"
  )
  
  #heatmaps#
  ggplot(df_long, aes(x = Sample, y= comparison, fill = GD)) + geom_tile() +
    labs(fill = "Genetic Distance:\n", x= "\nTissue", y = "Tissue\n") +
    scale_fill_viridis_c(option = "B", limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    scale_x_discrete(labels = c(S = "SPLEEN")) +
    scale_y_discrete(labels = c(S = "SPLEEN")) +
    theme_bw() + theme(legend.title = element_text(size = 20),
                       legend.position = "bottom",
                       plot.title = element_text(hjust = 0.5, size = 24),
                       panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.title = element_text(size=20), 
                       axis.line = 
                         element_line(colour = "black"), 
                       axis.text.x = element_text(size = 16, colour = "black", vjust = 0.5, angle = 90),
                       axis.text.y = element_text(size = 16, colour = "black"),
                       legend.text = element_text(size = 16, colour = "black"))
  ggsave(paste0(file = "Results/GD_heatmaps/GD_dist_heatmap_",paste0("mouse_", mouse),".png"), height = 10, width=10)
  ggsave(paste0(file = "Results/GD_heatmaps/GD_dist_heatmap_",paste0("mouse_", mouse),".pdf"), height = 10, width=10)
  

  #save GD data#
  write.csv(chord_df_data, file = paste0("Results/GD_data/GD_dist_","mouse_", mouse,".csv"))
  
  df_long <- mutate(df_long, 
                    mouse_ID = mouse)

  #add GD data to full dataframe#
  all_GD_data <- bind_rows(all_GD_data, df_long)
  
  #network#
  
  chord_df_data[is.na(chord_df_data)] <- 1
  #Create a graph adjacency based on GD between cell types in  pairwise fashion#
  g <- graph_from_adjacency_matrix(
    as.matrix(subset(chord_df_data, select = -c(comparison))),
    mode = c("undirected"),
    weighted = TRUE,
    diag = FALSE,
    add.colnames = NULL,
    add.rownames = NA
  )
  #Remove any vertices remaining that have no edges#
  g <- delete_vertices(g, igraph::degree(g)==0)
  V(g)$shape <- "sphere"
  
  #delete edges above 0.5 GD#
  g <- delete_edges(g, E(g)[which(E(g)$weight>0.5)])
  g <- delete_vertices(g, igraph::degree(g)==0)
  
  edgeweights <- 1-E(g)$weight #so we are plotting 'relatedness' on the graph#
  E(g)$weight <- edgeweights
  edgeweights2 <- 1-edgeweights
  
  #plot network#
  ggraph(g, layout = "circle") +
    geom_edge_fan( aes(edge_colour = weight, width  = weight, alpha = weight),
                   show.legend = TRUE) +
    geom_node_point(colour = "grey", size = 3) +
    geom_node_text(aes(label = name),  colour = 'black', size=8,
                   show.legend = FALSE) +
    scale_edge_color_viridis(
      begin = 0,
      end = 1,
      option = "B",
      direction = -1,
      space = "Lab",
      guide = "edge_colourbar",
      aesthetics = "edge_colour",
      limits = c(0.5, 1)
    ) +
    theme_bw() +
    labs(title = "Full distance network\n", edge_colour = "Genetic Relatedness:") +
    theme(legend.title = element_text(size = 16, hjust = 0.5),
          legend.position = "bottom",
          legend.title.position = "top",
          plot.title = element_text(hjust = 0.5, size = 24),
          panel.border = element_blank(),
          plot.caption = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 12, hjust = 1),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(),  
          axis.ticks.y=element_blank(),
          axis.title = element_blank()) +
    guides (edge_width = "none", edge_alpha = "none")
  ggsave(paste0("Results/GD_networks/full_network",paste0("mouse",mouse),".png"), width = 12, height = 10)
  ggsave(paste0("Results/GD_networks/full_network",paste0("mouse",mouse),".pdf"), width = 12, height = 10)
  
  #mst#
  
  #Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm#
  mst <- mst(g, algorithm="prim", weights = edgeweights2) #so we use distance again here for mst calculation#
  
  #find communities#
  mst.communities <- cluster_louvain(g)
  #plot the graph#
  mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
  V(mst)$color <- as.factor(mst.communities$membership + 1)

#plot MST#
  ggraph(mst, layout = "dh") +
    geom_edge_link( aes(edge_colour = weight, width  = weight, alpha = weight),
                    show.legend = TRUE) +
    geom_node_point(aes(colour = as.factor(mst.communities$membership)), size = 8) +
    geom_node_text(aes(label = name),  colour = 'black', size=8,
                   show.legend = FALSE) +
    scale_edge_color_viridis(
      begin = 0,
      end = 1,
      option = "B",
      direction = -1,
      space = "Lab",
      guide = "edge_colourbar",
      aesthetics = "edge_colour", 
      breaks = c(0.5, 0.75, 1), 
      limits = c(0.5, 1)
    ) +
    scale_colour_viridis_d(option = "H", alpha = 0.7) +
    theme_bw() +
    labs(title = "MST network\n", edge_colour = "Genetic Relatedness:") +
    theme(legend.title = element_text(size = 16, hjust = 0.5),
          legend.position = "bottom",
          legend.title.position = "top",
          plot.title = element_text(hjust = 0.5, size = 24),
          panel.border = element_blank(),
          plot.caption = element_text(size = 14),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 12, hjust = 1),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(),  
          axis.ticks.y=element_blank(),
          axis.title = element_blank()) +
    guides (edge_width = "none", edge_alpha = "none", colour = "none")
  
  ggsave(paste0("Results/GD_networks/mst_network",paste0("mouse",mouse),".png"), width = 12, height = 10)
  ggsave(paste0("Results/GD_networks/mst_network",paste0("mouse",mouse),".pdf"), width = 12, height = 10)
  
  #metrics to save as data files#
  
  bet <- betweenness(g, weights = edgeweights2)
  harm <- harmonic_centrality(g, weights = edgeweights2)
  bet <- bet[order(names(bet))]
  harm <- harm[order(names(harm))]
  eigen <- eigen_centrality(g, weights = NULL)
  degree <- degree(g)
  
  #save betweenness and harmonic and eigen centrality scores#
  bet_data <- bind_rows(bet_data, bet)
  harm_data <- bind_rows(harm_data, harm)
  eigen_data <- bind_rows(eigen_data, eigen$vector)
  degree_data <- bind_rows(degree_data, degree)
}

#save data#
write.csv(eigen_data, "Results/GD_metrics/eigen_data.csv")

write.csv(bet_data, "Results/GD_metrics/bet_data.csv")

write.csv(harm_data, "Results/GD_metrics/harm_data.csv")

write.csv(degree_data, "Results/GD_metrics/degree_data.csv")

write.csv(all_GD_data, "Results/GD_data/Reinf3_all_GD_data.csv")

#################################################################################
