#validation of TSFP and FFP#
#campbell hotinger data used for this re-analysis#

#set working directory#
setwd("") #add working directory details here#

#package loading#
library(tidyverse)
library(ggsci)
library(ggpubr)
library(vegan)
library(ggraph)
library(igraph)
library(naniar)
library(Seurat)
library(fitdistrplus)

#founder pop function#
founder_pop <- function(output_count, input_count){
  output_count = as.vector(output_count)
  input_count = as.vector(input_count)
  output_freq = output_count/sum(output_count)
  input_freq = input_count/sum(input_count)
  diffsq = sum(((output_freq - input_freq)^2)/(input_freq*(1-input_freq)))
  fhat = (1/(length(input_count)))*diffsq
  input_seq_total = sum(input_count)
  output_seq_total = sum(output_count)
  founderpop = 1/(fhat - (1/input_seq_total) - (1/output_seq_total))
  founderpop
}



#############################################################################################################
#############################################################################################################
#comparing input for migration/FFP validation#
#############################################################################################################
#############################################################################################################

#load the dataframes#
df1 <- read.csv("jh3-7_freqNohop.csv")
df2 <- read.csv("jh08_freqNohop.csv")
df3 <- read.csv("jh09_freqNohop.csv")
df4 <- read.csv("jh10_freqNohop.csv")

#remove columns, unclear what they are#
remove_this_col <- colnames(df1[2:12])
df1[colnames(df1) %in% remove_this_col] <- NA 
df1 <- df1[,setdiff(names(df1), remove_this_col)]

remove_this_col2 <- colnames(df2[2:12])
df2[colnames(df2) %in% remove_this_col2] <- NA 
df2 <- df2[,setdiff(names(df2), remove_this_col2)]

remove_this_col3 <- colnames(df3[2:12])
df3[colnames(df3) %in% remove_this_col3] <- NA 
df3 <- df3[,setdiff(names(df3), remove_this_col3)]

remove_this_col4 <- colnames(df4[2:12])
df4[colnames(df4) %in% remove_this_col4] <- NA 
df4 <- df4[,setdiff(names(df4), remove_this_col4)]

#merge together by barcode#
barcode_counts5 <- merge(df1, df2, by = "barcode")
barcode_counts5 <- merge(barcode_counts5, df3, by = "barcode")
barcode_counts5 <- merge(barcode_counts5, df4, by = "barcode")
#832 samples total#

#remove low read samples as before#
remove_columns_x <- sapply(subset(barcode_counts5, select = -c(barcode)), sum)
remove_columns <- names(remove_columns_x[remove_columns_x < 50000])
barcode_counts5[colnames(barcode_counts5) %in% remove_columns] <- NA 

#137 columns identified as too low and need removing#
barcode_counts5 <- barcode_counts5[,setdiff(names(barcode_counts5), remove_columns)]

#total samples 745#
barcode_counts6 <- barcode_counts5[-c(1)] #remove barcode list column#

#change AO samples to correctly identify individual mice#
colnames(barcode_counts6) <- sub("AO1_", "AO1", colnames(barcode_counts6))
colnames(barcode_counts6) <- sub("AO2_", "AO2", colnames(barcode_counts6))
colnames(barcode_counts6) <- sub("AO3_", "AO3", colnames(barcode_counts6))

#new empty dataframe#
FP_df5 <- data.frame()

#loop TSFP through all pairs including the inputs#
for (tissue in barcode_counts6){
  
  df <- barcode_counts6
  input_names <- colnames(df) #list of column names#
  #subset dataframe so all columns contain only barcodes which are over 0 in the tissue used for input#
  df <- subset(df, tissue > 0)
  #new input tissue vector with 0 removed#
  tissue2 <- tissue[tissue != 0]
  #name input as tissue2 vector in the FP calculation#
  FP_df4 <- apply(df, 2, function(col_in) {
    return(founder_pop(col_in, tissue2))
  })
  FP_df5 <- rbind(FP_df5, FP_df4)
  
}
FP_df6 <- data.frame(FP_df5) #all pairwise comparisons in dataframe#


colnames(FP_df6) <- names(df) #column names add#
colnames(FP_df6) <- paste("OUTPUT_", colnames(FP_df6), sep = "")
rownames(FP_df6) <- names(df) #rownames add#
rownames(FP_df6) <- paste("INPUT_", rownames(FP_df6), sep = "")

FP_df6 <- mutate(FP_df6, "INPUT_POOL" = rownames(FP_df6)) #add new column with input sample details#

#pivot data longer#
df_long2 <- pivot_longer(
  FP_df6,
  cols = starts_with("OUTPUT"),
  names_to = "OUTPUT_POOL",
  values_to = "FP"
)

df_long2$INPUT2 <- gsub("INPUT_", "", df_long2$INPUT_POOL) #add input prefix#
df_long2$OUTPUT2 <- gsub("OUTPUT_", "", df_long2$OUTPUT_POOL)#add output prefix#

#change NA or negative to 0#
df_long2$new_FP <- df_long2$FP
df_long2$new_FP[df_long2$FP <0] <- 0
df_long2$new_FP[df_long2$FP == "NA"] <- 0

#vector for input sample names#
INPUTS_VEC <- c("Oral_Inoc_A", "Oral_Inoc_B", "Oral_Inoc_C", "pqr_inoc_a", "pqr_inoc_b", "pqr_inoc_c",
                "s_inoc_a", "s_inoc_b", "s_inoc_c", "v_inoc_a", "v_inoc_b", "v_inoc_c", "f_inoc_a", "f_inoc_b2",
                "f_inoc_c", "g_inoc_a", "g_inoc_b2", "g_inoc_c", "tu_inoc_a", "tu_inoc_b", "tu_inoc_c", "no_inoc_a",
                "IP_Inoc_A", "IP_Inoc_B", "IP_Inoc_C", "iv_inoc_a", "iv_inoc_b", "iv_inoc_c")

INPUTS_VEC2 <- paste("INPUT_", INPUTS_VEC, sep = "")

OUTPUTS_VEC <- paste("OUTPUT_", INPUTS_VEC, sep = "")

#FFP#

#calculate fractional TSFP for INPUT POOL TSFPs in both directions#
#minus one from other#
FP_mat <- as.matrix(subset(FP_df6, select = -c(INPUT_POOL)))
FP_trans <- data.frame(t(subset(FP_df6, select = -c(INPUT_POOL))))
FP_mat2 <- as.matrix(FP_trans)

FP_residual <- FP_mat - FP_mat2 #matrix minus its transposed matrix#
FP_df_resid <- data.frame(FP_residual)
rownames(FP_df_resid) <- rownames(FP_df6)
FP_df_resid$INPUT_POOL <- rownames(FP_df_resid)
#gives a signed value for FP for each combo#

#remove values below 1 (so all negatives plus those which are negligible) to show only positive seeding#
FFP_df_final <- FP_df_resid %>% replace_with_na_all(condition = ~.x < 1)

#long format#
#pivot data longer#


#the new FFP after correction#
#pivot data longer#
df_long3 <- pivot_longer(
  FFP_df_final,
  cols = starts_with("OUTPUT"),
  names_to = "OUTPUT_POOL",
  values_to = "FFP"
)

#change the NAs to 0s for plots#
df_long3$FFP[is.na(df_long3$FFP)] <- 0


df_long3$INPUT2 <- gsub("INPUT_", "", df_long3$INPUT_POOL) #add input prefix#
df_long3$OUTPUT2 <- gsub("OUTPUT_", "", df_long3$OUTPUT_POOL)#add output prefix#

df_long2$INPUT2 <- gsub("INPUT_", "", df_long2$INPUT_POOL) #add input prefix#
df_long2$OUTPUT2 <- gsub("OUTPUT_", "", df_long2$OUTPUT_POOL)#add output prefix#

#remove . from names#
#df_long2$INPUT3 <- gsub("\\.", "", df_long2$INPUT2)
#df_long2$OUTPUT3 <- gsub("\\.", "", df_long2$OUTPUT2)
#df_long3$INPUT3 <- gsub("\\.", "", df_long3$INPUT2)
#df_long3$OUTPUT3 <- gsub("\\.", "", df_long3$OUTPUT2)

#mouse IDs#

df_long2$mouse_ID_in <- parse_number(df_long2$INPUT2)
df_long2$mouse_ID_in[df_long2$INPUT2 %in% INPUTS_VEC] <- NA
df_long3$mouse_ID_in <- parse_number(df_long3$INPUT2)
df_long3$mouse_ID_in[df_long3$INPUT2 %in% INPUTS_VEC] <- NA

df_long2$mouse_ID_out <- parse_number(df_long2$OUTPUT2)
df_long2$mouse_ID_out[df_long2$OUTPUT2 %in% INPUTS_VEC] <- NA
df_long3$mouse_ID_out <- parse_number(df_long3$OUTPUT2)
df_long3$mouse_ID_out[df_long3$OUTPUT2 %in% INPUTS_VEC] <- NA

#group IDs#
df_long2$group_in <- str_extract(df_long2$INPUT2, "^[^[1-9]_]+")
df_long2$group_out <- str_extract(df_long2$OUTPUT2, "^[^[1-9]_]+")
df_long3$group_in <- str_extract(df_long3$INPUT2, "^[^[1-9]_]+")
df_long3$group_out <- str_extract(df_long3$OUTPUT2, "^[^[1-9]_]+")

df_long2$group_in2 <- str_to_upper(df_long2$group_in)
df_long2$group_out2 <- str_to_upper(df_long2$group_out)
df_long3$group_in2 <- str_to_upper(df_long3$group_in)
df_long3$group_out2 <- str_to_upper(df_long3$group_out)

#input groups#
df_long2$INOC <- NULL
df_long2$INOC[df_long2$group_in2 %in% c("P", "Q", "R")] <- "PQR"
df_long2$INOC[df_long2$INPUT2 %in% c("pqr_inoc_a", "pqr_inoc_b", "pqr_inoc_c")] <- "PQR"

df_long2$INOC[df_long2$group_in2 %in% c("S")] <- "S"
df_long2$INOC[df_long2$INPUT2 %in% c("s_inoc_a", "s_inoc_b", "s_inoc_c")] <- "S"

df_long2$INOC[df_long2$group_in2 %in% c("F")] <- "F"
df_long2$INOC[df_long2$INPUT2 %in% c("f_inoc_a", "f_inoc_b2", "f_inoc_c")] <- "F"

df_long2$INOC[df_long2$group_in2 %in% c("T", "U")] <- "TU"
df_long2$INOC[df_long2$INPUT2 %in% c("tu_inoc_a", "tu_inoc_b", "tu_inoc_c")] <- "TU"

df_long2$INOC[df_long2$group_in2 %in% c("G")] <- "G"
df_long2$INOC[df_long2$INPUT2 %in% c("g_inoc_a", "g_inoc_b", "g_inoc_c")] <- "G"

df_long2$INOC[df_long2$group_in2 %in% c("N", "O")] <- "NO"
df_long2$INOC[df_long2$INPUT2 %in% c("no_inoc_a")] <- "NO"

df_long2$INOC[df_long2$group_in2 %in% c("V")] <- "V"
df_long2$INOC[df_long2$INPUT2 %in% c("v_inoc_a", "v_inoc_b", "v_inoc_c")] <- "V"

df_long2$INOC[df_long2$group_in2 %in% c("L", "K", "C")] <- "IP"
df_long2$INOC[df_long2$INPUT2 %in% c("IP_Inoc_A", "IP_Inoc_B", "IP_Inoc_C")] <- "IP"

df_long2$INOC[df_long2$group_in2 %in% c("J", "H")] <- "IV"
df_long2$INOC[df_long2$INPUT2 %in% c("iv_inoc_a", "iv_inoc_b", "iv_inoc_c")] <- "IV"

df_long2$INOC[df_long2$group_in2 %in% c("D", "AO")] <- "ORAL"
df_long2$INOC[df_long2$INPUT2 %in% c("Oral_Inoc_A", "Oral_Inoc_B", "Oral_Inoc_C")] <- "ORAL"

#input groups#
df_long3$INOC <- NULL
df_long3$INOC[df_long3$group_in2 %in% c("P", "Q", "R")] <- "PQR"
df_long3$INOC[df_long3$INPUT2 %in% c("pqr_inoc_a", "pqr_inoc_b", "pqr_inoc_c")] <- "PQR"

df_long3$INOC[df_long3$group_in2 %in% c("S")] <- "S"
df_long3$INOC[df_long3$INPUT2 %in% c("s_inoc_a", "s_inoc_b", "s_inoc_c")] <- "S"

df_long3$INOC[df_long3$group_in2 %in% c("F")] <- "F"
df_long3$INOC[df_long3$INPUT2 %in% c("f_inoc_a", "f_inoc_b2", "f_inoc_c")] <- "F"

df_long3$INOC[df_long3$group_in2 %in% c("T", "U")] <- "TU"
df_long3$INOC[df_long3$INPUT2 %in% c("tu_inoc_a", "tu_inoc_b", "tu_inoc_c")] <- "TU"

df_long3$INOC[df_long3$group_in2 %in% c("G")] <- "G"
df_long3$INOC[df_long3$INPUT2 %in% c("g_inoc_a", "g_inoc_b", "g_inoc_c")] <- "G"

df_long3$INOC[df_long3$group_in2 %in% c("N", "O")] <- "NO"
df_long3$INOC[df_long3$INPUT2 %in% c("no_inoc_a")] <- "NO"

df_long3$INOC[df_long3$group_in2 %in% c("V")] <- "V"
df_long3$INOC[df_long3$INPUT2 %in% c("v_inoc_a", "v_inoc_b", "v_inoc_c")] <- "V"

df_long3$INOC[df_long3$group_in2 %in% c("L", "K", "C")] <- "IP"
df_long3$INOC[df_long3$INPUT2 %in% c("IP_Inoc_A", "IP_Inoc_B", "IP_Inoc_C")] <- "IP"

df_long3$INOC[df_long3$group_in2 %in% c("J", "H")] <- "IV"
df_long3$INOC[df_long3$INPUT2 %in% c("iv_inoc_a", "iv_inoc_b", "iv_inoc_c")] <- "IV"

df_long3$INOC[df_long3$group_in2 %in% c("D", "AO")] <- "ORAL"
df_long3$INOC[df_long3$INPUT2 %in% c("Oral_Inoc_A", "Oral_Inoc_B", "Oral_Inoc_C")] <- "ORAL"

#make a new FFP just for log based plotting which does not have 0 or negatives in it#
df_long3$new_FFP <- df_long3$FFP
df_long3$new_FFP[df_long3$new_FFP == 0] <- 0.9


#plots#
#paper supp figures#
#histograms for thresholding#

ggplot(subset(df_long3, INPUT2 != OUTPUT2 & INPUT2 %in% INPUTS_VEC & INOC == INOC & !(OUTPUT2 %in% INPUTS_VEC)), aes(x = log(new_FFP))) + 
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = 5, size = 1, linetype = "dashed", colour = "red") +
  labs(title = "True migration A to B", x = "\nlog(FFP)", y = "Count\n") +
  theme_bw() + theme(legend.title = element_text(size = 20),
                     legend.position = "bottom",
                     plot.title = element_text(hjust = 0.5, size = 26),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size=24), 
                     axis.line = 
                       element_line(colour = "black"), 
                     axis.text.x = element_text(size = 16, colour = "black"),
                     axis.text.y = element_text(size = 16, colour = "black"),
                     legend.text = element_text(size = 16, colour = "black"))
ggsave("FFP_HISTOGRAM_hotinger.png", height = 6, width = 6)
ggsave("FFP_HISTOGRAM_hotinger.pdf", height = 6, width = 6)


ggplot(subset(df_long3, INPUT2 != OUTPUT2 & !(INPUT2 %in% INPUTS_VEC) & OUTPUT2 %in% INPUTS_VEC), aes(x = log(new_FFP))) + 
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = 5, size = 1, linetype = "dashed", colour = "red") +
  labs(title = "Back migration B to A", x = "\nlog(FFP)", y = "Count\n") +
  theme_bw() + theme(legend.title = element_text(size = 20),
                     legend.position = "bottom",
                     plot.title = element_text(hjust = 0.5, size = 26),
                     panel.border = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.title = element_text(size=24), 
                     axis.line = 
                       element_line(colour = "black"), 
                     axis.text.x = element_text(size = 16, colour = "black"),
                     axis.text.y = element_text(size = 16, colour = "black"),
                     legend.text = element_text(size = 16, colour = "black"))
ggsave("FFP_HISTOGRAM_hotinger2.png", width = 6, height = 6)
ggsave("FFP_HISTOGRAM_hotinger2.pdf", width = 6, height = 6) 


ggplot(subset(df_long3, direction != "NA"), aes(x = log(new_FFP))) + 
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = 5, size = 1, linetype = "dashed", colour = "red") +
  labs(title = "Migration", x = "\nlog(FFP)", y = "Count\n") +
  theme_bw() + theme(legend.title = element_text(size = 20),
                     legend.position = "bottom",
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
                     strip.text = element_text(size = 20)) +
  facet_wrap(~direction, nrow = 2, labeller = as_labeller(c("1TRUE" = "True: A to B", "2BACK" = "Back: B to A")))
ggsave("FFP_HISTOGRAM_hotinger_facet.png", width = 6, height = 10)
ggsave("FFP_HISTOGRAM_hotinger_facet.pdf", width = 6, height = 10) 

ggplot(subset(df_long2, direction != "NA"), aes(x = log(new_FP))) + 
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = 5, size = 1, linetype = "dashed", colour = "red") +
  labs(title = "Migration", x = "\nlog(TSFP)", y = "Count\n") +
  theme_bw() + theme(legend.title = element_text(size = 20),
                     legend.position = "bottom",
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
                     strip.text = element_text(size = 20)) +
  facet_wrap(~direction, nrow = 2, labeller = as_labeller(c("1TRUE" = "True: A to B", "2BACK" = "Back: B to A")))
ggsave("TSFP_HISTOGRAM_hotinger_facet.png", width = 6, height = 10)
ggsave("TSFP_HISTOGRAM_hotinger_facet.pdf", width = 6, height = 10) 

######################################################################################
########################################################################################
#########################################################################################


