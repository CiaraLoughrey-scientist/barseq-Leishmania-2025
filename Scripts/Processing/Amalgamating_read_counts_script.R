#this script creates the dataframe from the read counts in the RUST script. You must use the RUST script on the FASTQ data before running this script#

dataset_label <- "2025-JCM-002"
io <- list(
    "all_barcodes" = file.path("..", "metadata", "all-barcodes.txt"),
    "dataset_dir" = file.path("..", "data", dataset_label)
)
# Load all the barcodes:
all_barcodes <- readLines(io$all_barcodes)
# Find all of the sample input files:
sample_names <- gsub("-barcodes.txt", "", list.files(file.path(io$dataset_dir, "barcodes"), pattern = "*-barcodes.txt"), fixed = TRUE)
# Load all the data:
message("loading barcode counts...")
d_raw <- sapply(sample_names, function(sample_name){
    barcode_path <- file.path(io$dataset_dir, "barcodes", sprintf("%s-barcodes.txt", sample_name))
    res <- read.table(barcode_path, sep = "\t", header = FALSE, colClasses=c("numeric", "character"), col.names=c("count", "barcode"), row.names=2)
    return(res)
}, simplify=FALSE)
# Count the number of unmatched reads and save these to file:
message("loading unmatched counts...")
n_total <- sapply(d_raw, function(i){sum(i[,1])})
n_unmatched <- sapply(d_raw, function(i){i["no_barcode", 1]})
counts <- data.frame(
    matched = n_total - n_unmatched,
    unmatched = n_unmatched,
    row.names = names(d_raw)
)

# Build a matrix of all our data:
message("merging barcode data...")
barcode_counts <- data.frame(
    row.names = all_barcodes
)
for(i in names(d_raw)){
    x <- d_raw[[i]]
    barcode_counts[[i]] <- d_raw[[i]][rownames(barcode_counts), ]
}
barcode_counts[is.na(barcode_counts)] <- 0
counts$canonical <- colSums(barcode_counts) # NEW: Add the number of "canonical barcodes"
# Write the output counts to file:
write.csv(barcode_counts, file.path(io$dataset_dir, sprintf("%s-barcodes.csv", dataset_label)))

write.csv(counts, file=file.path(io$dataset_dir, sprintf("%s-counts.csv", dataset_label)))
