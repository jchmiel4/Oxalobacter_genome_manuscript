# This script filters out low-abundance and rare ASVs

# Load in counts table from dada2
dada2_counts <- read.table("~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/demultiplex_dada2/dada2_nochim_tax.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Make a taxa table with the SVs
taxa <- data.frame(dada2_counts$tax_vector) 
rownames(taxa) <- rownames(dada2_counts)

# Remove tax_vector from dada2_counts
dada2_counts <- dada2_counts[, !names(dada2_counts) %in% ("tax_vector")]


# Find initial samples, ASVs, and reads
sum(dada2_counts) # Output: 500129
dim(dada2_counts) # Output: 441  24

# Generate a relative abundance table to remove ASVs accounting for less than 0.1% of reads in every sample
props <- apply(dada2_counts, 2, function(x) {x/sum(x)})
filtered_by_props <- dada2_counts[apply(props, 1, max) >= 0.001, ]

# Re-check dimensions
sum(filtered_by_props) # Output: 489013
dim(filtered_by_props) # Output: 160  24

any(colSums(filtered_by_props) == 0) # Output: FALSE

# Paste on tax table
filtered_counts <- merge(filtered_by_props, taxa, by = 0)
rownames(filtered_counts) <- filtered_counts$Row.names
filtered_counts$Row.names <- NULL

# Write output table
write.table(filtered_counts, file = "~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/filtered_counts.txt", sep = "\t", col.names = NA, quote = FALSE)
