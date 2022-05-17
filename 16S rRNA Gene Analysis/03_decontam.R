# This scripts executes the decontam pipeline --> looks for ASVs that are contaminants (based on control samples)

# Load libraries
library(phyloseq)
library(decontam) 

# Load in the filtered counts table
filtered_counts <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/filtered_counts.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Separate taxa 
filtered_taxa <- as.matrix(filtered_counts$dada2_counts.tax_vector)
rownames(filtered_taxa) <- rownames(filtered_counts)

# Remove taxa from counts table
filtered_counts <- filtered_counts[, !names(filtered_counts) %in% ("dada2_counts.tax_vector")]

# Load in meta data table
metadata <- read.table("~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/meta_data.txt", header = T, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE, row.names = 1, comment.char = "")

s# Make phyloseq parts
count_phylo <- otu_table(filtered_counts, taxa_are_rows = TRUE)
tax_phylo <- tax_table(filtered_taxa)
meta_phylo <- sample_data(metadata)
# Merge them all together
physeq <- phyloseq(count_phylo, tax_phylo, meta_phylo)

# Inspect library size
df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize), ]
df$Index <- seq(nrow(df))
pdf("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/library_size.pdf")
ggplot(data = df, aes(x = Index, y = LibrarySize, color = Sample_or_Control)) + geom_point()
dev.off()


# Identify contaminants by prevalence ---------------------------

# Add on TRUE (control) and FALSE (sample) into $is.neg
sample_data(physeq)$is.neg <- sample_data(physeq)$Sample_or_Control == "Control"

# Run decontam, batch is based on PCR plate
contamdf_prev <- isContaminant(physeq, method = "prevalence", neg = "is.neg", batch = physeq@sam_data[["PCR Plate"]], threshold = 0.5)

# Overview table
table(contamdf_prev$contaminant)
# Output: 
# FALSE  TRUE 
# 112    48 

# Create vector of non-contam ASVs
not_contam <- contamdf_prev[contamdf_prev$contaminant == "FALSE", ]
not_contam <- rownames(not_contam)

# Remove the contaminating samples
decontam_counts <- filtered_counts[(rownames(filtered_counts) %in% not_contam), ]

# Add on tax table
decontam_tax <- filtered_taxa[(rownames(filtered_taxa) %in% not_contam), ]

# Merge into one
decontam_counts_table <- merge(decontam_counts, decontam_tax, by = 0)
rownames(decontam_counts_table) <- decontam_counts_table$Row.names
decontam_counts_table$Row.names <- NULL
colnames(decontam_counts_table)[ncol(decontam_counts_table)] <-"tax_vector"

# Write output table
write.table(decontam_counts_table, file = "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/decontam_counts_table.txt", sep = "\t", col.names = NA, quote = FALSE)