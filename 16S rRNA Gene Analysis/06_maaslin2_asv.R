# This script run MaAsLin2 for ASV abundance (P values only)

# Load libraries
library(Maaslin2)
library(zCompositions)


# Load in the data table
final_counts_table <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/decontam_counts_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Only keep samples from tax table
final_counts <- as.matrix(final_counts_table[,c(grep("RV_no", colnames(final_counts_table)), grep("RV_ox", colnames(final_counts_table)))])


# Load in the metadata table
metadata <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/meta_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = F, row.names = 1, comment.char = "")

# Only keep samples from metadata table
final_metadata <-  as.data.frame(metadata[c(grep("RV_no", rownames(metadata)), grep("RV_ox", rownames(metadata))), ])



# Perform Bayesian-Multiplicative replacement of count zeros.
final_czm <- cmultRepl(t(final_counts), label = 0, method = "CZM")
# CLR transform the data
final_clr <- t(apply(final_czm, 1, function(x){log(x) - mean(log(x))}))

maaslin2_results <- Maaslin2(
  input_data = final_czm,
  input_metadata = final_metadata,
  output = "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/maaslin2_czm",
  transform = "NONE",
  normalization = "CLR",
  standardize = FALSE,
  fixed_effects = c("Group"),
  reference = c("no"),
  random_effects = c("Time"),
  correction = "BH",
  min_prevalence = 0,
  analysis_method = "LM")
