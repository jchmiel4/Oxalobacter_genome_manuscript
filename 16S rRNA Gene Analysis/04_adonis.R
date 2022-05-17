# This script runs adonis (PERANOVA) to determine if the centroids and/or dispersion of the data are different


# Setup ------------------------------------------------

# Load libraries
library(zCompositions)
library(vegan)

# Load in the data table
final_counts_table <- read.table("~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/decontam_counts_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Load in the metadata table
metadata <- read.table("~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/meta_data.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", check.names = F, row.names = 1, comment.char = "")



# Check samples vs control (PERMANOVA) ----------------------------

# Remove taxa column
all_counts <- as.matrix(final_counts_table[, c(grep("RV_no", colnames(final_counts_table)), grep("RV_ox", colnames(final_counts_table)), grep("SW", colnames(final_counts_table)))])

# Can't have zeroes. Bayesian-Multiplicative replacement of count zeros transpose so samples are as rows for this function (original table has samples as cols)
all_czm <- cmultRepl(t(all_counts), label = 0, method = "CZM")
# The table needs to be transposed again (samples as COLUMNS) clr transform the transposed data 
all_clr <- t(apply(all_czm, 1, function(x){log(x) - mean(log(x))}))

# Check the homogeneity condition
all_dist <- vegdist(all_clr, method = "euclidean")
anova(betadisper(all_dist, metadata$Sample_or_Control))
# Output:
# Pr(>F) = 0.7929

# Perform PERMANOVA
set.seed(1996)
all_permanova <- adonis2(all_clr ~ Sample_or_Control, data = metadata, permuations = 999, by = NULL, method = "euclidean")
all_permanova
# Output:
# Pr(>F) = 0.001
# Save permanova results
write.table(all_permanova, file = "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/samples_control_permanova.txt", sep = "\t", col.names = NA, quote = F)



# Check oxalate vs control while blocking for time (PERMANOVA) ----------------------------

# Only keep samples from tax table
group_counts <- as.matrix(final_counts_table[, c(grep("RV_no", colnames(final_counts_table)), grep("RV_ox", colnames(final_counts_table)))])

# Only keep samples from metadata table
group_metadata <-  as.data.frame(metadata[c(grep("RV_no", rownames(metadata)), grep("RV_ox", rownames(metadata))), ])

# Perform Bayesian-Multiplicative replacement of count zeros.
group_czm <- cmultRepl(t(group_counts),  label = 0, method = "CZM")
# CLR transform the data
group_clr <- t(apply(group_czm, 1, function(x){log(x) - mean(log(x))}))

# Check the homogeneity condition
group_dist <- vegdist(group_clr, method = "euclidean")
anova(betadisper(group_dist, group_metadata$Group))
# Output:
# Pr(>F) = 0.1945
anova(betadisper(group_dist, group_metadata$Time))
# Output:
# Pr(>F) = 0.9715

# Perform PERMOANOVA with blocking by time
perm <- how(nperm = 999)
# This is the blocking design
setBlocks(perm) <- with(group_metadata, Time)
set.seed(1996)
group_permanova <- adonis2 (group_clr ~ Group, data = group_metadata, permutations = perm, method = "euclidean", by = "term")
group_permanova
# Output:
# Group     Pr(>F) = 0.001


# Save PERMANOVA results
write.table(group_permanova, file = "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/group_permanova.txt", sep = "\t", col.names = NA, quote = F)