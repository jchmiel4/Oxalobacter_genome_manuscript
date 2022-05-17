# This script runs generates a PCA biplot


# Setup ------------------------------------------------

# Load libraries
library(zCompositions)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)

# Load in the data table
final_counts_table <- read.table("~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/decontam_counts_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Only keep samples from tax table
group_counts <- as.matrix(final_counts_table[, c(grep("RV_no", colnames(final_counts_table)), grep("RV_ox", colnames(final_counts_table)))])


# Load in the metadata table
metadata <- read.table("~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/meta_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE, row.names = 1, comment.char = "")

# Only keep samples from metadata table
group_metadata <- as.data.frame(metadata[c(grep("RV_no", rownames(metadata)), grep("RV_ox", rownames(metadata))), ])


# Define the taxa
tax <- final_counts_table$tax_vector

# Split to the 6th taxonomic level -> genus (separated by :)
split6 <- sapply(strsplit(as.character(tax), ":"), "[", 6)
split6 <- as.data.frame(split6)

# Tack on SV number
rownames(split6) <- rownames(group_counts)



# Setup and perform PCA ----------------------------

# Perform Bayesian-Multiplicative replacement of count zeros.
group_czm <- cmultRepl(t(group_counts), label = 0, method = "CZM")
# CLR transform the data
group_clr <- t(apply(group_czm, 1, function(x){log(x) - mean(log(x))}))
# perform PCA
group_pcx <- prcomp(group_clr)


# Generate figure ----------------------------

# Define parameters for PCA
sample_positions <- data.frame(group_pcx[["x"]])
sv_positions <- data.frame(group_pcx[["rotation"]])
sample_names <- data.frame(rownames(group_pcx[["x"]])) 
colnames(sample_names) <- "samples"

sample_names$group <- ifelse(grepl("RV_no", sample_names$samples), "no", "ox")
sample_names$time <- ifelse(grepl("-t0-", sample_names$samples), "t0", ifelse(grepl("-t24-", sample_names$samples), "t24", "t48"))
sample_names$grouping <- ifelse(grepl("no-t0-", sample_names$samples), "t0_no", ifelse(grepl("no-t24-", sample_names$samples), "t24_no", ifelse(grepl("no-t48-", sample_names$samples), "t48_no", ifelse(grepl("ox-t0-", sample_names$samples), "t0_ox", ifelse(grepl("ox-t24-", sample_names$samples), "t24_ox", "t48_ox")))))

# Merge on genus table
sv_pos <- merge(sv_positions, split6, by = 0)

# Calculate Euclidean distance
arrow_len <- function(x, y) {
  sqrt((x - 0)^2 + (y - 0)^2)
}

sv_pos_dist <- mapply(arrow_len, sv_pos$PC1, sv_pos$PC2)
sv_pos_dist <- as.data.frame(cbind(sv_pos_dist, sv_pos$Row.names, sv_pos$split6, sv_pos$PC1, sv_pos$PC2))
colnames(sv_pos_dist) <- c("Distance", "SV", "Genus", "PC1", "PC2")


# Filter arrow based on distance
filter <- sv_pos_dist %>% filter(Distance >= 0.15) %>% filter(Genus != "NA") %>% bind_rows(sv_pos_dist["67", ])


# Plot PCA
pdf("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/figures/PCA_genus_test.pdf")
ggplot(sample_positions, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.9, aes(fill = sample_names$group, shape = sample_names$time), size = 6) + # Plot samples
  geom_segment(data = sv_pos, aes(x = 0, y = 0, xend = 18 * PC1, yend = 18 * PC2),
               arrow = arrow(length = unit(1/2, 'picas')),
               color = "grey69", alpha = 0.8, size = 0.15) + # Plot features
  geom_text_repel(data = filter, aes(x = 18 * as.numeric(PC1), y = 18 * as.numeric(PC2)), nudge_x = 0.25, nudge_y = 0.25,
                  segment.color = 'transparent',
                  force_pull = 50,
                  force = 7,
            label = filter$Genus, color = "black", size = 3, fontface = "italic") +
  scale_fill_manual(name = "Group", labels = c("No oxalate", "Oxalate"), values = c("#7470B1", "#F49A7D")) +
  guides(fill = guide_legend(override.aes = list(shape = 16, size = 4, colour = c("#7470B1", "#F49A7D")))) +
  guides(shape = guide_legend(override.aes = list(size=4))) +
  scale_shape_manual(name = "Time point", labels = c("0h", "24h", "48h"), values = c(21, 24, 22)) +
  labs(
    x = paste("PC1: ", sprintf("%.1f", group_pcx$sdev[1]^2/sum(group_pcx$sdev^2) * 100), "%"),
    y = paste("PC2: ", sprintf("%.1f", group_pcx$sdev[2]^2/sum(group_pcx$sdev^2) * 100), "%"),
    shape = "Time point") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 1, fill = NA),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = c(0.1, 0.825),
        legend.text = element_text(size = 11),
        axis.text=element_text(size=12)
        )
dev.off()
