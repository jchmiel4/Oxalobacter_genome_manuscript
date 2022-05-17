# This script makes the differential abundance figure (using ANCOM as base and confirming results with MaAsLin2)

# Load libraries
library(ggplot2)



# Load in manually curated table of ANCOM and MaAsLin2 results
da_table <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/differential_abundance_output.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)


# Plot results of ANCOM and MaAsLin2 comparison of no oxalate vs oxalate
ggsave("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/figures/da_strip_small.pdf", width = 9, height = 3)
ggplot(da_table, aes(x = com_id, y = ANCOM_CLR)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey77", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "grey77", size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey77", size = 0.5) +
  geom_point(alpha = 0.8) +
  coord_flip(clip = "off") +
  scale_x_discrete(limits = rev) +
  annotate("segment", y = 0.1, yend = 0.9, x = 0.6, xend = 0.6, arrow = arrow(length = unit(0.03, "npc")), lineend = "butt", linejoin = "round", col = "gray75") +
  annotate("text", x = 0.77, y = 0.5, label = "More Abundant with Oxalate", size = 2, col="gray69") +
  annotate("segment", y = -0.1, yend = -0.9, x = 0.6, xend = 0.6, arrow = arrow(length = unit(0.03, "npc")), lineend = "butt", linejoin = "round", col = "gray75") +
  annotate("text", x = 0.77, y = -0.5, label = "More Abundant without Oxalate", size = 2, col="gray69") +
  labs(x = "Taxonmic Identity", y = "ANCOM Effect Size") +
  theme_bw() +
  theme(panel.spacing.x = unit(2, "lines"),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(face = "italic", colour = "gray21"),
        axis.text.x = element_text(face = NULL, colour = "gray21"))
dev.off()




# Relative abundance of Oxalobacter only ----------------------------

# Load in counts table
final_counts_table <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/decontam_counts_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Only take samples
final_counts <- as.matrix(final_counts_table[, c(grep("RV_no", colnames(final_counts_table)), grep("RV_ox", colnames(final_counts_table)))])

# Determine relative abundance
prop_counts <- apply(final_counts, 2, function(x) {x/sum(x)})

# Only take Oxalobacter numbers
oxalobacter_counts <- as.matrix(prop_counts[grep("SV_51", rownames(prop_counts)), c(grep("RV_no", colnames(prop_counts)), grep("RV_ox", colnames(prop_counts)))])