# This script generates genus level bar plots

library(ggplot2)
library(tidyr)
library(RColorBrewer)

# Taxa below 1% abundance (total across ALL samples) will be groups as "remainder"
# WITHIN EACH SAMPLE taxa below 1% will be grouped into "remainder"


# Setup ----------------------------

# Load in the data table
final_counts_table <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/decontam_counts_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Define the taxa
tax <- final_counts_table$tax_vector

# Only keep samples from tax table
group_counts <- as.matrix(final_counts_table[, c(grep("RV_no", colnames(final_counts_table)), grep("RV_ox", colnames(final_counts_table)))])


# Data handling ----------------------------

# Sum all the reads by the 6th taxonomic level -> genus (separated by :)
split6 <- sapply(strsplit(as.character(tax), ":"), "[", 6)

# Aggregate data based on genus identifier
agg6 <- aggregate(group_counts, by = list(split6), FUN = sum)
# Move labels to rownames so the data are numeric
rownames(agg6) <- agg6$Group.1
agg6$Group.1 <- NULL

# Convert to percent abundances
agg6p <- apply(agg6, 2, function(x){x/sum(x)})

# 1% abundance
abund <- 0.01

# Define samples with respect to cutoff
agg6.a <- agg6[apply(agg6p, 1, max) >= abund,]	#everything above cutoff
agg6.b <- agg6[apply(agg6p, 1, max) < abund,]	#everything below cutoff

# Sum the remainder (below the cutoff) into one group
x <- colSums(agg6.b, dims=1)

# Get percent abundance to order the taxa
agg6.ap <- apply(agg6.a, 2, function(x){x/sum(x)})

# Order the taxa from most to least abundant. This is so the largest fraction are plotted first
agg6.ao <- agg6.a[order(rowSums(agg6.ap), decreasing = TRUE), ]

# Merge the above cutoff and the "remainder" at the counts level
p <- rbind(agg6.ao, x)
rownames(p)[nrow(p)] <- "Remainder"

# Sum "Remaineder" and "NA"
p2 <- aggregate(p, list(Group = replace(rownames(p),rownames(p) %in% c("NA", "Remainder"), "Remainder")), sum)
row.names(p2) <- p2$Group
p2$Group <- NULL

# Get percent abundances for plotting
y <- as.data.frame(apply(p2, 2, function(x){x/sum(x)}))
y <- y[order(rowSums(y), decreasing = TRUE), ]


# Check sample columns sums to 1 (100%)
colSums(y)
dim(y)



# Graphing ----------------------------


# Convert the data into "long" format for ggplot and 
y$Taxon <- rownames(y)
rownames(y) <- NULL

long_y <- y %>%
  pivot_longer(!Taxon, names_to = "Sample", values_to = "Rel_Abun")

long_y$Time <- ifelse(grepl("t0", long_y$Sample), "0h", ifelse(grepl("t24", long_y$Sample), "24h", "48h"))
long_y$Group <- ifelse(grepl("ox", long_y$Sample), "Oxalate", "No Oxalate")


long_y <- long_y[order(long_y$Rel_Abun, decreasing = TRUE), ]


pal <- c("#F0F0F0","#FF7979","#8BABD3","#FFD757","#DAE3F3","#FFC0CB","#A2CD5A","#CD950C","#9A32CD","#00FF00","#9999CC","#663366","#999966","#9999FF","#FFF5EE","#87CEFF","#FFFF00","#FF0000","#9ACD32","#FA8072","#FFED6F","#B452CD","#7F7F7F","#FFA54F","#66CDAA","#C0C0C0","#27408B","#C71585","#999933","#8B0A50","#FFE7BA")
set.seed(1996)
pal2 <- sample(pal, 22, replace = FALSE)



# RUN THIS GGPLOT
pdf("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/figures/barplot.pdf")
ggplot(long_y, aes(x = Sample, y = Rel_Abun, fill = Taxon)) +
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  geom_col(colour = "black", alpha = 0.7) + # puts outlines on bar plot
  scale_fill_manual(values = pal2) + # 22 colours needed
  labs(y = "Relative Abundance", x = NULL) +
  theme_minimal() +
  theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(legend.text = element_text(face = "italic"), legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.75,"cm")) +
  guides(fill = guide_legend(ncol = 1)) + # Make it one column on the legend
  facet_wrap(~ Group + Time, scales = "free_x", labeller = labeller(Group = label_context, Time = label_both), strip.position = "top")

dev.off()