### This script generates the plots for Figure 2

# Load libraries

library(ggtree)
library(tidyverse)
library(ggmsa)
library(Biostrings)
library(treeio)

### Script for phylogenetic tree in Figure 2A

# Load data
rRNA_tree_file <- read.tree("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/oxf_16S_seqs/16S_rRNA_new_tree/RAxML_bestTree.subtree")
rRNA_meta <- read.csv("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/oxf_16S_seqs/rRNA_meta.csv")

# Load and format bootstrapping values
# Load bootstrap file
bootstrap_file <- read.raxml("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/oxf_16S_seqs/16S_rRNA_new_tree/RAxML_bipartitionsBranchLabels.subtree")

# Combine the bootstrap data
bootstrap <- data.frame(bootstrap_file@data[["node"]], bootstrap_file@data[["bootstrap"]])
colnames(bootstrap) <- c("node", "bs_val")

# Prep the filler/dummy data
node_fill <- 1:26
bs_fill <- numeric(26)
fill <- data.frame(node_fill, bs_fill)
colnames(fill) <- c("node", "bs_val")

# Merge them all
bootstrap_all <- rbind(fill, bootstrap)
bootstrap_all[is.na(bootstrap_all)] <- 0

# Check nodes
#ggtree(rRNA_tree_file, size = 0.25) %<+% rRNA_meta +
#geom_text(aes(label = node))

# Generate 16S rRNA tree
rRNA_tree <- ggtree(rRNA_tree_file, size = 0.25) %<+% rRNA_meta +
  guides(color = guide_legend(override.aes = list(size = 1.5, shape = 15))) +
  scale_colour_manual(name="New Species Designation", values = c("#9999ff", "#bf546a", "#CD950C", "#8BABD3")) +
  geom_tippoint(aes(color = group), size = -100) +
  geom_tiplab(aes(subset = !grepl("WoOx3|B_cepacia_ATCC_25416", label), color = group), hjust = 0, offset = 0.001, show.legend = F, size = 2.5, linesize = 0.25) +
  geom_tiplab(aes(subset=grepl("WoOx3", label), label = paste0('italic("O. vibrioformis")~', "WoOx3",'')), parse = TRUE, hjust = 0, offset = 0.001, show.legend = F, size = 2.5, linesize = 0.25) +
  xlim(0, 0.151) +
  geom_tiplab(aes(subset=grepl("B_cepacia_ATCC_25416", label), label = paste0('italic("B. cepacia")~', '"ATCC 25416"','')), parse = TRUE, hjust = 0, offset = 0.001, show.legend = F, size = 2.5, linesize = 0.25) +
  theme(legend.position=c(0.2,0.9), 
        legend.text = element_text(size = 7.25, face = "italic"), 
        legend.title = element_text(size=7.25, face = "bold"), 
        legend.key.size = unit(3, "mm", 'lines'), 
        legend.background=element_blank()) +
  geom_nodepoint(aes(label = bootstrap_all$bs_val, subset = as.numeric(bootstrap_all$bs_val) > 50), size = 1.5, colour ="steelblue", alpha = 0.75, pch = 16) +
  geom_treescale(x = -0.0, y = 19, fontsize = 2.5, linesize = 0.25)


# Save tree
ggsave(rRNA_tree, filename = "~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/rRNA_tree.pdf", height = 80, width = 80, units = "mm")


### Script for multiple scequence alignment in Figure 2B

# Load in data
DNA_seqs <- readDNAMultipleAlignment("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/oxf_16S_seqs/consensus_species_all.fna")

# Generate ggmsa plot
pdf("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/16S_align2.pdf")
ggmsa(DNA_seqs, start = 445, end = 472, char_width = 0.5, seq_name = TRUE, use_dot = TRUE, disagreement = TRUE, color = "Chemistry_NT", consensus_views = TRUE)
dev.off()