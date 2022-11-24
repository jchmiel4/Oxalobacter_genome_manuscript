### This script generates the plots for Figure 1

# Load libraries

library(ggtree)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(UpSetR)
library(treeio)

### Script for phylogenetic tree in Figure 1A

# Load files
core_tree_file <- read.tree("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/Orthofinder_RAxML_new/RAxML_bestTree.core")
core_meta <- read.csv("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/Orthofinder_RAxML_new/tree_meta.csv")

# Load bootstrap file
bootstrap_file <- read.raxml("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/Orthofinder_RAxML_new/RAxML_bipartitionsBranchLabels.core")

# Combine the bootstrap data
bootstrap <- data.frame(bootstrap_file@data[["node"]], bootstrap_file@data[["bootstrap"]])
colnames(bootstrap) <- c("node", "bs_val")

# Prep the filler/dummy data
node_fill <- 1:23
bs_fill <- numeric(23)
fill <- data.frame(node_fill, bs_fill)
colnames(fill) <- c("node", "bs_val")

# Merge them all
bootstrap_all <- rbind(fill, bootstrap)
bootstrap_all[is.na(bootstrap_all)] <- 0

# Generate tree
core_tree <- ggtree(core_tree_file, size = 0.25) %<+% core_meta +
  theme(legend.position=c(0.12,0.9), legend.text = element_text(size = 7.5), legend.title = element_text(size=7.5), legend.key.size = unit(1, "mm", 'lines')) +
  guides(color = guide_legend(override.aes = list(size = 1.5, shape = 15))) +
  scale_colour_manual(name = "Original Grouping", values = c("gray5", "gray40", "gray65")) +
  geom_treescale(x = 0, y = -0.5,linesize = 0.35, offset=0.15, width = 0.05, fontsize = 2.25) +
  geom_tippoint(aes(color=group), size=-100) +
  geom_tiplab(aes(subset=!grepl("WoOx3|B_cepacia", label), color=group), hjust=0, offset = 0.001, show.legend = F, align = T, size = 2.5, linesize = 0.25) +
  geom_tiplab(aes(subset = grepl("WoOx3", label), label = paste0('italic("O. vibrioformis")~', "WoOx3")), parse = TRUE, hjust = 0, offset = 0.001, show.legend = F, align = T, size = 2.5, linesize = 0.25) +
  geom_tiplab(aes(subset=grepl("B_cepacia", label), label = paste0('italic("B. cepacia")~', '"ATCC 25416"','')), parse = TRUE, hjust = 0, offset = 0.001, show.legend = F, align = T, size = 2.5, linesize = 0.25) +
  # Species strips
  geom_strip("HC-1", "ERR2013569", barsize = 0.35, offset = 0.265, color = "#bf546a", label = "O. formigenes", offset.text = 0.01, fontsize = 2.5, extend = 0.4) + 
  geom_strip("HOxNP-1", "BA1", barsize = 0.35, offset = 0.265, color = "#9999ff", label = "O. aliiformigenes", offset.text = 0.01, fontsize = 2.5, extend = 0.4) +
  geom_strip("OxGP1", "OxGP1", barsize = 0.35, offset = 0.265, color = "#CD950C", label = "O. paeniformigenes", offset.text = 0.01, fontsize = 2.5, extend = 0.4) +
  geom_strip("MGYG-HGUT-02505", "HOxBLS", barsize = 0.35, offset = 0.265, color = "#8BABD3", label = "O. paraformigenes", offset.text = 0.01, fontsize = 2.5, extend = 0.4) +
  xlim(0, 1.3) +
  geom_nodepoint(aes(label = bootstrap_all$bs_val, subset = as.numeric(bootstrap_all$bs_val) > 80), size = 1.5, color="steelblue", alpha = 0.75, pch = 16)

# Save tree
ggsave(core_tree, filename = "~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/core_phylo_tree.pdf", height = 80, width = 140, units = "mm")


### Script for ANI heatmap tree in Figure 1B

# Load files
ani <- read.table("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/fastani/oxf_ani.txt.matrix_clean.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)

# Save file
# 80 mm = 3.15 inches
pdf("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/ani_heatmap.pdf", width = 2.5, height = 2.5)

cn = colnames(ani)
ComplexHeatmap::pheatmap(
  mat = as.matrix(ani), 
  #color = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(400),
  color = colorRampPalette(brewer.pal(9,"Reds"))(400),
  #color = colorRampPalette(brewer.pal(n = 3, name = "YlOrRd"))(400),
  border_color = "grey60",
  fontsize_row = 13,
  cellwidth=12, cellheight=12,
  show_column_dend = FALSE,
  show_colnames = FALSE,
  heatmap_legend_param = list(title = "ANI", legend_height = unit(7, "cm")),
  bottom_annotation = HeatmapAnnotation(
    text = anno_text(cn, rot = 55, location = unit(1, "npc"), just = "right", gp = gpar(fontsize = 12)),
    annotation_height = max_text_width(cn)),
)
dev.off()



### Script for UpSetR plot of pangenome tree in Figure 1C

# Load presence absence from Roary
genes <- read.csv("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/upsetR/gene_presence_absence_upsetR.csv", header = TRUE, sep = ",", check.names = FALSE)

# Plot UpSetR graph
pdf("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/upsetR.pdf", height = 10, width = 6)

upset(genes, keep.order = F, order.by = c("freq"), nsets = 23, sets.x.label = "Number of CDS", mainbar.y.label = "CDS Intersections", queries = list(
  list(query = intersects, params = list("ERR2013569", "HOxHM18", "MGYG_HGUT_01331", "OxB", "OxCC13", "HOxSLD_1", "HOxSLD_2", "OxWR1", "HC_1", "SSYG_15"), color = "#bf546a", active = T, query.name = "Species 1"),	
  list(query = intersects, params = list("HOxNP_1", "HOxNP_2", "HOxNP_3", "BA1", "BA2", "Va3", "OxK"), color = "#9999ff", active = T, query.name = "Species 2"),
  list(query = intersects, params = list("OxGP1"), color = "#CD950C", active = T, query.name = "Species 3"),
  list(query = intersects, params = list("HOxBLS", "MGYG_HGUT_02505"), color = "#8BABD3", active = T, query.name = "Species 4")
), 
query.legend = "none", intersections = list(
  list("WoOx3"),
  list("HOxBLS", "MGYG_HGUT_02505"),
  list("OxGP1"),
  list("HOxNP_1", "HOxNP_2", "HOxNP_3", "BA1", "BA2", "Va3", "OxK"),
  list("ERR2013569", "HOxHM18", "HRGM_Genome_1325", "MGYG_HGUT_01331", "OxB", "OxCC13", "HOxSLD_1", "HOxSLD_2", "OxWR1", "HC_1", "SSYG_15"),
  list("HRGM_Genome_1325"),
  list("BA1", "BA2"),
  list("ERR2013569", "HOxHM18", "MGYG_HGUT_01331", "OxB", "OxCC13", "HOxSLD_1", "HOxSLD_2", "OxWR1", "HC_1", "SSYG_15"),
  list("ERR2013569", "HOxHM18", "HRGM_Genome_1325", "MGYG_HGUT_01331", "OxB", "OxCC13", "HOxSLD_1", "HOxSLD_2", "OxWR1", "HC_1", "SSYG_15", "HOxNP_1", "HOxNP_2", "HOxNP_3", "BA1", "BA2", "Va3", "OxK", "OxGP1", "HOxBLS", "MGYG_HGUT_02505"),
  list("HOxNP_1", "HOxNP_2", "HOxNP_3",  "Va3", "OxK"),
  list("ERR2013569", "HOxHM18", "HRGM_Genome_1325", "MGYG_HGUT_01331", "OxB", "OxCC13", "HOxSLD_1", "HOxSLD_2", "OxWR1", "HC_1", "SSYG_15", "WoOx3", "HOxBLS", "MGYG_HGUT_02505", "OxGP1", "HOxNP_1", "HOxNP_2", "HOxNP_3", "BA1", "BA2", "Va3", "OxK"),
  list("ERR2013569", "HOxHM18", "HRGM_Genome_1325", "MGYG_HGUT_01331", "OxCC13", "HOxSLD_1", "HOxSLD_2", "HC_1", "SSYG_15", "HOxNP_1", "HOxNP_2", "HOxNP_3", "BA1", "BA2", "Va3", "OxK", "HOxBLS", "MGYG_HGUT_02505")
), 
text.scale =  c(1.75, 1.8, 1.6, 1.6, 1.7, 1.75),
# axis title / axis/ small axis title / small axis / names / bar number
mb.ratio = c(0.59, 0.41)
) 
  
legend("topright", legend = c("Species_1", "Species_2", "Species_3", "Species_4"), fill = c("#bf546a", "#9999ff", "#CD950C", "#8BABD3"), box.lty = 0, cex = 1, xjust = 0, title.adj = 0)

dev.off()
