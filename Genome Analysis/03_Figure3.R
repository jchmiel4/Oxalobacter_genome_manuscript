### This script generates the plots for Figure 3

# Load libraries

library(circlize)
library(stringr)
library(vegan)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(ggforce)
library(gggenes)
library(pairwiseAdonis)
library(viridis)
library(ComplexHeatmap)
library(RColorBrewer)

### Script for PCoA in Figure 3A

# Setup names
strain <- c("BA1", "BA2", "ERR2013569", "HC-1", "HOxHM18", "MGYG-HGUT-01331", "MGYG-HGUT-02505", "HOxNP-1", "HOxNP-2", "HOxNP-3", "OxB", "OxCC13", "OxGP1", "OxK", "OxWR1", "HOxSLD-1", "HOxSLD-2", "SSYG-15", "Va3", "HOxBLS")

species <- c("O. aliformigenes" ,"O. aliformigenes", "O. formigenes", "O. formigenes", "O. formigenes", "O. formigenes", "O. paraformigenes","O. aliformigenes", "O. aliformigenes", "O. aliformigenes", "O. formigenes", "O. formigenes", "O. paeniformigenes", "O. aliformigenes", "O. formigenes", "O. formigenes", "O. formigenes", "O. formigenes", "O. aliformigenes", "O. paraformigenes")


# Setup function for splitting cogs
count_cogs <- function(strain, species) {
  raw_nog <- read.delim(paste0("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/eggnog_mapper/", strain, ".emapper.annotations.tsv"), header = FALSE, comment.char = "#")%>%
    # replace all blank values with "S" in column 7 (COG column)
    mutate(V7 = ifelse(V7 == "-", "S", V7))
  # split doubles up, eg. KL -> K L
  clean_nog <- strsplit(raw_nog$V7, split="") %>% unlist %>%
    table() %>%
    as.data.frame() %>%
    # add genome name and species
    mutate(Genome = strain, Species = species)#%>%
  # Commented out to keep variable
  #write.csv(file = paste0("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/eggnog_mapper/clean_nog/", strain, ".csv"), row.names=FALSE)
}

# Run function for each strain
for (i in 1:20) {
  assign(paste0(strain[i], "_annotation"), count_cogs(strain[i], species[i]))
}

# Merge all cog files
all_nog <- do.call("rbind", mget(paste0(strain, "_annotation"), envir = as.environment(-1)))

# Rename column header
colnames(all_nog) <- c("Function", "Frequency", "Genome", "Taxon")  
write.table(all_nog, file = "~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/eggnog_mapper/clean_nog/all_nog.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# PCoA function
full_pcoa <- all_nog %>%
  #choose columns to work with
  dplyr::select(Function, Frequency, Genome) %>%
  #make pairs across columns with key, value is the value relative to the key
  spread(key = Function, value = Frequency, fill = 0) %>%
  #remove row names
  remove_rownames() %>%
  #convert column to rownames
  column_to_rownames(var = "Genome") %>%
  #converts data table to matrix
  as.matrix() %>%
  #calculate dissimilarity via Bray-Curtis
  vegdist(method = "bray") %>%
  # Run to here for variance
  # max dimension of space that the principal coordinates are represented in
  cmdscale(k = 2) %>%
  #makes data to data frame
  as.data.frame() %>%
  #convert rownames to column
  rownames_to_column(var = "Genome") %>%
  mutate(Taxon = ifelse(Genome %in% c("ERR2013569", "HC-1", "HOxHM18", "MGYG-HGUT-01331", "OxB", "OxCC13", "OxWR1", "HOxSLD-1", "HOxSLD-2", "SSYG-15"), "O. formigenes",
                        ifelse(Genome %in% c("HOxNP-1", "HOxNP-2", "HOxNP-3", "BA1", "BA2", "Va3", "OxK"), "O. aliformigenes",
                               ifelse(Genome %in% "OxGP1", "O. paeniformigenes", "O. paraformigenes"))))

# Plot PCoA
oxf_pcoa <-  ggplot(full_pcoa, aes(x = V1, y = V2, fill = Taxon, color = Taxon)) +
  # plot data points
  geom_point(size = 2) +
  geom_text_repel(data = full_pcoa,
                  segment.color = "transparent",
                  label = full_pcoa$Genome, color = "black", size = 2.65) +
  scale_fill_manual(values = c("#9999FF", "#bf546a", "#CD950C", "#8BABD3")) +
  scale_color_manual(values = c("#9999FF", "#bf546a", "#CD950C", "#8BABD3")) +
  theme_classic() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 1, fill = NA),
        legend.title = element_text(size = 9,  face = "bold"),
        legend.text = element_text(size = 8, face = "italic"),
        legend.key.size = unit(4, "mm", 'lines'),
        legend.key = element_rect(colour = "grey95", fill = "grey95"),
        legend.position = c(0.85, 0.875),
        legend.background = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) 


# Save PCoA
ggsave(oxf_pcoa, filename = "~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/oxf_pcoa.pdf", width = 130, height = 100, units = "mm")  

## PERMANOVA
# Load data
all_nog <- read.table("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/eggnog_mapper/clean_nog/all_nog.txt", header = TRUE, check.names = FALSE, row.names = 1, sep = "\t")

# Convert missing values to zeros and make a matrix
rownames(all_nog) <- NULL
all_nog_matrix <- pivot_wider(all_nog, names_from = Function, values_from = Frequency)
all_nog_matrix[is.na(all_nog_matrix)] <- 0

# Build a metadata table
taxa <- as.data.frame(all_nog_matrix[, c("Genome", "Taxon")])

# Clean up the matrix
all_nog_matrix_clean <- as.data.frame(all_nog_matrix)
all_nog_matrix_clean$Taxon <- NULL
rownames(all_nog_matrix_clean) <- all_nog_matrix_clean$Genome
all_nog_matrix_clean$Genome <- NULL

# Check the homogeneity condition (PERMDISP)
nog_dist <- vegdist(all_nog_matrix_clean, method = "bray")
anova(betadisper(nog_dist, taxa$Taxon))
# Output:
# Pr(>F) = 0.3219

# Run PERMANOVA
set.seed(1996)
nog_permanova <- adonis2(all_nog_matrix_clean ~ Taxon, data = taxa, permuations = 999, by = "terms", method = "bray")
nog_permanova
# Output:
# Pr(>F) = 0.001 ***
# Save PERMANOVA results
write.table(nog_permanova, file="~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/eggnog_mapper/clean_nog/nog_permanova.txt", sep="\t", col.names=NA, quote=F)

## Pairwise PERMANOVA
nog_pwpermanova <- pairwiseAdonis::pairwise.adonis2(all_nog_matrix_clean ~ Taxon, data = taxa, permuations = 999, method = "bray", p.adjust.m = "BH")
nog_pwpermanova
# Write output
cat(capture.output(print(nog_pwpermanova), file="~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/eggnog_mapper/clean_nog/nog_pwpermanova.txt"))


### Script for heatmap in Figure 3B

cogs <- read.table("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/eggnog_mapper/clean_nog/cogs3.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)

# Save file
pdf("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/cog_heatmap.pdf", width = 3, height = 2.5)

cn = colnames(cogs)
ComplexHeatmap::pheatmap(
  mat = as.matrix(cogs), 
  gaps_row = c(10, 17, 18),
  color = (colorRampPalette(brewer.pal(9,"RdYlBu"))(400)),
  border_color = "grey60",
  fontsize_row = 13,
  cellwidth=12, cellheight=12,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = F,
  heatmap_legend_param = list(title = "Number \n\ of COGS", legend_height = unit(7, "cm")),
  bottom_annotation = HeatmapAnnotation(
    text = anno_text(cn, rot = 0, location = unit(0.5, "npc"), just = "centre", gp = gpar(fontsize = 12)),
    annotation_height = max_text_width(cn)),
)
dev.off()


### Script for Cas gene ribbon structures in Figure 3C
# Load data
cas <- read.table("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/CRISPRCasFinder/cas_regions3.txt", header = TRUE, check.names = FALSE, row.names = NULL, sep = "\t")

# Plot cas sequences
cas_plot <- ggplot(cas, aes(xmin = start, xmax = end, y = genome, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), arrowhead_width = unit(1, "mm"), arrow_body_height = unit(5, "mm")) +
  #ggplot2::layer(params = list(max.size = 6)) +
  geom_gene_label(fontface = "italic", size = 12, padding.x = unit(0.75, "mm"), padding.y = unit(0, "mm"), grow = FALSE) +
  facet_wrap(~ genome, scales = "free_y", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes() +
  theme(legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 10, face = "italic", margin = margin(b=5)))

# Save plot
ggsave(cas_plot, filename = "~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/cas_plot2.pdf", width = 250, height = 90, units = "mm")

### Script for circular heatmap in Figure 3D
# Load in data
mm <- read.table("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/genome_facts/genome_facts.txt", header = TRUE, check.names = FALSE, sep = "\t", row.names = 1)

group <- as.data.frame(mm[, "Group"])
rownames(group) <- rownames(mm)
colnames(group) <- "Group"

bacteriocins <- as.data.frame(mm[, c(grep("Sactipeptide", colnames(mm)), grep("Lantipeptide", colnames(mm)))])
abx <- as.data.frame(mm[, c(grep("bla", colnames(mm)), grep("tetW", colnames(mm)))])


# Plot bacteriocins and antibiotic resistance
pdf("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/circ_features.pdf")

circos.par(start.degree = 0, gap.degree = c(7, 7, 7, 7, 90))
col_bacteriocins = colorRamp2(c(0, 1), c("#78A2CC", "#3548BA"))
circos.heatmap(rev(bacteriocins), col = col_bacteriocins, cluster = FALSE, rownames.side = "outside", track.height = 0.14, rownames.cex = 1, split = group$Group)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 5) {
    cn = colnames(bacteriocins)
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(0.01, "mm"), 
                1:n - 0.5, cn, cex = 1, adj = c(-0.01, 0.5), facing = "inside")
  }
}, bg.border = NA)
col_abx = colorRamp2(c(0, 1, 2), c("#E48290", "#ED4750", "#F60C10"))
circos.heatmap(abx, col = col_abx, cluster = FALSE, track.height = 0.14, split = group$Group)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 5) { 
    cn = rev(colnames(abx))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(0.01, "mm"),
                1:n - 0.5, cn, cex = 1, adj = c(-0.01, 0.5), facing = "inside")
  }
}, bg.border = NA)
legend(-0.97, 1, legend = c("0", "1"), title = "Bacteriocins", fill = c("#78A2CC", "#3548BA"), box.lty = 0, cex = 1.15, xjust = 0, title.adj = 4)
legend(-1.07, 0.735, legend = c("0", "1", "2"), title = "Antibiotic reistance", fill = c("#E48290", "#ED4750", "#F60C10"), box.lty = 0, cex = 1.15, xjust = 0, title.adj = 7.5)

circos.clear()
dev.off()

### Script for presence/absence plot in Figure 3E

# Load in data
prophage <- read.csv("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/genome_facts/prophage_summary_clean.csv", check.names = FALSE, header = TRUE)

# Keep prophages with completeness of Intact or Questionable
prophage_filt <- prophage[c(grep("Intact", prophage$completeness), grep("Questionable", prophage$completeness)), ]

# Plot data
prophage_plot <- ggplot(prophage_filt, aes(x = prophage, y = factor(strain, levels = rev(levels(factor(strain)))), fill = completeness)) +
  geom_tile(color = "black") +
  scale_fill_grey(name = "Completeness") +
  labs(y = "", x = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 50, hjust = 1, size = 11, vjust = 1),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(face = "italic", size = 10)
  ) +
  ggforce::facet_col(vars(group), scales = 'free_y', space = 'free')


# Save plot
ggsave(prophage_plot, filename = "~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/raw_figures/prophages.pdf", height = 20, width = 12, units = "cm")
