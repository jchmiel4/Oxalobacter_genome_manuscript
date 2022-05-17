# This script runs the preparation part of PICRUSt2

# Setup ------------------------------------------------

# Load library
library(seqinr) # Pre-step
library(ggplot2)
library(dplyr)
library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
source("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/src/ANCOM_scripts/ancom_v2.1.R")

# Pre-PICRUSt2 -----------------------------------------

# Load in counts table
final_counts_table <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/decontam_counts_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Drop tax column
final_counts <- as.matrix(final_counts_table[,c(grep("RV_no", colnames(final_counts_table)), grep("RV_ox", colnames(final_counts_table)))])

# Load in sv_seqs file
sv_seqs <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/sv_seqs.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE, row.names = NULL, comment.char = "")

# Only keep ASVs that are in the final counts table
ASVs <- sv_seqs[sv_seqs$V1 %in% rownames(final_counts_table), ]

# Prep. PICRUSt2 input data
write.fasta(as.list(ASVs$V2), ASVs$V1,
            "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/picrust2_svs.fasta",
            as.string = TRUE)

write.table(final_counts,
            file = "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/picrust2_counts.txt", 
            sep = "\t", col.names = NA, quote = F)








# ANCOM on PICRUSt 2 ----------------------------------

# Load in the metacyc athways
picrust_pathways <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE, row.names = 1, comment.char = "")

# Load in meta data
metadata <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/meta_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE, row.names = 1, comment.char = "")

# Only keep samples from metadata table
picrust_metadata <-  as.data.frame(metadata[c(grep("RV_no", rownames(metadata)), grep("RV_ox", rownames(metadata))), ])
# Make SampleID column
picrust_metadata$SampleID <- rownames(picrust_metadata)


# Step 1: Data pre-processing
prepro <- feature_table_pre_process(feature_table = picrust_pathways, 
                                    meta_data = picrust_metadata, 
                                    sample_var = "SampleID", 
                                    group_var = NULL, 
                                    out_cut = 0.05, 
                                    zero_cut = 0.90, 
                                    lib_cut = 0, 
                                    neg_lb = FALSE)


# Run ANCOM
t_start <- Sys.time()
res <- ANCOM(feature_table = prepro$feature_table, 
             meta_data = prepro$meta_data, 
             struc_zero = prepro$structure_zeros, 
             main_var = "Group", 
             p_adj_method = "BH", 
             alpha = 0.05, 
             adj_formula = NULL, 
             rand_formula = "~1 | Time")
t_end = Sys.time()
t_run = t_end - t_start
t_run
# Output:
# Time difference of 29.78661 secs

res$out[res[["out"]][["detected_0.6"]], "taxa_id"]




res <- ANCOM(feature_table = prepro$feature_table, 
             meta_data = prepro$meta_data, 
             struc_zero = prepro$structure_zeros, 
             main_var = "Group", 
             p_adj_method = "BH", 
             alpha = 0.05, 
             adj_formula = NULL)



# Extract EC numbers of interest ----------------------------------

# Load in EC numbers
EC_numbers <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/picrust2_out_pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, comment.char="")

# Calculate proportions
EC_props <- apply(EC_numbers, 2, function(x) {x/sum(x)})

# Pull out genes associated with oxalate degradation
ox_genes <- as.data.frame(EC_props[rownames(EC_props) == "EC:1.2.7.10" | rownames(EC_props) == "EC:4.1.1.8" | rownames(EC_props) == "EC:2.8.3.16" | rownames(EC_props) == "EC:4.1.1.2", ])

ox_genes$EC_number <- rownames(ox_genes)
rownames(ox_genes) <- NULL
  
long_ox_genes <- ox_genes %>%
  pivot_longer(!EC_number, names_to = "Sample", values_to = "Rel_Abun")  

long_ox_genes$Time <- ifelse(grepl("t0", long_ox_genes$Sample), "0h", ifelse(grepl("t24", long_ox_genes$Sample), "24h", "48h"))
long_ox_genes$Group <- ifelse(grepl("ox", long_ox_genes$Sample), "Oxalate", "No Oxalate")





KO_numbers <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE, row.names = 1, comment.char = "")

# Calculate proportions
KO_props <- apply(KO_numbers, 2, function(x) {x/sum(x)})


ox_genes <- as.data.frame(KO_props[rownames(KO_props) == "K01569" | rownames(KO_props) == "K08177" | rownames(KO_props) == "K18702" | rownames(KO_props) == "K19070" | rownames(KO_props) == "K19071" | rownames(KO_props) == "K19072" | rownames(KO_props) == "K07749" | rownames(KO_props) == "K01577", ])

sum <- as.data.frame(colSums(ox_genes))
sum$Sample <- row.names(sum)
sum$Time <- ifelse(grepl("t0", sum$Sample), "0h", ifelse(grepl("t24", sum$Sample), "24h", "48h"))
sum$Group <- ifelse(grepl("ox", sum$Sample), "Oxalate", "No Oxalate")
sum$Grouping <- ifelse(grepl("no-t0", sum$Sample), "no_0h", ifelse(grepl("no-t24", sum$Sample), "no_24h", ifelse(grepl("no-t48", sum$Sample), "no_48h", ifelse(grepl("ox-t0", sum$Sample), "ox_0h", ifelse(grepl("ox-t24", sum$Sample), "ox_24h","ox_48h")))))



ggplot(sum, aes(x = Time, y = colSums(ox_genes), fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(width = 0.75)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0), alpha = 0.8) +
  scale_fill_manual(values = brewer.pal(3,"Set3"), name = "") +
  labs(title = NULL, x = NULL, y = "colSums(ox_genes)") + 
  scale_y_continuous(trans='log10') +
  theme_classic() +
  theme(legend.position = "right")



#K01569	oxdD; oxalate decarboxylase [EC:4.1.1.2]
#K08177	oxlT; MFS transporter, OFA family, oxalate/formate antiporter
#K18702	uctC; CoA:oxalate CoA-transferase [EC:2.8.3.19]
#K19070	oorA; oxalate oxidoreductase subunit alpha [EC:1.2.7.10]
#K19071	oorB; oxalate oxidoreductase subunit beta [EC:1.2.7.10]
#K19072	oorD; oxalate oxidoreductase subunit delta [EC:1.2.7.10]
#K07749	frc; formyl-CoA transferase [EC:2.8.3.16]
#K01577	oxc; oxalyl-CoA decarboxylase [EC:4.1.1.8]



########## Percent contribution ############

# Load in KO percent contributions
KO_contrib <- read.table("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", check.names = FALSE, row.names = NULL, comment.char = "")

KO_contrib_ox <- as.data.frame(KO_contrib[KO_contrib$KO_number == "K01569" | KO_contrib$KO_number == "K08177" | KO_contrib$KO_number == "K18702" | KO_contrib$KO_number == "K19070" | KO_contrib$KO_number == "K19071" | KO_contrib$KO_number == "K19072" | KO_contrib$KO_number == "K07749" | KO_contrib$KO_number == "K01577", ])