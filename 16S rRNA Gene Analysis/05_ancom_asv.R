# This script run ANCOM for ASV abundance

# Load libraries
library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
source("~/SCIENCE/Project_oxalobacter/SRK-16S_new/src/ANCOM_scripts/ancom_v2.1.R")


# Load in the data table
final_counts_table <- read.table("~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/decontam_counts_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = "", stringsAsFactors = FALSE)

# Only keep samples from tax table
final_counts <- as.matrix(final_counts_table[, c(grep("RV_no", colnames(final_counts_table)), grep("RV_ox", colnames(final_counts_table)))])


# Load in the metadata table
metadata <- read.table("~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/meta_data.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", check.names = F, row.names = 1, comment.char = "")

# Only keep samples from metadata table
final_metadata <-  as.data.frame(metadata[c(grep("RV_no", rownames(metadata)), grep("RV_ox", rownames(metadata))), ])
# Make SampleID column
final_metadata$SampleID <- rownames(final_metadata)


# Step 1: Data pre-processing
prepro <- feature_table_pre_process(feature_table = final_counts, 
                                    meta_data = final_metadata, 
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
# Time difference of 1.928861 mins


res$out[res[["out"]][["detected_0.7"]], "taxa_id"]

# Output:
#  [1] "SV_107" "SV_131" "SV_17"  "SV_24"  "SV_27"  "SV_36"  "SV_40"  "SV_45"  "SV_71"  "SV_8"   "SV_98" 

# Write out which ASVs are significant (TRUE) or not (FALSE)
write.table(res$out, file = "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/ANCOM_results.txt", sep = "\t", col.names = NA, quote = F)
# Write out the effect sizes of all ASVs
write.table(res[["fig"]][["data"]], file = "~/SCIENCE/Project_oxalobacter/SRK-16S_new/data/ANCOM_results_CLR.txt", sep = "\t", col.names = NA, quote = F)
