# DADA2 pipeline for reads resulting from the custom demultiplexing script.
# The code is very minimally modified from the tutorial at 
# https://benjjneb.github.io/dada2/tutorial.html

# Setup ------------------------------------------------

# Load packages
library(dada2)
library(dplyr)

#Dump the R sessioninfo for later
writeLines(capture.output(sessionInfo()), "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/RsessionInfo_dada2.txt")

# Establish the file path to the directory containing the demultiplexed samples

reads <- "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_reads_ox"
list.files(reads) # Check that all the files are there 


# Forward and reverse fastq filenames have the format:  
# SAMPLENAME-R1.fastq and SAMPLENAME-R2.fastq
# Use this fact to sort them into two groups
fnFs <- sort(list.files(reads, pattern = "-R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(reads, pattern = "-R2.fastq", full.names = TRUE))

# Extract sample names (i.e., exclude the forward/reverse identifier), 
# assuming filenames have format SAMPLENAME-Rn.fastq and SAMPLENAME
# does not include any "-"
sample_names <- sapply(strsplit(basename(fnFs), "-R"), `[`, 1)
any(duplicated(sample_names))



# Check read quality -----------------------------------
# Grab four samples, the reads of which will be examined in terms of their quality profiles
ids <- round(runif(4,1,length(sample_names)))

pdf("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/quality_profiles.pdf")
plotQualityProfile(fnFs[ids])
plotQualityProfile(fnRs[ids])
dev.off()



# Filter reads based on QC -----------------------------

# Make filenames for the filtered fastq files
filtFs <- paste0(reads, "/", sample_names, "-F-filt.fastq.gz")
filtRs <- paste0(reads, "/", sample_names, "-R-filt.fastq.gz")

# Perform qualtiy filtering
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                   truncLen = c(175,150), #set based off of quality profiles
                   truncQ = 2,
                   maxN = 0,
                   maxEE = c(2,2), rm.phix = TRUE,
                   compress = TRUE, verbose = TRUE, multithread = FALSE)

write.table(out, file = "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/after_filter.txt", sep = "\t", col.names = NA, quote = F)


# Learn the error rates --------------------------------

errF <- learnErrors(filtFs, multithread = TRUE, randomize = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE, randomize = TRUE)

# Plot the error rates
pdf("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/error_profiles.pdf")
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()



# Dereplication ----------------------------------------

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample_names
names(derepRs) <- sample_names



# Sample inference, merge paired reads, remove chimeras ------------------------------------------------

# Use the filtered files and error rates to perform
# sample inference (without pooling)
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = TRUE)

# Merge the forward and reverse paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)


# Make the sequence table, samples are by rows
seqtab_raw <- makeSequenceTable(mergers)

# Summarize the output by sequence length
table(nchar(getSequences(seqtab_raw)))
# Output:
# 223  224  252  253  254  270 
#  1    1  169 1762   22    3  

# Most ASVs are acceptable lengths, but the others must be removed (keep ASVs with length of 252 - 254 bp)
seqtab_trimmed <- seqtab_raw[, nchar(colnames(seqtab_raw)) %in% seq(252, 254)]


# Filter chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab_trimmed, method = "consensus", verbose = TRUE, multithread = FALSE)
# Output :
# Identified 1512 bimeras out of 1953 input sequences.

sum(seqtab_nochim)/sum(seqtab_trimmed)
# Output:
# 0.9002185

# Write the table as traditional read count table with ASVs
#samples are rows
write.table(seqtab_nochim, file = "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/temp_dada2_nochim.txt", sep = "\t", col.names = NA, quote = F)



# Sanity check -----------------------------------------

# Check how many reads made it through the pipeline
# This is good to report in your methods/results
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab_trimmed), rowSums(seqtab_nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample_names
track <- cbind(track, percent_pass = mapply(function(x,y) x/y*100, track[ , "nonchim"], track[ ,"input"]))

write.table(track, file="~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/readsout.txt", sep="\t", col.names=NA, quote=F)


# Assign taxonomy --------------------------------------

tax_nochim <- assignTaxonomy(seqtab_nochim, "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

tax_nochim <- addSpecies(tax_nochim, "~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/silva_species_assignment_v138.1.fa.gz")


# Get the taxonomy string by merging columns 1 to 7 to get the full taxonomy (to species)
untax <- unname(tax_nochim)
tax_vector <- apply(untax, 1, function(x){paste(x[1:7], collapse=":")})

# Add taxonomy to the table
seqtab_nochim_tax <- rbind(seqtab_nochim, tax_vector)

# Transpose the table so samples are columns
t_seqtab_nochim_tax <-t (seqtab_nochim_tax)

# Remove the rownames (SV sequences) to a separate table and replace with arbitrary SV (sequence variants) numbers
sv_seqs <- rownames(t_seqtab_nochim_tax)
sv_num <- paste("SV", seq(from = 0, to = nrow(t_seqtab_nochim_tax)-1), sep="_")

# Replace reads with ASVs
rownames(t_seqtab_nochim_tax)<-sv_num

# Write output tables
write.table(t_seqtab_nochim_tax, file="~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/dada2_nochim_tax.txt", sep="\t", col.names=NA, quote=F)
write.table(sv_seqs, file="~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/sv_seqs.txt", sep="\t", row.names=sv_num, col.names=F,  quote=F)



save.image("~/SCIENCE/Project_Oxalate_venema/SRK-16S_new/data/demultiplex_dada2/dada2.RData")

