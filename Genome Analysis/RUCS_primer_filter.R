### This script filters the RUCS primer output


# Load library
library(stringr)

# Function for primer filtering
primer_filter <- function(x) {
  # Filter based on GC (must be >= 45C and <= 60C in both)
  primers1 <- x[x$`forward_gc%` >= 45 & x$`reverse_gc%` >= 45 & 
                        x$`forward_gc%` <= 55 & x$`reverse_gc%` <= 55, ]
  
  # Remove primers with AA or TT at 3' end
  primers2 <- primers1[substr(primers1$forward_primer, 20, 21) != "AA" & 
                         substr(primers1$forward_primer, 20, 21) != "TT" & 
                         substr(primers1$reverse_primer, 20, 21) != "AA" & 
                         substr(primers1$reverse_primer, 20, 21) != "TT", ]
  
  # Remove primers with A or T at 3' end
  primers3 <- primers2[substr(primers2$forward_primer, 21, 21) != "A" & 
                         substr(primers2$forward_primer, 21, 21) != "T" & 
                         substr(primers2$reverse_primer, 21, 21) != "A" & 
                         substr(primers2$reverse_primer, 21, 21) != "T", ]
  
  # Make a column with the last 5 nt's
  primers3$fwd_sub <- substr(primers3$forward_primer, 17, 21)
  primers3$rev_sub <- substr(primers3$reverse_primer, 17, 21)
  
  # Keep primers with 3/5 GC at 3' end nt's
  primers4 <- primers3[str_count(primers3$fwd_sub, "[GC]") == 3 & str_count(primers3$rev_sub, "[GC]") == 3, ]
  
  # Remove primers with >3 consecutive G/C residues at the 5'
  primers5 <- primers4[substr(primers4$forward_primer, 19, 21) != "CCC" & 
                         substr(primers4$forward_primer, 19, 21) != "GCC" & 
                         substr(primers4$forward_primer, 19, 21) != "CGC" & 
                         substr(primers4$forward_primer, 19, 21) != "GGC" &
                         substr(primers4$forward_primer, 19, 21) != "CCG" & 
                         substr(primers4$forward_primer, 19, 21) != "GCG" & 
                         substr(primers4$forward_primer, 19, 21) != "CGG" &
                         substr(primers4$forward_primer, 19, 21) != "GGG" & 
                         substr(primers4$reverse_primer, 19, 21) != "CCC" & 
                         substr(primers4$reverse_primer, 19, 21) != "GCC" & 
                         substr(primers4$reverse_primer, 19, 21) != "CGC" & 
                         substr(primers4$reverse_primer, 19, 21) != "GGC" &
                         substr(primers4$reverse_primer, 19, 21) != "CCG" & 
                         substr(primers4$reverse_primer, 19, 21) != "GCG" & 
                         substr(primers4$reverse_primer, 19, 21) != "CGG" &
                         substr(primers4$reverse_primer, 19, 21) != "GGG"
                       , ]
  # Remove primers with >3 consecutive G/C residues at the 3'
  primers6 <- primers5[substr(primers5$forward_primer, 0, 4) != "CCCC" & 
                         substr(primers5$forward_primer, 0, 4) != "GCCC" & 
                         substr(primers5$forward_primer, 0, 4) != "CGCC" & 
                         substr(primers5$forward_primer, 0, 4) != "GGCC" &
                         substr(primers5$forward_primer, 0, 4) != "CCGC" & 
                         substr(primers5$forward_primer, 0, 4) != "GCGC" & 
                         substr(primers5$forward_primer, 0, 4) != "CGGC" &
                         substr(primers5$forward_primer, 0, 4) != "GGGC" & 
                         substr(primers5$forward_primer, 0, 4) != "CCCG" & 
                         substr(primers5$forward_primer, 0, 4) != "GCCG" &
                         substr(primers5$forward_primer, 0, 4) != "CGCG" & 
                         substr(primers5$forward_primer, 0, 4) != "GGCG" & 
                         substr(primers5$forward_primer, 0, 4) != "CCGG" &
                         substr(primers5$forward_primer, 0, 4) != "GCGG" & 
                         substr(primers5$forward_primer, 0, 4) != "CGGG" & 
                         substr(primers5$forward_primer, 0, 4) != "GGGG" &
                         substr(primers5$reverse_primer, 0, 4) != "CCCC" & 
                         substr(primers5$reverse_primer, 0, 4) != "GCCC" & 
                         substr(primers5$reverse_primer, 0, 4) != "CGCC" & 
                         substr(primers5$reverse_primer, 0, 4) != "GGCC" &
                         substr(primers5$reverse_primer, 0, 4) != "CCGC" & 
                         substr(primers5$reverse_primer, 0, 4) != "GCGC" & 
                         substr(primers5$reverse_primer, 0, 4) != "CGGC" &
                         substr(primers5$reverse_primer, 0, 4) != "GGGC" & 
                         substr(primers5$reverse_primer, 0, 4) != "CCCG" & 
                         substr(primers5$reverse_primer, 0, 4) != "GCCG" &
                         substr(primers5$reverse_primer, 0, 4) != "CGCG" & 
                         substr(primers5$reverse_primer, 0, 4) != "GGCG" & 
                         substr(primers5$reverse_primer, 0, 4) != "CCGG" &
                         substr(primers5$reverse_primer, 0, 4) != "GCGG" & 
                         substr(primers5$reverse_primer, 0, 4) != "CGGG" & 
                         substr(primers5$reverse_primer, 0, 4) != "GGGG" 
                       , ]
  
  # Keep forward primers that terminate (5') with either C or G
  primers7 <- primers6[substr(primers6$forward_primer, 20, 21) == "GG" | 
                         substr(primers6$forward_primer, 20, 21) == "GC" | 
                         substr(primers6$forward_primer, 20, 21) == "CG" | 
                         substr(primers6$forward_primer, 20, 21) == "CC", ]
  
  # Keep reverse primers that terminate (5') with either C or G
  filtered_primers <<- primers7[substr(primers7$reverse_primer, 20, 21) == "GG" | 
                          substr(primers7$reverse_primer, 20, 21) == "GC" | 
                          substr(primers7$reverse_primer, 20, 21) == "CG" | 
                          substr(primers7$reverse_primer, 20, 21) == "CC", ]
  
}


# Load input list of primers
primers <- read.table("~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/qPCR_primers/Trial_1_primers/group1/results.tsv", header = TRUE, check.names = FALSE, sep = "\t", row.names = NULL)

# Run primer filter function
primer_filter(primers)

# Save output
write.table(filtered_primers, file="~/SCIENCE/Project_oxalobacter/genome_analysis/!key_files/qPCR_primers/group4/filtered_primers.tsv", sep="\t", col.names=NA, quote=F)
