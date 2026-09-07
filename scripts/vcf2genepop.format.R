################################################################################
### VCF to Genepop (Fixed Mapping for NeEstimator)
###
### Description: Maps VCF genotypes to Genepop 6-digit format:
###   VCF 0/0  -> Genepop 001001 (Homozygous Major)
###   VCF 0/1  -> Genepop 001002 (Heterozygote)
###   VCF 1/1  -> Genepop 002002 (Homozygous Minor)
###   VCF ./.  -> Genepop 000000 (Missing Data)
################################################################################

library(vcfR)
library(tidyverse)
library(beepr)

# Set working directory
setwd("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis")

# 1. Load the Filtered VCF
message("--- Loading Filtered VCF ---")
vcf <- read.vcfR("Output/populations.filtered.final.vcf.gz")

# ==============================================================================
# Robust Genepop Export (File Connection + writeLines)
# ==============================================================================

write_ne_genepop_fast <- function(mat, filename, pop_labels, description) {
  # 1. Open a file connection for writing ("w")
  con <- file(filename, "w")
  
  # 2. Write the description (Line 1)
  writeLines(description, con)
  
  # 3. Write the locus names (Line 2) - separated by commas
  writeLines(paste(loci_names, collapse = ", "), con)
  
  # 4. Get unique populations
  unique_pops <- unique(pop_labels)
  
  for(p in unique_pops) {
    # Write the Pop delimiter
    writeLines("Pop", con)
    
    # Identify samples in this pop
    pop_indices <- which(pop_labels == p)
    
    for(i in pop_indices) {
      # Extract genotypes for this individual
      gts <- mat[, i]
      
      # Construct the line: ID , GTY1 GTY2 ...
      # writeLines automatically handles the newline (\n) at the end
      line_data <- paste0(sample_names[i], " , ", paste(gts, collapse = " "))
      writeLines(line_data, con)
    }
  }
  
  # 5. Close the connection and flush the buffer
  close(con)
  message(paste("Successfully wrote:", filename))
}

# ==============================================================================
# Run the recoding and export
# ==============================================================================

# Ensure VCF is loaded and GT matrix is extracted (from previous steps)
gt <- extract.gt(vcf, element = "GT")
gt_recode <- gt

message("Recoding genotypes to 6-digit format...")
gt_recode[gt_recode == "0/0"] <- "001001"
gt_recode[gt_recode == "0/1" | gt_recode == "1/0"] <- "001002"
gt_recode[gt_recode == "1/1"] <- "002002"
gt_recode[is.na(gt_recode) | gt_recode == "./."] <- "000000"

# Metadata
loci_names <- rownames(gt_recode)
sample_names <- colnames(gt_recode)

# 3. Use an absolute path for the export to be safe
out_path_year <- file.path(getwd(), "Output/round_goby_by_year.gen")
out_path_global <- file.path(getwd(), "Output/round_goby_global.gen")

# 4. Run the export
# Final Exports
write_ne_genepop_fast(mat = gt_recode, filename = "Output/round_goby_by_year.gen", pop_labels = substr(sample_names, 1, 4),description =  "Round_Goby_By_Year")

write_ne_genepop_fast(gt_recode, "Output/round_goby_global.gen", 
                      rep("Global", length(sample_names)), "Round_Goby_Global")

beep(sound = 8)
