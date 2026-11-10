# ==============================================================================
#
# R SCRIPT FOR COHORT & LOCATION GENEPOP (.gen) EXPORT (NeEstimator)
#
# Description:
# Integrates biological data to calculate Birth Year (Cohort) and groups by 
# Region. Subsets the VCF to only include individuals from a given cohort/region
# combination that meets a minimum sample size of 30. Recodes genotypes to 
# Genepop 6-digit format, and writes a .gen file for each valid group.
#
# Genepop Mapping:
#   VCF 0/0  -> Genepop 001001 (Homozygous Major)
#   VCF 0/1  -> Genepop 001002 (Heterozygote)
#   VCF 1/1  -> Genepop 002002 (Homozygous Minor)
#   VCF ./.  -> Genepop 000000 (Missing Data)
#
# ==============================================================================

# --- Required libraries ---
library(vcfR)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(readxl)
library(beepr)

# --- Set Working Directory ---
working_directory <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/pedigree reconstruction"
setwd(working_directory)

# --- Input & Output Paths ---
vcf_file_path <- "../post_RAD analysis/Output/populations.filtered.final.vcf.gz"
bio_data_path <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/plate_barcodes/Input/GobyRad_tabularPlateWell with Bio Data.xlsx"
output_dir <- "Output/genepop_cohorts_regions"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# SECTION 1: ROBUST GENEPOP WRITER FUNCTION
# ==============================================================================
write_ne_genepop_fast <- function(mat, filename, pop_labels, description, loci_names, sample_names) {
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
      line_data <- paste0(sample_names[i], " , ", paste(gts, collapse = " "))
      writeLines(line_data, con)
    }
  }
  
  # 5. Close the connection and flush the buffer
  close(con)
  message(paste("Successfully wrote:", filename))
}

# ==============================================================================
# SECTION 2: SETUP, BIO DATA & INFER GROUPS
# ==============================================================================
message("\n--- Loading Biological Data & Inferring Cohorts/Regions ---")
bio_data <- read_xlsx(path = bio_data_path, sheet = 1)

# Clean column names
names(bio_data) <- trimws(names(bio_data))
names(bio_data) <- gsub(" ", "_", names(bio_data))
names(bio_data) <- gsub("\\(|\\)|#|`|'", "", names(bio_data)) 

# Filter out NTCs and dummy samples
bio_data <- bio_data[!grepl(pattern = "NTC", x = bio_data$ID), ]
bio_data <- bio_data[!grepl(pattern = "999", x = bio_data$ID), ]

# Clean up Region spacing issues (e.g. "Albany ")
bio_data$Region <- trimws(as.character(bio_data$Region))

# Parse Sample Year & Age
bio_data$Sample_Year <- as.numeric(substr(x = bio_data$ID, start = 1, stop = 4))
bio_data$Age_num <- suppressWarnings(as.numeric(bio_data$Age))

# Infer Birth Year
bio_data$Birth_Year <- bio_data$Sample_Year - bio_data$Age_num
bio_data_clean <- bio_data[!is.na(bio_data$Birth_Year) & !is.na(bio_data$Fish_ID) & !is.na(bio_data$Region) & bio_data$Region != "", ]

# Extract unique combinations of Birth Year and Region
group_combos <- unique(bio_data_clean[, c("Birth_Year", "Region")])
message(paste("Successfully inferred", nrow(group_combos), "distinct Cohort + Region combinations."))


# ==============================================================================
# SECTION 3: LOAD VCF & GENERATE GENEPOPS
# ==============================================================================
message("\n--- Loading Global VCF file ---")
vcf <- read.vcfR(vcf_file_path, verbose = FALSE)
vcf_samples <- colnames(vcf@gt)[-1]

message("\n--- Generating Genepop files (Min N = 30) ---")

for (i in 1:nrow(group_combos)) {
  
  cohort <- group_combos$Birth_Year[i]
  region <- group_combos$Region[i]
  
  # Format a clean string for filenames and pop labels (removes spaces/special chars)
  clean_region <- gsub("[^A-Za-z0-9]", "", region)
  group_name <- paste(cohort, clean_region, sep = "_")
  
  # 1. Identify individuals in this specific cohort AND region
  group_fish <- as.character(bio_data_clean$Fish_ID[
    bio_data_clean$Birth_Year == cohort & 
      bio_data_clean$Region == region
  ])
  
  # Match with VCF
  valid_samples <- intersect(vcf_samples, group_fish)
  
  message(sprintf("\nProcessing %s: Retained %d matched samples.", group_name, length(valid_samples)))
  
  # APPLY MINIMUM SAMPLE SIZE FILTER
  if (length(valid_samples) < 30) {
    message("  -> Skipping: Sample size is below the minimum threshold of 30.")
    next
  }
  
  # 2. Subset VCF
  vcf_cohort <- vcf[, c("FORMAT", valid_samples)]
  
  # 3. Clean Loci (Drop 100% missing)
  gt_matrix <- extract.gt(vcf_cohort, element = "GT")
  missing_rates <- rowMeans(is.na(gt_matrix) | gt_matrix == "./.")
  valid_loci <- which(missing_rates < 1.0) 
  vcf_cohort <- vcf_cohort[valid_loci, ]
  
  message(sprintf("  -> Retaining %d valid loci.", nrow(vcf_cohort)))
  
  # 4. Extract and Recode GT Matrix
  message("  -> Recoding genotypes to 6-digit format...")
  gt_recode <- extract.gt(vcf_cohort, element = "GT")
  
  gt_recode[gt_recode == "0/0"] <- "001001"
  gt_recode[gt_recode == "0/1" | gt_recode == "1/0"] <- "001002"
  gt_recode[gt_recode == "1/1"] <- "002002"
  gt_recode[is.na(gt_recode) | gt_recode == "./."] <- "000000"
  
  # 5. Metadata for Genepop Writer
  loci_names <- rownames(gt_recode)
  sample_names <- colnames(gt_recode)
  
  # Using the clean Cohort_Region name as the single 'Pop' group for the entire file
  pop_labels <- rep(group_name, length(sample_names))
  
  # 6. Execute Fast Write
  out_path <- file.path(output_dir, paste0("round_goby_", group_name, ".gen"))
  desc <- paste("Round_Goby", group_name)
  
  write_ne_genepop_fast(
    mat = gt_recode, 
    filename = out_path, 
    pop_labels = pop_labels, 
    description = desc,
    loci_names = loci_names,
    sample_names = sample_names
  )
}

message("\n--- Genepop Pipeline Complete ---")
beep(sound = 8)