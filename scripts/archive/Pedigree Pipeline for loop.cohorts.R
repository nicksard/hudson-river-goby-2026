# ==============================================================================
#
# R SCRIPT FOR COHORT-SPECIFIC COLONY INPUT FILE GENERATION
#
# Author: Gemini & Nick Sard
#
# Description:
# This script integrates biological data to calculate Birth Year (Cohort),
# cross-references against the genotyped VCF, diagnoses missing samples,
# subsets the VCF to only include individuals from a given cohort, and 
# generates a single COLONY .dat file per cohort to estimate Nb.
#
# SETTINGS APPLIED:
#   - Markers: Highest number of valid markers possible per cohort
#   - Run Length: Medium (2)
#   - Likelihood Precision: Medium (1)
#
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1A: CORE CONFIGURATION
# ------------------------------------------------------------------------------

# --- Set Working Directory ---
working_directory <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/pedigree reconstruction"
setwd(working_directory)

# --- Required libraries ---
library(vcfR)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(readxl)

# Stop automatic scientific notations in COLONY dat file
options(scipen=999)

# --- Input & Output Paths ---
vcf_file_path <- "../post_RAD analysis/Output/populations.filtered.final.vcf.gz"
output_dir <- "Output/colony_runs_cohorts"
bio_data_path <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/plate_barcodes/Input/GobyRad_tabularPlateWell with Bio Data.xlsx"

# --- Source the FLEX Function ---
flex_function_path <- "colony.dat.create_FLEX_v2.R"
if (!file.exists(flex_function_path)) {
  stop(paste("FATAL ERROR: The function file", flex_function_path, "was not found."))
}
source(flex_function_path)

# ------------------------------------------------------------------------------
# SECTION 1B: COLONY RUN PARAMETERS (FIXED)
# ------------------------------------------------------------------------------
COLONY_PARAMETERS <- list(
  project_name = "rgSibship", 
  analysis_type = "sibship",   
  
  # Core Run Parameters
  spp.type = 2,             # 2 = Dioecious
  inbreeding = 0,           # 0 = No
  ploidy = 0,               # 0 = Diploid
  fem.gamy = 0,             # 0 = Polygamy
  mal.gamy = 0,             # 0 = Polygamy
  
  # Likelihood & Alleles (FIXED BASED ON REQUEST)
  full.likelihood = 1,      # 1 = Full Likelihood
  likelihood.precision = 1, # 1 = Medium Precision
  run.length = 2,           # 2 = Medium Run Length
  update.alfs = 0,          # 0 = No (Standard for subset runs)
  known.alfs = 0,
  
  # Error Rates
  marker_type = 0,          
  error_rate_geno = 0.02,   
  error_rate_other = 0.001, 
  
  # Other params
  clone = 0,
  sib.scale = 1,
  sib.prior = 0,
  run.number = 1,
  monitor = 0,
  windows.version = 0,
  prob.mom = 0.0,
  prob.dad = 0.0
)

# ==============================================================================
# SECTION 2: SETUP, BIO DATA & INFER BIRTH YEARS
# ==============================================================================
message("\n--- Starting COHORT COLONY File Generation ---")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("\n--- Loading Biological Data & Inferring Cohorts ---")

# Read Excel
bio_data <- read_xlsx(path = bio_data_path, sheet = 1)

# Clean column names programmatically (removing spaces and special characters)
names(bio_data) <- trimws(names(bio_data))
names(bio_data) <- gsub(" ", "_", names(bio_data))
names(bio_data) <- gsub("\\(|\\)|#|`|'", "", names(bio_data)) 

# Filter out NTCs and dummy samples
bio_data <- bio_data[!grepl(pattern = "NTC", x = bio_data$ID), ]
bio_data <- bio_data[!grepl(pattern = "999", x = bio_data$ID), ]

# Parse Sample Year directly from the first 4 digits of the ID string
bio_data$Sample_Year <- as.numeric(substr(x = bio_data$ID, start = 1, stop = 4))

# Force Age to numeric (safely converting "NA" text to true NAs)
bio_data$Age_num <- suppressWarnings(as.numeric(bio_data$Age))

# Infer Birth Year
bio_data$Birth_Year <- bio_data$Sample_Year - bio_data$Age_num

# Filter out individuals lacking age data or IDs
bio_data_clean <- bio_data[!is.na(bio_data$Birth_Year) & !is.na(bio_data$Fish_ID), ]

cohorts <- sort(unique(bio_data_clean$Birth_Year))
message(paste("Successfully inferred", length(cohorts), "distinct birth cohorts:", paste(cohorts, collapse=", ")))


# ==============================================================================
# SECTION 3: LOAD VCF & DIAGNOSE MISSING SAMPLES
# ==============================================================================
message("\n--- Loading Global VCF file ---")
vcf <- read.vcfR(vcf_file_path, verbose = FALSE)
vcf_samples <- colnames(vcf@gt)[-1]
message(paste("Global VCF contains:", nrow(vcf), "loci and", length(vcf_samples), "samples."))

# --- DIAGNOSTICS: Identify Dropped Fish ---
message("\n--- Running Missing Sample Diagnostics ---")
bio_sample_ids <- as.character(bio_data_clean$Fish_ID) 
missing_from_vcf <- setdiff(bio_sample_ids, vcf_samples)

if(length(missing_from_vcf) > 0) {
  message(paste("IDENTIFIED:", length(missing_from_vcf), "biologically aged samples are MISSING from the final VCF."))
  dropped_fish_metadata <- bio_data_clean[bio_data_clean$Fish_ID %in% missing_from_vcf, ]
  
  message("-> Missing Fish by Plate:")
  print(table(dropped_fish_metadata$Plate))
  
  message("-> Missing Fish by Sample Year:")
  print(table(dropped_fish_metadata$Sample_Year))
} else {
  message("All biological samples were successfully found in the VCF.")
}
write.table(x = dropped_fish_metadata,file = "Output/dropped.fish.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

# ==============================================================================
# SECTION 4: HELPER FUNCTION (The "Translator")
# ==============================================================================
convert_vcf_to_colony_format <- function(vcf_obj, marker_type, error_rate_1, error_rate_2) {
  n_loci <- nrow(vcf_obj)
  n_inds <- ncol(vcf_obj@gt) - 1
  
  original_locus_names <- vcf_obj@fix[, "ID"]
  if (all(is.na(original_locus_names)) || all(original_locus_names == ".")) {
    original_locus_names <- paste0(vcf_obj@fix[, "CHROM"], ":", vcf_obj@fix[, "POS"])
  }
  
  generic_locus_names <- paste0("LOCUS", 1:n_loci)
  locus_map_df <- data.frame(generic_name = generic_locus_names, original_name = original_locus_names, stringsAsFactors = FALSE)
  
  markers_df <- data.frame(matrix(data = 0, nrow = 3, ncol = n_loci))
  colnames(markers_df) <- generic_locus_names
  markers_df[1, ] <- marker_type
  markers_df[2, ] <- error_rate_1
  markers_df[3, ] <- error_rate_2
  
  gt_alleles <- extract.gt(vcf_obj, element = "GT", return.alleles = TRUE)
  gt_t <- as.data.frame(t(gt_alleles), stringsAsFactors = FALSE)
  gt_t <- tibble::rownames_to_column(gt_t, var = "Sample_ID")
  colnames(gt_t) <- c("Sample_ID", generic_locus_names)
  
  genotypes_df <- gt_t[, "Sample_ID", drop = FALSE]
  
  for (locus in generic_locus_names) {
    new_col_A <- paste0(locus, "_A")
    new_col_B <- paste0(locus, "_B")
    raw_gts <- gt_t[[locus]]
    
    hom_idx <- !grepl("/", raw_gts, fixed = TRUE) & !grepl(".", raw_gts, fixed = TRUE) & raw_gts != "NA"
    hom_idx[is.na(hom_idx)] <- FALSE
    if(any(hom_idx)) raw_gts[hom_idx] <- paste(raw_gts[hom_idx], raw_gts[hom_idx], sep = "/")
    
    raw_gts <- gsub(".", "0", raw_gts, fixed = TRUE)
    
    split_cols <- tidyr::separate(data.frame(gts = raw_gts), col = "gts", into = c(new_col_A, new_col_B), sep = "/", remove = TRUE)
    split_cols[is.na(split_cols)] <- "0"
    
    genotypes_df[[new_col_A]] <- split_cols[[new_col_A]]
    genotypes_df[[new_col_B]] <- split_cols[[new_col_B]]
  }
  
  genotypes_df[is.na(genotypes_df)] <- "0"
  all_alleles <- unique(unlist(genotypes_df[, -1]))
  alleles_to_map <- all_alleles[all_alleles != "0"]
  allele_map <- setNames(as.character(1:length(alleles_to_map)), alleles_to_map)
  allele_map["0"] <- "0"
  
  genotypes_df[, -1] <- lapply(genotypes_df[, -1], function(col) { unname(allele_map[as.character(col)]) })
  
  return(list(genotypes_df = genotypes_df, markers_df = markers_df, locus_map_df = locus_map_df))
}


# ==============================================================================
# SECTION 5: GENERATE COHORT-SPECIFIC FILES
# ==============================================================================
message("\n--- Generating COLONY .dat files per Cohort ---")

generated_dat_files <- c()

for (cohort in cohorts) {
  
  # 1. Identify individuals in this cohort
  cohort_fish <- as.character(bio_data_clean$Fish_ID[bio_data_clean$Birth_Year == cohort])
  
  # 2. Match with VCF Sample IDs
  valid_samples <- intersect(vcf_samples, cohort_fish)
  
  message(sprintf("\nProcessing Cohort %s: Found %d biological records matching %d samples in the VCF.", 
                  cohort, length(cohort_fish), length(valid_samples)))
  
  if (length(valid_samples) < 2) {
    message("  -> Skipping: Not enough individuals to run sibship analysis.")
    next
  }
  
  # 3. Subset VCF to ONLY this cohort (FORMAT is col 1, then the valid samples)
  vcf_cohort <- vcf[, c("FORMAT", valid_samples)]
  
  # 4. Clean Loci for this specific cohort (Drop loci that are 100% missing in these specific fish)
  gt_matrix <- extract.gt(vcf_cohort, element = "GT")
  missing_rates <- rowMeans(is.na(gt_matrix) | gt_matrix == "./.")
  valid_loci <- which(missing_rates < 1.0) 
  vcf_cohort <- vcf_cohort[valid_loci, ]
  
  message(sprintf("  -> Retaining %d valid loci for Cohort %s", nrow(vcf_cohort), cohort))
  
  # 5. Convert to COLONY Format
  colony_data <- convert_vcf_to_colony_format(
    vcf_obj = vcf_cohort,
    marker_type = COLONY_PARAMETERS$marker_type,
    error_rate_1 = COLONY_PARAMETERS$error_rate_geno,
    error_rate_2 = COLONY_PARAMETERS$error_rate_other
  )
  
  # 6. Write the .dat file
  project_name_full <- paste0(COLONY_PARAMETERS$project_name, "_Cohort_", cohort)
  dat_filename <- file.path(output_dir, paste0(project_name_full, ".dat"))
  
  colony.dat.create_FLEX.R(
    kids = colony_data$genotypes_df,
    markers = colony_data$markers_df,
    moms = NA, 
    dads = NA, 
    output_filename = dat_filename,
    project_name = project_name_full,
    analysis_type = COLONY_PARAMETERS$analysis_type,
    
    update.alfs = COLONY_PARAMETERS$update.alfs,
    run.length = COLONY_PARAMETERS$run.length,
    likelihood.precision = COLONY_PARAMETERS$likelihood.precision,
    spp.type = COLONY_PARAMETERS$spp.type,
    inbreeding = COLONY_PARAMETERS$inbreeding,
    ploidy = COLONY_PARAMETERS$ploidy,
    fem.gamy = COLONY_PARAMETERS$fem.gamy,
    mal.gamy = COLONY_PARAMETERS$mal.gamy,
    clone = COLONY_PARAMETERS$clone,
    sib.scale = COLONY_PARAMETERS$sib.scale,
    sib.prior = COLONY_PARAMETERS$sib.prior,
    known.alfs = COLONY_PARAMETERS$known.alfs,
    run.number = COLONY_PARAMETERS$run.number,
    monitor = COLONY_PARAMETERS$monitor,
    windows.version = COLONY_PARAMETERS$windows.version,
    full.likelihood = COLONY_PARAMETERS$full.likelihood,
    prob.mom = COLONY_PARAMETERS$prob.mom,
    prob.dad = COLONY_PARAMETERS$prob.dad
  )
  
  message(paste("  -> Saved:", basename(dat_filename)))
  generated_dat_files <- c(generated_dat_files, basename(dat_filename))
  
  # Write mapping file
  map_filename <- file.path(output_dir, paste0(project_name_full, "_LOCUS_MAP.csv"))
  write.csv(colony_data$locus_map_df, file = map_filename, row.names = FALSE, quote = TRUE)
}

# ==============================================================================
# SECTION 6: GENERATE FILE LIST FOR HPC
# ==============================================================================
job_list_filename <- file.path(output_dir, "sbatch_cohort_job_list.txt")
writeLines(generated_dat_files, job_list_filename)
message("\n--- Pipeline Complete ---")
message(paste("Master Job List written to:", job_list_filename))