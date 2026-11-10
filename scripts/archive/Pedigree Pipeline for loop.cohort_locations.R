# ==============================================================================
# R SCRIPT FOR COHORT & LOCATION SPECIFIC COLONY FILE GENERATION (MAX 30)
# ==============================================================================

# --- Set Working Directory ---
working_directory <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/pedigree reconstruction"
setwd(working_directory)

library(vcfR)
library(dplyr)
library(readxl)
library(tidyr)
library(tibble)

# Stop automatic scientific notations in COLONY dat file
options(scipen=999)

# --- Input & Output Paths ---
vcf_file_path <- "../post_RAD analysis/Output/populations.filtered.final.vcf.gz"
output_dir <- "Output/colony_runs_cohort_location_max30"
bio_data_path <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/plate_barcodes/Input/GobyRad_tabularPlateWell with Bio Data.xlsx"

# --- Source the FLEX Function ---
flex_function_path <- "colony.dat.create_FLEX_v2.R"
if (!file.exists(flex_function_path)) {
  stop(paste("FATAL ERROR: The function file", flex_function_path, "was not found."))
}
source(flex_function_path)

# ------------------------------------------------------------------------------
# SECTION 1A: HELPER FUNCTION
# ------------------------------------------------------------------------------
convert_vcf_to_colony_format <- function(vcf_obj, marker_type, error_rate_1, error_rate_2) { 
  n_loci <- nrow(vcf_obj) 
  n_inds <- ncol(vcf_obj@gt) - 1 
  
  # --- 1. Get Locus Names & Create Map --- 
  original_locus_names <- vcf_obj@fix[, "ID"] 
  if (all(is.na(original_locus_names)) || all(original_locus_names == ".")) { 
    original_locus_names <- paste0(vcf_obj@fix[, "CHROM"], ":", vcf_obj@fix[, "POS"]) 
  } 
  generic_locus_names <- paste0("LOCUS", 1:n_loci) 
  
  locus_map_df <- data.frame( 
    generic_name = generic_locus_names, 
    original_name = original_locus_names, 
    stringsAsFactors = FALSE 
  ) 
  
  # --- 2. Create the 'markers_df' --- 
  markers_df <- data.frame(matrix(data = 0, nrow = 3, ncol = n_loci)) 
  colnames(markers_df) <- generic_locus_names 
  markers_df[1, ] <- marker_type 
  markers_df[2, ] <- error_rate_1 
  markers_df[3, ] <- error_rate_2 
  
  # --- 3. Create the 'genotypes_df' --- 
  gt_alleles <- extract.gt(vcf_obj, element = "GT", return.alleles = TRUE) 
  gt_t <- as.data.frame(t(gt_alleles), stringsAsFactors = FALSE) 
  gt_t <- tibble::rownames_to_column(gt_t, var = "Sample_ID") 
  colnames(gt_t) <- c("Sample_ID", generic_locus_names) 
  
  genotypes_df <- gt_t[, "Sample_ID", drop = FALSE] 
  
  for (locus in generic_locus_names) { 
    new_col_A <- paste0(locus, "_A") 
    new_col_B <- paste0(locus, "_B") 
    raw_gts <- gt_t[[locus]] 
    
    # Handle Homozygotes (e.g. "T" -> "T/T")
    hom_idx <- !grepl("/", raw_gts, fixed = TRUE) & !grepl(".", raw_gts, fixed = TRUE) & raw_gts != "NA" 
    hom_idx[is.na(hom_idx)] <- FALSE 
    if(any(hom_idx)) { raw_gts[hom_idx] <- paste(raw_gts[hom_idx], raw_gts[hom_idx], sep = "/") } 
    
    # Handle Missing Data
    raw_gts <- gsub(".", "0", raw_gts, fixed = TRUE) 
    
    split_cols <- tidyr::separate(data.frame(gts = raw_gts), col = "gts", into = c(new_col_A, new_col_B), sep = "/", remove = TRUE) 
    split_cols[is.na(split_cols)] <- "0" 
    
    genotypes_df[[new_col_A]] <- split_cols[[new_col_A]] 
    genotypes_df[[new_col_B]] <- split_cols[[new_col_B]] 
  } 
  
  genotypes_df[is.na(genotypes_df)] <- "0"
  
  # --- 5. Map Character Alleles to Integers --- 
  all_alleles <- unique(unlist(genotypes_df[, -1])) 
  alleles_to_map <- all_alleles[all_alleles != "0"] 
  allele_map <- setNames(as.character(1:length(alleles_to_map)), alleles_to_map) 
  allele_map["0"] <- "0" 
  
  genotypes_df[, -1] <- lapply(genotypes_df[, -1], function(col) { 
    unname(allele_map[as.character(col)]) 
  }) 
  
  return(list(genotypes_df = genotypes_df, markers_df = markers_df, locus_map_df = locus_map_df)) 
}

# ------------------------------------------------------------------------------
# SECTION 1B: COLONY RUN PARAMETERS 
# ------------------------------------------------------------------------------
COLONY_PARAMETERS <- list(
  project_name = "rgSibship", 
  analysis_type = "sibship",   
  
  # Core Run Parameters
  spp.type = 2,              # 2 = Dioecious
  inbreeding = 0,            # 0 = No
  ploidy = 0,                # 0 = Diploid
  fem.gamy = 0,              # 0 = Polygamy
  mal.gamy = 0,              # 0 = Polygamy
  
  # Likelihood & Alleles
  full.likelihood = 1,       # 1 = Full Likelihood
  likelihood.precision = 1,  # 1 = Medium Precision
  run.length = 2,            # 2 = Medium Run Length
  update.alfs = 0,           # 0 = No 
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
message("\n--- Starting COHORT + LOCATION COLONY File Generation ---")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("\n--- Loading Biological Data & Inferring Cohorts ---")

# Read CSV and clean names
bio_data <- read_xlsx(path = bio_data_path, sheet = 1)
names(bio_data) <- trimws(names(bio_data)) 
names(bio_data) <- trimws(gsub(pattern = " ", replacement = "_", x = names(bio_data)))
names(bio_data) <- trimws(gsub(pattern = "\\(|\\)", replacement = "", x = names(bio_data)))
names(bio_data) <- trimws(gsub(pattern = "_\\[|\\]_", replacement = "", x = names(bio_data)))

# Clean bad IDs
bio_data <- bio_data[!grepl(pattern = "NTC", x = bio_data$ID), ]
bio_data <- bio_data[!grepl(pattern = "999", x = bio_data$ID), ]

# Infer Cohorts
bio_data$Sample_Year <- as.numeric(substr(x = bio_data$ID, start = 1, stop = 4))
bio_data$Age_num <- suppressWarnings(as.numeric(bio_data$Age))
bio_data$Birth_Year <- bio_data$Sample_Year - bio_data$Age_num

# Filter out NA birth years and clean the Region column
bio_data_clean <- bio_data[!is.na(bio_data$Birth_Year) & !is.na(bio_data$Fish_ID), ]
bio_data_clean$Region <- trimws(bio_data_clean$Region)

# ------------------------------------------------------------------------------
# SECTION 2B: GROUP BY COHORT & REGION, FILTER MINIMUM 30 SAMPLES
# ------------------------------------------------------------------------------
message("\n--- Filtering Cohort & Region groups (Minimum 30 samples required) ---")

# Group by cohort and region, and only keep groups with 30 or more individuals
bio_data_limited <- bio_data_clean %>%
  group_by(Birth_Year, Region) %>%
  filter(n() >= 30) %>%
  ungroup()

# Create a master list of unique populations we need to process
unique_pops <- bio_data_limited %>%
  select(Birth_Year, Region) %>%
  distinct() %>%
  arrange(Birth_Year, Region)
unique_pops
message(paste("Generated", nrow(unique_pops), "unique Cohort-Location populations to process."))

# ==============================================================================
# SECTION 3: LOAD GLOBAL VCF
# ==============================================================================
message("\n--- Loading Global VCF file ---")
vcf <- read.vcfR(vcf_file_path, verbose = FALSE)

vcf_samples <- colnames(vcf@gt)[-1]
message(paste("Global VCF contains:", nrow(vcf), "loci and", length(vcf_samples), "samples."))

# ==============================================================================
# SECTION 4: MAIN PROCESSING LOOP (SUBSET, FILTER, & WRITE)
# ==============================================================================

# Start loop over unique populations
for (i in 1:nrow(unique_pops)) {
  
  current_cohort <- unique_pops$Birth_Year[i]
  current_region <- unique_pops$Region[i]
  
  # Format name without spaces (e.g., "rgSibship_Albany_2022")
  pop_name_clean <- paste0(COLONY_PARAMETERS$project_name, "_", current_region, "_", current_cohort)
  
  message(paste0("\n========================================================"))
  message(paste0("Processing ", i, "/", nrow(unique_pops), ": ", pop_name_clean))
  
  # Get the sample IDs for this specific cohort + region
  pop_target_samples <- bio_data_limited %>%
    filter(Birth_Year == current_cohort, Region == current_region) %>%
    pull(Fish_ID)
  
  # Convert numeric IDs to character to match VCF column names
  pop_target_samples <- as.character(pop_target_samples)
  
  # Find which of these target samples actually exist in the VCF
  valid_samples_in_vcf <- intersect(pop_target_samples, vcf_samples)
  missing_samples <- setdiff(pop_target_samples, vcf_samples)
  
  message(paste("  -> Target bio samples:", length(pop_target_samples)))
  message(paste("  -> Samples successfully found in VCF:", length(valid_samples_in_vcf)))
  
  if (length(valid_samples_in_vcf) < 2) {
    message("  -> SKIPPING: Not enough genotyped samples to run COLONY.")
    next
  }
  
  # --- Step 4.1: Subset the VCF to only include this group's samples ---
  vcf_subset <- vcf[, c("FORMAT", valid_samples_in_vcf)]
  
  # --- Step 4.2: Drop Loci that are 100% missing in this specific subset ---
  gt_subset <- extract.gt(vcf_subset, element = "GT")
  
  # Calculate missingness per locus (rows)
  missing_rates <- rowMeans(is.na(gt_subset) | gt_subset == "./.")
  
  # Identify loci that are NOT 100% missing
  loci_to_keep_indices <- which(missing_rates < 1.0)
  loci_dropped <- nrow(vcf_subset) - length(loci_to_keep_indices)
  
  message(paste("  -> Filtering Loci: Dropping", loci_dropped, "loci that are 100% missing in this population."))
  message(paste("  -> Retaining", length(loci_to_keep_indices), "polymorphic loci."))
  
  # Apply the loci filter
  vcf_final_filtered <- vcf_subset[loci_to_keep_indices, ]
  
  # --- Step 4.3: Convert to COLONY Format ---
  message("  -> Converting to COLONY format...")
  colony_data <- convert_vcf_to_colony_format(
    vcf_obj = vcf_final_filtered,
    marker_type = COLONY_PARAMETERS$marker_type,
    error_rate_1 = COLONY_PARAMETERS$error_rate_geno,
    error_rate_2 = COLONY_PARAMETERS$error_rate_other
  )
  
  # --- Step 4.4: Write the .dat file ---
  dat_filename <- file.path(output_dir, paste0(pop_name_clean, ".dat"))
  
  message(paste("  -> Writing .dat file:", basename(dat_filename)))
  
  colony.dat.create_FLEX.R(
    kids = colony_data$genotypes_df,
    markers = colony_data$markers_df,
    moms = NA, 
    dads = NA, 
    output_filename = dat_filename,
    project_name = pop_name_clean,
    
    # Pass parameters
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
}

message("\n========================================================")
message("--- Pipeline Finished! All location-cohort files written. ---")
