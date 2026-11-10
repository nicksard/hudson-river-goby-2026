# ==============================================================================
#
# R SCRIPT FOR AUTOMATED COLONY INPUT FILE GENERATION
#
# Author: Gemini & Nick Sard
# Date: 20 October 2025
# Version: 4.2 (Bug Fixes: Integer Mapping & Homozygotes)
#
# Description:
# This script integrates the VCF-subsetting pipeline (V3.3) with the
# robust 'colony.dat.create_FLEX.R' function.
#
# IT PERFORMS 5 KEY TASKS:
#   1. Loads a VCF and calculates MAF.
#   2. Runs a data quality check for partially missing data.
#   3. Subsets the VCF based on an analysis plan (MAF or random).
#   4. FOR EACH SUBSET:
#      a. Converts VCF genotypes (e.g., "A/G", "T") into the 2-column,
#         INTEGER format (e.g., "1 2", "3 3") required by COLONY.
#      b. Auto-generates the required "markers" data frame.
#      c. Calls 'colony.dat.create_FLEX.R' to write the final,
#         all-in-one .dat file for sibship reconstruction.
#      d. Writes a .csv "map" file to link generic LOCUS names
#         (LOCUS1) back to their original VCF names (scaffold_1:12345).
#
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1A: CORE CONFIGURATION
# ------------------------------------------------------------------------------

# --- Set Working Directory ---
working_directory <- "G:/Shared drives/Sard Lab/Round Goby/RADseq/R/pedigree reconstruction"
setwd(working_directory)
list.files()

# --- Required libraries ---
library(vcfR)
library(dplyr)
library(readr)
library(tidyr)  # Needed for tidyr::separate()
library(tibble) # Needed for tibble::rownames_to_column()

#this is to stop automatic scientific notations in COLONY dat file
options(scipen=999)

# --- Input & Output Paths ---
vcf_file_path <- "../post_RAD analysis/Output/populations.filtered.final.vcf.gz"
output_dir <- "colony_runs_output_all_years_V1"

# --- Source the FLEX Function ---
# (Assumes the file is in your working directory)
flex_function_path <- "colony.dat.create_FLEX.R"
if (!file.exists(flex_function_path)) {
  stop(paste("FATAL ERROR: The function file", flex_function_path, "was not found."))
}
message(paste("Loading flexible function from:", flex_function_path))
source(flex_function_path)


# ------------------------------------------------------------------------------
# SECTION 1B: COLONY RUN PARAMETERS
# ------------------------------------------------------------------------------
# These values will be passed to the colony.dat.create_FLEX.R function
# for every file generated.
COLONY_PARAMETERS <- list(
  project_name = "rgSibship", # Base name for project
  analysis_type = "sibship",    # CRITICAL: Ensures 0 parents are set
  
  # Core Run Parameters
  run.length = 2,           # 2 = Medium
  spp.type = 2,             # 2 = Dioecious
  inbreeding = 0,           # 0 = No
  ploidy = 0,               # 0 = Diploid
  fem.gamy = 0,             # 0 = Polygamy
  mal.gamy = 0,             # 0 = Polygamy
  
  # Likelihood & Alleles
  full.likelihood = 1,      # 1 = Full Likelihood
  likelihood.precision = 2, # 2 = Medium
  update.alfs = 0,          # 0 = No, use empirical
  known.alfs = 0,
  
  # Error Rates (for Markers DF)
  marker_type = 0,          # 0 = Codominant
  error_rate_geno = 0.02,   # 2% error (e.g., allelic dropout)
  error_rate_other = 0.001, # 0.1% error (e.g., mutation)
  
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

# ------------------------------------------------------------------------------
# SECTION 1C: ANALYSIS PLAN ("Kitty Data Set" Test)
# ------------------------------------------------------------------------------
# Per your request, we are starting with a simple 1-run test.
#analysis_plan <- expand.grid(
#  num_loci = c(250),
#  subset_type = c("random"),
#  stringsAsFactors = FALSE
#)
#message("--- Using 'Kitty Data Set' Analysis Plan ---")
#print(analysis_plan)

# To run the full analysis later, uncomment this line:
analysis_plan <- expand.grid(
   num_loci = c(250, 500, 1000, 2000, 5000, 8000, 11664),
   subset_type = c("maf", "random"),
   stringsAsFactors = FALSE
)
print(analysis_plan)
# ==============================================================================
# --- DO NOT EDIT BELOW THIS LINE UNLESS YOU ARE SURE ---
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 2: SETUP & PRE-CHECKS
# ------------------------------------------------------------------------------
message("\n--- Starting COLONY File Generation Pipeline (V4.1) ---")
message(paste("Current working directory:", getwd()))

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste("Output directory created:", output_dir))
} else {
  message(paste("Output directory already exists:", output_dir))
}


# ------------------------------------------------------------------------------
# SECTION 3: LOAD AND VALIDATE VCF FILE
# ------------------------------------------------------------------------------
message("\n--- SECTION 3: Loading VCF file ---")
if (!file.exists(vcf_file_path)) {
  stop(paste("FATAL ERROR: VCF file not found at:", vcf_file_path))
}
vcf <- read.vcfR(vcf_file_path, verbose = FALSE)
message("VCF file loaded into memory.")
message(paste(
  "VCF contains:",
  nrow(vcf), "variants and",
  (ncol(vcf@gt) - 1), "samples."
))
if (max(analysis_plan$num_loci) > nrow(vcf)) {
  stop(paste0(
    "ERROR: The analysis plan requests ", max(analysis_plan$num_loci), " loci, ",
    "but the VCF file only contains ", nrow(vcf), " loci."
  ))
}

# ------------------------------------------------------------------------------
# SECTION 3B: DATA QUALITY CHECK (NEW)
# ------------------------------------------------------------------------------
message("\n--- SECTION 3B: Running Data Quality Check ---")

# Extract all GT (genotype) data
gt_matrix <- extract.gt(vcf, element = "GT")

# Check for partially missing genotypes (e.g., "./A" or "A/.")
# These are genotypes that contain a '.' but are not *only* './.'
partial_missing_any <- any(grepl("\\.", gt_matrix) & gt_matrix != "./.")

if (partial_missing_any) {
  # Find the specific loci and samples
  partial_indices <- which(grepl("\\.", gt_matrix) & gt_matrix != "./.", arr.ind = TRUE)
  num_partial <- nrow(partial_indices)
  
  # Get locus names (CHROM:POS) for the problematic rows
  locus_names <- paste0(vcf@fix[partial_indices[, "row"], "CHROM"], ":", vcf@fix[partial_indices[, "row"], "POS"])
  # Get sample names for the problematic columns
  sample_names <- colnames(gt_matrix)[partial_indices[, "col"]]
  
  problem_report <- data.frame(
    Sample = sample_names,
    Locus = locus_names,
    Genotype = gt_matrix[partial_indices]
  )
  
  message(paste("WARNING: Found", num_partial, "partially missing genotypes (e.g., 'A/.' or './T')."))
  message("This may indicate an 'upstream bug' in VCF filtering.")
  message("These will be converted to 'Allele 0' (e.g., '1 0').")
  message("First 10 instances:")
  print(head(problem_report, 10))
  
  # Optional: Write a full report
  # report_filename <- file.path(output_dir, "partial_missing_data_report.csv")
  # write.csv(problem_report, report_filename, row.names = FALSE)
  # message(paste("Full report saved to:", report_filename))
  
} else {
  message("Data quality check passed: No partially missing genotypes (like 'A/.') found.")
}


# ------------------------------------------------------------------------------
# SECTION 4: PRE-CALCULATE AND VALIDATE MINOR ALLELE FREQUENCY (MAF)
# ------------------------------------------------------------------------------
message("\n--- SECTION 4: Calculating MAF ---")

# We only need to do this if 'maf' is in the plan
if ("maf" %in% analysis_plan$subset_type) {
  maf_matrix <- maf(vcf, element = 2)
  maf_vector <- maf_matrix[, 2]
  
  num_na_maf <- sum(is.na(maf_vector))
  if (num_na_maf > 0) {
    message(paste(
      "Warning:", num_na_maf, "loci have NA MAF values and will be excluded."
    ))
  }
  valid_loci_indices <- which(!is.na(maf_vector))
  
  if (max(analysis_plan$num_loci[analysis_plan$subset_type == "maf"]) > length(valid_loci_indices)) {
    stop(paste0(
      "ERROR: After excluding NA MAF loci, only ", length(valid_loci_indices),
      " loci remain, but the 'maf' plan requests up to ", 
      max(analysis_plan$num_loci[analysis_plan$subset_type == "maf"]),
      " loci. Cannot proceed."
    ))
  }
  
} else {
  message("Skipping MAF calculation: 'maf' method not in analysis plan.")
  # We can just use all loci for 'random'
  valid_loci_indices <- 1:nrow(vcf)
  maf_vector <- NULL # Not needed
}
message(paste("Total number of loci available for subsetting:", length(valid_loci_indices)))


# ------------------------------------------------------------------------------
# SECTION 5: NEW HELPER FUNCTION (The "Translator")
# ------------------------------------------------------------------------------

#' Converts a vcfR object into the 2-column format required
#' by colony.dat.create_FLEX.R
#'
#' @param vcf_obj A vcfR object (already subsetted).
#' @param marker_type The value for the marker type row (e.g., 0).
#' @param error_rate_1 The value for the first error row (e.g., 0.02).
#' @param error_rate_2 The value for the second error row (e.g., 0.001).
#'
#' @return A list containing three data frames:
#'   $genotypes_df (for 'kids' arg)
#'   $markers_df (for 'markers' arg)
#'   $locus_map_df (for traceability)

convert_vcf_to_colony_format <- function(vcf_obj, marker_type, error_rate_1, error_rate_2) {
  
  
  n_loci <- nrow(vcf_obj)
  n_inds <- ncol(vcf_obj@gt) - 1
  
  # --- 1. Get Locus Names & Create Map ---
  message("Creating locus map for traceability...")
  
  # Get original names from VCF 'fix' slot
  original_locus_names <- vcf_obj@fix[, "ID"]
  
  # Fallback if ID column is missing or empty (e.g., ".")
  if (all(is.na(original_locus_names)) || all(original_locus_names == ".")) {
    message("VCF 'ID' column is empty. Using 'CHROM:POS' as original names.")
    original_locus_names <- paste0(vcf_obj@fix[, "CHROM"], ":", vcf_obj@fix[, "POS"])
  }
  
  # Create generic names (LOCUS1, LOCUS2, ...)
  generic_locus_names <- paste0("LOCUS", 1:n_loci)
  
  # Create the map data frame
  locus_map_df <- data.frame(
    generic_name = generic_locus_names,
    original_name = original_locus_names,
    stringsAsFactors = FALSE
  )
  
  # --- 2. Create the 'markers_df' ---
  message("Creating 'markers' data frame...")
  markers_df <- data.frame(matrix(
    data = 0,
    nrow = 3,
    ncol = n_loci
  ))
  
  # Set column names (e.g., LOCUS1, LOCUS2)
  colnames(markers_df) <- generic_locus_names
  
  # Fill in the error rate rows
  markers_df[1, ] <- marker_type
  markers_df[2, ] <- error_rate_1
  markers_df[3, ] <- error_rate_2
  
  
  # --- 3. Create the 'genotypes_df' (for 'kids') ---
  message("Extracting and reformatting genotypes...")
  
  # Extract as "A/G", "C/T", etc.
  gt_alleles <- extract.gt(vcf_obj, element = "GT", return.alleles = TRUE)
  gt_alleles
  
  # Transpose and format as data frame
  gt_t <- as.data.frame(t(gt_alleles), stringsAsFactors = FALSE)
  gt_t <- tibble::rownames_to_column(gt_t, var = "Sample_ID")
  
  # Rename VCF columns (e.g., "scaffold_1:12345") to our standard names
  # (Sample_ID, LOCUS1, LOCUS2, ...)
  colnames(gt_t) <- c("Sample_ID", generic_locus_names)
  
  # This is the final data frame we will build
  genotypes_df <- gt_t[, "Sample_ID", drop = FALSE]
  
  message("Splitting 'A/G' columns into 'A' 'G' format...")

  for (locus in generic_locus_names) {
    
    new_col_A <- paste0(locus, "_A")
    new_col_B <- paste0(locus, "_B")
    
    # Get the raw genotype column (e.g., "A/G", "T/T", "A", "A/.", "./.")
    raw_gts <- gt_t[[locus]]
    raw_gts
    
    # --- BUG FIX 1: Handle Homozygotes ---
    # Find genotypes that are NOT missing and do NOT have a "/"
    # (e.g., "A", "T"). These are simplified homozygotes.
    hom_idx <- !grepl("/", raw_gts, fixed = TRUE) & 
      !grepl(".", raw_gts, fixed = TRUE) & 
      raw_gts != "NA"
    hom_idx[is.na(hom_idx)] <- FALSE

    # Duplicate them (e.g., "T" -> "T/T")
    if(any(hom_idx)) {
      raw_gts[hom_idx] <- paste(raw_gts[hom_idx], raw_gts[hom_idx], sep = "/")
    }
    
    # --- BUG FIX 2: Handle Missing Data ---
    # Globally replace the VCF missing character "." with "0"
    # "A/." becomes "A/0"
    # "./." becomes "0/0"
    raw_gts <- gsub(".", "0", raw_gts, fixed = TRUE)
    
    # --- Now, all data is "A/G", "T/T", "A/0", "0/0" ---
    # We can safely separate it.
    split_cols <- tidyr::separate(
      data.frame(gts = raw_gts),
      col = "gts",
      into = c(new_col_A, new_col_B),
      sep = "/",
      remove = TRUE
    )
    
    # Final safety net: replace any lingering NAs with "0"
    split_cols[is.na(split_cols)] <- "0"
    
    # Add the two new, correct columns
    genotypes_df[[new_col_A]] <- split_cols[[new_col_A]]
    genotypes_df[[new_col_B]] <- split_cols[[new_col_B]]
  }
  
  # --- 4. Clean Missing Data (Redundant but safe) ---
  # This section is now largely redundant but acts as a final safety check.
  message("Cleaning missing data (NA -> 0)...")
  genotypes_df[is.na(genotypes_df)] <- "0"
  
  
  # --- 5. Map Character Alleles to Integers ---
  # COLONY requires integer alleles (e.g., 1, 2, 3) not characters (A, C, G)
  message("Mapping character alleles (A,C,T,G) to integers (1,2,3,4)...")
  
  # Get all unique allele characters from the genotype columns
  all_alleles <- unique(unlist(genotypes_df[, -1]))
  
  # Remove the missing data character "0" from the list to be mapped
  alleles_to_map <- all_alleles[all_alleles != "0"]
  
  # Create the map (e.g., A -> "1", C -> "2", ...)
  # We use characters for the map values for safe substitution
  allele_map <- setNames(as.character(1:length(alleles_to_map)), alleles_to_map)
  
  # Add the "0" -> "0" mapping back in to preserve missing data
  allele_map["0"] <- "0"
  
  # Apply the map to all genotype columns (all columns except Sample_ID)
  # Use as.character() inside to protect against factors
  genotypes_df[, -1] <- lapply(genotypes_df[, -1], function(col) {
    unname(allele_map[as.character(col)])
  })
  
  message("Data reformatting complete.")
  return(list(
    genotypes_df = genotypes_df,
    markers_df = markers_df,
    locus_map_df = locus_map_df
  ))
}


# ------------------------------------------------------------------------------
# SECTION 6: MAIN PROCESSING LOOP
# ------------------------------------------------------------------------------
message("\n--- SECTION 6: Starting main processing loop to generate files ---")

# Initialize vectors to store the names of the generated files
generated_dat_files <- c()
generated_map_files <- c()
i <- 1
for (i in 1:nrow(analysis_plan)) {
  
  # --- Get parameters for the current run ---
  run_params <- analysis_plan[i, ]
  num_loci_to_select <- run_params$num_loci
  subset_method <- run_params$subset_type
  
  # --- NEW: Create shorter run name for filenames ---
  # Gsub "random" to "rand"
  short_method <- gsub("random", "rand", subset_method)
  # Format as R01_rand_L0250
  run_name_short <- sprintf("R%02d_%s_L%04d", i, short_method, num_loci_to_select)
  
  # Create the full project name (respects 40-char limit)
  project_name_full <- paste(COLONY_PARAMETERS$project_name, run_name_short, sep = "_")
  
  message(paste0(
    "\n--- Processing Run ", i, "/", nrow(analysis_plan),
    ": (", project_name_full, ") ---"
  ))
  
  # --- Step 6.1: Select Loci Indices based on Method ---
  message("Selecting loci based on method '", subset_method, "'...")
  
  if (subset_method == "maf") {
    ordered_valid_indices <- valid_loci_indices[order(maf_vector[valid_loci_indices], decreasing = TRUE)]
    indices_to_keep <- ordered_valid_indices[1:num_loci_to_select]
    
  } else if (subset_method == "random") {
    indices_to_keep <- sample(valid_loci_indices, size = num_loci_to_select, replace = FALSE)
    
  } else {
    warning(paste("Subset type '", subset_method, "' is not recognized. Skipping run.", i))
    next # Skip to the next iteration
  }
  
  # --- Step 6.2: Subset the vcfR Object ---
  message("Subsetting main vcfR object...")
  vcf_subset <- vcf[indices_to_keep, ]
  
  # --- Step 6.3: Convert VCF subset to COLONY format ---
  message("Calling 'convert_vcf_to_colony_format'...")
  colony_data <- convert_vcf_to_colony_format(
    vcf_obj = vcf_subset,
    marker_type = COLONY_PARAMETERS$marker_type,
    error_rate_1 = COLONY_PARAMETERS$error_rate_geno,
    error_rate_2 = COLONY_PARAMETERS$error_rate_other
  )
  
  # --- Step 6.4: Generate COLONY .dat file using FLEX function ---
  # Use the new full project name for the .dat file
  dat_filename <- file.path(output_dir, paste0(project_name_full, ".dat"))
  
  message(paste("Calling 'colony.dat.create_FLEX.R' to write:", dat_filename))
  
  colony.dat.create_FLEX.R(
    # Data Inputs
    kids = colony_data$genotypes_df,
    markers = colony_data$markers_df,
    moms = NA, # Ignored by "sibship" mode
    dads = NA, # Ignored by "sibship" mode
    
    # File/Project Names
    output_filename = dat_filename,
    project_name = project_name_full,
    
    # Pass all other parameters from our list
    analysis_type = COLONY_PARAMETERS$analysis_type,
    update.alfs = COLONY_PARAMETERS$update.alfs,
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
    run.length = COLONY_PARAMETERS$run.length,
    monitor = COLONY_PARAMETERS$monitor,
    windows.version = COLONY_PARAMETERS$windows.version,
    full.likelihood = COLONY_PARAMETERS$full.likelihood,
    likelihood.precision = COLONY_PARAMETERS$likelihood.precision,
    prob.mom = COLONY_PARAMETERS$prob.mom,
    prob.dad = COLONY_PARAMETERS$prob.dad
  )
  
  message(".dat file written successfully.")
  generated_dat_files <- c(generated_dat_files, basename(dat_filename))
  
  # --- Step 6.5: Write Locus Map for Traceability ---
  # Use the new full project name for the map file
  map_filename <- file.path(output_dir, paste0(project_name_full, "_LOCUS_MAP.csv"))
  message(paste("Writing locus map for traceability:", map_filename))
  write.csv(colony_data$locus_map_df, file = map_filename, row.names = FALSE, quote = TRUE)
  generated_map_files <- c(generated_map_files, basename(map_filename))
}

# ------------------------------------------------------------------------------
# SECTION 7: GENERATE FILE LIST FOR SBATCH ARRAY JOB
# ------------------------------------------------------------------------------
message("\n--- SECTION 7: Generating file lists for HPC job array ---")
job_list_filename <- file.path(output_dir, "sbatch_job_list.txt")
writeLines(generated_dat_files, job_list_filename)
message(paste("-> A file containing the names of all .dat files has been saved to:"))
message(paste("   ", job_list_filename))
message("-> Use this file on your HPC for your sbatch array job.")

message(paste("->", length(generated_map_files), "locus map CSV files also written to the output directory."))

message("\n--- All files generated successfully. ---")
message("--- Pipeline Finished ---")