# ==============================================================================
#
# R SCRIPT FOR AUTOMATED COLONY INPUT FILE GENERATION
#
# Author: Gemini
# Date: 20 October 2025
# Version: 3.3 (Robust MAF Handling)
#
# Description:
# This script automates the creation of input files for the pedigree
# analysis software, COLONY. It takes a VCF file as input, calculates
# the Minor Allele Frequency (MAF) for each SNP, and then generates
# multiple subsets of the data based on a predefined analysis plan.
#
# The analysis plan specifies different methods for subsetting:
#   1. 'maf': Selects loci with the highest MAF values.
#   2. 'random': Selects a random set of loci.
#
# This version focuses on generating the .dat files locally and also creates
# a list of all generated filenames, suitable for use with an HPC
# sbatch array job.
#
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: CONFIGURATION
# --- Please edit the variables in this section for your analysis ---
# ------------------------------------------------------------------------------

# --- Set Working Directory ---
# This path is specific to your machine.
working_directory <- "G:/Shared drives/Sard Lab/Round Goby/RADseq/R/pedigree reconstruction"
setwd(working_directory)
list.files()
# Required libraries. Ensure they are installed.
# install.packages(c("vcfR", "dplyr", "readr"))
library(vcfR)
library(dplyr)
library(readr)

# --- Input & Output Paths ---
# IMPORTANT: Paths should be relative to your working directory or absolute.
vcf_file_path <- "../post_RAD analysis/Output/populations.filtered.final.vcf.gz" # <-- Set path to your VCF
output_dir <- "colony_runs_output"           # Directory to save generated files

# --- Analysis Plan ---
# Define the combinations of subset types and number of loci to generate.
analysis_plan <- expand.grid(
  num_loci = c(250, 500, 1000, 2000, 5000, 8000, 11664),
  subset_type = c("maf", "random"),
  stringsAsFactors = FALSE
)

# ==============================================================================
# --- DO NOT EDIT BELOW THIS LINE UNLESS YOU ARE SURE ---
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 2: SETUP & PRE-CHECKS
# ------------------------------------------------------------------------------
message("--- Starting COLONY File Generation Pipeline ---")
message(paste("Current working directory:", getwd()))

# Create the output directory if it doesn't exist
message(paste("Checking for output directory:", output_dir))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Output directory created.")
} else {
  message("Output directory already exists.")
}

# Validate that the analysis plan is a non-empty data frame
stopifnot(
  is.data.frame(analysis_plan),
  nrow(analysis_plan) > 0,
  "num_loci" %in% names(analysis_plan),
  "subset_type" %in% names(analysis_plan)
)
message("Analysis plan validated successfully. Number of runs to process: ", nrow(analysis_plan))
print(analysis_plan)


# ------------------------------------------------------------------------------
# SECTION 3: LOAD AND VALIDATE VCF FILE
# ------------------------------------------------------------------------------
message("\n--- SECTION 3: Loading VCF file ---")
message(paste("Reading VCF from:", vcf_file_path))

# Check if file exists before trying to read
if (!file.exists(vcf_file_path)) {
  stop(paste("FATAL ERROR: VCF file not found at:", vcf_file_path))
}

vcf <- read.vcfR(vcf_file_path, verbose = FALSE)
message("VCF file loaded into memory.")

# --- Validation of vcfR object ---
message("Validating the loaded vcfR object...")
stopifnot(
  inherits(vcf, "vcfR"),
  nrow(vcf) > 0,
  ncol(vcf) > 1 # Must have at least one sample
)
message(paste(
  "Validation successful. VCF contains:",
  nrow(vcf), "variants and",
  (ncol(vcf@gt) - 1), "samples."
))
# Ensure the maximum requested number of loci does not exceed available loci
if (max(analysis_plan$num_loci) > nrow(vcf)) {
  stop(paste0(
    "ERROR: The analysis plan requests ", max(analysis_plan$num_loci), " loci, ",
    "but the VCF file only contains ", nrow(vcf), " loci. ",
    "Please reduce num_loci in the analysis plan."
  ))
}
message("The number of requested loci is valid.")


# ------------------------------------------------------------------------------
# SECTION 4: PRE-CALCULATE AND VALIDATE MINOR ALLELE FREQUENCY (MAF)
# ------------------------------------------------------------------------------
message("\n--- SECTION 4: PRE-CALCULATE AND VALIDATE MINOR ALLELE FREQUENCY (MAF) ---")
message("Calculating MAF for all loci...")

# The maf() function returns a matrix with multiple columns of information.
maf_matrix <- maf(vcf, element = 2)
message("MAF calculation complete.")

# --- Validation of MAF results ---
message("Validating the structure of the MAF matrix...")
# The user correctly identified that we should check nrow, not length.
stopifnot(
  is.matrix(maf_matrix),
  nrow(maf_matrix) == nrow(vcf)
)
message("Validation successful: MAF matrix has the correct number of rows.")


# We are interested in the second column, which contains the MAF values.
maf_vector <- maf_matrix[, 2]
message("Extracted MAF vector from the second column of the matrix.")

# --- Temporary code for debugging NA MAF loci ---
# You can insert this right after the MAF validation in SECTION 4

message("\n--- DEBUGGING: Investigating loci with NA MAF values ---")

# 1. Get the row indices of the loci where MAF is NA
na_maf_indices <- which(is.na(maf_vector))
message(paste("Found", length(na_maf_indices), "loci with NA MAF."))

if (length(na_maf_indices) > 0) {
  # 2. Extract the 'fix' data (like CHROM, POS, ID) for these loci
  problematic_loci_info <- vcf@fix[na_maf_indices, ]
  
  # 3. Print the first 10 problematic loci to the console
  message("--- First 10 loci with NA MAF: ---")
  print(head(problematic_loci_info, 10))
  
  # 4. (Optional) To see all of them, save to a CSV file
  problem_loci_filename <- file.path(output_dir, "loci_with_na_maf.csv")
  write.csv(problematic_loci_info, problem_loci_filename, row.names = FALSE)
  message(paste("-> Full list of problematic loci saved to:", problem_loci_filename))
}
message("--- End of DEBUGGING section ---\n")

# Now, we validate the extracted vector itself.
stopifnot(is.numeric(maf_vector))

# Identify and report loci with NA MAF values from the vector
num_na_maf <- sum(is.na(maf_vector))
if (num_na_maf > 0) {
  message(paste(
    "Warning:", num_na_maf, "loci have undefined (NA) MAF values.",
    "These loci will be excluded from selection."
  ))
}

# Get a vector of indices for loci that have a valid, non-NA MAF value.
# This is the pool of loci we can select from.
valid_loci_indices <- which(!is.na(maf_vector))
message(paste("Total number of loci available for subsetting:", length(valid_loci_indices)))

# Final check to ensure we have enough valid loci for the analysis plan
if (max(analysis_plan$num_loci) > length(valid_loci_indices)) {
  stop(paste0(
    "ERROR: After excluding loci with NA MAF, only ", length(valid_loci_indices),
    " loci remain, but the analysis plan requests up to ", max(analysis_plan$num_loci),
    " loci. Cannot proceed."
  ))
}


# ------------------------------------------------------------------------------
# SECTION 5: MAIN PROCESSING LOOP
# ------------------------------------------------------------------------------
message("\n--- SECTION 5: Starting main processing loop to generate files ---")

# Initialize a vector to store the names of the generated .dat files
generated_dat_files <- c()

for (i in 1:nrow(analysis_plan)) {
  
  # --- Get parameters for the current run ---
  run_params <- analysis_plan[i, ]
  num_loci_to_select <- run_params$num_loci
  subset_method <- run_params$subset_type
  
  run_name <- sprintf("run_%02d_%s_%04d-loci", i, subset_method, num_loci_to_select)
  
  message(paste0(
    "\n--- Processing Run ", i, "/", nrow(analysis_plan),
    ": (", run_name, ") ---"
  ))
  
  # --- Step 5.1: Select Loci Indices based on Method ---
  message("Selecting loci based on method '", subset_method, "'...")
  
  # ...
  
  if (subset_method == "maf") {
    # Order the VALID loci by their MAF values (descending)
    ordered_valid_indices <- valid_loci_indices[order(maf_vector[valid_loci_indices], decreasing = TRUE)]
    
    # Select the top N indices from this ordered list
    indices_to_keep <- ordered_valid_indices[1:num_loci_to_select]
    
  } else if (subset_method == "random") {
    # Randomly sample N indices from the pool of valid loci
    indices_to_keep <- sample(valid_loci_indices, size = num_loci_to_select, replace = FALSE)
    
  } else {
    warning(paste("Subset type '", subset_method, "' is not recognized. Skipping run.", i))
    next # Skip to the next iteration
  }
  
  # --- Validation of selected indices ---
  message("Loci selected. Validating index vector...")
  stopifnot(
    length(indices_to_keep) == num_loci_to_select,
    !any(is.na(indices_to_keep)),
    is.numeric(indices_to_keep)
  )
  message("Index vector validation successful.")
  
  # --- Step 5.2: Subset the vcfR Object ---
  message("Subsetting main vcfR object...")
  
  # Use the validated numeric index vector to subset the vcfR object.
  # This is the critical step that previously failed.
  vcf_subset <- vcf[indices_to_keep, ]
  
  # --- Validation of the subsetted object ---
  message("Subsetting complete. Validating new vcfR object...")
  stopifnot(
    nrow(vcf_subset) == num_loci_to_select,
    ncol(vcf_subset) == ncol(vcf)
  )
  message("Subsetted vcfR object is valid.")
  
  # --- Step 5.3: Generate COLONY .dat file ---
  dat_filename <- file.path(output_dir, paste0(run_name, ".dat"))
  message(paste("Generating COLONY .dat file:", dat_filename))
  
  # Extract genotypes as alleles (e.g., "A/G", "C/C") for robustness
  gt_alleles <- extract.gt(vcf_subset, element = "GT", return.alleles = TRUE)
  
  # Replace the '/' separator with a space
  gt_formatted <- gsub("/", " ", gt_alleles)
  
  # Transpose the matrix so that rows are samples and columns are loci
  gt_transposed <- t(gt_formatted)
  
  # Create one line per individual, with SampleID followed by genotypes
  colony_lines <- apply(gt_transposed, 1, function(row) paste(row, collapse = " "))
  colony_lines <- paste(rownames(gt_transposed), colony_lines, sep = " ")
  
  # Prepare the final content for the .dat file
  file_content <- c(
    paste("NUMLOCI =", num_loci_to_select),
    paste("NUMINDS =", (ncol(vcf@gt) - 1)),
    "LENGTHS = 0", "PLOIDY = 0", "GENOTYPE = 1", "PHASED = 0",
    "MARKERTYPE = 0", "MISSINGALLELE = '0'",
    "", # Blank line
    paste0("LOCUS", 1:num_loci_to_select, collapse = " "), # Locus names
    "", # Blank line
    colony_lines # Individual genotypes
  )
  
  # Write the file
  writeLines(file_content, dat_filename)
  message(".dat file written successfully.")
  
  # Add the generated filename to our list for the sbatch array
  generated_dat_files <- c(generated_dat_files, basename(dat_filename))
}

# ------------------------------------------------------------------------------
# SECTION 6: GENERATE FILE LIST FOR SBATCH ARRAY JOB
# ------------------------------------------------------------------------------
message("\n--- SECTION 6: Generating file list for HPC job array ---")
job_list_filename <- file.path(output_dir, "sbatch_job_list.txt")
writeLines(generated_dat_files, job_list_filename)
message(paste("-> A file containing the names of all generated .dat files has been saved to:"))
message(paste("   ", job_list_filename))
message("-> Use this file on your HPC for your sbatch array job.")


message("\n--- All .dat files generated successfully. ---")
message("--- Pipeline Finished ---")


