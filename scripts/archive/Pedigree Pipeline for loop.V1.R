# ==============================================================================
# Master Generator for COLONY Pedigree Analysis on HPC (v2.5 - Boolean Subsetting)
#
# AUTHOR: Nicholas Sard, Ph.D.
# ASSISTANT: Gemini
# DATE: October 17, 2025
#
# DESCRIPTION:
# This single, self-contained script automates the generation of input files
# and corresponding shell scripts for a large-scale COLONY analysis.
# This version now uses a more robust boolean vector method for subsetting
# the vcfR object, avoiding the need to manually set rownames and preventing
# potential errors. MAF calculation remains outside the loop for efficiency.
#
# ==============================================================================


# ==============================================================================
# Part 1: Load Libraries & Set Global Options
# ==============================================================================
message("--- Loading required libraries ---")
library(vcfR)
library(tidyverse)

# --- Prevent scientific notation in output files ---
options(scipen = 999)

# --- VCF File Path ---
# Path to your final, filtered VCF file.
VCF_FILE_PATH <- "../post_RAD analysis/Output/populations.filtered.final.vcf.gz"
message(paste(" -> Checking for VCF file at:", VCF_FILE_PATH))
if (!file.exists(VCF_FILE_PATH)) {
  stop("VCF file not found. Please check the 'VCF_FILE_PATH' variable.")
}
message(" -> VCF file found.")


# ==============================================================================
# Part 2: Define the Analysis Plan
# ==============================================================================
# This is the control center for your entire analysis.
# We use expand_grid to create all combinations of locus counts and subset methods.
message("--- Defining the analysis plan ---")

analysis_plan <- expand_grid(
  locus_subset = c("250", "2000", "4000", "5000", "8000", "all"),
  subset_method = c("MAF", "Random"),
  analysis_method = "FL",
  run_length = "medium",
  likelihood_precision = "medium",
  error_rate = 0.001,
  dropout_rate = 0.001
) %>%
  # The "all" subset is the same regardless of method, so we remove the duplicate
  # to avoid running it twice. We'll keep the "MAF" designated one.
  filter(!(locus_subset == "all" & subset_method == "Random"))

message(paste(" -> Analysis plan created with", nrow(analysis_plan), "unique runs."))


# ==============================================================================
# Part 3: Pre-computation and Data Loading
# ==============================================================================
message("\n--- Pre-computation Step ---")

# --- Pre-load the VCF file once to avoid repeated reading ---
message(paste(" -> Pre-loading VCF from:", VCF_FILE_PATH))
full_vcf <- read.vcfR(VCF_FILE_PATH)
message(" -> VCF loaded successfully.")
stopifnot("vcfR" %in% class(full_vcf)) # Check if it's a valid vcfR object
message(paste(" -> Full VCF dimensions:", nrow(full_vcf), "loci by", ncol(full_vcf@gt) -1, "individuals."))

# --- Pre-calculate MAF once for efficiency ---
message(" -> Pre-calculating MAF for all loci...")
maf_info <- vcfR::maf(full_vcf, element = 2)
message(paste("  -> Initial MAF calculation complete for", nrow(maf_info), "loci."))
maf_info_clean <- maf_info[complete.cases(maf_info), ]
message(paste("  -> Removed", nrow(maf_info) - nrow(maf_info_clean), "loci with NA MAF values."))

# --- Create a sorted list of locus IDs by MAF ---
sorted_maf_loci <- names(sort(maf_info_clean[,1], decreasing = TRUE))
message(paste(" -> Created a sorted list of", length(sorted_maf_loci), "valid loci for MAF subsetting."))


# ==============================================================================
# Part 4: Execute the Batch Generation Loop
# ==============================================================================
message(paste("\n--- Starting batch generation for", nrow(analysis_plan), "runs ---"))

for (i in 1:nrow(analysis_plan)) {
  
  # --- Get parameters for the current run ---
  params <- analysis_plan[i, ]
  locus_param <- if (params$locus_subset == "all") "all" else as.numeric(params$locus_subset)
  
  # --- Define a descriptive label for directories and job names ---
  locus_label <- if (is.numeric(locus_param)) {
    paste0("Loci_", locus_param, "_", params$subset_method)
  } else {
    "Loci_All" # The 'all' run is unique
  }
  
  # --- Define a structured output directory for this specific run ---
  run_dir <- file.path("Output",
                       paste0(params$analysis_method, "_", params$run_length, "_", params$likelihood_precision),
                       locus_label)
  if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE)
  
  message(paste("\n[", i, "/", nrow(analysis_plan), "] Processing Run:", locus_label))
  message(paste(" -> Output directory:", run_dir))
  
  
  # ============================================================================
  # <<< START: INLINED FILE GENERATION LOGIC (for easy debugging) >>>
  # ============================================================================
  
  ### --- 4.1: Data Preparation & Subsetting ---
  message(" -> Subsetting VCF...")
  vcf_subset <- NULL # Initialize to be safe
  
  if (!is.numeric(locus_param) || locus_param >= nrow(full_vcf)) {
    message(paste(" -> Using all", nrow(full_vcf), "available loci."))
    vcf_subset <- full_vcf
  } else if (params$subset_method == "MAF") {
    message(paste(" -> Selecting top", locus_param, "loci using boolean vector..."))
    
    # --- LOGIC CHECK: Ensure enough valid loci remain for subsetting ---
    if (length(sorted_maf_loci) < locus_param) {
      stop(paste("ERROR: Not enough valid loci to subset. Requested", locus_param, "but only", length(sorted_maf_loci), "remain after removing NAs."))
    }
    
    top_loci_ids <- sorted_maf_loci[1:locus_param]
    loci_to_keep_vector <- full_vcf@fix[, "ID"] %in% top_loci_ids
    vcf_subset <- full_vcf[loci_to_keep_vector, ]
    
  } else if (params$subset_method == "Random") {
    message(paste(" -> Randomly selecting", locus_param, "loci using boolean vector..."))
    all_loci_ids <- full_vcf@fix[, "ID"]
    random_loci_ids <- sample(all_loci_ids, size = locus_param)
    loci_to_keep_vector <- full_vcf@fix[, "ID"] %in% random_loci_ids
    vcf_subset <- full_vcf[loci_to_keep_vector, ]
  }
  
  # --- LOGIC CHECK ---
  stopifnot(!is.null(vcf_subset)) # Ensure the subset object was created.
  expected_loci <- if(is.numeric(locus_param)) min(locus_param, nrow(full_vcf)) else nrow(full_vcf)
  stopifnot(nrow(vcf_subset) == expected_loci) # Ensure it has the correct number of loci.
  message(paste(" -> VCF successfully subsetted. Dimensions:", nrow(vcf_subset), "loci by", ncol(vcf_subset@gt) - 1, "inds."))
  
  
  n_ind <- ncol(vcf_subset@gt) - 1
  n_loc <- nrow(vcf_subset)
  ind_names <- colnames(vcf_subset@gt)[-1]
  project_name <- paste0("Goby_", n_loc, "Loci_", params$analysis_method, "_", params$run_length)
  output_filename <- file.path(run_dir, paste0(project_name, ".dat"))
  
  message(paste(" -> Calculated variables: n_ind =", n_ind, "| n_loc =", n_loc))
  # --- LOGIC CHECK ---
  stopifnot(n_ind > 0 & n_loc > 0)
  
  ### --- 4.2: Convert User Inputs to COLONY Integer Codes ---
  method_code <- switch(toupper(params$analysis_method), "PLS" = 0, "FL" = 1, "FPLS" = 2, 1)
  length_code <- switch(tolower(params$run_length), "short" = 1, "medium" = 2, "long" = 3, "verylong" = 4, 2)
  precision_code <- switch(tolower(params$likelihood_precision), "low" = 0, "medium" = 1, "high" = 2, "veryhigh" = 3, 1)
  random_seed <- sample(1:9999, 1)
  
  ### --- 4.3: Prepare all file content in R memory first ---
  message(" -> Preparing file content in memory...")
  
  # A. Create the header content
  header_content <- c(
    paste0(project_name, "\t! Project name"), paste0(project_name, "\t! Output file name"),
    paste(n_ind, "\t! Number of offspring in the sample"), paste(n_loc, "\t! Number of loci"),
    paste(random_seed, "\t! Seed for random number generator"), "0\t! 0/1=Not updating/updating allele frequency",
    "2\t! 2/1=Dioecious/Monoecious species", "0\t! 0/1=no inbreeding/inbreeding",
    "0\t! 0/1=Diploid species/HaploDiploid species", "0 0\t! 0/1=Polygamy/Monogamy for males & females",
    "0\t! 0/1=Clone inference =No/Yes", "1\t! 0/1=Scale full sibship=No/Yes",
    "0\t! 0/1/2/3=No/Weak/Medium/Strong sibship prior", "0\t! 0/1=Unknown/Known population allele frequency",
    "1\t! Number of runs", paste(length_code, "\t! 1/2/3/4=short/medium/long/very long run"),
    "0\t! 0/1=Monitor method by Iterate#/Time in second", "10000\t! Monitor interval",
    "0\t! 0/1=No/Yes for run with Windows GUI", paste(method_code, "\t! 0/1/2=PairLikelihood score/Fulllikelihood/FPLS"),
    paste(precision_code, "\t! 0/1/2/3=Low/Medium/High/Very high precision")
  )
  
  # B. Create the marker data content
  marker_content <- c(
    paste0(paste0("Locus_", 1:n_loc), collapse = "\t"),
    paste(rep("0", n_loc), collapse = "\t"), # Marker Type
    paste(rep(params$dropout_rate, n_loc), collapse = "\t"),
    paste(rep(params$error_rate, n_loc), collapse = "\t")
  )
  
  # C. Format the genotype data
  message(" -> Formatting genotype matrix...")
  gt_matrix <- extract.gt(vcf_subset, element = "GT")
  # --- LOGIC CHECK ---
  stopifnot(nrow(gt_matrix) == n_loc, ncol(gt_matrix) == n_ind)
  message(paste("  -> Extracted genotype matrix with dimensions:", nrow(gt_matrix), "x", ncol(gt_matrix)))
  
  gt_matrix[gt_matrix == "0/0"] <- "1 1"; gt_matrix[gt_matrix == "0/1"] <- "1 2"
  gt_matrix[gt_matrix == "1/1"] <- "2 2"; gt_matrix[is.na(gt_matrix) | gt_matrix == "./."] <- "0 0"
  genotype_content <- apply(cbind(ind_names, gt_matrix), 1, paste, collapse = " ")
  # --- LOGIC CHECK ---
  stopifnot(length(genotype_content) == n_ind)
  message("  -> Genotype formatting complete.")
  
  
  # D. Create the footer content (all defaults)
  footer_content <- c("0.5 0.5\t! Probabilities father/mother in candidates", "0 0\t! Numbers of candidate males and females",
                      "0 0\t! Number of known paternity, exclusion threshold", "0 0\t! Number of known maternity, exclusion threshold",
                      "0\t! Number of known paternal sibships", "0\t! Number of known maternal sibships",
                      "0\t! Number of offspring with known excluded paternity", "0\t! Number of offspring with known excluded maternity",
                      "0\t! Number of offspring with known excluded paternal sibships", "0\t! Number of offspring with known excluded maternal sibships")
  
  # E. Combine all content into a single character vector
  full_file_content <- c(header_content, marker_content, genotype_content, footer_content)
  message(paste(" -> Final file content assembled with", length(full_file_content), "lines."))
  
  ### --- 4.4: Write the entire file in a single operation ---
  message(paste(" -> Writing COLONY input file to:", basename(output_filename)))
  writeLines(full_file_content, con = output_filename)
  message(" -> Write operation successful.")
  
  # ============================================================================
  # <<< END: INLINED FILE GENERATION LOGIC >>>
  # ============================================================================
  
  
  ### --- 4.5: Generate the Cluster Execution Script ---
  shell_script_path <- file.path(run_dir, "run_colony.sh")
  message(paste(" -> Generating execution script:", basename(shell_script_path)))
  
  # Define the content of the shell script
  shell_content <- c(
    "#!/bin/bash", "#", "# --- Example SLURM Directives ---",
    "# This is a template. You MUST modify these for your cluster's scheduler.",
    paste0("#SBATCH --job-name=colony_", locus_label), "#SBATCH --time=48:00:00",
    "#SBATCH --nodes=1", "#SBATCH --ntasks-per-node=1", "#SBATCH --mem=4G",
    "#SBATCH --output=colony_run.log", "# --------------------------------", "",
    "echo \"Starting COLONY run at $(date)\"", "",
    "# --- IMPORTANT: UPDATE THIS PATH ---",
    "PATH_TO_COLONY_EXECUTABLE=\"/path/on/server/to/colony2s.gfortran.x\"", "",
    "# COLONY requires the input file to be named colony2.dat",
    "# We rename our generated file to meet this requirement.",
    paste("mv", basename(output_filename), "colony2.dat"), "",
    "# Execute the COLONY program", "# The executable should be in your PATH or specified directly.",
    "$PATH_TO_COLONY_EXECUTABLE", "", "echo \"COLONY run finished at $(date)\""
  )
  
  # Write the shell script to the run directory
  writeLines(shell_content, con = shell_script_path)
  # Make the script executable (for Linux/macOS systems)
  system(paste("chmod +x", shell_script_path))
  message(" -> Execution script created successfully.")
  
} # --- End of for loop ---

message("\n--- Batch generation complete! All files and scripts have been created. ---")
message("Next steps: Upload the 'Output' directory to your cluster and submit each 'run_colony.sh' script.")

