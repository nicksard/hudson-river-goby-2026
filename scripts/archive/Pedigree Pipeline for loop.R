# ==============================================================================
# Master Generator for COLONY Pedigree Analysis on HPC (v2.0)
#
# AUTHOR: Nicholas Sard, Ph.D.
# ASSISTANT: Gemini
# DATE: October 16, 2025
#
# DESCRIPTION:
# This single, self-contained script automates the generation of input files
# and corresponding shell scripts for a large-scale COLONY analysis.
#
# The core logic for file creation is placed directly within the main loop
# to facilitate debugging with print statements. It is designed for use with
# a High-Performance Computing (HPC) cluster by generating a separate
# 'run_colony.sh' script for each analysis, which can then be submitted to a
# job scheduler (e.g., Slurm, SGE) for parallel execution.
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
if (!file.exists(VCF_FILE_PATH)) {
  stop("VCF file not found. Please check the 'VCF_FILE_PATH' variable.")
}


# ==============================================================================
# Part 2: Define the Analysis Plan
# ==============================================================================
# This is the control center for your entire analysis.
# Each row represents a unique COLONY run to be generated.
message("--- Defining the analysis plan (6 runs) ---")
analysis_plan <- tibble::tribble(
  ~locus_subset, ~analysis_method, ~run_length, ~likelihood_precision, ~error_rate, ~dropout_rate,
  "250",         "FL",             "medium",    "medium",                0.001,       0.001,
  "2000",        "FL",             "medium",    "medium",                0.001,       0.001,
  "4000",        "FL",             "medium",    "medium",                0.001,       0.001,
  "5000",        "FL",             "medium",    "medium",                0.001,       0.001,
  "8000",        "FL",             "medium",    "medium",                0.001,       0.001,
  "all",         "FL",             "medium",    "medium",                0.001,       0.001
)


# ==============================================================================
# Part 3: Execute the Batch Generation Loop
# ==============================================================================
message(paste("\n--- Starting batch generation for", nrow(analysis_plan), "runs ---"))

# --- Pre-load the VCF file once to avoid repeated reading ---
message(paste("Pre-loading VCF from:", VCF_FILE_PATH))
full_vcf <- read.vcfR(VCF_FILE_PATH)
message("VCF loaded successfully.")

for (i in 1:nrow(analysis_plan)) {
  
  # --- Get parameters for the current run ---
  params <- analysis_plan[i, ]
  locus_param <- if (params$locus_subset == "all") "all" else as.numeric(params$locus_subset)
  locus_label <- if (is.numeric(locus_param)) paste0("Loci_", locus_param) else "Loci_All"
  
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
  
  ### --- 3.1: Data Preparation & Subsetting ---
  vcf_subset <- NULL # Initialize to be safe
  if (!is.numeric(locus_param) || locus_param >= nrow(full_vcf)) {
    message(paste(" -> Using all", nrow(full_vcf), "available loci."))
    vcf_subset <- full_vcf
  } else {
    message(paste(" -> Calculating MAF to select top", locus_param, "loci..."))
    maf_info <- vcfR::maf(full_vcf, element = 2)
    top_loci_ids <- names(sort(maf_info[,1], decreasing = TRUE)[1:locus_param])
    vcf_subset <- full_vcf[top_loci_ids, ]
  }
  message(paste(" -> VCF successfully subsetted to", nrow(vcf_subset), "loci."))
  # DEBUG: print(vcf_subset)
  
  n_ind <- ncol(vcf_subset@gt) - 1
  n_loc <- nrow(vcf_subset)
  ind_names <- colnames(vcf_subset@gt)[-1]
  project_name <- paste0("Goby_", n_loc, "Loci_", params$analysis_method, "_", params$run_length)
  output_filename <- file.path(run_dir, paste0(project_name, ".dat"))
  
  ### --- 3.2: Convert User Inputs to COLONY Integer Codes ---
  method_code <- switch(toupper(params$analysis_method), "PLS" = 0, "FL" = 1, "FPLS" = 2, 1)
  length_code <- switch(tolower(params$run_length), "short" = 1, "medium" = 2, "long" = 3, "verylong" = 4, 2)
  precision_code <- switch(tolower(params$likelihood_precision), "low" = 0, "medium" = 1, "high" = 2, "veryhigh" = 3, 1)
  random_seed <- sample(1:9999, 1) # Generate a new seed for each run
  
  ### --- 3.3: Prepare all file content in R memory first ---
  message(" -> Preparing file content in memory...")
  
  # A. Create the header content
  header_content <- c(
    paste0(project_name, "\t! Project name"),
    paste0(project_name, "\t! Output file name"),
    paste(n_ind, "\t! Number of offspring in the sample"),
    paste(n_loc, "\t! Number of loci"),
    paste(random_seed, "\t! Seed for random number generator"),
    "0\t! 0/1=Not updating/updating allele frequency",
    "2\t! 2/1=Dioecious/Monoecious species",
    "0\t! 0/1=no inbreeding/inbreeding",
    "0\t! 0/1=Diploid species/HaploDiploid species",
    "0 0\t! 0/1=Polygamy/Monogamy for males & females",
    "0\t! 0/1=Clone inference =No/Yes", "1\t! 0/1=Scale full sibship=No/Yes",
    "0\t! 0/1/2/3=No/Weak/Medium/Strong sibship prior",
    "0\t! 0/1=Unknown/Known population allele frequency", "1\t! Number of runs",
    paste(length_code, "\t! 1/2/3/4=short/medium/long/very long run"),
    "0\t! 0/1=Monitor method by Iterate#/Time in second", "10000\t! Monitor interval",
    "0\t! 0/1=No/Yes for run with Windows GUI",
    paste(method_code, "\t! 0/1/2=PairLikelihood score/Fulllikelihood/FPLS"),
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
  gt_matrix <- extract.gt(vcf_subset, element = "GT")
  gt_matrix[gt_matrix == "0/0"] <- "1 1"
  gt_matrix[gt_matrix == "0/1"] <- "1 2"
  gt_matrix[gt_matrix == "1/1"] <- "2 2"
  gt_matrix[is.na(gt_matrix) | gt_matrix == "./."] <- "0 0"
  genotype_content <- apply(cbind(ind_names, gt_matrix), 1, paste, collapse = " ")
  
  # D. Create the footer content (all defaults)
  footer_content <- c("0.5 0.5\t! Probabilities father/mother in candidates",
                      "0 0\t! Numbers of candidate males and females",
                      "0 0\t! Number of known paternity, exclusion threshold",
                      "0 0\t! Number of known maternity, exclusion threshold",
                      "0\t! Number of known paternal sibships",
                      "0\t! Number of known maternal sibships",
                      "0\t! Number of offspring with known excluded paternity",
                      "0\t! Number of offspring with known excluded maternity",
                      "0\t! Number of offspring with known excluded paternal sibships",
                      "0\t! Number of offspring with known excluded maternal sibships")
  
  # E. Combine all content into a single character vector
  full_file_content <- c(header_content, marker_content, genotype_content, footer_content)
  
  ### --- 3.4: Write the entire file in a single operation ---
  message(paste(" -> Writing all content to:", basename(output_filename)))
  writeLines(full_file_content, con = output_filename)
  
  # ============================================================================
  # <<< END: INLINED FILE GENERATION LOGIC >>>
  # ============================================================================
  
  
  ### --- 3.5: Generate the Cluster Execution Script ---
  shell_script_path <- file.path(run_dir, "run_colony.sh")
  message(paste(" -> Generating execution script:", basename(shell_script_path)))
  
  # Define the content of the shell script
  shell_content <- c(
    "#!/bin/bash",
    "#",
    "# --- Example SLURM Directives ---",
    "# This is a template. You MUST modify these for your cluster's scheduler.",
    paste0("#SBATCH --job-name=colony_", locus_label),
    "#SBATCH --time=48:00:00",
    "#SBATCH --nodes=1",
    "#SBATCH --ntasks-per-node=1",
    "#SBATCH --mem=4G",
    "#SBATCH --output=colony_run.log",
    "# --------------------------------",
    "",
    "echo \"Starting COLONY run at $(date)\"",
    "",
    "# --- IMPORTANT: UPDATE THIS PATH ---",
    "PATH_TO_COLONY_EXECUTABLE=\"/path/on/server/to/colony2s.gfortran.x\"",
    "",
    "# COLONY requires the input file to be named colony2.dat",
    "# We rename our generated file to meet this requirement.",
    paste("mv", basename(output_filename), "colony2.dat"),
    "",
    "# Execute the COLONY program",
    "# The executable should be in your PATH or specified directly.",
    "$PATH_TO_COLONY_EXECUTABLE",
    "",
    "echo \"COLONY run finished at $(date)\""
  )
  
  # Write the shell script to the run directory
  writeLines(shell_content, con = shell_script_path)
  # Make the script executable (for Linux/macOS systems)
  system(paste("chmod +x", shell_script_path))
  
} # --- End of for loop ---

message("\n--- Batch generation complete! All files and scripts have been created. ---")
message("Next steps: Upload the 'Output' directory to your cluster and submit each 'run_colony.sh' script.")

