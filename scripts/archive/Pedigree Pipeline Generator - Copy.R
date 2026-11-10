################################################################################
### Batch Runner for COLONY Pedigree Analysis
###
### Author: Nicholas Sard, Ph.D.
### Assistant: Gemini
### Date: October 15, 2025
###
### Description: This script automates the generation of COLONY input files
### by sourcing the "Pedigree Pipeline Generator" script and running it with
### different parameter combinations. It creates a structured directory system
### within the "Output" folder, making it easy to manage and upload multiple
### parallel runs to a high-performance computing cluster.
###
################################################################################


# ==============================================================================
# Part 1: Setup
# ==============================================================================

# --- Set the Working Directory ---
# The script assumes this R file is in the main project directory.
# Uncomment and set this to your specific path if needed.
setwd("G:/Shared drives/Sard Lab/Round Goby/RADseq/R/pedigree reconstruction")

# --- Source the Generator Script ---
# This loads the create_colony_input() and create_ckmrsim_input() functions.
# It assumes 'Pedigree Pipeline Generator.R' is in the same directory.
message("--- Loading the Pedigree Pipeline Generator script ---")
source("Pedigree Pipeline Generator.R")


# ==============================================================================
# Part 2: Define the Analysis Plan
# ==============================================================================
# This is the control center for your entire analysis.
# Create a data frame where each row represents a unique COLONY run.

# --- SHORT-TERM PLAN: Expanded locus subset analysis ---
# As requested, runs for 250, 2000, 4000, 5000, 8000, and all loci.
# All other parameters are held constant at medium/medium.
# NOTE: locus_subset is now character to handle the "all" value.
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


# --- LONG-TERM PLAN: Example of a more comprehensive analysis ---
# When you have more time, you can comment out the short-term plan above
# and uncomment this section to generate a full suite of analyses.
# This example would create 32 different input files (2 precisions * 2 lengths * 2 errors * 4 locus sets).
#
# message("--- Defining the comprehensive long-term analysis plan ---")
# analysis_plan <- expand.grid(
#   locus_subset = c("1000", "2000", "5000", "all"), # Must be character
#   analysis_method = "FL", # Sticking with Full Likelihood
#   run_length = c("medium", "long"),
#   likelihood_precision = c("medium", "high"),
#   error_rate = c(0.001, 0.005),
#   dropout_rate = 0.001 # Keeping dropout rate constant for this example
# ) %>% as_tibble()


# ==============================================================================
# Part 3: Execute the Batch Generation
# ==============================================================================
message(paste("\n--- Starting batch generation for", nrow(analysis_plan), "runs ---"))

for (i in 1:nrow(analysis_plan)) {
  
  # Get parameters for the current run from the plan
  params <- analysis_plan[i, ]
  
  # --- FIX for mixed data types ---
  # Convert the locus_subset back to a numeric type if it's not "all",
  # as the create_colony_input function expects either a number or "all".
  locus_param <- params$locus_subset
  if (locus_param != "all") {
    locus_param <- as.numeric(locus_param)
  }
  
  # Sanitize locus_subset for directory naming
  locus_label <- ifelse(is.numeric(locus_param),
                        paste0("Loci_", locus_param),
                        "Loci_All")
  
  # Define a structured output directory for this specific run
  run_dir <- file.path("Output",
                       paste0(params$analysis_method, "_", params$run_length, "_", params$likelihood_precision),
                       locus_label)
  
  # Create the directory if it doesn't already exist
  if (!dir.exists(run_dir)) {
    dir.create(run_dir, recursive = TRUE)
  }
  
  # Print a status message to the console
  message(paste("\n[", i, "/", nrow(analysis_plan), "] Generating files for:", locus_label,
                "| Method:", params$analysis_method,
                "| Length:", params$run_length))
  message(paste(" -> Output directory:", run_dir))
  
  # Call the function from the sourced script
  create_colony_input(
    vcf_path = "Input/populations.filtered.final.vcf.gz",
    locus_subset = locus_param, # Use the corrected parameter here
    analysis_method = params$analysis_method,
    run_length = params$run_length,
    likelihood_precision = params$likelihood_precision,
    dropout_rate = params$dropout_rate,
    error_rate = params$error_rate,
    output_dir = run_dir # IMPORTANT: Use the specific directory for this run
  )
}

message("\n--- Batch generation complete! All files have been created. ---")
