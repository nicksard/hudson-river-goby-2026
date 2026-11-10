################################################################################
### R Pipeline for Pedigree Reconstruction Input File Generation (Version 2.3)
###
### Author: Nicholas Sard, Ph.D.
### Assistant: Gemini
### Date: October 15, 2025
###
### Description: This script takes a filtered VCF file as input and generates
### formatted files for COLONY and CKMRsim. This version disables scientific
### notation to ensure compatibility with external programs.
###
################################################################################


# ==============================================================================
# Part 1: Load Libraries & Set Global Options
# ==============================================================================
message("--- Loading required libraries ---")
library(vcfR)
library(tidyverse)

# --- Prevent scientific notation in output files ---
# This is critical for compatibility with programs like COLONY that expect
# plain numbers (e.g., 1000000 instead of 1e+06).
options(scipen = 999)


# ==============================================================================
# Part 2: Helper Function for Locus Subsetting
# ==============================================================================
.subset_vcf_by_maf <- function(vcf_path, locus_subset) {
  message(paste("Reading VCF from:", vcf_path))
  vcf <- read.vcfR(vcf_path)
  
  if (!is.numeric(locus_subset) || locus_subset >= nrow(vcf)) {
    message(paste("Using all", nrow(vcf), "available loci."))
    return(vcf)
  }
  
  message(paste("Calculating MAF to select the top", locus_subset, "most informative loci..."))
  maf_info <- vcfR::maf(vcf, element = 2)
  top_loci_ids <- names(sort(maf_info[,1], decreasing = TRUE)[1:locus_subset])
  vcf_subset <- vcf[top_loci_ids, ]
  message(paste("VCF successfully subsetted to", nrow(vcf_subset), "loci."))
  return(vcf_subset)
}


# ==============================================================================
# Part 3: COLONY Input File Generation
# ==============================================================================
create_colony_input <- function(vcf_path,
                                locus_subset,
                                analysis_method = "FL",
                                run_length = "medium",
                                likelihood_precision = "medium",
                                dropout_rate = 0.001,
                                error_rate = 0.001,
                                random_seed = NULL,
                                output_dir = "Output") {
  
  ### --- 3.1: Data Preparation ---
  vcf_subset <- .subset_vcf_by_maf(vcf_path, locus_subset)
  n_ind <- ncol(vcf_subset@gt) - 1
  n_loc <- nrow(vcf_subset)
  ind_names <- colnames(vcf_subset@gt)[-1]
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  project_name <- paste0("Goby_", n_loc, "Loci_", analysis_method, "_", run_length)
  output_filename <- file.path(output_dir, paste0(project_name, ".dat"))
  
  ### --- 3.2: Convert User Inputs to COLONY Integer Codes ---
  message(paste("Setting COLONY parameters: Method=", analysis_method,
                "| Length=", run_length, "| Precision=", likelihood_precision))
  method_code <- switch(toupper(analysis_method), "PLS" = 0, "FL" = 1, "FPLS" = 2, 1)
  length_code <- switch(tolower(run_length), "short" = 1, "medium" = 2, "long" = 3, "verylong" = 4, 2)
  precision_code <- switch(tolower(likelihood_precision), "low" = 0, "medium" = 1, "high" = 2, "veryhigh" = 3, 1)
  if (is.null(random_seed)) random_seed <- sample(1:9999, 1)
  
  ### --- 3.3: Prepare all file content in R memory first ---
  message("Preparing file content in memory...")
  
  # A. Create the header content as a vector of strings
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
    "0\t! 0/1=Clone inference =No/Yes",
    "1\t! 0/1=Scale full sibship=No/Yes",
    "0\t! 0/1/2/3=No/Weak/Medium/Strong sibship prior",
    "0\t! 0/1=Unknown/Known population allele frequency",
    "1\t! Number of runs",
    paste(length_code, "\t! 1/2/3/4=short/medium/long/very long run"),
    "0\t! 0/1=Monitor method by Iterate#/Time in second",
    "10000\t! Monitor interval in Iterate# / in seconds",
    "0\t! 0/1=No/Yes for run with Windows GUI",
    paste(method_code, "\t! 0/1/2=PairLikelihood score/Fulllikelihood/FPLS"),
    paste(precision_code, "\t! 0/1/2/3=Low/Medium/High/Very high precision with Fulllikelihood")
  )
  
  # B. Create the marker data content
  marker_content <- c(
    paste0(paste0("Locus_", 1:n_loc), collapse = ","),
    paste(rep("0", n_loc), collapse = ","),
    paste(rep(dropout_rate, n_loc), collapse = ","),
    paste(rep(error_rate, n_loc), collapse = ",")
  )
  
  # C. Format the genotype data
  gt_matrix <- extract.gt(vcf_subset, element = "GT")
  gt_matrix[gt_matrix == "0/0"] <- "1 1"
  gt_matrix[gt_matrix == "0/1"] <- "1 2"
  gt_matrix[gt_matrix == "1/1"] <- "2 2"
  gt_matrix[is.na(gt_matrix) | gt_matrix == "./."] <- "0 0"
  
  # Use paste to create a single string for each individual's row
  genotype_content <- apply(cbind(ind_names, gt_matrix), 1, paste, collapse = " ")
  
  # D. Create the footer content
  footer_content <- c(
    "0.5 0.5\t! Probabilities that the father and mother are included in candidates",
    "0 0\t! Numbers of candidate males and females",
    "0 0\t! Number of known paternity, exclusion threshold",
    "0 0\t! Number of known maternity, exclusion threshold",
    "0\t! Number of known paternal sibships",
    "0\t! Number of known maternal sibships",
    "0\t! Number of offspring with known excluded paternity",
    "0\t! Number of offspring with known excluded maternity",
    "0\t! Number of offspring with known excluded paternal sibships",
    "0\t! Number of offspring with known excluded maternal sibships"
  )
  
  # E. Combine all content into a single character vector
  full_file_content <- c(header_content, marker_content, genotype_content, footer_content)
  
  ### --- 3.4: Write the entire file in a single operation ---
  message(paste("Writing all content to:", output_filename))
  writeLines(full_file_content, con = output_filename)
  
  message("--- COLONY file generation complete! ---")
  return(invisible(TRUE))
}


# ==============================================================================
# Part 4: CKMRsim Input File Generation (Unchanged)
# ==============================================================================
create_ckmrsim_input <- function(vcf_path, locus_subset, output_dir = "Output") {
  
  vcf_subset <- .subset_vcf_by_maf(vcf_path, locus_subset)
  n_loc <- nrow(vcf_subset)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Allele Frequency File
  ref_freqs <- vcfR::maf(vcf_subset, element = 1)
  freq_df <- as_tibble(ref_freqs, rownames = "Locus") %>%
    select(Locus, Freq_Ref = `Frequency of 1st allele`)
  freq_filename <- file.path(output_dir, paste0("ckmrsim_allele_freqs_", n_loc, "Loci.txt"))
  write_tsv(freq_df, freq_filename)
  message(paste("CKMRsim allele frequency file written to:", freq_filename))
  
  # Genotypes File
  gt_vcf <- extract.gt(vcf_subset, element = "GT")
  gt_numeric <- matrix(NA, nrow = nrow(gt_vcf), ncol = ncol(gt_vcf), dimnames = dimnames(gt_vcf))
  gt_numeric[gt_vcf == "0/0"] <- 0
  gt_numeric[gt_vcf == "0/1"] <- 1
  gt_numeric[gt_vcf == "1/1"] <- 2
  genotypes_df <- as_tibble(t(gt_numeric), rownames = "Sample")
  geno_filename <- file.path(output_dir, paste0("ckmrsim_genotypes_", n_loc, "Loci.txt"))
  write_tsv(genotypes_df, geno_filename)
  message(paste("CKMRsim genotypes file written to:", geno_filename))
  
  message("--- CKMRsim file generation complete! ---")
  return(invisible(TRUE))
}

