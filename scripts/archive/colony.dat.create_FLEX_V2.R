#' Create a Flexible COLONY .dat file (V2 - Robust Writing)
#'
#' This function creates a COLONY .dat file.
#' V2 uses a single, persistent file connection (file_conn)
#' for all write operations, which is more efficient and
#' prevents I/O buffering errors (which caused random newlines).
#'
#' @param kids A data frame of offspring genotypes. Sample_ID MUST
#'   be the first column.
#' @param ... (All other parameters are the same as V1)
#'
colony.dat.create_FLEX.R <- function(moms = NA, dads = NA, kids, markers,
                                     analysis_type = "parentage",
                                     output_filename = "colony_FLEX.dat",
                                     project_name = "COLONY_Run",
                                     update.alfs = 0,
                                     spp.type = 2, inbreeding = 0, ploidy = 0,
                                     fem.gamy = 0, # Defaulting to Polygamy
                                     mal.gamy = 0, # Defaulting to Polygamy
                                     clone = 0, sib.scale = 1, sib.prior = 0,
                                     known.alfs = 0, run.number = 1,
                                     run.length = 2, monitor = 0, windows.version = 0,
                                     full.likelihood = 1, likelihood.precision = 2,
                                     prob.mom = 1.0, prob.dad = 1.0) {
  
  message(paste("Opening file connection for:", output_filename))
  # Open the file connection once (in "write" mode, which truncates)
  file_conn <- file(output_filename, "w")
  
  # Ensure the connection is closed when the function exits,
  # even if it errors out.
  on.exit(close(file_conn))
  
  my.wd <- paste0("'", project_name, "_Input'")
  my.wd2 <- paste0("'", project_name, "_Output'")
  
  # --- 1. VALIDATE AND SET COUNTS BASED ON ANALYSIS TYPE ---
  analysis_type <- tolower(analysis_type)
  
  if (is.null(kids) || !is.data.frame(kids) || nrow(kids) == 0) {
    stop("Offspring/Individuals genotypes ('kids') are required for all analysis types.")
  }
  noff <- nrow(kids) # Number of offspring/individuals
  
  if (analysis_type == "sibship") {
    message("Mode: 'sibship'. Ignoring 'moms' and 'dads' arguments.")
    nmoms <- 0
    ndads <- 0
    moms <- NULL # Force NULL to prevent writing
    dads <- NULL # Force NULL to prevent writing
    
  } else if (analysis_type == "parentage") {
    message("Mode: 'parentage'. Checking for candidate 'moms' and 'dads' data frames.")
    
    # Check Moms
    if (is.data.frame(moms) && nrow(moms) > 0) {
      nmoms <- nrow(moms)
      message(paste("Found", nmoms, "candidate mothers."))
    } else {
      message("No valid 'moms' data frame provided. Setting candidate mothers = 0.")
      nmoms <- 0
      moms <- NULL # Ensure it's NULL
    }
    
    # Check Dads
    if (is.data.frame(dads) && nrow(dads) > 0) {
      ndads <- nrow(dads)
      message(paste("Found", ndads, "candidate fathers."))
    } else {
      message("No valid 'dads' data frame provided. Setting candidate fathers = 0.")
      ndads <- 0
      dads <- NULL # Ensure it's NULL
    }
    
    if (nmoms == 0 && ndads == 0) {
      warning("Analysis type is 'parentage' but no candidate parents were found. Did you mean 'sibship'?")
    }
    
  } else {
    stop(paste("Unknown 'analysis_type':", analysis_type, ". Must be 'sibship' or 'parentage'."))
  }
  
  # --- 2. FORMAT PARAMETERS FOR FILE ---
  noff_str <- paste(noff, "! Number of offspring in the sample", sep = "\t")
  
  # Get number of loci
  nloci <- ncol(markers)
  nloci_str <- paste(nloci, "! Number of loci", sep = "\t")
  
  # Set random number seed
  random.seed <- round(runif(n = 1, min = 1, max = 9999))
  random.seed_str <- paste(random.seed, "! Seed for random number generator", sep = "\t")
  
  # Add comments to input parameters
  update.alfs_str <- paste(update.alfs, "! 0/1=Not update/Update allele frequency", sep = "\t")
  spp.type_str <- paste(spp.type, "! 2/1=Dioecious/Monoecious", sep = "\t")
  inbreeding_str <- paste(inbreeding, "! 0/1=No inbreeding/Inbreeding", sep = "\t")
  ploidy_str <- paste(ploidy, "! 0/1=Diploid species/HaploDiploid species", sep = "\t")
  gamy_str <- paste(paste(mal.gamy, fem.gamy), "! 0/1=Polygamy/Monogamy for males & females", sep = "\t")
  clone_str <- paste(clone, "! 0/1=Clone inference =No/Yes", sep = "\t")
  sib.scale_str <- paste(sib.scale, "! 0/1=Scale full sibship=No/Yes", sep = "\t")
  sib.prior_str <- paste(paste(sib.prior, 1, 1, sep = " "), "! 0/1/2/3=No sibship prior/Weak/Medium/Strong", sep = "\t")
  known.alfs_str <- paste(known.alfs, "! 0/1=Unknown/Known population allele frequency", sep = "\t")
  run.number_str <- paste(run.number, "! Number of runs", sep = "\t")
  run.length_str <- paste(run.length, "! Length of run", sep = "\t")
  monitor_str <- paste(monitor, "! 0/1=Monitor method by Iterate#/Time in second", sep = "\t")
  monitor.interval_str <- paste(1000000, "! Monitor interval in Iterate# / in seconds", sep = "\t")
  windows.version_str <- paste(windows.version, "! Windows version", sep = "\t")
  full.likelihood_str <- paste(full.likelihood, "! 0/1/2 = PL,FL,FPL", sep = "\t")
  likelihood.precision_str <- paste(likelihood.precision, "! 1/2/3=low/medium/high Precision for Fulllikelihood", sep = "\t")
  
  # Collate info for parents
  if (prob.dad == 0 & prob.mom == 0) {
    probs_str <- "0.01 0.01" # Use a small non-zero probability
  } else {
    probs_str <- as.character(paste(prob.dad, prob.mom))
  }
  
  # This is now SAFE because nmoms/ndads are guaranteed numbers
  npars_str <- as.character(paste(ndads, nmoms))
  
  my.value <- "0" # Just the string, writeLines will add the \n
  
  # --- 3. MAKING THE ACTUAL FILE (using the single connection) ---
  message("Writing parameters to connection...")
  
  # Collate all parameter strings
  param_lines <- c(
    my.wd, my.wd2, noff_str, nloci_str, random.seed_str, update.alfs_str, spp.type_str,
    inbreeding_str, ploidy_str, gamy_str, clone_str, sib.scale_str, sib.prior_str,
    known.alfs_str, run.number_str, run.length_str, monitor_str, monitor.interval_str,
    windows.version_str, full.likelihood_str, likelihood.precision_str,
    "" # Extra newline
  )
  
  # Write all parameter lines at once
  writeLines(param_lines, con = file_conn)
  
  # --- Write Marker Block ---
  # write.table handles its own newlines
  write.table(
    x = markers, file = file_conn, append = TRUE,
    quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE # Use space separator
  )
  
  # --- Write Offspring (All Individuals) Block ---
  writeLines("\n", con = file_conn) # Add a blank line
  
  # This is the critical write step.
  # kids df now has Sample_ID as column 1, which is correct.
  # col.names=FALSE and row.names=FALSE ensures it's written as:
  # [Sample_ID] [Locus1_A] [Locus1_B] ...
  write.table(
    x = kids, file = file_conn, append = TRUE,
    quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE
  )
  
  # --- Write Parent Blocks (Data-dependent) ---
  writeLines(c(
    "\n", # Blank line
    probs_str,
    npars_str,
    "\n"  # Blank line
  ), con = file_conn)
  
  # --- This is the new, ROBUST logic for writing parents ---
  if (ndads > 0 && !is.null(dads)) {
    message("Writing candidate fathers to file...")
    write.table(x = dads, file = file_conn, append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  }
  
  if (nmoms > 0 && !is.null(moms)) {
    message("Writing candidate mothers to file...")
    write.table(x = moms, file = file_conn, append = TRUE, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
  }
  
  # --- Write Final Footer Lines ---
  # Collate all footer lines
  footer_lines <- c(
    "\n",
    "0 0",
    "\n",
    "0 0",
    "\n",
    my.value,
    "\n",
    my.value,
    "\n",
    my.value,
    "\n",
    my.value,
    "\n",
    my.value,
    "\n",
    my.value
  )
  
  # Write all footer lines at once
  writeLines(footer_lines, con = file_conn)
  
  message(paste("--- File", output_filename, "written successfully. ---"))
  # The on.exit() command will automatically close file_conn now
}
