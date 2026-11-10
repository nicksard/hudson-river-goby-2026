# R Script to Analyze process_radtags.log and Barcode Distances
#
# This script reads the log file from the Stacks 'process_radtags' tool,
# parses read counts, quantifies specific barcode patterns (e.g., GG, GN, NG prefixes),
# and compares unassigned barcodes to an official list using string distance.

# --- 0. Package Installation ---
# Install and load necessary packages
# If you don't have these packages, uncomment the lines below to install them
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")
# install.packages("stringdist")

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringdist)

#setting working directory
setwd("G:/Shared drives/Sard Lab/Round Goby/RADseq/R/process.rad.tags.analysis")

# --- 1. Configuration ---
log_file_path <- "Input/process_radtags.log" # Make sure this file is in your R working directory
barcode_file_path <- "Input/Plate1Barcodes.txt" # IMPORTANT: Add your barcode file name here

# --- 2. Read and Parse Log File (Per-Sample Counts) ---

# Check if the log file exists
if (!file.exists(log_file_path)) {
  stop("Error: The file '", log_file_path, "' was not found in the working directory.")
}

log_lines <- readLines(log_file_path)

# Find the start and end of the 'per_barcode_raw_read_counts' section
start_index_counts <- grep("BEGIN per_barcode_raw_read_counts", log_lines)
end_index_counts <- grep("END per_barcode_raw_read_counts", log_lines)

rad_data <- data.frame()
if (length(start_index_counts) > 0 && length(end_index_counts) > 0) {
  table_lines <- log_lines[(start_index_counts + 2):(end_index_counts - 1)]
  data_list <- lapply(table_lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    if(length(parts) < 8) parts <- c(parts, rep(NA, 8 - length(parts)))
    return(parts)
  })
  rad_data <- as.data.frame(do.call(rbind, data_list), stringsAsFactors = FALSE)
  colnames(rad_data) <- c("Barcode", "Filename", "Retained", "LowQuality", "RAD_Cutsite_Not_Found", "Properly_Paired", "Retained_Percentage", "Total")
  
  # Clean and process data
  cols_to_convert <- c("Retained", "LowQuality", "RAD_Cutsite_Not_Found", "Properly_Paired", "Total")
  rad_data[cols_to_convert] <- lapply(rad_data[cols_to_convert], function(x) as.numeric(gsub(",", "", x)))
  rad_data <- rad_data %>% mutate(Discarded = Total - Retained)
  
  print("Successfully parsed per-sample read counts. Head of the data:")
  print(head(rad_data))
} else {
  warning("Could not find the 'per_barcode_raw_read_counts' section in the log file.")
}

# --- 3. Visualization of Read Counts ---
if (nrow(rad_data) > 0) {
  plot_data <- rad_data %>%
    select(Barcode, Retained, Discarded) %>%
    pivot_longer(cols = c(Retained, Discarded), names_to = "Read_Status", values_to = "Count")
  
  read_plot <- ggplot(plot_data, aes(x = reorder(Barcode, -Count), y = Count, fill = Read_Status)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Retained" = "#4CAF50", "Discarded" = "#F44336")) +
    labs(title = "Retained vs. Discarded Reads per Sample", x = "Sample Barcode", y = "Number of Reads") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave("radtags_read_counts.png", plot = read_plot, width = 12, height = 7, dpi = 300)
  cat("\n--- Visualization ---\n")
  cat("A plot named 'radtags_read_counts.png' has been saved.\n")
}

# --- 4. Barcode Analysis ---

# Check if barcode file exists
if (!file.exists(barcode_file_path)) {
  stop("Error: Barcode file '", barcode_file_path, "' not found. Please check the `barcode_file_path` variable.")
}

# Read the official barcodes
official_barcodes_df <- read.table(barcode_file_path, sep = "\t", stringsAsFactors = FALSE)
official_barcodes <- official_barcodes_df[[1]]

# Find and parse the 'barcodes_not_recorded' section
start_index_unassigned <- grep("BEGIN barcodes_not_recorded", log_lines)
end_index_unassigned <- grep("END barcodes_not_recorded", log_lines)

if (length(start_index_unassigned) > 0 && length(end_index_unassigned) > 0) {
  unassigned_lines <- log_lines[(start_index_unassigned + 1):(end_index_unassigned - 1)]
  unassigned_df <- read.table(text = unassigned_lines, sep = "\t", col.names = c("Unassigned_Barcode", "Count"))
  
  # --- 4a. Barcode Pattern Analysis ---
  cat("\n--- Barcode Pattern Analysis ---\n")
  
  # Quantify reads for specific patterns in unassigned barcodes
  gg_prefix_reads <- unassigned_df %>% filter(startsWith(Unassigned_Barcode, "GG")) %>% summarise(total = sum(Count)) %>% pull(total)
  gn_prefix_reads <- unassigned_df %>% filter(startsWith(Unassigned_Barcode, "GN")) %>% summarise(total = sum(Count)) %>% pull(total)
  ng_prefix_reads <- unassigned_df %>% filter(startsWith(Unassigned_Barcode, "NG")) %>% summarise(total = sum(Count)) %>% pull(total)
  
  cat(paste("Total reads in unassigned barcodes starting with 'GG':", format(gg_prefix_reads, big.mark = ","), "\n"))
  cat(paste("Total reads in unassigned barcodes starting with 'GN':", format(gn_prefix_reads, big.mark = ","), "(Possible error in 2nd position)\n"))
  cat(paste("Total reads in unassigned barcodes starting with 'NG':", format(ng_prefix_reads, big.mark = ","), "(Possible error in 1st position)\n"))
  
  # --- 4b. Barcode Distance Analysis ('MaxDiff') ---
  # Calculate the closest match for each unassigned barcode
  results <- lapply(1:nrow(unassigned_df), function(i) {
    unassigned_bc <- unassigned_df$Unassigned_Barcode[i]
    count <- unassigned_df$Count[i]
    
    # Calculate Levenshtein distance to all official barcodes
    distances <- stringdist(unassigned_bc, official_barcodes, method = "lv")
    
    # Find the minimum distance and the index of the best match
    min_dist <- min(distances)
    best_match_index <- which.min(distances)
    closest_official_bc <- official_barcodes[best_match_index]
    
    return(data.frame(
      Unassigned_Barcode = unassigned_bc,
      Count = count,
      Closest_Correct_Barcode = closest_official_bc,
      Edit_Distance = min_dist
    ))
  })
  
  distance_analysis_results <- do.call(rbind, results)
  
  cat("\n--- Barcode Distance Analysis ---\n")
  print("Top unassigned barcodes that are 1-2 edits away from an official barcode:")
  
  # Filter for barcodes that are likely sequencing errors (distance 1 or 2) and print
  potential_errors <- distance_analysis_results %>%
    filter(Edit_Distance > 0, Edit_Distance <= 2) %>%
    arrange(Edit_Distance, desc(Count))
  
  print(head(potential_errors, 20)) # Print the top 20 potential errors
  
} else {
  warning("Could not find the 'barcodes_not_recorded' section in the log file.")
}

