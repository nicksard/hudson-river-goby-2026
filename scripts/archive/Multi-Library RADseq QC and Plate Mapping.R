# ==============================================================================
# RADseq Pipeline Bias and Quality Control Analysis Script (V2)
# Developed for the NYSDEC Round Goby Hudson River Expansion Study
#
# This script processes multiple sequencing libraries (Plates 1-4) by:
#   1. Parsing demultiplexing results from Stacks 'process_radtags.log' files.
#   2. Mapping inline barcode counts to a 96-well plate layout.
#   3. Plotting plate heatmaps (individual and unified facet comparisons).
#   4. Plotting comparative boxplots of reads per library.
#   5. Compiling library-level summary statistics.
# ==============================================================================

# --- 0. Package Setup ---
# Load required libraries (installing first if missing)
required_packages <- c("ggplot2", "dplyr", "tidyr", "scales")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# --- 1. Working Directory & Configuration ---
# Update this path to your local working directory if needed.
# If commented out, R will use your active project directory.
# setwd("G:/Shared drives/Sard Lab/Round Goby/RADseq/R/process.rad.tags.analysis")

# Define plate layout filling direction:
# "row-wise"    = A1, A2, A3... A12, B1, B2... (Standard for consecutive sample lists)
# "column-wise" = A1, B1, C1... H1, A2, B2...
PLATE_LAYOUT_DIRECTION <- "row-wise" 

# Define input files for the 4 libraries
libraries_config <- data.frame(
  Lib_ID = 1:4,
  Lib_Name = c("Library 1 (Plate 1)", "Library 2 (Plate 2)", "Library 3 (Plate 3)", "Library 4 (Plate 4)"),
  Log_File = c("Input/lib1_process_radtags.log", 
               "Input/lib2_process_radtags.log", 
               "Input/lib3_process_radtags.log", 
               "Input/lib4_process_radtags.log"),
  Barcode_File = c("Input/Plate_1_goby.barcodes.txt", 
                   "Input/Plate_2_goby.barcodes.txt", 
                   "Input/Plate_3_goby.barcodes.txt", 
                   "Input/Plate_4_goby.barcodes.txt"),
  stringsAsFactors = FALSE
)

# Create an output directory for plots and tables if it doesn't exist
dir.create("Output", showWarnings = FALSE)


# --- 2. Run Pipeline Loop ---

all_libraries_data <- list()

cat("=========================================\n")
cat("Starting RADseq Multi-Library Analysis Pipeline\n")
cat("=========================================\n\n")

i <- 1
for (i in 1:nrow(libraries_config)) {
  lib_name <- libraries_config$Lib_Name[i]
  log_file <- libraries_config$Log_File[i]
  barcode_file <- libraries_config$Barcode_File[i]
  
  cat(sprintf("Processing %s...\n", lib_name))
  
  # Check file existence before proceeding
  if (!file.exists(log_file)) stop(paste("Error: Log file not found:", log_file))
  if (!file.exists(barcode_file)) stop(paste("Error: Barcodes file not found:", barcode_file))
  
  # ==========================================
  # Section A: Load Barcodes & Identify Coordinates
  # ==========================================
  barcodes_df <- read.table(barcode_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(barcodes_df) <- c("Barcode", "Sample_ID")
  barcodes_df$Barcode <- trimws(barcodes_df$Barcode)
  
  # Assign 96-well plate coordinates based on layout direction
  n_wells <- nrow(barcodes_df)
  rows_letters <- LETTERS[1:8]
  cols_numbers <- 1:12
  
  if (PLATE_LAYOUT_DIRECTION == "row-wise") {
    barcodes_df$Row <- rep(rows_letters, each = 12, length.out = n_wells)
    barcodes_df$Column <- rep(cols_numbers, times = 8, length.out = n_wells)
  } else {
    barcodes_df$Row <- rep(rows_letters, times = 12, length.out = n_wells)
    barcodes_df$Column <- rep(cols_numbers, each = 8, length.out = n_wells)
  }
  barcodes_df$Well <- paste0(barcodes_df$Row, barcodes_df$Column)
  barcodes_df
  
  # ==========================================
  # Section B: Parse the process_radtags Log
  # ==========================================
  log_lines <- readLines(log_file, warn = FALSE)
  
  # Locate the start of the barcode counts table
  header_idx <- which(grepl("^BEGIN per_barcode_raw_read_counts", log_lines))
  table_end_idx <- which(grepl("^END per_barcode_raw_read_counts", log_lines))
  
  header_line <- log_lines[header_idx+1]
  headers <- unlist(strsplit(header_line, "\t"))
  
  # Extract data lines belonging to the barcode table
  data_lines <- log_lines[(header_idx+2):(table_end_idx-1)]
  head(data_lines)
  tail(data_lines)
  
  # Convert parsed text to structured dataframe
  counts_df <- read.table(text = paste(data_lines, collapse = "\n"), 
                          header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(counts_df) <- headers
  counts_df$Barcode <- trimws(counts_df$Barcode)
  head(counts_df)
  head(barcodes_df)
  
  # ==========================================
  # Section C: Merge and Clean
  # ==========================================
  plate_data <- left_join(barcodes_df, counts_df, by = "Barcode")
  head(plate_data)
  
  # Clean up NA values and enforce types
  plate_data[is.na(plate_data)] <- 0
  plate_data$Retained <- as.numeric(plate_data$Retained)
  plate_data$Total <- as.numeric(plate_data$Total)
  plate_data$Library <- lib_name
  
  # Store to list for binding later
  all_libraries_data[[i]] <- plate_data
}

# Bind lists into master dataframes
master_plate_df <- do.call(rbind, all_libraries_data)
head(master_plate_df)

# Export processed spatial count table
write.csv(master_plate_df, "Output/master_parsed_plate_counts.csv", row.names = FALSE)

# --- 3. Summary Statistics Table for R1 (Retained Reads) ---

summary_stats_table <- master_plate_df %>%
  filter(!grepl(pattern = "NTC|BLA|blank",x = Sample_ID)) %>%
  group_by(Library) %>%
  summarise(
    Min = min(Retained, na.rm = TRUE),
    Max = max(Retained, na.rm = TRUE),
    Median = median(Retained, na.rm = TRUE),
    Mean = mean(Retained, na.rm = TRUE),
    SD = sd(Retained, na.rm = TRUE),
    Total_Retained_Reads = sum(Retained, na.rm = TRUE),
    Active_Wells = sum(Retained > 0, na.rm = TRUE)
  )
summary_stats_table

cat("========================================================================\n")
cat("SUMMARY STATISTICS FOR R1 (RETAINED READS) PER LIBRARY\n")
cat("========================================================================\n")
print(as.data.frame(summary_stats_table))
cat("========================================================================\n\n")

# Save table to CSV
write.csv(summary_stats_table, "Output/libraries_summary_statistics.csv", row.names = FALSE)

summary_stats_table_negs <- master_plate_df %>%
  filter(grepl(pattern = "NTC|BLA|blank",x = Sample_ID)) %>%
  group_by(Library) %>%
  summarise(
    Min = min(Retained, na.rm = TRUE),
    Max = max(Retained, na.rm = TRUE),
    Median = median(Retained, na.rm = TRUE),
    Mean = mean(Retained, na.rm = TRUE),
    SD = sd(Retained, na.rm = TRUE),
    Total_Retained_Reads = sum(Retained, na.rm = TRUE),
    Active_Wells = sum(Retained > 0, na.rm = TRUE)
  )
summary_stats_table_negs

cat("========================================================================\n")
cat("SUMMARY STATISTICS FOR R1 (RETAINED READS) NEGs PER LIBRARY\n")
cat("========================================================================\n")
print(as.data.frame(summary_stats_table_negs))
cat("========================================================================\n\n")

# Save table to CSV
write.csv(summary_stats_table, "Output/libraries_summary_statistics_negs.csv", row.names = FALSE)

# --- 4. Data Visualizations ---

# --- 4a. Plot 1: Boxplot Distribution of Reads per Library ---
boxplot_plot <- ggplot(master_plate_df, aes(x = factor(Library), y = Retained, fill = factor(Library))) +
  geom_boxplot(alpha = 0.65, outlier.shape = 21, outlier.size = 2.5, outlier.fill = "white", color = "#2c3e50", size = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.2, color = "#34495e") +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Distribution of Retained Reads (R1) per Library",
    subtitle = "Analysis of well-to-well yield variations and outliers across all 4 plates",
    x = "Library Number",
    y = "Retained Reads per Well (R1)",
    fill = "Library ID"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "#2c3e50"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic", color = "#7f8c8d"),
    axis.title = element_text(face = "bold", color = "#2c3e50"),
    axis.text.x = element_text(face = "bold", size = 10, color = "#2c3e50"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
boxplot_plot

# Save boxplot in dual formats (PDF and JPEG)
ggsave("Output/libraries_retained_reads_boxplot.pdf", plot = boxplot_plot, width = 8, height = 6)
ggsave("Output/libraries_retained_reads_boxplot.jpg", plot = boxplot_plot, width = 8, height = 6, dpi = 300)

# --- 4b. Plot 2: Individual and Faceted Plate Heatmaps ---

# Set up spatial labeling classes and assign distinctive popping colors
master_plate_df <- master_plate_df %>%
  mutate(
    # Identify control well types
    Well_Type = case_when(
      grepl("NTC", Sample_ID, ignore.case = TRUE) ~ "N",      # Negative Control (NTC)
      grepl("BLA|blank", Sample_ID, ignore.case = TRUE) ~ "B", # Blank well
      TRUE ~ "S"                                               # Biological Sample
    ),
    # Setup heatmap text labels (No raw read counts inside the tiles as requested)
    Label_Text = case_when(
      Well_Type == "N" ~ "N",
      Well_Type == "B" ~ "B",
      TRUE ~ "S"
    )
  ) %>%
  group_by(Library) %>%
  ungroup()

# Plot individual heatmaps for each plate
for (i in 1:nrow(libraries_config)) {
  lib_name <- libraries_config$Lib_Name[i]
  single_plate_df <- master_plate_df %>% filter(Library == lib_name)
  
  single_heatmap <- ggplot(single_plate_df, aes(x = factor(Column), y = factor(Row, levels = rev(LETTERS[1:8])))) +
    geom_tile(aes(fill = Retained), color = "#ffffff", size = 0.4) +
    # Smooth Red-to-Yellow-to-Green continuous color ramp
    scale_fill_gradientn(colors = c("red", "yellow", "green"), labels = scales::comma_format()) +
    # Text labels (S + ID, N for negative, B for blank)
    geom_text(aes(label = Label_Text, color = Text_Color), 
              size = 1.8, fontface = "bold") +
    scale_color_identity() +
    labs(
      title = paste("Plate Layout Heatmap:", lib_name),
      subtitle = paste("R1 reads mapped across 96-well grid (Controls: N = Orange Negatives, B = Black Blanks)"),
      x = "Plate Column",
      y = "Plate Row",
      fill = "Retained Reads"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, color = "#2c3e50"),
      plot.subtitle = element_text(size = 9, hjust = 0.5, face = "italic", color = "#7f8c8d"),
      axis.title = element_text(face = "bold", color = "#2c3e50"),
      axis.text = element_text(face = "bold", size = 10, color = "#2c3e50"),
      panel.grid = element_blank(),
      aspect.ratio = 8/12,
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 9)
    )
  
  # Save individual heatmaps in dual formats (PDF and JPEG)
  ggsave(sprintf("Output/lib%d_plate_heatmap.pdf", i), plot = single_heatmap, width = 11, height = 8)
  ggsave(sprintf("Output/lib%d_plate_heatmap.jpg", i), plot = single_heatmap, width = 11, height = 8, dpi = 300)
}

head(master_plate_df)
master_plate_df$Retained2 <- factor(x = as.character(str_length(master_plate_df$Retained)),
                                    levels = as.character(1:7))
levels(master_plate_df$Retained2)

# --- 4c. Plot 3: Faceted Multi-Plate Heatmap Comparison (All Plates in One Plot) ---
faceted_heatmap <- ggplot(master_plate_df, aes(x = factor(Column),
                                               y = factor(Row, levels = rev(LETTERS[1:8])))) +
  geom_tile(aes(fill = Retained2), color = "#ffffff", size = 0.2) +
  # Smooth Red-to-Yellow-to-Green continuous color ramp
  scale_fill_manual(values = c(
    "1" = "#990000",
    "2" = "#d7301f",
    "3" = "#ef6548",
    "4" = "#fc8d59",
    "5" = "#99d8c9",
    "6" = "#41ae76",
    "7" = "#238b45"
  )) +
  #scale_fill_gradient2(low ="#990000",mid = "#99d8c9",high =  "#238b45")
  # Simplified text labeling for faceted comparison
  geom_text(aes(label = Label_Text), 
            size = 1.5, fontface = "bold") +
  facet_wrap(~ Library, ncol = 2) +
  labs(
    title = "Comparison of RADseq Plate Layout Yields",
    subtitle = "Faceted comparison of R1 reads across 4 plates (Controls: N = Orange Negatives, B = Black Blanks)",
    x = "Plate Column",
    y = "Plate Row",
    fill = "Retained Reads"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#2c3e50"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, face = "italic", color = "#7f8c8d"),
    axis.title = element_text(face = "bold", size = 11, color = "#2c3e50"),
    axis.text = element_text(face = "bold", size = 8, color = "#2c3e50"),
    strip.text = element_text(face = "bold", size = 12, color = "#2c3e50"),
    panel.grid = element_blank(),
    aspect.ratio = 8/12,
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10)
  )
faceted_heatmap

# Save faceted heatmap in dual formats (PDF and JPEG)
ggsave("Output/faceted_all_plates_heatmaps.pdf", plot = faceted_heatmap, width = 14, height = 10)
ggsave("Output/faceted_all_plates_heatmaps.jpg", plot = faceted_heatmap, width = 14, height = 10, dpi = 300)

cat("=========================================\n")
cat("Pipeline Finished Successfully! All output files saved to current directory under 'Output/'.\n")
cat("=========================================\n")
