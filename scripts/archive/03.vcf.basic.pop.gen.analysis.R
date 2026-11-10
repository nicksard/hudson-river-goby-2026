################################################################################
### Population Genetics Analysis of Filtered VCF Data
###
### Author: Nicholas Sard, Ph.D.
### Assistant: Gemini
### Date: October 14, 2025
###
### Description: This script reads in a filtered VCF file and performs a
### standard suite of population genetics analyses, including calculations of
### heterozygosity, Fis, pairwise Fst, and Principal Components Analysis (PCA).
### All results are compiled into a summary PDF.
###
################################################################################


# ==============================================================================
# Part 1: Setup & Data Loading
# ==============================================================================

### Step 1.1: Load Libraries ----
# ------------------------------------------------------------------------------
message("--- Part 1: Setup & Data Loading ---")

library(vcfR)
library(tidyverse)
library(adegenet) # Essential for data conversion and PCA
library(hierfstat) # For Ho, Hs, Fis, and Fst calculations
library(patchwork) # For combining plots
library(gridExtra) # For creating table grobs
library(beepr)

# Set your working directory (should be the same as the filtering script)
setwd("G:/Shared drives/Sard Lab/Round Goby/RADseq/R/post_RAD analysis")


### Step 1.2: Load Filtered VCF Data ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.2: Loading the Filtered VCF File ---")

vcf_filtered <- read.vcfR("Output/populations.filtered.final.vcf.gz")
message(paste("Successfully loaded VCF with", nrow(vcf_filtered), "loci and", ncol(vcf_filtered@gt)-1, "individuals."))


### Step 1.3: Prepare Data for Analysis (CORRECTED) ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.3: Converting VCF to adegenet's genind format ---")

# --- THE FIX IS HERE ---
# Convert the vcfR object directly to a genind object.
genind_obj <- vcfR2genind(vcf_filtered)

# Create a population map from the sample names (first 4 characters = year)
pop_map <- tibble(
  Sample = indNames(genind_obj),
  Population = substr(Sample, 1, 4)
)

# Assign the population information to the genind object
pop(genind_obj) <- pop_map$Population
message("Population information has been assigned to the dataset.")
beep(sound = 2)

# ==============================================================================
# Part 2: Genetic Diversity & Inbreeding Analysis
# ==============================================================================
message("\n--- Part 2: Genetic Diversity & Inbreeding Analysis ---")

### Step 2.1: Calculate Basic Stats (Ho, Hs, Fis) (CORRECTED) ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.1: Calculating Heterozygosity and Fis per population ---")

# --- THE FIX IS HERE ---
# Convert the genind object to hierfstat's format. This is the correct function.
hierf_data <- genind2hierfstat(genind_obj)

# Calculate stats per locus for each population
basic_stats <- basic.stats(hierf_data, diploid = TRUE)
beep(sound = 2)

# --- Create Heterozygosity Plot ---
# Combine Ho and Hs into a single tidy data frame for plotting
het_long <- as.data.frame(basic_stats$Ho) %>%
  rownames_to_column("Locus") %>%
  pivot_longer(-Locus, names_to = "Population", values_to = "Ho") %>%
  left_join(
    as.data.frame(basic_stats$Hs) %>%
      rownames_to_column("Locus") %>%
      pivot_longer(-Locus, names_to = "Population", values_to = "Hs"),
    by = c("Locus", "Population")
  ) %>%
  pivot_longer(cols = c(Ho, Hs), names_to = "Stat", values_to = "Value")

p_het_boxplot <- ggplot(het_long, aes(x = Population, y = Value, fill = Stat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    title = "Heterozygosity per Year",
    y = "Heterozygosity"
  ) +
  scale_fill_manual(
    name = "",
    values = c("Ho" = "cornflowerblue", "Hs" = "lightgoldenrodyellow"),
    labels = c("Observed (Ho)", "Expected (Hs)")
  ) +
  theme_bw()

print(p_het_boxplot)
beep(sound = 2)

# --- Create Fis Plot ---
fis_long <- as.data.frame(basic_stats$Fis) %>%
  rownames_to_column("Locus") %>%
  pivot_longer(-Locus, names_to = "Population", values_to = "Fis")

p_fis_boxplot <- ggplot(fis_long, aes(x = Population, y = Fis, fill = Population)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = "Inbreeding Coefficient (Fis) per Year",
    y = "Fis"
  ) +
  theme_bw() +
  theme(legend.position = "none")

print(p_fis_boxplot)
beep(sound = 2)


# ==============================================================================
# Part 3: Population Structure Analysis
# ==============================================================================
message("\n--- Part 3: Population Structure Analysis ---")

### Step 3.1: Pairwise Fst ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.1: Calculating and plotting pairwise Fst ---")

# Calculate Weir and Cockerham's Fst
pairwise_fst_matrix <- pairwise.WCfst(hierf_data)

# Prepare the matrix for plotting with ggplot
fst_df <- as.matrix(pairwise_fst_matrix)
fst_df[upper.tri(fst_df)] <- NA # Set upper triangle to NA to avoid redundancy
fst_df <- as.data.frame(fst_df) %>%
  rownames_to_column("Pop1") %>%
  pivot_longer(-Pop1, names_to = "Pop2", values_to = "Fst")
fst_df <- fst_df[!is.na(fst_df$Fst),]

p_fst_heatmap <- ggplot(fst_df, aes(x = Pop1, y = Pop2, fill = Fst, label = round(Fst, 3))) +
  geom_tile(color = "black") +
  geom_text(size = 4, color = "black") + # Add Fst values as text
  scale_fill_gradient(low = "white", high = "firebrick", name = "Fst") +
  labs(
    title = "Pairwise Fst Matrix",
    x = "", y = ""
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_fst_heatmap)
beep(sound = 2)


### Step 3.2: Principal Components Analysis (PCA) ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.2: Performing Principal Components Analysis (PCA) ---")

# Scale the data (center but don't scale variance) - standard for genetic PCA
# This function works perfectly on a genind object.
scaled_genind <- scaleGen(genind_obj, center = TRUE, scale = FALSE, NA.method = "mean")

# Perform the PCA
pca_result <- dudi.pca(scaled_genind, cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)

# Calculate the variance explained by each axis
variance_explained <- (pca_result$eig / sum(pca_result$eig)) * 100

# Create a tidy data frame for plotting
pca_scores <- as_tibble(pca_result$li) %>%
  mutate(Sample = indNames(genind_obj)) %>%
  left_join(pop_map, by = "Sample")

# --- Create PCA Plots ---
# Create axis labels with variance explained
pc1_lab <- paste0("PC1 (", round(variance_explained[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(variance_explained[2], 1), "%)")
pc3_lab <- paste0("PC3 (", round(variance_explained[3], 1), "%)")

p_pca_1v2 <- ggplot(pca_scores, aes(x = Axis1, y = Axis2, color = Population)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of Genetic Variation", x = pc1_lab, y = pc2_lab) +
  theme_bw()

p_pca_1v3 <- ggplot(pca_scores, aes(x = Axis1, y = Axis3, color = Population)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of Genetic Variation", x = pc1_lab, y = pc3_lab) +
  theme_bw()

p_pca_2v3 <- ggplot(pca_scores, aes(x = Axis2, y = Axis3, color = Population)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of Genetic Variation", x = pc2_lab, y = pc3_lab) +
  theme_bw()

# Print the main PCA plot
print(p_pca_1v2)
print(p_pca_1v3)
print(p_pca_2v3)
beep(sound = 2)

# ==============================================================================
# Part 4: Final Report Generation
# ==============================================================================
message("\n--- Part 4: Generating Final PDF Report ---")

# Combine all plots into a single report using patchwork
final_report_layout <- (p_het_boxplot + p_fis_boxplot) /
  (p_fst_heatmap + p_pca_1v2) /
  (p_pca_1v3 + p_pca_2v3) +
  plot_annotation(
    title = 'Population Genetics Summary of Filtered VCF Data',
    subtitle = paste("Analysis of", nLoc(genind_obj), "loci in", nInd(genind_obj), "individuals")
  )

# Save the report to a PDF
ggsave(
  "Output/PopGen_Analysis_Report.pdf",
  plot = final_report_layout,
  width = 14,
  height = 16
)

message("Population genetics analysis report saved to: Output/PopGen_Analysis_Report.pdf")
message("\n--- SCRIPT COMPLETE ---")
beep(sound = 8)
