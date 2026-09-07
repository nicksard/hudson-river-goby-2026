################################################################################
### Population Genetics Analysis of Filtered VCF Data
###
### Author: Nicholas Sard, Ph.D.
### Assistant: Gemini
### Date: October 14, 2025
###
### Description: This script reads in a filtered VCF file and performs a
### standard suite of population genetics analyses. It includes a summary
### statistics table and compiles all plots into a multi-page PDF report.
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
library(adegenet)
library(hierfstat)
library(patchwork)
library(gridExtra)
library(grid)
library(beepr)

# Set your working directory
setwd("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis")


### Step 1.2: Load Filtered VCF Data ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.2: Loading the Filtered VCF File ---")

vcf_filtered <- read.vcfR("Output/populations.filtered.final.vcf.gz")
message(paste("Successfully loaded VCF with", nrow(vcf_filtered), "loci and", ncol(vcf_filtered@gt)-1, "individuals."))


### Step 1.3: Prepare Data for Analysis ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.3: Converting VCF to adegenet's genind format ---")

genind_obj <- vcfR2genind(vcf_filtered)
pop_map <- tibble(
  Sample = indNames(genind_obj),
  Population = substr(Sample, 1, 4)
)
pop(genind_obj) <- pop_map$Population
message("Population information has been assigned to the dataset.")


# ==============================================================================
# Part 2: Genetic Diversity & Inbreeding Analysis
# ==============================================================================
message("\n--- Part 2: Genetic Diversity & Inbreeding Analysis ---")

### Step 2.1: Calculate Basic Per-Locus Stats ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.1: Calculating Heterozygosity and Fis per locus ---")
hierf_data <- genind2hierfstat(genind_obj)
basic_stats <- basic.stats(hierf_data, diploid = TRUE)


### Step 2.2: Create Summary Statistics Table ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.2: Generating Summary Statistics Table ---")
pop_n <- pop_map %>% count(Population, name = "N")
het_long <- as.data.frame(basic_stats$Ho) %>%
  rownames_to_column("Locus") %>% pivot_longer(-Locus, names_to = "Population", values_to = "Ho") %>%
  left_join(as.data.frame(basic_stats$Hs) %>% rownames_to_column("Locus") %>%
              pivot_longer(-Locus, names_to = "Population", values_to = "Hs"), by = c("Locus", "Population"))
fis_long <- as.data.frame(basic_stats$Fis) %>%
  rownames_to_column("Locus") %>% pivot_longer(-Locus, names_to = "Population", values_to = "Fis")
het_summary <- het_long %>% group_by(Population) %>%
  summarise(Mean_Ho = mean(Ho, na.rm = TRUE), Mean_Hs = mean(Hs, na.rm = TRUE))
fis_summary <- fis_long %>% group_by(Population) %>%
  summarise(Median_Fis = median(Fis, na.rm = TRUE))
popgen_summary_table <- pop_n %>%
  left_join(het_summary, by = "Population") %>%
  left_join(fis_summary, by = "Population") %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

print("--- Population Genetics Summary Statistics ---")
print(as.data.frame(popgen_summary_table))
beep(sound = 2)


### Step 2.3: Create Boxplots ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.3: Generating visualization plots ---")
p_het_boxplot <- ggplot(het_long %>% pivot_longer(cols = c(Ho, Hs), names_to = "Stat", values_to = "Value"),
                        aes(x = Population, y = Value, fill = Stat)) +
  geom_boxplot(outlier.alpha = 0.3) + labs(title = "Heterozygosity per Year", y = "Heterozygosity") +
  scale_fill_manual(name="", values=c("Ho"="cornflowerblue", "Hs"="lightgoldenrodyellow"), labels=c("Observed (Ho)", "Expected (Hs)")) +
  theme_bw()
print(p_het_boxplot)

p_fis_boxplot <- ggplot(fis_long, aes(x = Population, y = Fis, fill = Population)) +
  geom_boxplot(outlier.alpha = 0.3) + geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Inbreeding Coefficient (Fis) per Year", y = "Fis") + theme_bw() + theme(legend.position = "none")
print(p_fis_boxplot)


# ==============================================================================
# Part 3: Population Structure Analysis
# ==============================================================================
message("\n--- Part 3: Population Structure Analysis ---")

### Step 3.1: Pairwise Fst ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.1: Calculating and plotting pairwise Fst ---")
pairwise_fst_matrix <- pairwise.WCfst(hierf_data)

# --- THE FIX IS HERE: Your more robust method for handling NAs ---
fst_df <- as.matrix(pairwise_fst_matrix)
fst_df[upper.tri(fst_df)] <- NA
fst_df <- as.data.frame(fst_df) %>%
  rownames_to_column("Pop1") %>%
  pivot_longer(-Pop1, names_to = "Pop2", values_to = "Fst")
# Explicitly remove rows where Fst is NA
fst_df <- fst_df[!is.na(fst_df$Fst),]

p_fst_heatmap <- ggplot(fst_df, aes(x = Pop1, y = Pop2, fill = Fst, label = sprintf("%.4f", Fst))) +
  geom_tile(color = "black") + geom_text(size = 4, color = "black") +
  scale_fill_gradient(low = "white", high = "firebrick", name = "Fst") +
  labs(title = "Pairwise Fst Matrix", x = "", y = "") + theme_bw()
print(p_fst_heatmap)

p_fst_heatmap1 <- ggplot(fst_df, aes(x = Pop1, y = Pop2, fill = Fst, label = sprintf("%.4f", Fst))) +
  geom_tile(color = "black") + geom_text(size = 4, color = "black") +
  scale_fill_gradient(low = "white", high = "firebrick", name = "Fst") +
  labs(x = "", y = "") + theme_bw() +
  theme(axis.text = element_text(size=14,color="black"),
        legend.text  = element_text(size=14,color="black"),
        axis.title = element_text(size=14,color="black"),
        legend.title  = element_text(size=14,color="black"))
ggsave(filename = "Output/pairwise_Fst.png",plot = p_fst_heatmap1,device = "png",width = 11,height = 8,units = "in",dpi = 300)
### Step 3.2: Principal Components Analysis (PCA) ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.2: Performing Principal Components Analysis (PCA) ---")
scaled_genind <- scaleGen(genind_obj, center = TRUE, scale = FALSE, NA.method = "mean")
pca_result <- dudi.pca(scaled_genind, cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)
variance_explained <- (pca_result$eig / sum(pca_result$eig)) * 100
pca_scores <- as_tibble(pca_result$li) %>% mutate(Sample = indNames(genind_obj)) %>% left_join(pop_map, by = "Sample")

pc1_lab <- paste0("PC1 (", round(variance_explained[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(variance_explained[2], 1), "%)")
pc3_lab <- paste0("PC3 (", round(variance_explained[3], 1), "%)")

# 1. PCA: PC1 vs PC2
p_pca_1v2 <- ggplot(pca_scores, aes(x = Axis1, y = Axis2, color = Population)) + 
  geom_point(size = 3, alpha = 0.3) + 
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed", size = 0.7, show.legend = FALSE) +
  # Bivariate Centroid Fix
  geom_point(data = pca_scores %>% group_by(Population) %>% summarise(Axis1 = mean(Axis1), Axis2 = mean(Axis2)),
             aes(x = Axis1, y = Axis2), size = 6, shape = 18, show.legend = FALSE) +
  labs(title = "PCA of Genetic Variation", x = pc1_lab, y = pc2_lab) + 
  theme_bw()

# 2. PCA: PC1 vs PC3
p_pca_1v3 <- ggplot(pca_scores, aes(x = Axis1, y = Axis3, color = Population)) + 
  geom_point(size = 3, alpha = 0.3) + 
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed", size = 0.7, show.legend = FALSE) +
  # Bivariate Centroid Fix
  geom_point(data = pca_scores %>% group_by(Population) %>% summarise(Axis1 = mean(Axis1), Axis3 = mean(Axis3)),
             aes(x = Axis1, y = Axis3), size = 6, shape = 18, show.legend = FALSE) +
  labs(title = "PCA of Genetic Variation", x = pc1_lab, y = pc3_lab) + 
  theme_bw()

# 3. PCA: PC2 vs PC3
p_pca_2v3 <- ggplot(pca_scores, aes(x = Axis2, y = Axis3, color = Population)) + 
  geom_point(size = 3, alpha = 0.3) + 
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed", size = 0.7, show.legend = FALSE) +
  # Bivariate Centroid Fix
  geom_point(data = pca_scores %>% group_by(Population) %>% summarise(Axis2 = mean(Axis2), Axis3 = mean(Axis3)),
             aes(x = Axis2, y = Axis3), size = 6, shape = 18, show.legend = FALSE) +
  labs(title = "PCA of Genetic Variation", x = pc2_lab, y = pc3_lab) + 
  theme_bw()

# --- Print all three PCA plots to the console ---
message("\nDisplaying all three PCA comparisons with correct spatial centroids:")
print(p_pca_1v2)
print(p_pca_1v3)
print(p_pca_2v3)
beep(sound = 2)

# ==============================================================================
# Part 4: Final Report Generation
# ==============================================================================
message("\n--- Part 4: Generating Final PDF Report ---")

# Convert the summary table into a graphical object
p_popgen_summary_table <- tableGrob(popgen_summary_table, rows = NULL)

# Save the multi-page report to a PDF
pdf("Output/PopGen_Analysis_Report.pdf", width = 14, height = 11)

# --- THE FIX IS HERE: Using grid.arrange for PDF layout ---
# Page 1: Main results
grid.arrange(
  p_het_boxplot, p_fis_boxplot,
  p_fst_heatmap, p_pca_1v2,
  ncol = 2,
  top = grid::textGrob(
    paste("Population Genetics Summary of Filtered VCF Data\nAnalysis of", nLoc(genind_obj), "loci in", nInd(genind_obj), "individuals"),
    gp = grid::gpar(fontsize = 14)
  )
)

# Page 2: Supplementary PCA axes
grid.arrange(
  p_pca_1v3, p_pca_2v3,
  ncol = 2,
  top = grid::textGrob("Supplementary PCA Plots", gp = grid::gpar(fontsize = 14))
)

# Page 3: The quantitative summary table
grid.arrange(
  p_popgen_summary_table,
  top = grid::textGrob("Quantitative Summary Statistics by Year", gp = grid::gpar(fontsize = 14))
)

dev.off()

message("Population genetics analysis report saved to: Output/PopGen_Analysis_Report.pdf")
message("\n--- SCRIPT COMPLETE ---")
beep(sound = 8)