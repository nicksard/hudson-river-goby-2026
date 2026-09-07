################################################################################
### Population Genetics Analysis of Filtered VCF Data -- BIRTH COHORT VERSION
###
### Author: Nicholas Sard, Ph.D.
### Assistant: Claude
### Date: August 2026
###
### Description: This script reads in a filtered VCF file and performs a
### standard suite of population genetics analyses. It includes a summary
### statistics table and compiles all plots into a multi-page PDF report.
###
### ---------------------------------------------------------------------------
### WHAT CHANGED FROM 02_vcf_basic_pop_gen_analysis.R
### ---------------------------------------------------------------------------
### The original grouped individuals by COLLECTION YEAR:
###
###     pop_map <- tibble(Sample = indNames(genind_obj),
###                       Population = substr(Sample, 1, 4))
###
### A collection-year sample is a MIXTURE of two birth cohorts (Age-0 + Age-1),
### so adjacent years share fish spawned in the same year. That shared drift
### history deflates pairwise Fst, which inflates temporal Ne. It also means
### Figure 2 was grouped differently from Table 3, which is already cohort-based.
###
### This version assigns populations by BIRTH COHORT, using the SAME derivation
### already in Pedigree_Pipeline_for_loop_cohort_locations_combo.R and
### vcf2genepop_format_cohort_Albany_Troy.R:
###
###     Birth_Year <- Sample_Year - Age
###
### so all four scripts now agree on what a cohort is.
###
### Changes are marked  ### NEW  or  ### CHANGED  so you can diff against 02.
###   Step 1.3  NEW      load bio data, derive Birth_Year, length QC
###   Step 1.4  CHANGED  pop_map now uses Cohort, not substr(Sample, 1, 4)
###   Step 3.2  NEW      temporal Ne from the cohort Fst matrix
### Everything else is structurally identical to 02.
################################################################################


# ==============================================================================
# Part 0: Options  ### NEW
# ==============================================================================

bio_data_path <- paste0("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/",
                        "plate_barcodes/Input/",
                        "GobyRad_tabularPlateWell with Bio Data.xlsx")

MIN_YEARS_APART   <- 2      # Waples 1990; Waples & Yokota 2007
GENERATION_LENGTH <- 1      # years per generation. 1 = reproduce the draft's
# year-based numbers. Set to ~1.5-2 for the
# biologically correct value -- it scales EVERY
# Ne linearly, so state your choice in Methods.
DROP_COHORT_BELOW <- 5      # drop cohorts this small (e.g. the single 2020 fish)
RUN_BOOTSTRAP     <- TRUE   # bootstrap CIs on pairwise Fst (slow)
N_BOOT            <- 200

# Expected cohort N, back-calculated from Table 3. Used only as a sanity check.
EXPECTED_COHORT_N <- c("2021" = 44, "2022" = 16, "2023" = 166, "2024" = 94)


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
library(readxl)      ### NEW -- to read the bio data workbook
library(beepr)

# Set your working directory
setwd("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis")


### Step 1.2: Load Filtered VCF Data ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.2: Loading the Filtered VCF File ---")

vcf_filtered <- read.vcfR("Output/populations.filtered.final.vcf.gz")
message(paste("Successfully loaded VCF with", nrow(vcf_filtered), "loci and", ncol(vcf_filtered@gt)-1, "individuals."))


### Step 1.3: Load Bio Data & Infer Birth Cohorts ----   ### NEW
# ------------------------------------------------------------------------------
# Identical name-cleaning and cohort derivation to the COLONY / Genepop scripts.
message("\n--- Step 1.3: Loading Biological Data & Inferring Cohorts ---")

bio_data <- read_xlsx(path = bio_data_path, sheet = 1)
names(bio_data) <- trimws(names(bio_data))
names(bio_data) <- trimws(gsub(pattern = " ", replacement = "_", x = names(bio_data)))
names(bio_data) <- trimws(gsub(pattern = "\\(|\\)", replacement = "", x = names(bio_data)))
names(bio_data) <- trimws(gsub(pattern = "_\\[|\\]_", replacement = "", x = names(bio_data)))

# Clean bad IDs
bio_data <- bio_data[!grepl(pattern = "NTC", x = bio_data$ID), ]
bio_data <- bio_data[!grepl(pattern = "999", x = bio_data$ID), ]

# Infer Cohorts
bio_data$Sample_Year <- as.numeric(substr(x = bio_data$ID, start = 1, stop = 4))
bio_data$Age_num <- suppressWarnings(as.numeric(bio_data$Age))
bio_data$Birth_Year <- bio_data$Sample_Year - bio_data$Age_num

# --- Report age-data completeness. This is the load-bearing assumption of the
# --- whole script, so it gets printed rather than silently filtered.
message(sprintf("  Bio records: %d | with usable Age: %d | missing Age: %d",
                nrow(bio_data), sum(!is.na(bio_data$Age_num)),
                sum(is.na(bio_data$Age_num))))
if (any(is.na(bio_data$Age_num))) {
  message("  !! Fish with no Age get no cohort and will be dropped downstream.")
}

# --- Length QC. NOTE: length is NOT used to assign cohorts here (Age is).
# --- This is a data-integrity check for the dispersal figure.
len_col <- intersect(c("Length_mm","Total_Length_mm","TL_mm","Length","TL","Length_MM"),
                     names(bio_data))
if (length(len_col) > 0) {
  bio_data$Length_QC <- suppressWarnings(as.numeric(bio_data[[len_col[1]]]))
  bad_len <- bio_data[!is.na(bio_data$Length_QC) &
                        (bio_data$Length_QC < 12 | bio_data$Length_QC > 250), ]
  if (nrow(bad_len) > 0) {
    message(sprintf("  !! %d implausible length value(s) flagged:", nrow(bad_len)))
    print(as.data.frame(bad_len[, intersect(c("ID","Fish_ID","Region","Sample_Year",
                                              "Age_num","Birth_Year","Length_QC"),
                                            names(bad_len))]))
    write_csv(bad_len, "Output/length_QC_flags.csv")
  }
  if (any(bio_data$ID == "2024221", na.rm = TRUE)) {
    message(sprintf("  >> Sample 2024221 length = %s mm  (Lock C1)",
                    paste(bio_data$Length_QC[bio_data$ID == "2024221"], collapse = ", ")))
    message("  >> 2 mm is below hatch size (~4-5 mm) and cannot be taken by a")
    message("  >> backpack electrofisher. Data-entry error. Anchors Cluster 70.")
  }
} else {
  message("  No length column found; skipping length QC.")
  bio_data$Length_QC <- NA_real_
}

# Region cleaning, same as the COLONY pipeline
bio_data$Region <- trimws(bio_data$Region)
bio_data$Region[grepl(pattern = "Albany|Troy", x = bio_data$Region)] <- "Albany_Troy"


### Step 1.4: Prepare Data for Analysis ----   ### CHANGED
# ------------------------------------------------------------------------------
message("\n--- Step 1.4: Converting VCF to adegenet's genind format ---")

genind_obj <- vcfR2genind(vcf_filtered)

# --- WAS: Population = substr(Sample, 1, 4)   [collection year]
# --- NOW: Population = Birth_Year             [cohort]
# Fish_ID is the key that matches VCF sample names (see the genepop script).
pop_map <- tibble(Sample = indNames(genind_obj)) %>%
  left_join(
    bio_data %>%
      transmute(Sample       = as.character(Fish_ID),
                Collect_Year = as.character(Sample_Year),
                Age          = Age_num,
                Population   = as.character(Birth_Year),   # <- the cohort
                Region       = Region,
                Length_mm    = Length_QC),
    by = "Sample"
  )

# Fall back to the ID prefix for collection year if the join missed anything
pop_map$Collect_Year[is.na(pop_map$Collect_Year)] <-
  substr(pop_map$Sample[is.na(pop_map$Collect_Year)], 1, 4)

n_unmatched <- sum(is.na(pop_map$Population))
if (n_unmatched > 0) {
  message(sprintf("  !! %d genotyped fish have no cohort and will be dropped:", n_unmatched))
  print(pop_map$Sample[is.na(pop_map$Population)])
}

# --- Sanity check against Table 3 ---
cohort_n <- pop_map %>% filter(!is.na(Population)) %>% count(Population, name = "N")
message("\n  Cohort sample sizes:")
print(as.data.frame(cohort_n))
for (ch in names(EXPECTED_COHORT_N)) {
  got <- cohort_n$N[cohort_n$Population == ch]; got <- if (length(got) == 0) 0 else got
  message(sprintf("  Cohort %s: observed %3d | Table 3 expects %3d  [%s]",
                  ch, got, EXPECTED_COHORT_N[[ch]],
                  ifelse(got == EXPECTED_COHORT_N[[ch]], "OK", "MISMATCH")))
}
message("  A MISMATCH means this script and Table 3 disagree on cohort membership.")
message("  Resolve that before reporting either one.")

# --- Collection year x cohort. This cross-tab IS the age structure, and it
# --- shows exactly how much the old grouping was mixing cohorts.
message("\n  Collection year x birth cohort:")
print(table(pop_map$Collect_Year, pop_map$Population, useNA = "ifany"))

write_csv(pop_map, "Output/cohort_assignment_table.csv")

# --- Drop unassigned and tiny cohorts, then subset ---
small <- cohort_n$Population[cohort_n$N < DROP_COHORT_BELOW]
if (length(small) > 0) {
  message(sprintf("\n  Dropping cohort(s) with N < %d: %s",
                  DROP_COHORT_BELOW, paste(small, collapse = ", ")))
}
keep <- which(!is.na(pop_map$Population) & !(pop_map$Population %in% small))

genind_obj <- genind_obj[keep, ]
pop_map    <- pop_map[keep, ]
pop(genind_obj) <- pop_map$Population
message(sprintf("Cohort information has been assigned: %d individuals, %d cohorts.",
                nInd(genind_obj), nlevels(pop(genind_obj))))


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
# Ratio of sums, not median of ratios. Equals WC's f to 4 decimals for a
# single population, and unlike the median it is unbiased at small n.
fis_summary <- het_long %>% group_by(Population) %>%
  summarise(Fis = 1 - sum(Ho, na.rm = TRUE) / sum(Hs, na.rm = TRUE))

# WC per-cohort Fis with bootstrap CIs (fast -- no pairwise loop)
boot_fis <- boot.ppfis(dat = hierf_data, nboot = 200, quant = c(0.025, 0.975))
fis_ci <- as.data.frame(boot_fis$fis.ci)
fis_ci$Population <- levels(pop(genind_obj))
fis_summary <- left_join(fis_summary, fis_ci, by = "Population")
popgen_summary_table <- pop_n %>%
  left_join(het_summary, by = "Population") %>%
  left_join(fis_summary, by = "Population") %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

print("--- Population Genetics Summary Statistics (by birth cohort) ---")
print(as.data.frame(popgen_summary_table))
write_csv(popgen_summary_table, "Output/diversity_by_cohort.csv")   ### NEW
beep(sound = 2)


### Step 2.3: Create Boxplots ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.3: Generating visualization plots ---")
p_het_boxplot <- ggplot(het_long %>% pivot_longer(cols = c(Ho, Hs), names_to = "Stat", values_to = "Value"),
                        aes(x = Population, y = Value, fill = Stat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(title = "Heterozygosity per Birth Cohort", y = "Heterozygosity", x = "Birth cohort") +
  scale_fill_manual(name="", values=c("Ho"="cornflowerblue", "Hs"="lightgoldenrodyellow"), labels=c("Observed (Ho)", "Expected (Hs)")) +
  theme_bw()
print(p_het_boxplot)

p_fis_boxplot <- ggplot(fis_long, aes(x = Population, y = Fis, fill = Population)) +
  geom_boxplot(outlier.alpha = 0.3) + geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Inbreeding Coefficient (Fis) per Birth Cohort", y = "Fis", x = "Birth cohort") +
  theme_bw() + theme(legend.position = "none")
print(p_fis_boxplot)


# ==============================================================================
# Part 3: Population Structure Analysis
# ==============================================================================
message("\n--- Part 3: Population Structure Analysis ---")

### Step 3.1: Pairwise Fst ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.1: Calculating and plotting pairwise Fst ---")
pairwise_fst_matrix <- pairwise.WCfst(hierf_data)

# --- Your more robust method for handling NAs ---
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
  labs(title = "Pairwise Fst Matrix (birth cohort)", x = "", y = "") + theme_bw()
print(p_fst_heatmap)

p_fst_heatmap1 <- ggplot(fst_df, aes(x = Pop1, y = Pop2, fill = Fst, label = sprintf("%.4f", Fst))) +
  geom_tile(color = "black") + geom_text(size = 4, color = "black") +
  scale_fill_gradient(low = "white", high = "firebrick", name = "Fst") +
  labs(x = "", y = "") + theme_bw() +
  theme(axis.text = element_text(size=14,color="black"),
        legend.text  = element_text(size=14,color="black"),
        axis.title = element_text(size=14,color="black"),
        legend.title  = element_text(size=14,color="black"))
ggsave(filename = "Output/pairwise_Fst_cohort.png",plot = p_fst_heatmap1,device = "png",width = 11,height = 8,units = "in",dpi = 300)

write_csv(fst_df, "Output/pairwise_Fst_cohort.csv")   ### NEW
beep(sound = 2)

### Step 3.2: Temporal Ne from cohort Fst ----   ### NEW
# ------------------------------------------------------------------------------
# Waples (1989/1990) temporal method. For two temporal samples F ~= 2*Fst, and
# Ne = t / (2F), which gives:
#
#       Ne = t / (4 * Fst)
#
# t is in GENERATIONS, hence GENERATION_LENGTH above. Only comparisons spanning
# >= MIN_YEARS_APART are retained (Waples 1990; Waples & Yokota 2007); one-year
# comparisons are excluded on methodological grounds.
message("\n--- Step 3.2: Estimating temporal Ne by cohort ---")

temporal_ne <- fst_df %>%
  mutate(Yr1 = as.numeric(Pop1),
         Yr2 = as.numeric(Pop2),
         Years_Apart = abs(Yr2 - Yr1),
         t_gen = Years_Apart / GENERATION_LENGTH) %>%
  filter(Years_Apart >= MIN_YEARS_APART) %>%
  left_join(pop_n %>% rename(Pop1 = Population, N1 = N), by = "Pop1") %>%
  left_join(pop_n %>% rename(Pop2 = Population, N2 = N), by = "Pop2") %>%
  mutate(Ne = t_gen / (4 * Fst),
         Min_N = pmin(N1, N2),
         Flag = ifelse(Min_N < 30, "THIN - low precision", "")) %>%
  select(Pop1, Pop2, N1, N2, Min_N, Years_Apart, t_gen, Fst, Ne, Flag) %>%
  arrange(Pop1, Pop2)

message("\n--- Temporal Ne by birth cohort ---")
print(as.data.frame(temporal_ne %>%
                      mutate(Fst = round(Fst, 6), Ne = round(Ne, 0))))

if (nrow(temporal_ne) > 0) {
  message(sprintf("\n  Cohort temporal Ne range: %.0f - %.0f",
                  min(temporal_ne$Ne), max(temporal_ne$Ne)))
  message("  Collection-year version in the current draft: 536 - 625")
  message("  Predicted from the cohort-mixture algebra:    401 - 514")
}
write_csv(temporal_ne, "Output/temporal_Ne_cohort.csv")

### Step 3.2b: Bootstrap CIs on pairwise Fst (optional) ----
# ------------------------------------------------------------------------------
# CAUTION: boot.ppfst resamples LOCI independently. Over t = 2-3 generations
# recombination has barely broken up the 23 chromosomes, so 10,915 loci are
# nowhere near independent and these CIs are optimistically narrow. Treat them
# as a lower bound on the true uncertainty, and say so if you report them.
if (RUN_BOOTSTRAP) {
  message("\n--- Step 3.2b: Bootstrapping pairwise Fst (this is slow) ---")
  boot_fst <- boot.ppfst(dat = hierf_data, nboot = N_BOOT, quant = c(0.025, 0.975))
  
  ci_df <- expand.grid(Pop1 = rownames(boot_fst$ll),
                       Pop2 = colnames(boot_fst$ll),
                       stringsAsFactors = FALSE) %>%
    mutate(Fst_L95 = as.vector(boot_fst$ll),
           Fst_U95 = as.vector(boot_fst$ul)) %>%
    filter(!is.na(Fst_L95))
  
  temporal_ne_ci <- temporal_ne %>%
    left_join(ci_df, by = c("Pop1", "Pop2")) %>%
    # bind both orientations in case the matrix triangle is flipped
    left_join(ci_df %>% rename(Pop1 = Pop2, Pop2 = Pop1,
                               Fst_L95b = Fst_L95, Fst_U95b = Fst_U95),
              by = c("Pop1", "Pop2")) %>%
    mutate(Fst_L95 = coalesce(Fst_L95, Fst_L95b),
           Fst_U95 = coalesce(Fst_U95, Fst_U95b),
           # NOTE the inversion: low Fst -> high Ne
           Ne_L95 = t_gen / (4 * Fst_U95),
           Ne_U95 = ifelse(Fst_L95 <= 0, Inf, t_gen / (4 * Fst_L95))) %>%
    select(Pop1, Pop2, N1, N2, Years_Apart, Fst, Fst_L95, Fst_U95,
           Ne, Ne_L95, Ne_U95, Flag)
  
  message("\n--- Temporal Ne with bootstrap CIs (narrow -- see caution above) ---")
  print(as.data.frame(temporal_ne_ci %>%
                        mutate(across(c(Fst, Fst_L95, Fst_U95), ~round(.x, 6)),
                               across(c(Ne, Ne_L95, Ne_U95), ~round(.x, 0)))))
  write_csv(temporal_ne_ci, "Output/temporal_Ne_cohort_withCI.csv")
  beep(sound = 2)
}


### Step 3.3: Principal Components Analysis (PCA) ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.3: Performing Principal Components Analysis (PCA) ---")
scaled_genind <- scaleGen(genind_obj, center = TRUE, scale = FALSE, NA.method = "mean")
pca_result <- dudi.pca(scaled_genind, cent = FALSE, scale = FALSE, nf = 3, scannf = FALSE)
variance_explained <- (pca_result$eig / sum(pca_result$eig)) * 100
pca_scores <- as_tibble(pca_result$li) %>% mutate(Sample = indNames(genind_obj)) %>%
  left_join(pop_map %>% select(Sample, Population), by = "Sample")

pc1_lab <- paste0("PC1 (", round(variance_explained[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(variance_explained[2], 1), "%)")
pc3_lab <- paste0("PC3 (", round(variance_explained[3], 1), "%)")

# 1. PCA: PC1 vs PC2
p_pca_1v2 <- ggplot(pca_scores, aes(x = Axis1, y = Axis2, color = Population)) +
  geom_point(size = 3, alpha = 0.3) +
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_point(data = pca_scores %>% group_by(Population) %>% summarise(Axis1 = mean(Axis1), Axis2 = mean(Axis2)),
             aes(x = Axis1, y = Axis2), size = 6, shape = 18, show.legend = FALSE) +
  labs(title = "PCA of Genetic Variation (birth cohort)", x = pc1_lab, y = pc2_lab) +
  theme_bw()

# 2. PCA: PC1 vs PC3
p_pca_1v3 <- ggplot(pca_scores, aes(x = Axis1, y = Axis3, color = Population)) +
  geom_point(size = 3, alpha = 0.3) +
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_point(data = pca_scores %>% group_by(Population) %>% summarise(Axis1 = mean(Axis1), Axis3 = mean(Axis3)),
             aes(x = Axis1, y = Axis3), size = 6, shape = 18, show.legend = FALSE) +
  labs(title = "PCA of Genetic Variation (birth cohort)", x = pc1_lab, y = pc3_lab) +
  theme_bw()

# 3. PCA: PC2 vs PC3
p_pca_2v3 <- ggplot(pca_scores, aes(x = Axis2, y = Axis3, color = Population)) +
  geom_point(size = 3, alpha = 0.3) +
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed", size = 0.7, show.legend = FALSE) +
  geom_point(data = pca_scores %>% group_by(Population) %>% summarise(Axis2 = mean(Axis2), Axis3 = mean(Axis3)),
             aes(x = Axis2, y = Axis3), size = 6, shape = 18, show.legend = FALSE) +
  labs(title = "PCA of Genetic Variation (birth cohort)", x = pc2_lab, y = pc3_lab) +
  theme_bw()

message("\nDisplaying all three PCA comparisons with correct spatial centroids:")
print(p_pca_1v2)
print(p_pca_1v3)
print(p_pca_2v3)
beep(sound = 2)


# ==============================================================================
# Part 4: Final Report Generation
# ==============================================================================
message("\n--- Part 4: Generating Final PDF Report ---")

p_popgen_summary_table <- tableGrob(popgen_summary_table, rows = NULL)

### NEW -- temporal Ne table as a report page
p_ne_table <- tableGrob(
  as.data.frame(temporal_ne %>%
                  mutate(Fst = round(Fst, 6), Ne = round(Ne, 0), t_gen = round(t_gen, 2))),
  rows = NULL)

pdf("Output/PopGen_Analysis_Report_COHORT.pdf", width = 14, height = 11)

# Page 1: Main results
grid.arrange(
  p_het_boxplot, p_fis_boxplot,
  p_fst_heatmap, p_pca_1v2,
  ncol = 2,
  top = grid::textGrob(
    paste("Population Genetics Summary of Filtered VCF Data -- BIRTH COHORT grouping\nAnalysis of",
          nLoc(genind_obj), "loci in", nInd(genind_obj), "individuals"),
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
  top = grid::textGrob("Quantitative Summary Statistics by Birth Cohort", gp = grid::gpar(fontsize = 14))
)

# Page 4: Temporal Ne   ### NEW
grid.arrange(
  p_ne_table,
  top = grid::textGrob(
    sprintf("Temporal Ne by Birth Cohort (>= %d years apart; generation length %.2f y)\nNe = t / (4 * Fst)",
            MIN_YEARS_APART, GENERATION_LENGTH),
    gp = grid::gpar(fontsize = 13))
)

dev.off()

message("Population genetics analysis report saved to: Output/PopGen_Analysis_Report_COHORT.pdf")
message("\n--- SCRIPT COMPLETE ---")
beep(sound = 8)
