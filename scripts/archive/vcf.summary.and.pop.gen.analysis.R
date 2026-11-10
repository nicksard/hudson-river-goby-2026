################################################################################
### VCF Filtering & Analysis Pipeline (Optimized Version)
###
### Description: This script follows a loci-first filtering approach, optimized
### to remove low-quality and uninformative markers as early as possible to
### improve computational efficiency.
###
### New Workflow:
### 1. Setup & Pre-filtering: Load data, remove NTCs, visualize raw depth.
### 2. Locus Culling: Apply MAF filter, then thin to one SNP per RAD locus
###    using an intelligent selection process that prioritizes low-coverage samples.
### 3. Locus Quality Control: On the reduced SNP set, perform HWE and
###    paralog filtering.
### 4. Genotype & Individual Filtering: Apply quality filters at the genotype
###    level, then remove low-call-rate individuals and loci.
### 5. Summary & Analysis: Generate final summary plots and proceed with
###    population genetics.
################################################################################


# ==============================================================================
# Part 1: Initial Setup & Pre-filtering
# ==============================================================================

### Step 1.1: Load Libraries and Data ----
# ------------------------------------------------------------------------------
# Loading necessary libraries
library(vcfR)
library(tidyverse)
library(adegenet)
library(hierfstat)
library(beepr)
library(patchwork)
library(gridExtra)
library(grid)

# Set yogrid# Set your working directory
setwd("G:/Shared drives/Sard Lab/Round Goby/RADseq/R/post_RAD analysis")

# Read in the VCF file
vcf <- read.vcfR("Input//populations.snps.vcf")

# Rename the problematic NTC sample
message("Renaming sample 'NTC4a' to '2024999'.")
colnames(vcf@gt)[colnames(vcf@gt) == "NTC4a"] <- "2024999"

# Initialize a clean filtering summary table
filtering_summary <- tibble(
  Step = "Initial State",
  Filter_Applied = "N/A",
  Loci_Removed = 0,
  Loci_Remaining = nrow(vcf),
  Individuals_Removed = 0,
  Individuals_Remaining = ncol(vcf@gt) - 1
)

### Step 1.2: Negative Control (NTC) Analysis & Removal ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.2: Negative Control (NTC) Analysis ---")

# Identify all NTC samples
ntc_samples <- colnames(vcf@gt)[grep("NTC", colnames(vcf@gt))]
message(paste("Found", length(ntc_samples), "NTC samples to remove."))

# Remove NTC samples from the VCF object
all_samples <- colnames(vcf@gt)
samples_to_keep <- setdiff(all_samples, ntc_samples)
vcf <- vcf[, samples_to_keep]

# Log this step
filtering_summary <- filtering_summary %>% add_row(
  Step = "NTC Removal",
  Filter_Applied = "Removed all NTC samples",
  Loci_Removed = 0,
  Loci_Remaining = nrow(vcf),
  Individuals_Removed = length(ntc_samples),
  Individuals_Remaining = ncol(vcf@gt) - 1
)
message(paste("Removed NTCs. Kept", ncol(vcf@gt)-1, "samples for analysis."))
beep(sound = 2)


### Step 1.3: Visualize Raw Read Depth ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.3: Visualizing Raw Read Depth By Year ---")
raw_dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
raw_pop_map <- tibble(Sample = colnames(raw_dp_matrix), Population = substr(colnames(raw_dp_matrix), 1, 4))

raw_dp_long <- as.data.frame(raw_dp_matrix) %>%
  rownames_to_column("Locus") %>%
  pivot_longer(-Locus, names_to="Sample", values_to="Depth") %>%
  left_join(raw_pop_map, by="Sample")

p_raw_depth_by_year <- ggplot(raw_dp_long, aes(x=Population, y=Depth, fill=Population)) +
  geom_boxplot(outlier.shape=NA) +
  coord_cartesian(ylim=c(0, 100)) +
  labs(title="Raw Read Depth by Collection Year (Before Filtering)",
       subtitle="Confirms that 2022 samples have lower coverage.",
       x="Year", y="Read Depth (DP)") +
  theme_bw() +
  theme(legend.position="none")

print(p_raw_depth_by_year)
rm(raw_dp_matrix, raw_pop_map, raw_dp_long) # Clean up
beep(sound = 2)


# ==============================================================================
# Part 2: Locus Culling (Major Efficiency Steps)
# ==============================================================================

### Step 2.1: Minor Allele Frequency (MAF) Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.1: Minor Allele Frequency (MAF) Filter ---")
loci_before <- nrow(vcf)
maf_cutoff <- 0.05

# Extract Allele Frequency (AF) from the INFO field. This is fast.
maf_data <- as_tibble(vcf@fix) %>%
  # The regular expression robustly extracts the AF value
  mutate(Info_AF = as.numeric(str_match(INFO, "AF=([0-9.]+)")[,2])) %>%
  mutate(MAF = ifelse(Info_AF > 0.5, 1 - Info_AF, Info_AF))

# Create a logical vector for which loci to keep
loci_to_keep <- maf_data$MAF >= maf_cutoff & !is.na(maf_data$MAF)

vcf <- vcf[loci_to_keep, ]

# Log this step
filtering_summary <- filtering_summary %>% add_row(
  Step = "MAF Filter",
  Filter_Applied = paste("MAF >=", maf_cutoff),
  Loci_Removed = loci_before - nrow(vcf),
  Loci_Remaining = nrow(vcf),
  Individuals_Removed = 0,
  Individuals_Remaining = ncol(vcf@gt) - 1
)
message(paste("Removed loci with MAF <", maf_cutoff, ". Loci remaining:", nrow(vcf)))
rm(maf_data) # Clean up
beep(sound = 2)


### Step 2.2: Thin to One SNP per RAD Locus ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.2: Thinning to One SNP per RAD Locus ---")
message("Strategy: Prioritize SNPs with the best data in low-coverage 2022 samples.")
loci_before <- nrow(vcf)

# Create a tidy data frame of genotypes to calculate call rates
tidy_gt <- vcfR2tidy(vcf, format_fields = c("GT"))$gt

# Calculate call rate (completeness) within the 2022 samples
snp_completeness_2022 <- tidy_gt %>%
  filter(substr(Indiv, 1, 4) == "2022") %>%
  mutate(is_genotyped = ifelse(!is.na(gt_GT) & gt_GT != "./.", 1, 0)) %>%
  group_by(ChromKey) %>%
  summarise(n_genotyped_2022 = sum(is_genotyped))

# Calculate overall call rate
snp_completeness_overall <- tidy_gt %>%
  mutate(is_genotyped = ifelse(!is.na(gt_GT) & gt_GT != "./.", 1, 0)) %>%
  group_by(ChromKey) %>%
  summarise(n_genotyped_overall = sum(is_genotyped))
head(snp_completeness_overall)
# Combine metrics and rank SNPs within each RAD locus
thinned_indices <- as_tibble(vcf@fix) %>%
  mutate(row_index = 1:n()) %>%
  # Extract the RAD locus ID (the part before the first colon)
  separate(ID, into = c("locus_id"), sep = ":", remove = FALSE, extra = "drop") %>%
  # Join our call rate stats
  left_join(snp_completeness_overall, by = c("row_index" = "ChromKey")) %>%
  left_join(snp_completeness_2022, by = c("row_index" = "ChromKey")) %>%
  # Replace NAs with 0 for loci not genotyped at all in 2022 samples
  mutate(n_genotyped_2022 = ifelse(is.na(n_genotyped_2022), 0, n_genotyped_2022)) %>%
  # Group by the RAD locus and arrange to find the best SNP
  group_by(locus_id) %>%
  arrange(
    desc(n_genotyped_2022),    # Priority 1: Keep SNP with most genotypes in 2022 samples
    desc(n_genotyped_overall), # Tie-breaker 1: Keep SNP with most genotypes overall
    .by_group = TRUE
  ) %>%
  # Keep only the top-ranked SNP for each locus
  slice(1) %>%
  ungroup() %>%
  pull(row_index)
thinned_indices
vcf <- vcf[thinned_indices, ]

# Log this step
filtering_summary <- filtering_summary %>% add_row(
  Step = "Thinning",
  Filter_Applied = "One SNP per locus (2022 priority)",
  Loci_Removed = loci_before - nrow(vcf),
  Loci_Remaining = nrow(vcf),
  Individuals_Removed = 0,
  Individuals_Remaining = ncol(vcf@gt) - 1
)
message(paste("Thinned to one SNP per locus. Loci remaining:", nrow(vcf)))
rm(tidy_gt, snp_completeness_2022, snp_completeness_overall, thinned_indices) # Clean up
beep(sound = 2)


# ==============================================================================
# Part 3: Locus Quality Control
# ==============================================================================

### Step 3.1: Hardy-Weinberg Equilibrium Filter (via Fis) ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.1: Hardy-Weinberg Equilibrium (HWE) Analysis ---")
# This is our fully debugged code from our previous session.
# It's now running on a much smaller, more relevant set of SNPs.
loci_before <- nrow(vcf)
final_fis_cutoff <- 0.25

# --- Calculations ---
message("Calculating Fis and HWE p-values...")
temp_genind <- vcfR2genind(vcf)
pop(temp_genind) <- rep("Overall", nInd(temp_genind))
# Fis from hierfstat
hierf_data <- genind2hierfstat(temp_genind)
fis_results <- hierfstat::basic.stats(hierf_data, diploid = TRUE)$perloc %>%
  rownames_to_column("Locus") %>% as_tibble()
fis_results$Locus <- locNames(temp_genind) # Correct mangled names
# P-values from pegas
hwe_results <- pegas::hw.test(temp_genind, B = 0)
hwe_p_values <- as_tibble(hwe_results, rownames = "Locus") %>%
  select(Locus, p_val = `Pr(chi^2 >)`)
# Combine
hwe_summary <- left_join(fis_results, hwe_p_values, by = "Locus")

# --- Plotting ---
message("Generating HWE volcano plot...")
p_hwe_volcano <- ggplot(hwe_summary, aes(x = Fis, y = -log10(p_val))) +
  geom_point(alpha = 0.4, color = "midnightblue") +
  geom_vline(xintercept = c(-final_fis_cutoff, final_fis_cutoff), color = "red", linetype = "dashed") +
  labs(title = "HWE Deviation: Magnitude vs. Significance",
       x = "Inbreeding Coefficient (Fis)", y = "-log10(p-value)") +
  theme_bw()
print(p_hwe_volcano)
beep(sound = 2)

# --- Filtering ---
loci_to_remove_hwe <- hwe_summary %>%
  filter(abs(Fis) > final_fis_cutoff & !is.na(Fis)) %>%
  pull(Locus)
logical_vector_to_keep <- !(locNames(temp_genind) %in% loci_to_remove_hwe)
vcf <- vcf[logical_vector_to_keep, ]

# Log this step
filtering_summary <- filtering_summary %>% add_row(
  Step = "HWE Filter (Fis)",
  Filter_Applied = paste("|Fis| <=", final_fis_cutoff),
  Loci_Removed = loci_before - nrow(vcf),
  Loci_Remaining = nrow(vcf),
  Individuals_Removed = 0,
  Individuals_Remaining = ncol(vcf@gt) - 1
)
message(paste("Removed loci with |Fis| >", final_fis_cutoff, ". Loci remaining:", nrow(vcf)))
rm(temp_genind, hierf_data, fis_results, hwe_results, hwe_p_values, hwe_summary) # Clean up
beep(sound = 2)

### Step 3.2: Paralog Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.2: Paralog Filter ---")
loci_before <- nrow(vcf)

# --- Calculations ---
message("Calculating Hobs and read-ratio deviation to identify paralogs...")
paralog_tidy <- vcfR2tidy(vcf, format_fields = c("GT", "AD"))$gt %>%
  rename(gt = gt_GT, allele_dp = gt_AD, locus_id = ChromKey) %>%
  separate(col = allele_dp, into = c("ref_dp", "alt_dp"), sep = ",", fill="right")

hobs_calc <- paralog_tidy %>%
  filter(!is.na(gt) & gt != "./.") %>%
  group_by(locus_id, gt) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  pivot_wider(names_from = gt, values_from = n, values_fill = 0) %>%
  mutate(
    n_total = `0/0` + `0/1` + `1/1`,
    Hobs = `0/1` / n_total
  ) %>%
  select(locus_id, Hobs)

z_score_calc <- paralog_tidy %>%
  filter(gt == "0/1", !is.na(ref_dp), !is.na(alt_dp)) %>%
  mutate(ref_dp = as.numeric(ref_dp), alt_dp = as.numeric(alt_dp)) %>%
  group_by(locus_id) %>%
  summarise(
    sum_ref = sum(ref_dp),
    sum_alt = sum(alt_dp),
    total_reads = sum_ref + sum_alt,
    z = (sum_ref - (total_reads / 2)) / sqrt(total_reads * 0.25),
    .groups = "drop"
  ) %>%
  select(locus_id, z)

# --- Plotting & Filtering ---
paralog_summary <- left_join(hobs_calc, z_score_calc, by = "locus_id") %>%
  mutate(is_singleton = ifelse(Hobs <= 0.6 & abs(z) <= 7 & !is.na(z), "Yes", "No"))

message("Generating paralog diagnostic plot...")
p_paralog_plot <- ggplot(paralog_summary, aes(x=Hobs, y=z, color=is_singleton)) +
  geom_point(alpha=0.6) +
  labs(
    title="Paralog Identification",
    subtitle="Loci classified as 'No' (potential paralogs) will be removed",
    x="Observed Heterozygosity (Hobs)",
    y="Read-Ratio Deviation (z-score)",
    color="Putative Singleton"
  ) +
  theme_bw()
print(p_paralog_plot)

# --- Filtering (CORRECTED) ---
loci_to_remove_paralogs <- paralog_summary %>%
  filter(is_singleton == "No") %>%
  pull(locus_id)

# FIX: Create a logical vector using the integer index of the VCF rows.
# The 'locus_id' in paralog_summary is the row number (1, 2, 3, etc.).
all_locus_indices <- 1:nrow(vcf)
logical_vector_to_keep <- !(all_locus_indices %in% loci_to_remove_paralogs)

vcf <- vcf[logical_vector_to_keep, ]

# Log this step
filtering_summary <- filtering_summary %>% add_row(
  Step = "Paralog Filter",
  Filter_Applied = "Hobs > 0.6 or |Z| > 7",
  Loci_Removed = loci_before - nrow(vcf),
  Loci_Remaining = nrow(vcf),
  Individuals_Removed = 0,
  Individuals_Remaining = ncol(vcf@gt) - 1
)
message(paste("Removed", length(loci_to_remove_paralogs), "putative paralogs. Loci remaining:", nrow(vcf)))
rm(paralog_tidy, hobs_calc, z_score_calc, paralog_summary, all_locus_indices) # Clean up
beep(sound = 2)

# ==============================================================================
# Part 4: Genotype & Individual Filtering
# ==============================================================================

### Step 4.1: Justifying and Applying Genotype Filters ----
# ------------------------------------------------------------------------------
message("\n--- Step 4.1: Visualizing Filter Thresholds & Masking Genotypes ---")

# --- Extract DP and GQ matrices ---
message("Extracting DP and GQ values from VCF...")
dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
gq_matrix <- extract.gt(vcf, element = "GQ", as.numeric = TRUE)

# Calculate the total number of non-missing genotypes to use as a denominator
total_genotypes <- sum(!is.na(vcf@gt) & vcf@gt != "./.")
message(paste("Found", format(total_genotypes, big.mark=","), "non-missing genotypes to assess."))

# --- Generate DP diagnostic plot ---
message("Analyzing impact of various DP thresholds...")
dp_thresholds_to_test <- seq(from = 0, to = 50, by = 5)

dp_loss_summary <- map_dfr(dp_thresholds_to_test, function(thresh) {
  # Count how many genotypes fall BELOW the current threshold
  masked_count <- sum(dp_matrix < thresh, na.rm = TRUE)
  # Calculate the percentage of total genotypes lost
  percent_loss <- (masked_count / total_genotypes) * 100
  tibble(threshold = thresh, percent_masked = percent_loss)
})

p_dp_impact <- ggplot(dp_loss_summary, aes(x = threshold, y = percent_masked)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Impact of Minimum Read Depth (DP) Filter",
    subtitle = "Red line shows our chosen threshold of 10x",
    x = "Minimum DP Threshold",
    y = "Cumulative % of Genotypes Masked"
  ) +
  theme_bw()

print(p_dp_impact)
beep(sound = 2)

# --- Generate GQ diagnostic plot ---
message("Analyzing impact of various GQ thresholds...")
gq_thresholds_to_test <- seq(from = 0, to = 40, by = 5)

gq_loss_summary <- map_dfr(gq_thresholds_to_test, function(thresh) {
  masked_count <- sum(gq_matrix < thresh, na.rm = TRUE)
  percent_loss <- (masked_count / total_genotypes) * 100
  tibble(threshold = thresh, percent_masked = percent_loss)
})

p_gq_impact <- ggplot(gq_loss_summary, aes(x = threshold, y = percent_masked)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(xintercept = 20, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = "Impact of Minimum Genotype Quality (GQ) Filter",
    subtitle = "Red line shows our chosen threshold of 20",
    x = "Minimum GQ Threshold",
    y = "Cumulative % of Genotypes Masked"
  ) +
  theme_bw()

print(p_gq_impact)
beep(sound = 2)

# --- Apply the chosen filters ---
message("\nProceeding with masking based on the chosen thresholds...")
min_depth_threshold <- 10
min_gq_threshold    <- 20

mask <- (dp_matrix < min_depth_threshold | gq_matrix < min_gq_threshold)
mask[is.na(mask)] <- FALSE # Ensure we don't mask existing NAs

message(paste("Masking", sum(mask), "genotypes with DP <", min_depth_threshold, "OR GQ <", min_gq_threshold))

gt_matrix <- vcf@gt
gt_matrix[mask] <- "./."
vcf@gt <- gt_matrix

message("Successfully masked low-quality genotypes.")
rm(dp_matrix, gq_matrix, gt_matrix, mask) # Clean up memory

### Step 4.2: Individual Call Rate Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 4.2: Filtering Individuals by Missing Data ---")
inds_before <- ncol(vcf@gt) - 1

# --- Calculate Missingness ---
# vcfR2tidy is an efficient way to get per-sample stats
missing_data_summary <- vcfR2tidy(vcf, format_fields = c("GT"))$gt %>%
  group_by(Indiv) %>%
  summarise(
    n_loci = n(),
    n_missing = sum(is.na(gt_GT) | gt_GT == "./."),
    prop_missing = n_missing / n_loci,
    .groups = "drop"
  ) %>%
  mutate(Population = substr(Indiv, 1, 4))

# --- Plot Missingness ---
max_missing_individual <- 0.40 # Set the maximum allowable proportion of missing data

p_missing_data_indivs <- ggplot(missing_data_summary, aes(x=Population, y=prop_missing, fill=Population)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_hline(yintercept = max_missing_individual, color = "red", linetype = "dashed", linewidth=1) +
  labs(
    title="Individual Missing Data (Post-Genotype Filtering)",
    subtitle=paste("Red line shows cutoff for removal at", max_missing_individual * 100, "% missing data"),
    x="Year",
    y="Proportion of Missing Genotypes"
  ) +
  theme_bw() +
  theme(legend.position="none")

print(p_missing_data_indivs)
beep(sound = 2)

# --- Filter Individuals ---
inds_to_remove <- missing_data_summary %>%
  filter(prop_missing > max_missing_individual) %>%
  pull(Indiv)

if(length(inds_to_remove) > 0) {
  message(paste("\nFound", length(inds_to_remove), "individuals with >", max_missing_individual * 100, "% missing data to remove."))
  # Get the sample names that are in the VCF header
  vcf_samples <- colnames(vcf@gt)
  # Keep only the samples that are NOT in our removal list
  inds_to_keep <- vcf_samples[!(vcf_samples %in% inds_to_remove)]
  vcf <- vcf[, inds_to_keep]
} else {
  message("No individuals exceeded the missing data threshold.")
}

# --- Log this step ---
filtering_summary <- filtering_summary %>% add_row(
  Step="Individual Call Rate",
  Filter_Applied=paste0(">", max_missing_individual * 100, "% missing"),
  Loci_Removed=0,
  Loci_Remaining=nrow(vcf),
  Individuals_Removed=inds_before - (ncol(vcf@gt) - 1),
  Individuals_Remaining=ncol(vcf@gt) - 1
)
message(paste("Filtering complete. Final number of individuals:", ncol(vcf@gt) - 1))
rm(missing_data_summary) # Clean up

### Step 4.3: Final Locus Call Rate Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 4.3: Final Locus Call Rate Filter ---")
loci_before <- nrow(vcf)

# --- Define Locus Missingness Threshold ---
max_missing_locus <- 0.20 # Remove loci with more than 20% missing data

# --- Calculate Missingness per Locus (Corrected Method) ---
message(paste("Calculating missing data for each of the", loci_before, "loci..."))

# THE FIX IS HERE: Apply is.na() directly to the genotype matrix (vcf@gt).
# We also exclude the first column, which is just the 'FORMAT' field.
genotype_matrix <- vcf@gt[, -1]
missing_per_locus <- rowSums(is.na(genotype_matrix))
n_individuals <- ncol(genotype_matrix)
prop_missing_per_locus <- missing_per_locus / n_individuals

# --- Filter Loci ---
# Create a logical vector: TRUE if a locus should be kept
loci_to_keep <- prop_missing_per_locus <= max_missing_locus

vcf <- vcf[loci_to_keep, ]

# --- Log this final step ---
filtering_summary <- filtering_summary %>% add_row(
  Step="Final Locus Call Rate",
  Filter_Applied=paste0("< ", max_missing_locus * 100, "% missing"),
  Loci_Removed=loci_before - nrow(vcf),
  Loci_Remaining=nrow(vcf),
  Individuals_Removed=0,
  Individuals_Remaining=ncol(vcf@gt) - 1
)
message(paste("Removed", sum(!loci_to_keep), "loci with >", max_missing_locus*100, "% missing data."))
message("\n==================================================")
message("FILTERING PIPELINE COMPLETE!")
message(paste("Final dataset contains", nrow(vcf), "loci for", ncol(vcf@gt)-1, "individuals."))
message("==================================================")
beep(sound = 8)

# --- Display Final Summary Table ---
print("--- FINAL FILTERING SUMMARY ---")
print(as.data.frame(filtering_summary))

# --- Clean up final variables ---
rm(genotype_matrix, loci_to_keep)

# ==============================================================================
# Part 5: Finalization & Export
# ==============================================================================

### Step 5.1: Final Data Sanity Check ----
# ------------------------------------------------------------------------------
message("\n--- Step 5.1: Final Sanity Check on Missing Data ---")

# Check missing data per individual
gt_matrix_final <- vcf@gt[, -1]
missing_per_indiv <- colSums(is.na(gt_matrix_final)) / nrow(gt_matrix_final)
message("\nSummary of MISSING data per INDIVIDUAL (should all be <= 0.40):")
print(summary(missing_per_indiv))

# Check missing data per locus
missing_per_locus <- rowSums(is.na(gt_matrix_final)) / ncol(gt_matrix_final)
message("\nSummary of MISSING data per LOCUS (should all be <= 0.20):")
print(summary(missing_per_locus))
rm(gt_matrix_final) # Clean up

### Step 5.2: Correct the Filtering Summary Table ----
# ------------------------------------------------------------------------------
message("\n--- Step 5.2: Correcting Final Summary Table Logic ---")
# This loop recalculates the 'Remaining' columns sequentially based on the 'Removed' columns.
for (i in 2:nrow(filtering_summary)) {
  filtering_summary$Loci_Remaining[i] <- filtering_summary$Loci_Remaining[i-1] - filtering_summary$Loci_Removed[i]
  filtering_summary$Individuals_Remaining[i] <- filtering_summary$Individuals_Remaining[i-1] - filtering_summary$Individuals_Removed[i]
}
message("The filtering summary table has been corrected.")
print("--- CORRECTED FINAL FILTERING SUMMARY ---")
print(as.data.frame(filtering_summary))


### Step 5.3: Export the Filtered VCF File ----
# ------------------------------------------------------------------------------
message("\n--- Step 5.3: Writing Final VCF to File ---")

# Create an 'Output' directory if it doesn't exist
if (!dir.exists("Output")) {
  dir.create("Output")
}

# Write the final VCF object to a compressed file
write.vcf(vcf, file = "Output/populations.filtered.final.vcf.gz")
message("Final VCF file has been written to: Output/populations.filtered.final.vcf.gz")
beep(sound = 2)


### Step 5.4: Generate Combined PDF Report of All Plots ----
# ------------------------------------------------------------------------------
message("\n--- Step 5.4: Generating PDF Report of All Diagnostic Plots ---")


# Convert the final summary table into a graphical object
p_summary_table <- tableGrob(filtering_summary, rows = NULL)

# Define the layout for the plots using patchwork syntax
# Page 1: Locus-level filtering decisions
page1_layout <- (p_raw_depth_by_year) / (p_hwe_volcano + p_paralog_plot) +
  plot_annotation(title = 'VCF Filtering Report: Locus-Level Decisions')

# Page 2: Genotype and Individual-level filtering decisions
page2_layout <- (p_dp_impact + p_gq_impact) / (p_missing_data_indivs) +
  plot_annotation(title = 'VCF Filtering Report: Genotype & Individual Decisions')


# Save the layouts to a PDF file
pdf("Output/VCF_Filtering_Report.pdf", width = 14, height = 11)
print(page1_layout)
print(page2_layout)
grid.arrange(p_summary_table)
dev.off()

message("All diagnostic plots have been saved to: Output/VCF_Filtering_Report.pdf")
message("\n--- SCRIPT COMPLETE ---")
beep(sound = 8)
