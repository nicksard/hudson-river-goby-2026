################################################################################
### VCF Filtering & Analysis Pipeline (Optimized & Final Version)
###
### Author: Nicholas Sard
### Assistant: Gemini
### Date: October 14, 2025
###
### Description: This script follows a loci-first filtering approach, optimized
### to remove low-quality and uninformative markers as early as possible to
### improve computational efficiency. It includes robust data visualization to
### justify all filtering decisions.
###
### Workflow:
### 1. Setup & Pre-filtering: Load data, analyze NTCs for contamination,
###    remove NTCs, and visualize raw read depth across sample years.
### 2. Locus Culling: Apply a Minor Allele Frequency (MAF) filter, then thin
###    to one SNP per RAD locus using a strategy that prioritizes the low-
###    coverage 2022 samples.
### 3. Locus Quality Control: On the reduced SNP set, perform rigorous
###  -- Hardy-Weinberg Equilibrium and paralog filtering.
### 4. Genotype & Individual Filtering: Mask low-quality individual genotypes,
###    then remove individuals and loci with high rates of missing data.
### 5. Finalization & Reporting: Perform final sanity checks, correct the
###    summary table, export the clean VCF, and generate a multi-page PDF
###    report summarizing the entire process.
################################################################################


# ==============================================================================
# Part 1: Initial Setup & Pre-filtering
# ==============================================================================

### Step 1.1: Load Libraries and Data ----
# ------------------------------------------------------------------------------
message("--- Part 1: Initial Setup & Pre-filtering ---")

# Loading necessary libraries
library(vcfR)
library(tidyverse)
library(adegenet)
library(hierfstat)
library(patchwork) # For combining ggplot objects
library(gridExtra) # For creating table objects
library(grid)      # For advanced PDF layout control
library(beepr)     # For audible notifications

# Set your working directory
setwd("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis")

# Read in the full VCF file
vcf_raw <- read.vcfR("Input//populations.snps.vcf")

# Rename the problematic NTC sample to ensure it's treated as a sample
message("Renaming sample 'NTC4a' to '2024999'.")
colnames(vcf_raw@gt)[colnames(vcf_raw@gt) == "NTC4a"] <- "2024999"


### Step 1.2: Negative Control (NTC) Contamination Analysis ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.2: Analyzing Negative Controls for Contamination ---")

# Identify all NTC samples from the raw data
ntc_samples <- colnames(vcf_raw@gt)[grep("NTC", colnames(vcf_raw@gt))]

# Extract allele depth (AD) for spurious genotypes found in NTCs
ntc_ad_data <- vcfR2tidy(vcf_raw, format_fields = c("GT", "AD"))$gt %>%
  filter(Indiv %in% ntc_samples) %>%
  filter(!is.na(gt_GT) & gt_GT != "./." & !is.na(gt_AD)) %>%
  separate(gt_AD, into = c("ad_ref", "ad_alt"), sep = ",") %>%
  mutate(
    ad_ref = as.numeric(ad_ref),
    ad_alt = as.numeric(ad_alt),
    total_depth = ad_ref + ad_alt
  ) %>%
  # Allele balance = proportion of reads supporting the alternate allele
  mutate(allele_balance = ad_alt / total_depth)

# Calculate key summary statistics for the plot subtitle
max_ntc_depth <- max(ntc_ad_data$total_depth, na.rm = TRUE)
percent_low_balance <- round(mean(ntc_ad_data$allele_balance < 0.25, na.rm = TRUE) * 100, 1)

# Create a diagnostic plot
p_ntc_allele_balance <- ggplot(ntc_ad_data, aes(x = allele_balance)) +
  geom_histogram(binwidth = 0.05, fill = "tomato", color = "black", boundary=0) +
  labs(
    title = "Allele Balance in Negative Controls (NTCs)",
    subtitle = paste0(
      "Max read depth in any spurious NTC genotype was ", max_ntc_depth, "x.\n",
      percent_low_balance, "% of spurious genotypes have allele balance < 0.25, suggesting index hopping."
    ),
    x = "Alternate Allele Balance (Alt Reads / Total Reads)",
    y = "Count of Spurious Genotypes"
  ) +
  theme_bw()

print(p_ntc_allele_balance)


### Step 1.3: Initialize Filtering & Remove NTCs ----
# ------------------------------------------------------------------------------
# We now create our main 'vcf' object for filtering by removing the NTCs
all_samples <- colnames(vcf_raw@gt)
samples_to_keep <- setdiff(all_samples, ntc_samples)
vcf <- vcf_raw[ , samples_to_keep]
rm(vcf_raw) # Clean up the raw object to save memory

# Initialize the filtering summary table
filtering_summary <- tibble(
  Step = "Initial State",
  Filter_Applied = "N/A",
  Loci_Removed = 0,
  Loci_Remaining = nrow(vcf),
  Individuals_Removed = 0,
  Individuals_Remaining = ncol(vcf@gt) - 1 + length(ntc_samples) # Add NTCs back for the initial count
)

# Log the NTC removal
filtering_summary <- filtering_summary %>% add_row(
  Step = "NTC Removal",
  Filter_Applied = "Removed all NTC samples",
  Loci_Removed = 0,
  Loci_Remaining = nrow(vcf),
  Individuals_Removed = length(ntc_samples),
  Individuals_Remaining = ncol(vcf@gt) - 1
)
message(paste("\nRemoved", length(ntc_samples), "NTCs. Kept", ncol(vcf@gt)-1, "samples for analysis."))
beep(sound = 2)


### Step 1.4: Visualize Raw Read Depth ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.4: Visualizing Raw Read Depth By Year ---")
raw_dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
raw_pop_map <- tibble(Sample = colnames(raw_dp_matrix), Population = substr(colnames(raw_dp_matrix), 1, 4))
raw_dp_long <- as.data.frame(raw_dp_matrix) %>%
  rownames_to_column("Locus") %>%
  pivot_longer(-Locus, names_to="Sample", values_to="Depth") %>%
  left_join(raw_pop_map, by="Sample")

p_raw_depth_by_year <- ggplot(raw_dp_long, aes(x=Population, y=Depth, fill=Population)) +
  geom_boxplot(outlier.shape=NA) +
  coord_cartesian(ylim=c(0, 100)) +
  labs(
    title="Raw Read Depth by Collection Year",
    subtitle="Confirms that 2022 samples have lower coverage.",
    x="Year", y="Read Depth (DP)"
  ) +
  theme_bw() + theme(legend.position="none")

print(p_raw_depth_by_year)
rm(raw_dp_matrix, raw_pop_map, raw_dp_long) # Clean up
beep(sound = 2)


# ==============================================================================
# Part 2: Locus Culling (Major Efficiency Steps)
# ==============================================================================
message("\n--- Part 2: Culling Loci for Efficiency ---")

### Step 2.1: Minor Allele Frequency (MAF) Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.1: Minor Allele Frequency (MAF) Filter ---")
loci_before <- nrow(vcf)
maf_cutoff <- 0.05

maf_data <- as_tibble(vcf@fix) %>%
  mutate(Info_AF = as.numeric(str_match(INFO, "AF=([0-9.]+)")[,2])) %>%
  mutate(MAF = ifelse(Info_AF > 0.5, 1 - Info_AF, Info_AF))
head(maf_data)
hist(maf_data$MAF)
loci_to_keep <- maf_data$MAF >= maf_cutoff & !is.na(maf_data$MAF)
vcf <- vcf[loci_to_keep, ]

filtering_summary <- filtering_summary %>% add_row(
  Step = "MAF Filter",
  Filter_Applied = paste("MAF >=", maf_cutoff),
  Loci_Removed = loci_before - nrow(vcf),
  Loci_Remaining = nrow(vcf), Individuals_Removed = 0,
  Individuals_Remaining = ncol(vcf@gt) - 1
)
filtering_summary
message(paste("Removed loci with MAF <", maf_cutoff, ". Loci remaining:", nrow(vcf)))
rm(maf_data) # Clean up
beep(sound = 2)


### Step 2.2: Thin to One SNP per RAD Locus ----
# ------------------------------------------------------------------------------
message("\n--- Step 2.2: Thinning to One SNP per RAD Locus ---")
message("Strategy: Prioritize SNPs with the best data in low-coverage 2022 samples.")
loci_before <- nrow(vcf)

tidy_gt <- vcfR2tidy(vcf, format_fields = c("GT"))$gt
snp_completeness_2022 <- tidy_gt %>%
  filter(substr(Indiv, 1, 4) == "2022") %>%
  group_by(ChromKey) %>% summarise(n_genotyped_2022 = sum(!is.na(gt_GT) & gt_GT != "./."))
snp_completeness_overall <- tidy_gt %>%
  group_by(ChromKey) %>% summarise(n_genotyped_overall = sum(!is.na(gt_GT) & gt_GT != "./."))

thinned_indices <- as_tibble(vcf@fix) %>%
  mutate(row_index = 1:n()) %>%
  separate(ID, into = c("locus_id"), sep = ":", remove = FALSE, extra = "drop") %>%
  left_join(snp_completeness_overall, by = c("row_index" = "ChromKey")) %>%
  left_join(snp_completeness_2022, by = c("row_index" = "ChromKey")) %>%
  mutate(n_genotyped_2022 = ifelse(is.na(n_genotyped_2022), 0, n_genotyped_2022)) %>%
  group_by(locus_id) %>%
  arrange(desc(n_genotyped_2022), desc(n_genotyped_overall), .by_group = TRUE) %>%
  slice(1) %>% ungroup() %>% pull(row_index)
head(thinned_indices)
vcf <- vcf[thinned_indices, ]

filtering_summary <- filtering_summary %>% add_row(
  Step = "Thinning",
  Filter_Applied = "One SNP per locus (2022 priority)",
  Loci_Removed = loci_before - nrow(vcf),
  Loci_Remaining = nrow(vcf), Individuals_Removed = 0,
  Individuals_Remaining = ncol(vcf@gt) - 1
)
filtering_summary
message(paste("Thinned to one SNP per locus. Loci remaining:", nrow(vcf)))
rm(tidy_gt, snp_completeness_2022, snp_completeness_overall, thinned_indices) # Clean up
beep(sound = 2)


# ==============================================================================
# Part 3: Locus Quality Control
# ==============================================================================
message("\n--- Part 3: Locus Quality Control ---")

### Step 3.1: Hardy-Weinberg Equilibrium Filter (via Fis) ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.1: Hardy-Weinberg Equilibrium (HWE) Analysis ---")
loci_before <- nrow(vcf)
final_fis_cutoff <- 0.25

message("Calculating Fis and HWE p-values...")
temp_genind <- vcfR2genind(vcf)
pop(temp_genind) <- rep("Overall", nInd(temp_genind))
hierf_data <- genind2hierfstat(temp_genind)
fis_results <- hierfstat::basic.stats(hierf_data, diploid = TRUE)$perloc %>%
  rownames_to_column("Locus") %>% as_tibble()
fis_results$Locus <- locNames(temp_genind)
hwe_results <- pegas::hw.test(temp_genind, B = 0)
hwe_p_values <- as_tibble(hwe_results, rownames = "Locus") %>% select(Locus, p_val = `Pr(chi^2 >)`)
hwe_summary <- left_join(fis_results, hwe_p_values, by = "Locus")

p_hwe_volcano <- ggplot(hwe_summary, aes(x = Fis, y = -log10(p_val))) +
  geom_point(alpha = 0.4, color = "midnightblue") +
  geom_vline(xintercept = c(-final_fis_cutoff, final_fis_cutoff), color = "red", linetype = "dashed") +
  labs(title = "HWE Deviation", x = "Inbreeding Coefficient (Fis)", y = "-log10(p-value)") +
  theme_bw()
print(p_hwe_volcano)

hwe_summary %>%
  filter(abs(Fis) > final_fis_cutoff)  %>%
  count()

hwe_summary %>%
  filter(p_val <= (0.05/nrow(hwe_summary)))  %>%
  count()

table(hwe_summary$p_val <=(0.05/nrow(hwe_summary)))

hwe_summary %>%
  filter(abs(Fis) >= final_fis_cutoff) %>%
  pull(p_val) < (0.05/12032)) %>%
  count()
loci_to_remove_hwe <- hwe_summary %>% filter(abs(Fis) > final_fis_cutoff & !is.na(Fis)) %>% pull(Locus)
vcf <- vcf[!(locNames(temp_genind) %in% loci_to_remove_hwe), ]

filtering_summary <- filtering_summary %>% add_row(
  Step = "HWE Filter (Fis)",
  Filter_Applied = paste("|Fis| <=", final_fis_cutoff),
  Loci_Removed = loci_before - nrow(vcf), Loci_Remaining = nrow(vcf),
  Individuals_Removed = 0, Individuals_Remaining = ncol(vcf@gt) - 1
)
filtering_summary
-log10(0.05/12032) #bonferroni cut-off for not deviating HWP

message(paste("Removed loci with |Fis| >", final_fis_cutoff, ". Loci remaining:", nrow(vcf)))
rm(temp_genind, hierf_data, fis_results, hwe_results, hwe_p_values, hwe_summary) # Clean up
beep(sound = 2)


### Step 3.2: Paralog Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 3.2: Paralog Filter ---")
loci_before <- nrow(vcf)

message("Calculating Hobs and read-ratio deviation to identify paralogs...")
paralog_tidy <- vcfR2tidy(vcf, format_fields = c("GT", "AD"))$gt %>%
  rename(gt = gt_GT, allele_dp = gt_AD, locus_id = ChromKey) %>%
  separate(col = allele_dp, into = c("ref_dp", "alt_dp"), sep = ",", fill="right")
hobs_calc <- paralog_tidy %>% filter(!is.na(gt) & gt != "./.") %>% group_by(locus_id, gt) %>%
  summarise(n = n(), .groups = "drop_last") %>% pivot_wider(names_from = gt, values_from = n, values_fill = 0) %>%
  mutate(n_total = `0/0` + `0/1` + `1/1`, Hobs = `0/1` / n_total) %>% select(locus_id, Hobs)
z_score_calc <- paralog_tidy %>% filter(gt == "0/1", !is.na(ref_dp), !is.na(alt_dp)) %>%
  mutate(ref_dp = as.numeric(ref_dp), alt_dp = as.numeric(alt_dp)) %>% group_by(locus_id) %>%
  summarise(total_reads = sum(ref_dp) + sum(alt_dp), z = (sum(ref_dp) - (total_reads / 2)) / sqrt(total_reads * 0.25), .groups = "drop") %>%
  select(locus_id, z)
paralog_summary <- left_join(hobs_calc, z_score_calc, by = "locus_id") %>%
  mutate(is_singleton = ifelse(Hobs <= 0.6 & abs(z) <= 7 & !is.na(z), "Yes", "No"))

p_paralog_plot <- ggplot(paralog_summary, aes(x=Hobs, y=z, color=is_singleton)) +
  geom_point(alpha=0.6) + labs(title="Paralog Identification", x="Observed Heterozygosity (Hobs)", y="Read-Ratio Deviation (z-score)", color="Putative Singleton") +
  theme_bw()
print(p_paralog_plot)

loci_to_remove_paralogs <- paralog_summary %>% filter(is_singleton == "No") %>% pull(locus_id)
vcf <- vcf[!( (1:nrow(vcf)) %in% loci_to_remove_paralogs), ]

filtering_summary <- filtering_summary %>% add_row(
  Step = "Paralog Filter",
  Filter_Applied = "Hobs > 0.6 or |Z| > 7",
  Loci_Removed = loci_before - nrow(vcf), Loci_Remaining = nrow(vcf),
  Individuals_Removed = 0, Individuals_Remaining = ncol(vcf@gt) - 1
)
message(paste("Removed", length(loci_to_remove_paralogs), "putative paralogs. Loci remaining:", nrow(vcf)))
rm(paralog_tidy, hobs_calc, z_score_calc, paralog_summary) # Clean up
beep(sound = 2)


# ==============================================================================
# Part 4: Genotype & Individual Filtering
# ==============================================================================
message("\n--- Part 4: Genotype & Individual Filtering ---")

### Step 4.1: Justifying and Applying Genotype Filters ----
# ------------------------------------------------------------------------------
message("\n--- Step 4.1: Visualizing Filter Thresholds & Masking Genotypes ---")
dp_matrix <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
gq_matrix <- extract.gt(vcf, element = "GQ", as.numeric = TRUE)
total_genotypes <- sum(!is.na(vcf@gt) & vcf@gt != "./.")

dp_loss_summary <- map_dfr(seq(0, 50, 5), ~tibble(threshold=.x, percent_masked=(sum(dp_matrix<.x, na.rm=T)/total_genotypes)*100))
p_dp_impact <- ggplot(dp_loss_summary, aes(x = threshold, y = percent_masked)) + geom_line(linewidth=1) + geom_point(size=3) +
  geom_vline(xintercept=10, color="red", linetype="dashed", linewidth=1) + labs(title="Impact of Min. Read Depth (DP) Filter", x="Minimum DP Threshold", y="Cumulative % of Genotypes Masked") + theme_bw()
print(p_dp_impact)

gq_loss_summary <- map_dfr(seq(0, 40, 5), ~tibble(threshold=.x, percent_masked=(sum(gq_matrix<.x, na.rm=T)/total_genotypes)*100))
p_gq_impact <- ggplot(gq_loss_summary, aes(x = threshold, y = percent_masked)) + geom_line(linewidth=1) + geom_point(size=3) +
  geom_vline(xintercept=20, color="red", linetype="dashed", linewidth=1) + labs(title="Impact of Min. Genotype Quality (GQ) Filter", x="Minimum GQ Threshold", y="Cumulative % of Genotypes Masked") + theme_bw()
print(p_gq_impact)

min_depth_threshold <- 10; min_gq_threshold <- 20
mask <- (dp_matrix < min_depth_threshold | gq_matrix < min_gq_threshold); mask[is.na(mask)] <- FALSE
vcf@gt[mask] <- "./."
message(paste("\nMasking", sum(mask), "genotypes with DP <", min_depth_threshold, "OR GQ <", min_gq_threshold))
rm(dp_matrix, gq_matrix, mask) # Clean up


### Step 4.2: Individual Call Rate Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 4.2: Filtering Individuals by Missing Data ---")
inds_before <- ncol(vcf@gt) - 1
max_missing_individual <- 0.40 # Stricter cutoff for pedigree analysis

missing_data_summary <- vcfR2tidy(vcf, format_fields = c("GT"))$gt %>%
  group_by(Indiv) %>% summarise(prop_missing = sum(is.na(gt_GT) | gt_GT == "./.") / n(), .groups="drop") %>%
  mutate(Population = substr(Indiv, 1, 4))

p_missing_data_indivs <- ggplot(missing_data_summary, aes(x=Population, y=prop_missing, fill=Population)) +
  geom_boxplot(outlier.alpha=0.5) + geom_hline(yintercept=max_missing_individual, color="red", linetype="dashed", linewidth=1) +
  labs(title="Individual Missing Data", subtitle=paste("Cutoff at", max_missing_individual*100, "% missing data"), x="Year", y="Proportion of Missing Genotypes") +
  theme_bw() + theme(legend.position="none")
print(p_missing_data_indivs)

inds_to_remove <- missing_data_summary %>% filter(prop_missing > max_missing_individual) %>% pull(Indiv)
vcf <- vcf[, !(colnames(vcf@gt) %in% inds_to_remove)]

filtering_summary <- filtering_summary %>% add_row(
  Step="Individual Call Rate", Filter_Applied=paste0(">", max_missing_individual*100, "% missing"), Loci_Removed=0,
  Loci_Remaining=nrow(vcf), Individuals_Removed=inds_before - (ncol(vcf@gt)-1), Individuals_Remaining=ncol(vcf@gt)-1
)
message(paste("Filtering complete. Final number of individuals:", ncol(vcf@gt)-1))
rm(missing_data_summary) # Clean up


### Step 4.3: Final Locus Call Rate Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 4.3: Final Locus Call Rate Filter ---")
loci_before <- nrow(vcf)
max_missing_locus <- 0.20

genotype_matrix <- vcf@gt[, -1]
# LOGIC FIX: Must count both NA AND "./." as missing
prop_missing_per_locus <- rowSums(is.na(genotype_matrix) | genotype_matrix == "./.") / ncol(genotype_matrix)
vcf <- vcf[prop_missing_per_locus <= max_missing_locus, ]

filtering_summary <- filtering_summary %>% add_row(
  Step="Final Locus Call Rate", Filter_Applied=paste0("< ", max_missing_locus*100, "% missing"),
  Loci_Removed=loci_before - nrow(vcf), Loci_Remaining=nrow(vcf),
  Individuals_Removed=0, Individuals_Remaining=ncol(vcf@gt)-1
)
message(paste("Removed", loci_before - nrow(vcf), "loci with >", max_missing_locus*100, "% missing data."))
message("\n==================================================")
message("FILTERING PIPELINE COMPLETE!")
message(paste("Final dataset contains", nrow(vcf), "loci for", ncol(vcf@gt)-1, "individuals."))
message("==================================================")
beep(sound = 2)

### Step 4.4: Final MAF & Monomorphic Loci Filter ----
# ------------------------------------------------------------------------------
message("\n--- Step 4.4: Final Monomorphic Loci Filter ---")
message("Removing loci that became monomorphic (MAF=0) after individual filtering.")
loci_before <- nrow(vcf)

# Calculate MAF on the final filtered set of individuals
# This returns a matrix (as we saw in the debug output)
final_maf_stats <- vcfR::maf(vcf)

# --- THIS IS THE FIX ---
# We must first extract the actual MAF values from the "Frequency" column
# The object is a matrix, so we subset by [row, column]
# We want all rows, and the "Frequency" column.
maf_column_values <- final_maf_stats[, "Frequency"]

# Now, create a logical vector (TRUE/FALSE) just from that column
logical_indices_to_keep <- maf_column_values > 0

# Treat any locus where MAF could not be calculated (NA) as a locus to remove.
logical_indices_to_keep[is.na(logical_indices_to_keep)] <- FALSE

# Filter the VCF object using the clean logical vector (now the correct length)
vcf <- vcf[logical_indices_to_keep, ]

# Log this final step
filtering_summary <- filtering_summary %>% add_row(
  Step="Final MAF Cleanup", 
  Filter_Applied="Removed monomorphic (MAF=0) or NA loci",
  Loci_Removed=loci_before - nrow(vcf), 
  Loci_Remaining=nrow(vcf),
  Individuals_Removed=0, 
  Individuals_Remaining=ncol(vcf@gt)-1
)
filtering_summary
message(paste("Removed", loci_before - nrow(vcf), "loci that became monomorphic or had NA MAF."))
message(paste("Final loci remaining:", nrow(vcf)))

# ==============================================================================
# Part 5: Finalization & Reporting
# ==============================================================================
message("\n--- Part 5: Finalization & Reporting ---")

### Step 5.1: Final Data Sanity Check ----
# ------------------------------------------------------------------------------
message("\n--- Step 5.1: Final Sanity Check on Missing Data ---")
gt_matrix_final <- vcf@gt[, -1]

missing_per_indiv <- colSums(is.na(gt_matrix_final) | gt_matrix_final == "./.") / nrow(gt_matrix_final)
missing_per_locus <- rowSums(is.na(gt_matrix_final) | gt_matrix_final == "./.") / ncol(gt_matrix_final)

# 1. Calculate underlying missingness logical matrix (TRUE = missing)
is_missing_mat <- is.na(gt_matrix_final) | gt_matrix_final == "./."

# Determine the collection year for each individual column
sample_years <- substr(colnames(gt_matrix_final), 1, 4)

# 2. Print summaries to the console for live sanity checking
message("\nSummary of MISSING data per INDIVIDUAL (should all be <= 0.40):")
print(summary(missing_per_indiv))

message("\nSummary of MISSING data per LOCUS (should all be <= 0.20):")
print(summary(missing_per_locus))

# 3. Structure data into tables for file export
indiv_missing_df <- tibble(
  Individual = names(missing_per_indiv),
  Prop_Missing = missing_per_indiv,
  Year = substr(Individual, 1, 4)) %>%
  group_by(Year) %>%
  summarise(
    Type = "Individual",
    Mean = mean(Prop_Missing, na.rm = TRUE),
    Median = median(Prop_Missing, na.rm = TRUE),
    SD = sd(Prop_Missing, na.rm = TRUE),
    Min = min(Prop_Missing, na.rm = TRUE),
    Max = max(Prop_Missing, na.rm = TRUE),
    .groups = "drop"
  )
indiv_missing_df

# We calculate missingness per locus within each year group independently
unique_years <- unique(sample_years)

locus_summary_by_year <- map_dfr(unique_years, function(yr) {
  # Subset matrix columns to only include individuals from the target year
  year_sub_mat <- is_missing_mat[, sample_years == yr, drop = FALSE]
  
  # Calculate proportion missing per locus for this specific year
  locus_missing_yr <- rowSums(year_sub_mat) / ncol(year_sub_mat)
  
  # Return a single summary row for this year
  tibble(
    Year = yr,
    Type = "Locus",
    Mean = mean(locus_missing_yr, na.rm = TRUE),
    Median = median(locus_missing_yr, na.rm = TRUE),
    SD = sd(locus_missing_yr, na.rm = TRUE),
    Min = min(locus_missing_yr, na.rm = TRUE),
    Max = max(locus_missing_yr, na.rm = TRUE)
  )
})
locus_summary_by_year


# 4. Save information to files
write_csv(indiv_missing_df, "Output/final_missing_data_per_individual_summary.csv")
write_csv(locus_summary_by_year, "Output/final_missing_data_summary_per_locus.csv")

rm(gt_matrix_final) # Clean up


### Step 5.2: Correct the Filtering Summary Table ----
# ------------------------------------------------------------------------------
message("\n--- Step 5.2: Correcting Final Summary Table Logic ---")
for (i in 2:nrow(filtering_summary)) {
  filtering_summary$Loci_Remaining[i] <- filtering_summary$Loci_Remaining[i-1] - filtering_summary$Loci_Removed[i]
  filtering_summary$Individuals_Remaining[i] <- filtering_summary$Individuals_Remaining[i-1] - filtering_summary$Individuals_Removed[i]
}
message("The filtering summary table has been corrected.")
print("--- CORRECTED FINAL FILTERING SUMMARY ---"); print(as.data.frame(filtering_summary))

### Step 5.3: Export the Filtered VCF File ----
# ------------------------------------------------------------------------------
message("\n--- Step 5.3: Writing Final VCF to File ---")
if (!dir.exists("Output")) { dir.create("Output") }
write.vcf(vcf, file = "Output/populations.filtered.final.vcf.gz")
message("Final VCF file has been written to: Output/populations.filtered.final.vcf.gz")

### Step 5.4: Generate Combined PDF Report of All Plots ----
# ------------------------------------------------------------------------------
message("\n--- Step 5.4: Generating PDF Report of All Diagnostic Plots ---")
p_summary_table <- tableGrob(filtering_summary, rows = NULL)

page1_layout <- (p_raw_depth_by_year) / (p_hwe_volcano + p_paralog_plot) +
  plot_annotation(title = 'VCF Filtering Report: Locus-Level Decisions')

page2_layout <- (p_dp_impact + p_gq_impact) / (p_ntc_allele_balance + p_missing_data_indivs) +
  plot_annotation(title = 'VCF Filtering Report: Genotype & Individual Decisions')

pdf("Output/VCF_Filtering_Report.pdf", width = 14, height = 11)
print(page1_layout)
print(page2_layout)
grid.arrange(p_summary_table)
dev.off()

message("All diagnostic plots have been saved to: Output/VCF_Filtering_Report.pdf")
message("\n--- SCRIPT COMPLETE ---")
beep(sound = 8)

gt <- vcf@gt[, -1]
print(paste("Loci with 100% missing data:", sum(rowSums(is.na(gt) | gt == "./.") == ncol(gt))))
print(paste("Max proportion missing per locus:", max(rowSums(is.na(gt) | gt == "./.") / ncol(gt))))
print(paste("Max proportion missing per individual:", max(colSums(is.na(gt) | gt == "./.") / nrow(gt))))
print(paste("Are all thinned 'locus_ids' unique?", !anyDuplicated(separate(as_tibble(vcf@fix), ID, into = c("locus_id"), sep = ":", extra = "drop")$locus_id)))
maf_final <- vcfR::maf(vcf)
print(paste("Min MAF in final VCF:", min(maf_final[, "Frequency"], na.rm = TRUE)))

