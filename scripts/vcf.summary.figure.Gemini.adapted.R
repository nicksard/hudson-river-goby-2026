# originally created on: October 18th, 2017
# by: Nick Sard
# Adapted on: October 7, 2025 to use the vcfR package for robust VCF parsing.

# ABOUT: This script was written to evaluate GT data associated with a VCF

# loading necessary libraries
library(vcfR)
library(tidyverse)
library(gridExtra)
library(grid)

# setting working directory and reading in data
# Make sure to set this to your own working directory
# setwd("G:/Shared drives/Sard Lab/Round Goby/RADseq/R/post_RAD analysis")

# listing files
# list.files("Input/")

######
### NOTE: The manual step of deleting the '#' before CHROM is no longer needed. vcfR handles this automatically. ###
######

### ADAPTED: Reading in the VCF file using vcfR
# This creates a special vcfR object that contains meta, fix, and gt data
vcf <- read.vcfR("Input//populations.snps.vcf")

# ADAPTED: Rename 'NTC4a' to 'Sample_X' as it is a real sample
message("Renaming sample 'NTC4a' to 'Sample_X'.")
colnames(vcf@gt)[colnames(vcf@gt) == "NTC4a"] <- "Sample_X"

### --- SEQUENTIAL FILTERING PIPELINE --- ###
message(paste("Original number of loci:", nrow(vcf)))

# --- Step 1: Negative Control (NTC) Analysis and Filtering ---
message("\n--- Step 1: Negative Control (NTC) Analysis ---")

# ADAPTED: Identify true NTC samples. 'NTC4a' is no longer present.
ntc_samples_to_check <- colnames(vcf@gt)[grep("NTC", colnames(vcf@gt))]
message(paste("Identified", length(ntc_samples_to_check), "true NTCs for contamination check."))

# Get genotype AND allele depth data for all samples before any filtering
tidy_vcf_unfiltered <- vcfR2tidy(vcf, format_fields = c("GT", "AD"))
geno_unfiltered <- tidy_vcf_unfiltered$gt

# Calculate the proportion of genotyped loci for each NTC
ntc_geno_summary <- geno_unfiltered %>%
  filter(Indiv %in% ntc_samples_to_check) %>%
  mutate(genotyped = ifelse(is.na(gt_GT) | gt_GT == "./.", 0, 1)) %>%
  group_by(Indiv) %>%
  summarise(
    n_loci_total = n(),
    n_loci_genotyped = sum(genotyped),
    prop_genotyped = round(n_loci_genotyped / n_loci_total, 4)
  )

message("Summary of genotyping in TRUE Negative Controls:")
print(as.data.frame(ntc_geno_summary))

# Identify all unique loci (by their row index, which is ChromKey) that have a genotype in any TRUE NTC
loci_with_ntc_genotypes <- geno_unfiltered %>%
  filter(Indiv %in% ntc_samples_to_check) %>%
  filter(!is.na(gt_GT) & gt_GT != "./.")

loci_to_remove_ntc <- unique(loci_with_ntc_genotypes$ChromKey)

### --- NEW: Deeper Dive into NTC Genotype Quality --- ###
# To justify removing these loci, we'll check their allele balance.
# Barcode spillover should result in highly unbalanced allele depths (e.g., AD near 0 or 1).
ntc_ad_data <- loci_with_ntc_genotypes %>%
  filter(!is.na(gt_AD)) %>%
  separate(gt_AD, into = c("ad_ref", "ad_alt"), sep = ",") %>%
  mutate(ad_ref = as.numeric(ad_ref), ad_alt = as.numeric(ad_alt)) %>%
  mutate(total_depth = ad_ref + ad_alt) %>%
  # Calculate allele balance for the alternate allele
  mutate(allele_balance = ad_alt / total_depth)

# Plot the distribution of allele balances in the NTCs
p_ntc_allele_balance <- ggplot(ntc_ad_data, aes(x = allele_balance)) +
  geom_histogram(binwidth = 0.05, fill = "tomato", color = "black") +
  labs(title = "Allele Balance for Genotypes in Negative Controls",
       x = "Alternate Allele Balance (AD_alt / Total Depth)",
       y = "Count of Spurious Genotypes",
       subtitle = "Values near 0 or 1 suggest barcode spillover.") +
  theme_bw()
p_ntc_allele_balance
### --- END OF NTC DIVE --- ###


### --- NTC FILTERING (DEACTIVATED) --- ###
# Based on the allele balance plot, the NTC genotypes are low-quality noise.
# This noise is unlikely to affect high-coverage genotypes in real samples.
# Therefore, we are now proceeding WITHOUT removing these loci to be less stringent.
# To reactivate, simply uncomment the following lines.
#
# message(paste("\nFound", length(loci_to_remove_ntc), "loci with genotype calls in NTCs. These will be removed."))
# if(length(loci_to_remove_ntc) > 0){
#     original_indices <- 1:nrow(vcf)
#     ntc_indices_to_keep <- setdiff(original_indices, loci_to_remove_ntc)
#     vcf <- vcf[ntc_indices_to_keep, ]
# }
# message(paste("Number of loci remaining after NTC filter:", nrow(vcf)))


# Now, remove the true NTC samples themselves from the VCF object so they are not in downstream analyses
all_samples <- colnames(vcf@gt)
samples_to_keep <- setdiff(all_samples, ntc_samples_to_check)
vcf <- vcf[, samples_to_keep]
message(paste("\nRemoved", length(ntc_samples_to_check), "NTC samples. Kept", ncol(vcf@gt)-1, "samples for analysis."))


# --- Step 2: MAF Filtering ---
message("\n--- Step 2: MAF Filtering ---")

# Create a data frame from the NTC-filtered data to calculate MAF
maf_data_before_filter <- as_tibble(vcf@fix) %>%
  mutate(Info_AF = as.numeric(str_remove(str_extract(INFO, "AF=[0-9.,]+"), "AF="))) %>%
  mutate(MAF = ifelse(Info_AF > 0.5, 1 - Info_AF, Info_AF)) %>%
  mutate(Filter_Step = "Before MAF Filter")

# Identify loci to keep based on MAF >= 0.05
maf_indices_to_keep <- which(maf_data_before_filter$MAF >= 0.05 & !is.na(maf_data_before_filter$MAF))
vcf <- vcf[maf_indices_to_keep, ]
message(paste("Number of loci remaining after MAF >= 0.05 filter:", nrow(vcf)))

# Create a data frame for the data AFTER MAF filtering
maf_data_after_filter <- as_tibble(vcf@fix) %>%
  mutate(Info_AF = as.numeric(str_remove(str_extract(INFO, "AF=[0-9.,]+"), "AF="))) %>%
  mutate(MAF = ifelse(Info_AF > 0.5, 1 - Info_AF, Info_AF)) %>%
  mutate(Filter_Step = "After MAF Filter")

# Combine for plotting
combined_maf_data <- bind_rows(maf_data_before_filter, maf_data_after_filter)
p_maf_comparison <- ggplot(combined_maf_data, aes(x = Filter_Step, y = MAF, fill = Filter_Step)) +
  geom_boxplot() +
  labs(title = "MAF Distribution Before and After Filtering", x = "Filtering Step", y = "Minor Allele Frequency") +
  theme_bw() +
  theme(legend.position = "none")
p_maf_comparison

# --- Step 3: Thin to one SNP per RAD locus (Maximizing MAF) ---
message("\n--- Step 3: Thinning to one SNP per RAD locus ---")

# First, visualize the "before" state (SNPs per locus after MAF filter, but before thinning)
temp_fix_before_thinning <- as_tibble(vcf@fix) %>%
  separate(ID, into = c("locus_id", "snp_pos_in_read", "strand"), sep = ":", remove = FALSE)
snp_counts_before <- temp_fix_before_thinning %>% count(locus_id, name = "snps_per_locus")
p_snps_per_locus_before <- ggplot(snp_counts_before, aes(x = snps_per_locus)) +
  geom_histogram(binwidth = 1, fill = "darkseagreen", color = "black") +
  labs(x = "Number of SNPs per RAD Locus", y = "Count of Loci", title = "SNPs per Locus (Before Thinning)") +
  theme_bw()
p_snps_per_locus_before

# Now, perform the thinning by selecting the SNP with the highest MAF for each locus
temp_fix_maf_filtered <- as_tibble(vcf@fix) %>%
  mutate(row_index = 1:n()) %>%
  separate(ID, into = c("locus_id", "snp_pos_in_read", "strand"), sep = ":", remove = FALSE) %>%
  mutate(Info_AF = as.numeric(str_remove(str_extract(INFO, "AF=[0-9.,]+"), "AF="))) %>%
  mutate(MAF = ifelse(Info_AF > 0.5, 1 - Info_AF, Info_AF))

thinned_indices <- temp_fix_maf_filtered %>%
  group_by(locus_id) %>%
  slice_max(order_by = MAF, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pull(row_index)
vcf <- vcf[thinned_indices, ]
message(paste("Number of loci remaining after thinning to one SNP/locus:", nrow(vcf)))


# --- Step 4: Prune for Linkage Disequilibrium (Physical Distance) ---
message("\n--- Step 4: Pruning SNPs for linkage (within 1000 bp) ---")
temp_fix_thinned <- as_tibble(vcf@fix) %>%
  mutate(POS = as.numeric(POS)) %>%
  arrange(CHROM, POS)

indices_to_keep_ld <- c()
if (nrow(temp_fix_thinned) > 0) {
  indices_to_keep_ld <- c(1)
  last_kept_index <- 1
  for (i in 2:nrow(temp_fix_thinned)) {
    if (temp_fix_thinned$CHROM[i] != temp_fix_thinned$CHROM[last_kept_index]) {
      indices_to_keep_ld <- c(indices_to_keep_ld, i)
      last_kept_index <- i
    } else {
      distance <- temp_fix_thinned$POS[i] - temp_fix_thinned$POS[last_kept_index]
      if (distance >= 1000) {
        indices_to_keep_ld <- c(indices_to_keep_ld, i)
        last_kept_index <- i
      }
    }
  }
}
vcf <- vcf[indices_to_keep_ld, ]
message(paste("Final number of loci after LD pruning:", nrow(vcf)))
### --- END OF FILTERING --- ###


### --- FINAL DATA PREPARATION FOR PLOTTING --- ###
message("\n--- Generating final plots from filtered data ---")
tidy_vcf <- vcfR2tidy(vcf, format_fields = c("GT", "DP", "AD"))

fixed_data <- as_tibble(vcf@fix) %>%
  mutate(Info_AF = as.numeric(str_remove(str_extract(INFO, "AF=[0-9.,]+"), "AF="))) %>%
  mutate(MAF = ifelse(Info_AF > 0.5, 1 - Info_AF, Info_AF)) %>%
  separate(ID, into = c("locus_id", "snp_pos_in_read", "strand"), sep = ":", remove = FALSE) %>%
  mutate(snp_pos_in_read = as.numeric(snp_pos_in_read)) %>%
  mutate(ID = 1:n())

p_maf <- ggplot(fixed_data, aes(x=MAF)) +
  geom_histogram(binwidth = 0.01, fill="steelblue", color="black") +
  labs(x="Minor Allele Frequency", y="Number of loci", title="Distribution of MAF (Post-Filtering)") +
  theme_bw()
p_maf
p1 <- ggplot(fixed_data, aes(x=snp_pos_in_read)) +
  geom_histogram(binwidth = 1, fill = "cornflowerblue", color = "black") +
  labs(x="SNP Position Within RAD Tag", y="Number of Loci", title = "SNP Positions Along Reads (Post-Filtering)") +
  theme_bw()
p1

snp_counts_after <- fixed_data %>% count(locus_id, name = "snps_per_locus")
#p_snps_per_locus_after <- ggplot(snp_counts_after, aes(x = snps_per_locus)) +
#  geom_histogram(binwidth = 1, fill = "darkseagreen", color = "black") +
#  labs(x = "Number of SNPs per RAD Locus", y = "Count of Loci", title = "SNPs per Locus (After Thinning)") +
#  scale_x_continuous(breaks=c(1)) +
#  theme_bw()
#p_snps_per_locus_after
geno <- tidy_vcf$gt
geno <- geno %>%
  rename(id = Indiv, gt = gt_GT, dp = gt_DP, allele_dp = gt_AD, ID = ChromKey) %>%
  mutate(loc = substr(x = id, start = 3, stop = 5)) %>%
  separate(col = allele_dp, into = c("allele_a", "allele_b"), sep = ",", fill="right")

geno$gt <- gsub(pattern = "\\./\\.", replacement = NA, x = geno$gt)
geno$gt <- gsub(pattern = "0/0", replacement = "AA", x = geno$gt)
geno$gt <- gsub(pattern = "0/1", replacement = "AB", x = geno$gt)
geno$gt <- gsub(pattern = "1/1", replacement = "BB", x = geno$gt)
geno$allele_a <- gsub(pattern = "\\.", replacement = NA, x = geno$allele_a)
geno$allele_b <- gsub(pattern = "\\.", replacement = NA, x = geno$allele_b)
geno <- geno %>% mutate(het = ifelse(test = gt == "AB", 1, 0))
geno$test <- ifelse(is.na(geno$gt), "No", "Yes")

p2 <- ggplot(geno, aes(x=id, y=ID, fill=gt)) +
  facet_wrap(~loc, nrow = 1, scales = "free_x") +
  geom_tile() +
  labs(x= "Samples", y="Loci", fill = "Genotyped") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), strip.text = element_text(size=8), legend.position = "bottom")
p2

loc_count <- geno %>%
  group_by(loc, id, test) %>% count() %>% spread(key = test, value = n, fill = 0) %>%
  mutate(total = No + Yes, prop = No/total, fail = ifelse(prop >= .95, 1, 0)) %>%
  group_by(loc) %>% summarize(count = n(), n_failp95 = sum(fail))
loc_count
p2.1 <- tableGrob(loc_count, rows = NULL)

rem_ind <- geno %>%
  group_by(loc, id, test) %>% count() %>% spread(key = test, value = n, fill = 0) %>%
  mutate(total = No + Yes, prop = round(No/total, 2), fail = ifelse(prop >= .95, 1, 0)) %>%
  filter(fail == 1) %>% select(-fail)
rem_ind
p2.2 <- tableGrob(rem_ind, rows = NULL)

geno$dp <- as.numeric(geno$dp)
geno$allele_a <- as.numeric(geno$allele_a)
geno$allele_b <- as.numeric(geno$allele_b)
geno1 <- geno %>% filter(!is.na(gt)) %>% group_by(ID) %>% summarize(med_cov = median(dp, na.rm = TRUE))
med_med_coverage <- median(geno1$med_cov, na.rm=TRUE)
p3 <- ggplot(geno1, aes(x=med_cov)) + geom_histogram() +
  annotate(geom = "text", label = med_med_coverage, x = med_med_coverage, y = 50, color="white") +
  labs(x="Median coverage per locus", y="Number of loci") + theme_bw()
p3

geno2 <- geno %>% filter(!is.na(gt)) %>% select(ID, gt) %>% group_by(ID, gt) %>% count() %>%
  spread(key = gt, value = n, fill = 0) %>%
  mutate(total_gt = AA + AB + BB, Hobs = round(AB/ total_gt, 3), p = round( ((2*AA) + AB)/ (2*total_gt), 3), q = round(1-p,3)) %>%
  left_join(fixed_data %>% select(ID, snp_pos_in_read), by = "ID")

p3.1 <- ggplot(geno2, aes(x=snp_pos_in_read, y=q)) + geom_point(alpha=.5, size=3) + geom_smooth() + theme_bw() +
  labs(x="SNP Position Within RAD Tag", y= "Minor allele frequency")
p3.1
p3.2 <- ggplot(geno2, aes(x=snp_pos_in_read, y=Hobs)) + geom_point(alpha=.5, size=3) + geom_smooth() + theme_bw() +
  labs(x="SNP Position Within RAD Tag", y= "Observed heterozygosity")
p3.2

geno2$check <- "good"
geno2$check[geno2$AA != 0 & geno2$AB == 0 & geno2$BB ==0] <- "AA_only"
geno2$check[geno2$AA == 0 & geno2$AB != 0 & geno2$BB ==0] <- "AB_only"
geno2$check[geno2$AA == 0 & geno2$AB == 0 & geno2$BB !=0] <- "BB_only"
geno2$check[geno2$AA != 0 & geno2$AB == 0 & geno2$BB !=0] <- "No_hets"
g2_sum <- geno2 %>% group_by(check) %>% count()
g2_sum
p4 <- tableGrob(g2_sum, rows = NULL)

no_good <- geno2 %>% filter(check != "good") %>% arrange(check) %>% as.data.frame()
no_good

med_med_coverage <- round(median(geno2 %>% filter(check == "good") %>% pull(Hobs)), 2)
p5 <- ggplot(geno2 %>% filter(check == "good"), aes(x=Hobs)) + geom_histogram() +
  annotate(geom = "text", label = med_med_coverage, x = med_med_coverage, y = 50, color="white") +
  labs(x="Observed heterozygosity", y="Number of loci") + theme_bw()
p5

p6 <- geno %>% mutate(gt2 = ifelse(test = is.na(gt), 0, 1)) %>% group_by(loc, ID) %>%
  summarise(counts = n(), gt3 = sum(gt2), prop  = round(gt3/counts, 2)) %>%
  mutate(test  = ifelse(prop >= .5, 1, 0)) %>% select(ID, loc, test) %>% group_by(ID) %>%
  summarize(total = sum(test)) %>% ggplot(aes(x=total)) + geom_histogram(binwidth = 1) +
  theme_bw() + labs(x="Number of populations where\n>=50% of ind. are genotyped", y="Number of loci")
p6

p7 <- geno %>% mutate(gt2 = ifelse(test = is.na(gt), 0, 1)) %>% group_by(loc, ID) %>%
  summarise(counts = n(), gt3 = sum(gt2), prop  = round(gt3/counts, 2)) %>%
  mutate(test  = ifelse(prop >= .5, 1, 0)) %>% group_by(loc) %>% summarise(nloc = sum(test)) %>%
  ggplot(aes(x=loc, y=nloc)) + geom_bar(stat="identity", color="black", fill="grey") + theme_bw() +
  labs(x="Populations", y="Number of loci where\n>=50% of ind. are genotyped") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p7

geno3 <- geno %>% filter(gt == "AB") %>% select(ID, id:het, allele_a, allele_b) %>%
  filter(!is.na(allele_a), !is.na(allele_b)) %>% group_by(ID) %>%
  summarize(count = n(), a1 = sum(as.numeric(allele_a)), a2 = sum(as.numeric(allele_b)),
            tot_reads = sum(a1, a2), sig = sqrt(tot_reads*.5*.5), tr_2 = tot_reads/2, z = (tr_2 - a1)/sig)

g2 <- merge(x = geno2, y = geno3, by = "ID")
g2$singleton <- ifelse(test = (g2$Hobs <= 0.6 & abs(g2$z) <=  7), yes = "Yes", no = "No")
p8 <- ggplot(g2, aes(x=Hobs, y=z, color=singleton)) + geom_point(size=3) + theme_bw() +
  labs(x="Proportion of heterozygotes", y="Read-ratio deviation", color="Putative\nsingleton\nlocus")
p8

g2 %>% filter(singleton == "No")
if (!dir.exists("Output")) { dir.create("Output") }

pdf("Output/VCF_filtering_summary_plots.pdf", onefile = TRUE, width = 11, height = 8.5)
# Added the new NTC allele balance plot
grid.arrange(p_ntc_allele_balance, p_maf_comparison, p_snps_per_locus_before, p_snps_per_locus_after, nrow=2)
dev.off()

pdf("Output/VCF_analysis_plots_1.pdf", onefile = TRUE, width = 11, height = 8.5)
lay <- rbind(c(1,2), c(3,4))
grid.arrange(p1, p_maf, p3.1, p3.2, layout_matrix = lay)
dev.off()

pdf("Output/VCF_analysis_plots_2.pdf", onefile = TRUE, width = 11, height = 8.5)
grid.arrange(p3, p5, p6, p7, nrow=2)
dev.off()

pdf("Output/VCF_paralog_plot.pdf", onefile = TRUE, width = 8.5, height = 8.5)
grid.arrange(p8)
dev.off()

pdf("Output/VCF_genotype_graphic.pdf", onefile = TRUE, width=14, height=8.5)
grid.arrange(p2)
dev.off()

# fin!
