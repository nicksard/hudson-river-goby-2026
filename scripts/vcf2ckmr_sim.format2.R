# ==============================================================================
# Script: goby_ckmr_local.R
# Purpose: VCF-to-CKMR conversion, Simulation, and Multi-Hypothesis Kinship Inference
# Directory: G:/Shared drives/Sard Lab/Round Goby/RADseq/R/post_RAD analysis
# ==============================================================================

library(vcfR)
library(tidyverse)
library(CKMRsim)
library(beepr)

# --- 1. SETTINGS & PATHS ---
working_dir <- "G:/Shared drives/Sard Lab/Round Goby/RADseq/R/post_RAD analysis"
setwd(working_dir)

vcf_file <- "Output/populations.filtered.final.vcf.gz"
out_csv <- "Output/CKMRsim_Genotypes_Ready.csv"
plot_dir <- "Plots/"

if (!dir.exists("Output")) dir.create("Output")
if (!dir.exists(plot_dir)) dir.create(plot_dir)

# --- 2. VCF CONVERSION ---
message("--- Loading and Converting VCF ---")
vcf <- read.vcfR(vcf_file)
tidy_gt <- vcfR2tidy(vcf, format_fields = c("GT"))$gt

#set.seed(42) 
#subset_loci <- unique(paste(tidy_gt$ChromKey, tidy_gt$POS, sep = "_"))
#if(length(subset_loci) > 500) {
#  message("--- Subsampling to 500 random SNPs ---")
#  selected_loci <- sample(subset_loci, 500)
#  tidy_gt <- tidy_gt %>% 
#    filter(paste(ChromKey, POS, sep = "_") %in% selected_loci)
#}

genos_long <- tidy_gt %>%
  filter(!is.na(gt_GT) & gt_GT != "./.") %>%
  separate(gt_GT, into = c("A1", "A2"), sep = "[|/]") %>%
  pivot_longer(cols = c(A1, A2), names_to = "GeneCopy", values_to = "Allele") %>%
  mutate(
    GeneCopy = ifelse(GeneCopy == "A1", 1, 2),
    Allele = as.character(Allele),
    Locus = paste(ChromKey, POS, sep = "_")
  ) %>%
  select(Locus, Indiv, GeneCopy, Allele)

#write_csv(genos_long, out_csv)

# --- 3. ALLELE FREQUENCIES & RE-INDEXING ---
message("--- Estimating Allele Frequencies ---")
alfs <- genos_long %>%
  group_by(Locus, Allele) %>%
  summarize(count = n(), .groups = "drop_last") %>%
  mutate(Freq = count / sum(count)) %>%
  ungroup() %>%
  mutate(
    Chrom = "Unk", 
    Pos = as.numeric(as.factor(Locus)),
    AlleIdx = NA,
    LocIdx = NA
  ) %>%
  select(Chrom, Locus, Pos, Allele, Freq, AlleIdx, LocIdx)

re_indexed <- reindex_markers(alfs)

# --- 4. CKMR OBJECT & SIMULATIONS ---
message("--- Building CKMR Model and Simulating ---")
epsilon_val <- 0.005 

kappas_custom <- matrix(c(
  0.25, 0.5, 0.25,  # Full Sib (FS)
  0.5, 0.5, 0,      # Half Sib (HS)
  0.75, 0.25, 0,    # First Cousin (FC)
  1, 0, 0           # Unrelated (U)
), nrow = 4, byrow = TRUE, dimnames = list(c("FS", "HS", "FC", "U"), NULL))

ck_obj <- create_ckmr(
  re_indexed, 
  kappa_matrix = kappas_custom,
  ge_mod_assumed = ge_model_TGIE,
  ge_mod_true = ge_model_TGIE,
  ge_mod_assumed_pars_list = list(epsilon = epsilon_val),
  ge_mod_true_pars_list = list(epsilon = epsilon_val)
)

sim_kin <- simulate_Qij(ck_obj, reps = 10000,
                        sim_relats = c("FS", "HS", "FC", "U"),
                        calc_relats = c("FS", "HS", "FC", "U"))

sim_logls <- extract_logls(sim_kin, numer = c(FS = 1), denom = c(U = 1))

# --- 5. EMPIRICAL PAIRWISE ANALYSIS ---
message("--- Running Pairwise Likelihood Calculations ---")

clean_genos <- genos_long %>%
  filter(Locus %in% re_indexed$Locus) %>%
  select(Indiv, Locus, gene_copy = GeneCopy, Allele) %>%
  mutate(
    Allele = as.character(Allele), 
    gene_copy = as.integer(gene_copy)
  ) %>%
  distinct(Indiv, Locus, gene_copy, .keep_all = TRUE)

n_cores <- parallel::detectCores() - 1
if(.Platform$OS.type == "windows") n_cores <- 1

obs_ratios <- pairwise_kin_logl_ratios(
  D1 = clean_genos, 
  D2 = clean_genos,
  CK = ck_obj, 
  numer = "FS", 
  denom = "U",
  num_cores = n_cores
)
obs_ratios

# --- 6. ANCESTRY CORRECTION & THRESHOLDING (FPR/FNR FOCUS) ---
message("--- Applying Percentile-Based Thresholds ---")

# Step 1: Center everything on the Unrelated Median
sim_u_median <- median(sim_logls$logl_ratio[sim_logls$true_relat == "U"])
obs_median <- median(obs_ratios$logl_ratio[obs_ratios$D1_indiv != obs_ratios$D2_indiv])

sim_logls$logl_adj <- sim_logls$logl_ratio
obs_ratios$logl_adj <- obs_ratios$logl_ratio - (obs_median-sim_u_median)

# Step 2: Calculate Percentiles for Kinship Boundaries
# We define thresholds based on where 99% of a relationship type falls
thresholds <- sim_logls %>%
  group_by(true_relat) %>%
  summarize(
    q01 = quantile(logl_adj, 0.01),
    q05 = quantile(logl_adj, 0.05),
    q50 = quantile(logl_adj, 0.50),
    q95 = quantile(logl_adj, 0.95),
    q99 = quantile(logl_adj, 0.99),
    q999 = quantile(logl_adj, 0.999) # For duplicate detection
  )
thresholds

# Step 3: Classification based on thresholds (Upper tail control for FPR)
# Using the user's logic: thresholds defined by the upper tails of overlapping distributions
dup_thresh <- thresholds$q999[thresholds$true_relat == "FS"]
fs_thresh <- thresholds$q99[thresholds$true_relat == "HS"]
hs_thresh <- thresholds$q01[thresholds$true_relat == "HS"]
fc_thresh <- thresholds$q999[thresholds$true_relat == "U"]

# Output thresholds to console for verification
message(paste("Duplicate Threshold:", dup_thresh))
message(paste("Full-Sib Threshold:", fs_thresh))
message(paste("Half-Sib Threshold:", hs_thresh))
message(paste("Cousin Threshold:", fc_thresh))

results <- obs_ratios %>%
  filter(D1_indiv != D2_indiv) %>%
  mutate(
    Category = case_when(
      logl_adj > dup_thresh ~ "Duplicated",
      logl_adj > fs_thresh  ~ "Full-Sib",
      logl_adj > hs_thresh  ~ "Half-Sib",
      logl_adj > fc_thresh  ~ " Cousin",
      TRUE ~ "Unrelated"
    )
  ) %>%
  arrange(desc(logl_adj))
head(results)
head(rev(sort(results$logl_adj)))

# Assign Category to true_relat for density plotting and table generation
results$true_relat <- factor(results$Category,
                             levels = c("Duplicated",
                                        "Full-Sib",
                                        "Half-Sib",
                                        "Cousin",
                                        "Unrelated"))

print(table(results$true_relat))

write_csv(results, "Output/Goby_Kinship_Final_Assignments.csv")
write_csv(thresholds, "Output/Goby_Kinship_Final_Thresholds.csv")

# --- 7. VISUALIZATION ---
# Page 1: Theoretical vs Empirical Overlay with Threshold Lines
p1 <- ggplot(sim_logls, aes(x = logl_adj, fill = true_relat)) +
  geom_density(alpha = 0.3) +
  geom_density(data = results, aes(x = logl_adj, fill = true_relat), alpha = 0.6) +
  geom_vline(xintercept = dup_thresh, linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = fs_thresh, linetype = "dashed", color = "green", size = 1) +
  geom_vline(xintercept = hs_thresh, linetype = "dashed", color = "purple", size = 1) +
  geom_vline(xintercept = fc_thresh, linetype = "dashed", color = "blue", size = 1) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Hypothesis Testing: Percentile Thresholds",
       subtitle = "Colored dashed lines indicate kinship thresholds (FPR focused).",
       x = "Centered Log-Likelihood Ratio (FS/U)", y = "Density")

print(p1)

ggsave(filename = paste0(plot_dir, "Goby_Threshold_Density.png"), plot = p1, width = 10, height = 7)

message("Analysis Complete.")
beep(sound = 8)
