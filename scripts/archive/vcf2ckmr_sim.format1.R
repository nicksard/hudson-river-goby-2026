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
working_dir <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis"
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

set.seed(42) 
subset_loci <- unique(paste(tidy_gt$ChromKey, tidy_gt$POS, sep = "_"))
if(length(subset_loci) > 500) {
  message("--- Subsampling to 500 random SNPs ---")
  selected_loci <- sample(subset_loci, 500)
  tidy_gt <- tidy_gt %>% 
    filter(paste(ChromKey, POS, sep = "_") %in% selected_loci)
}

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

# Simulate for all potential comparisons
sim_kin <- simulate_Qij(ck_obj, reps = 10000,
                        sim_relats = c("FS", "HS", "FC", "U"),
                        calc_relats = c("FS", "HS", "FC", "U"))

# Define the ratios we care about for classification
sim_logls <- bind_rows(
  extract_logls(sim_kin, numer = c(FS = 1), denom = c(U = 1)) %>% mutate(Ratio = "FS/U"),
  extract_logls(sim_kin, numer = c(HS = 1), denom = c(U = 1)) %>% mutate(Ratio = "HS/U"),
  extract_logls(sim_kin, numer = c(FC = 1), denom = c(U = 1)) %>% mutate(Ratio = "FC/U"),
  extract_logls(sim_kin, numer = c(FS = 1), denom = c(HS = 1)) %>% mutate(Ratio = "FS/HS"),
  extract_logls(sim_kin, numer = c(HS = 1), denom = c(FC = 1)) %>% mutate(Ratio = "HS/FC")
)

# --- 5. EMPIRICAL PAIRWISE ANALYSIS ---
message("--- Running Multi-Ratio Pairwise Likelihood Calculations ---")

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

# We calculate all requested ratios to provide a full suite of support metrics
obs_ratios <- list(
  FS_U  = pairwise_kin_logl_ratios(clean_genos, clean_genos, ck_obj, "FS", "U",  num_cores = n_cores),
  HS_U  = pairwise_kin_logl_ratios(clean_genos, clean_genos, ck_obj, "HS", "U",  num_cores = n_cores),
  FC_U  = pairwise_kin_logl_ratios(clean_genos, clean_genos, ck_obj, "FC", "U",  num_cores = n_cores),
  FC_HS = pairwise_kin_logl_ratios(clean_genos, clean_genos, ck_obj, "FC", "HS", num_cores = n_cores),
  FC_FS = pairwise_kin_logl_ratios(clean_genos, clean_genos, ck_obj, "FC", "FS", num_cores = n_cores),
  HS_FS = pairwise_kin_logl_ratios(clean_genos, clean_genos, ck_obj, "HS", "FS", num_cores = n_cores),
  FS_HS = pairwise_kin_logl_ratios(clean_genos, clean_genos, ck_obj, "FS", "HS", num_cores = n_cores)
)
str(obs_)
# Merge empirical results into a single large table
obs_combined <- obs_ratios$FS_U %>%
  rename(FS_U = logl_ratio) %>%
  inner_join(obs_ratios$HS_U %>% select(D1_indiv, D2_indiv, HS_U = logl_ratio), by = c("D1_indiv", "D2_indiv")) %>%
  inner_join(obs_ratios$FC_U %>% select(D1_indiv, D2_indiv, FC_U = logl_ratio), by = c("D1_indiv", "D2_indiv")) %>%
  inner_join(obs_ratios$FC_HS %>% select(D1_indiv, D2_indiv, FC_HS = logl_ratio), by = c("D1_indiv", "D2_indiv")) %>%
  inner_join(obs_ratios$FC_FS %>% select(D1_indiv, D2_indiv, FC_FS = logl_ratio), by = c("D1_indiv", "D2_indiv")) %>%
  inner_join(obs_ratios$HS_FS %>% select(D1_indiv, D2_indiv, HS_FS = logl_ratio), by = c("D1_indiv", "D2_indiv")) %>%
  inner_join(obs_ratios$FS_HS %>% select(D1_indiv, D2_indiv, FS_HS = logl_ratio), by = c("D1_indiv", "D2_indiv")) %>%
  filter(D1_indiv < D2_indiv)

# --- 6. CLASSIFICATION & ODDS ---
message("--- Performing Multi-Hypothesis Classification ---")

# Centering on Unrelated for the FS/U ratio (standard baseline)
sim_u_median <- median(sim_logls$logl_ratio[sim_logls$true_relat == "U" & sim_logls$Ratio == "FS/U"])
obs_median <- median(obs_combined$FS_U)

# Adjust for ancestry
obs_centered <- obs_combined %>%
  mutate(FS_U_adj = FS_U - obs_median)

# Refined Posterior Probability Function
# This takes the FS/U score and determines which distribution it fits best
calc_multi_posteriors <- function(obs_val, sim_data) {
  relats <- c("FS", "HS", "FC", "U")
  probs <- sapply(relats, function(r) {
    d_vals <- sim_data$logl_ratio[sim_data$true_relat == r & sim_data$Ratio == "FS/U"] - sim_u_median
    dens <- density(d_vals)
    approx(dens$x, dens$y, xout = obs_val, rule = 2)$y
  })
  posteriors <- probs / sum(probs)
  return(posteriors)
}

# Run classification on candidates
# Only evaluating pairs that show some support for being related
candidates <- obs_centered %>% filter(FS_U_adj > 0)
candidates
if(nrow(candidates) > 0) {
  classification_results <- candidates %>%
    mutate(posteriors = map(FS_U_adj, ~calc_multi_posteriors(.x, sim_logls))) %>%
    unnest_wider(posteriors) %>%
    rowwise() %>%
    mutate(
      Best_Relationship = c("FS", "HS", "FC", "U")[which.max(c(FS, HS, FC, U))],
      Confidence = max(c(FS, HS, FC, U))
    ) %>%
    ungroup()
  
  write_csv(classification_results, "Output/Goby_Kinship_Classification_Summary.csv")
}
classification_results
# --- 7. VISUALIZATION ---
message("--- Generating Final Reports ---")

p_multi <- ggplot(sim_logls, aes(x = logl_ratio, fill = true_relat)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Ratio, scales = "free") +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Theoretical Power Across Multi-Hypothesis Ratios",
       x = "Log-Likelihood Ratio", y = "Density")
p_multi
ggsave(filename = paste0(plot_dir, "Goby_Multi_Ratio_Power.png"), plot = p_multi, width = 12, height = 8)

message("Analysis Complete.")
beep(sound = 8)