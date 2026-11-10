# ==============================================================================
# Script: goby_family_chao.R
# Purpose: Direct Scenario Testing for Np (FS vs HS and Adjusted vs Unadjusted)
# Directory: G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis
# ==============================================================================

library(tidyverse)
library(igraph) 
library(SpadeR) 

# --- 1. SETTINGS & DATA LOAD ---
working_dir <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis"
setwd(working_dir)

results <- read_csv("Output/Goby_Kinship_Final_Assignments.csv")
thresholds <- read_csv("Output/Goby_Kinship_Final_Thresholds.csv")

# Constants for Thresholds
sim_u_median <- thresholds$q50[thresholds$true_relat == "U"]
obs_median   <- median(results$logl_ratio[results$D1_indiv != results$D2_indiv])
delta_shift  <- obs_median - sim_u_median

# Threshold Values (FPR control)
fs_thresh_adj <- thresholds$q99[thresholds$true_relat == "HS"]
hs_thresh_adj <- thresholds$q01[thresholds$true_relat == "HS"]

# --- 2. DATA PREPARATION ---
# Create Adjusted and Unadjusted columns
results <- results %>%
  mutate(
    logl_unadj = logl_ratio,
    logl_adj   = logl_ratio - delta_shift
  )

all_indivs <- unique(c(results$D1_indiv, results$D2_indiv))
n_offspring <- length(all_indivs)

# --- 3. SCENARIO 1: FULL-SIB HYPOTHESIS (ADJUSTED) ---
# Treat top 2 pairs as Full-Sibs, assume they share 2 parents.
message("--- Running Scenario 1: Full-Sib (Adjusted) ---")
fs_pairs_1 <- results %>% filter(logl_adj > fs_thresh_adj)
fs_pairs_1
g1 <- graph_from_data_frame(fs_pairs_1[, c("D1_indiv", "D2_indiv")], directed = FALSE)
cl1 <- components(g1)
m1_indivs <- names(membership(cl1))
f1_sizes <- as.vector(table(membership(cl1)))
# FS Logic: 2 parents per family
success_v1 <- c(rep(f1_sizes, each = 2), rep(1, n_offspring - length(m1_indivs)))

chao_1 <- ChaoSpecies(data = success_v1, datatype = "abundance", conf = 0.95)
rep_1 <- as.data.frame(chao_1$Species_table)[4, ] # iChao1

# --- 4. SCENARIO 2: HALF-SIB HYPOTHESIS (ADJUSTED) ---
# Treat top 2 pairs as Half-Sibs, assume they share 1 parent.
message("--- Running Scenario 2: Half-Sib (Adjusted) ---")
hs_pairs_2 <- results %>% filter(logl_adj > hs_thresh_adj)
hs_pairs_2
g2 <- graph_from_data_frame(hs_pairs_2[, c("D1_indiv", "D2_indiv")], directed = FALSE)
cl2 <- components(g2)
m2_indivs <- names(membership(cl2))
f2_sizes <- as.vector(table(membership(cl2)))
# HS Logic: 1 parent shared per family
success_v2 <- c(f2_sizes, rep(1, n_offspring - length(m2_indivs)))

chao_2 <- ChaoSpecies(data = success_v2, datatype = "abundance", conf = 0.95)
rep_2 <- as.data.frame(chao_2$Species_table)[4, ]

# --- 5. SCENARIO 3: UNADJUSTED BIAS TEST (FULL-SIB) ---
# Use raw Log-Likelihood (Biased High). More pairs will likely cross the threshold.
message("--- Running Scenario 3: Full-Sib (Unadjusted) ---")
fs_pairs_3 <- results %>% filter(logl_unadj > fs_thresh_adj) # Threshold remains fixed
fs_pairs_3
g3 <- graph_from_data_frame(fs_pairs_3[, c("D1_indiv", "D2_indiv")], directed = FALSE)
cl3 <- components(g3)
m3_indivs <- names(membership(cl3))
f3_sizes <- as.vector(table(membership(cl3)))
success_v3 <- c(rep(f3_sizes, each = 2), rep(1, n_offspring - length(m3_indivs)))

chao_3 <- ChaoSpecies(data = success_v3, datatype = "abundance", conf = 0.95)
rep_3 <- as.data.frame(chao_3$Species_table)[4, ]

# --- 6. SCENARIO 4: UNADJUSTED BIAS TEST (HALF-SIB) ---
message("--- Running Scenario 4: Half-Sib (Unadjusted) ---")
hs_pairs_4 <- results %>% filter(logl_unadj > hs_thresh_adj)
hs_pairs_4
g4 <- graph_from_data_frame(hs_pairs_4[, c("D1_indiv", "D2_indiv")], directed = FALSE)
cl4 <- components(g4)
m4_indivs <- names(membership(cl4))
f4_sizes <- as.vector(table(membership(cl4)))
success_v4 <- c(f4_sizes, rep(1, n_offspring - length(m4_indivs)))

chao_4 <- ChaoSpecies(data = success_v4, datatype = "abundance", conf = 0.95)
rep_4 <- as.data.frame(chao_4$Species_table)[4, ]

# --- 7. FINAL COMPARISON REPORT ---
final_report <- data.frame(
  Scenario = c("FS_Adjusted", "HS_Adjusted", "FS_Unadjusted", "HS_Unadjusted"),
  Assumed_Kinship = c("Full-Sib", "Half-Sib", "Full-Sib", "Half-Sib"),
  Ancestry_Correction = c("Yes", "Yes", "No", "No"),
  n_offspring = n_offspring,
  n_parents_obs = c(length(success_v1), length(success_v2), length(success_v3), length(success_v4)),
  Np_Estimate = c(round(rep_1$Estimate), round(rep_2$Estimate), round(rep_3$Estimate), round(rep_4$Estimate)),
  Np_LCL = c(round(rep_1$`95%Lower`), round(rep_2$`95%Lower`), round(rep_3$`95%Lower`), round(rep_4$`95%Lower`)),
  Np_UCL = c(round(rep_1$`95%Upper`), round(rep_2$`95%Upper`), round(rep_3$`95%Upper`), round(rep_4$`95%Upper`)),
  f1_singletons = c(sum(success_v1 == 1), sum(success_v2 == 1), sum(success_v3 == 1), sum(success_v4 == 1)),
  f2_doubletons = c(sum(success_v1 == 2), sum(success_v2 == 2), sum(success_v3 == 2), sum(success_v4 == 2))
)

print(final_report)
write_csv(final_report, "Output/Goby_Np_Sensitivity_Analysis.csv")

message("Sensitivity analysis complete. See Goby_Np_Sensitivity_Analysis.csv for comparisons.")
beep(sound = 8)