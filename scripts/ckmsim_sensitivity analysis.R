# ==============================================================================
# Script: goby_family_chao.R
# Purpose: Threshold Refinement, Family Clustering, and Chao1 Parentage Estimation
# Directory: G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis
# ==============================================================================

library(tidyverse)
library(igraph) # For pedigree/family clustering
library(SpadeR) # For advanced richness and Chao1 estimation

# --- 1. SETTINGS & DATA LOAD ---
working_dir <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis"
setwd(working_dir)

# Loading the assignments and thresholds from the previous CKMR run
results <- read_csv("Output/Goby_Kinship_Final_Assignments.csv")
thresholds <- read_csv("Output/Goby_Kinship_Final_Thresholds.csv")

# --- 2. REFINED ANCESTRY CORRECTION ---
# Using the specific delta-shift adjustment to align empirical unrelated to simulated baseline
# Results in empirical FS/U ~3000 for top pairs and ~1100-1300 for simulations
sim_u_median <- thresholds$q50[thresholds$true_relat == "U"]
obs_median <- median(results$logl_ratio[results$D1_indiv != results$D2_indiv])

# Step 1: Align empirical data to simulated baseline
# This preserves the raw scale of the simulated distributions
results$logl_adj <- results$logl_ratio - (obs_median - sim_u_median)

# --- 3. DYNAMIC THRESHOLD DEFINITION ---
message("--- Applying Multi-Tiered Thresholds (FPR Focus) ---")

# Pulling specific values based on your query logic
dup_thresh <- thresholds$q999[thresholds$true_relat == "FS"] # ~1310
fs_thresh  <- thresholds$q99[thresholds$true_relat == "HS"]  # ~59.0
hs_thresh  <- thresholds$q01[thresholds$true_relat == "HS"]  # ~ -167
fc_thresh  <- thresholds$q999[thresholds$true_relat == "U"]  # ~ -989

# --- 4. FINAL CLASSIFICATION LOGIC ---
# Categorizing based on your specific threshold adjustments
results_final <- results %>%
  mutate(
    Final_Category = case_when(
      logl_adj > dup_thresh ~ "Duplicate/MZ",
      logl_adj > fs_thresh  ~ "Full-Sib",
      logl_adj > hs_thresh  ~ "Half-Sib",
      logl_adj > fc_thresh  ~ "Cousin/Distantly Related",
      TRUE ~ "Unrelated"
    )
  )

# --- 5. PARENTAGE SCENARIO MODELING ---
message("--- Reconstructing Families under Two Scenarios ---")

# Scenario A: Conservative (Only Full-Sibs share parents)
# Scenario B: Inclusive (Full-Sibs and Half-Sibs share parents)
# Given your results (only 2 pairs likely siblings), we build success vectors for both.

all_indivs <- unique(c(results$D1_indiv, results$D2_indiv))

get_np_report <- function(category_filter, label) {
  # Step 1: Cluster pairs that pass the filter
  kin_pairs <- results_final %>% filter(Final_Category %in% category_filter)
  
  if(nrow(kin_pairs) > 0) {
    g <- graph_from_data_frame(kin_pairs[, c("D1_indiv", "D2_indiv")], directed = FALSE)
    cl <- components(g)
    family_membership <- data.frame(
      Indiv = names(membership(cl)),
      Family_ID = as.vector(membership(cl))
    )
    
    # Each family of size N implies parents (simplified to abundance for richness)
    # Full-Sibs count twice (mother/father), Half-Sibs are more complex but for 
    # richness estimation, we look at the frequency of parental occurrences.
    success_vector <- family_membership %>%
      group_by(Family_ID) %>%
      summarize(n_offspring = n()) %>%
      pull(n_offspring)
    
    # Add singletons
    singletons <- length(all_indivs) - length(family_membership$Indiv)
    success_vector <- c(success_vector, rep(1, singletons))
  } else {
    success_vector <- rep(1, length(all_indivs))
  }
  
  # Step 2: SpadeR Richness
  m1 <- ChaoSpecies(data = success_vector, datatype = "abundance", conf = 0.95)
  res <- as.data.frame(m1$Species_table)
  chao_row <- res[grepl("iChao1", rownames(res)), ]
  
  data.frame(
    Scenario = label,
    n_offspring = length(all_indivs),
    n_parents_obs = length(success_vector),
    Np_Estimate = round(chao_row$Estimate),
    Np_LCL = round(chao_row$`95%Lower`),
    Np_UCL = round(chao_row$`95%Upper`),
    f1 = sum(success_vector == 1),
    f2 = sum(success_vector == 2)
  )
}

report_fs_only <- get_np_report(c("Full-Sib"), "FS_Only")
report_fs_hs   <- get_np_report(c("Full-Sib", "Half-Sib"), "FS_and_HS")

final_report <- rbind(report_fs_only, report_fs_hs)
print(final_report)

# --- 6. EXPORT & PLOT ---
write_csv(results_final, "Output/Goby_Final_Pedigree_Assignments.csv")
write_csv(final_report, "Output/Goby_Parentage_Scenario_Report.csv")

p_final <- ggplot(results_final, aes(x = logl_adj, fill = Final_Category)) +
  geom_density(alpha = 0.7) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Goby Kinship Classification: Refined Correction",
       subtitle = "Empirical data adjusted by delta (Obs_med - Sim_med).",
       x = "Centered LogL (FS/U)", y = "Density")
p_final
ggsave("Plots/Goby_Final_Family_Split.png", p_final, width = 11, height = 7)

message("Comparison report generated for FS-only and FS+HS scenarios.")
beep(sound = 8)
