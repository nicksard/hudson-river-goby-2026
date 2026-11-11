# =============================================================================
# 08b_ldne_direct.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Table 3: linkage-disequilibrium Nb, excluding locus pairs on the same contig
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
# Updated:  2026-11-11 (sample-size correction of r^2; fish-group jackknife CIs)
#
# Description:
#   Estimates LD effective size for every group 07 wrote as a GENEPOP file
#   (all fish, each birth cohort, age-unambiguous cohorts, near-complete
#   subset). For a single cohort the estimate is Nb, the effective number of
#   breeders that produced it. The estimator lives in scripts/lib/
#   ldne_direct.R (method and references there). Two versions per group:
#     Ne_all_pairs      every locus pair; comparable to NeEstimator
#     Ne_diff_contig    pairs on different contigs only (primary): pairs on
#                       one contig are physically linked (r^2 0.061 vs 0.0044
#                       between contigs) and roughly halve the estimate.
#   The Waples, Larson & Waples (2016) chromosome-number correction (divisor
#   0.098 + 0.219 ln 23 = 0.7847) is reported for the all-pairs version only;
#   report it or the different-contig version, not both.
#
# Validation (see the supplement):
#   * Point estimates equal NeEstimator 2.1 on the same files: simulated
#     samples of a population with Ne = 500 (n = 45, 94, 166: 381.3, 602.2,
#     481.3 in both), the goby near-complete set (476.3 vs 476.9) and the 2023
#     cohort (503.6 vs 516.3).
#   * Fish-group jackknife CIs contained the true Ne in 15 of 16 simulated
#     samples. CIs are the primary columns (L95/U95); contig-jackknife CIs
#     (*_contigjk) are kept for comparison and are too narrow.
#   * NeEstimator's own large runs on the Windows PC (Nov 10) returned
#     corrupted output; the same exe under Wine on Linux did not.
#
# Order of operations:
#   1. Read each results/07_ldne/genepop/<group>.gen.
#   2. ld_direct_fish_ci(): estimate on all fish, then with each of 10 fish
#      groups left out (about 11 times one estimate; ~1 h for all groups).
#   3. Write one row per group.
#
# Inputs:   results/07_ldne/genepop/<group>.gen (from 07_ldne_inputs.R)
#           scripts/lib/ldne_direct.R
# Outputs:  results/07_ldne/direct/ld_ne_direct.csv
#
# References:
#   Waples RS, Do C (2008) LDNE. Molecular Ecology Resources 8:753-756.
#   Waples RK, Larson WA, Waples RS (2016) Estimating contemporary effective
#     population size in non-model species using linkage disequilibrium
#     across thousands of loci. Heredity 117:233-240.
#   Do C, Waples RS, Peel D, Macbeth GM, Tillett BJ, Ovenden JR (2014)
#     NeEstimator v2. Molecular Ecology Resources 14:209-214.
# =============================================================================

source("scripts/00_config.R")
suppressPackageStartupMessages({ library(readr); library(dplyr) })

genepop_dir <- project_path("results/07_ldne/genepop")
out_dir     <- project_path("results/07_ldne/direct")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pcrit       <- ldne_settings$pcrit
block_size  <- 1500
waples_divisor <- 0.098 + 0.219 * log(ldne_settings$n_chromosomes)

groups <- sub("\\.gen$", "", list.files(genepop_dir, pattern = "\\.gen$"))
if (exists("ld_direct_groups")) groups <- intersect(groups, ld_direct_groups)
msg("08b: direct LD Ne for ", length(groups), " group(s): ", paste(groups, collapse = ", "))

source("scripts/lib/ldne_direct.R", local = TRUE)   # read_genepop(), ld_direct(), ld_direct_fish_ci()

results_by_group <- list()
for (group_name in groups) {
  started <- Sys.time()
  genepop <- read_genepop(file.path(genepop_dir, paste0(group_name, ".gen")))
  result  <- ld_direct_fish_ci(genepop$dosage, genepop$loci, n_groups = 10, seed = random_seed)
  all_pairs <- result$all; between_contigs <- result$diff
  results_by_group[[group_name]] <- data.frame(
    group = group_name, n_fish = result$n_fish, loci_used = result$loci_used,
    contigs = result$contigs, formula = if (result$small_sample) "S<30" else "S>=30",
    pairs_all = all_pairs[["pairs"]], r2_all = all_pairs[["r2"]], E_r2 = all_pairs[["E"]],
    Ne_all_pairs = all_pairs[["Ne"]],
    L95_all = all_pairs[["Ne_L95_fish"]], U95_all = all_pairs[["Ne_U95_fish"]],
    L95_all_contigjk = all_pairs[["Ne_L95"]], U95_all_contigjk = all_pairs[["Ne_U95"]],
    Ne_all_waples2016 = all_pairs[["Ne"]] / waples_divisor,
    L95_all_waples2016 = all_pairs[["Ne_L95_fish"]] / waples_divisor,
    U95_all_waples2016 = all_pairs[["Ne_U95_fish"]] / waples_divisor,
    pairs_diff = between_contigs[["pairs"]], r2_diff = between_contigs[["r2"]],
    Ne_diff_contig = between_contigs[["Ne"]],
    L95_diff = between_contigs[["Ne_L95_fish"]], U95_diff = between_contigs[["Ne_U95_fish"]],
    L95_diff_contigjk = between_contigs[["Ne_L95"]], U95_diff_contigjk = between_contigs[["Ne_U95"]])
  show <- function(x) ifelse(is.finite(x), sprintf("%.0f", x), "inf")
  msg(sprintf("  %-20s n=%3d loci=%5d  all pairs: %s (%s-%s)  different contigs: %s (%s-%s)  [%.1f min]",
              group_name, result$n_fish, result$loci_used,
              show(all_pairs[["Ne"]]), show(all_pairs[["Ne_L95_fish"]]), show(all_pairs[["Ne_U95_fish"]]),
              show(between_contigs[["Ne"]]), show(between_contigs[["Ne_L95_fish"]]),
              show(between_contigs[["Ne_U95_fish"]]),
              as.numeric(difftime(Sys.time(), started, units = "mins"))))
}
out <- bind_rows(results_by_group)
write_csv(out, file.path(out_dir, "ld_ne_direct.csv"))
msg_table(out %>% select(group, n_fish, loci_used, Ne_all_pairs, L95_all, U95_all,
                         Ne_diff_contig, L95_diff, U95_diff))
msg("  Methods decision: report Ne_diff_contig (linked pairs removed) or",
    " Ne_all_waples2016 (all pairs, Waples et al. 2016 correction), not both corrections.")
