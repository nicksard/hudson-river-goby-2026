# =============================================================================
# 08b_ldne_direct.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# LD Ne computed directly, excluding locus pairs on the same contig
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
#
# Why this exists (findings of 2026-11-10):
#   * NeEstimator 2.1 itself is sound. On simulated data (true Ne = 100,
#     1,200-10,000 loci, random or clustered missingness) it returns 94-100
#     with CIs that contain the point. On the goby near-complete file the same
#     exe returns 476.9 (461.2-493.7) when run under Wine on Linux. The Windows
#     runs of the morning of Nov 10 (306.1 with CI ~111, etc.) were corrupted
#     in the run itself, not in the data.
#   * Locus pairs on the same contig are physically linked: r^2 = 0.061 vs
#     0.0044 for pairs on different contigs (near-complete set). Including
#     them roughly halves Ne (503 vs 974). NeEstimator has no option to drop
#     them, so this script does the estimate itself.
#
# Method:
#   r^2 per locus pair = squared Pearson correlation of allele dosages (the
#   composite r^2 for biallelic loci), loci with MAF >= pcrit in the group,
#   pairs with S >= 30 (S < 30 for groups with fewer fish: see below).
#   E(r^2) per pair from Waples & Do (2008): S >= 30: 1/S + 3.19/S^2;
#   S < 30: 0.0018 + 0.907/S + 4.44/S^2. r^2' = S^2-weighted mean of
#   (r^2 - E), then Ne from the matching Waples & Do (2008) formula.
#   CI: delete-one-contig jackknife on r^2' (pairs are not independent, and
#   contigs are the natural blocks), converted to Ne; a lower r^2' limit <= 0
#   gives an infinite upper Ne limit.
#   Computed in blocks of loci, accumulating contig x contig sums, so memory
#   stays bounded for the full 10,943-locus panel.
#
#   Two estimates per group:
#     Ne_all_pairs       comparable to NeEstimator (all pairs)
#     Ne_diff_contig     primary: pairs on different contigs only
#   The Waples et al. (2016) chromosome correction (divisor 0.7847) is reported
#   for the all-pairs estimate only. It corrects the bias from linked pairs on
#   the same chromosome, which the different-contig estimate already removes
#   in part (contigs on one chromosome are still loosely linked). Which to
#   report is a Methods decision: see the log.
#
# Inputs:  results/07_ldne/genepop/<group>.gen (written by 07)
# Outputs: results/07_ldne/direct/ld_ne_direct.csv
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

source("scripts/lib/ldne_direct.R", local = TRUE)      # read_genepop(), ne_from_r2p(), ld_direct()

rows <- list()
for (g in groups) {
  t0 <- Sys.time()
  dat <- read_genepop(file.path(genepop_dir, paste0(g, ".gen")))
  res <- ld_direct(dat$dosage, dat$loci)
  rows[[g]] <- data.frame(
    group = g, n_fish = res$n_fish, loci_used = res$loci_used, contigs = res$contigs,
    formula = if (res$small_sample) "S<30" else "S>=30",
    pairs_all = res$all[["pairs"]], r2_all = res$all[["r2"]], E_r2 = res$all[["E"]],
    Ne_all_pairs = res$all[["Ne"]], L95_all = res$all[["Ne_L95"]], U95_all = res$all[["Ne_U95"]],
    Ne_all_waples2016 = res$all[["Ne"]] / waples_divisor,
    L95_all_waples2016 = res$all[["Ne_L95"]] / waples_divisor,
    U95_all_waples2016 = res$all[["Ne_U95"]] / waples_divisor,
    pairs_diff = res$diff[["pairs"]], r2_diff = res$diff[["r2"]],
    Ne_diff_contig = res$diff[["Ne"]], L95_diff = res$diff[["Ne_L95"]], U95_diff = res$diff[["Ne_U95"]])
  msg(sprintf("  %-22s n=%3d loci=%5d  all pairs: %7.1f (%s-%s)  different contigs: %7.1f (%s-%s)  [%.1f min]",
              g, res$n_fish, res$loci_used,
              res$all[["Ne"]], format(round(res$all[["Ne_L95"]], 1)), format(round(res$all[["Ne_U95"]], 1)),
              res$diff[["Ne"]], format(round(res$diff[["Ne_L95"]], 1)), format(round(res$diff[["Ne_U95"]], 1)),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
out <- bind_rows(rows)
write_csv(out, file.path(out_dir, "ld_ne_direct.csv"))
msg_table(out %>% select(group, n_fish, loci_used, Ne_all_pairs, L95_all, U95_all,
                         Ne_diff_contig, L95_diff, U95_diff))
msg("  Methods decision: report Ne_diff_contig (linked pairs removed) or",
    " Ne_all_waples2016 (all pairs, Waples et al. 2016 correction), not both corrections.")
