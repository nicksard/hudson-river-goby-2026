# =============================================================================
# 08c_ldne_stratag.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Cross-check of 08b against an independent R implementation: strataG::ldNe
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-11
# Updated:  2026-11-11 (08b now applies the (S/(S-1))^2 correction; reframed)
#
# Description:
#   The Table 3 LD estimates come from our own code (08b, scripts/lib/
#   ldne_direct.R). NeEstimator 2.1 is the usual tool, but it is a
#   closed-source Windows/Java program that is hard to obtain and whose large
#   runs on our Windows PC produced corrupted output. This script records how
#   08b compares with strataG::ldNe (Archer et al. 2017), the only published R
#   implementation of the Waples & Do (2008) estimator.
#
#   What the comparison showed (logged 2026-11-11):
#     * Locus-pair r^2: 08b's uncorrected Pearson r^2 matches strataG's mean r^2
#       closely (same pairs, same fish). The r^2 calculation is not in question.
#     * strataG multiplies r^2 by (S/(S-1))^2 before comparing it with E(r^2),
#       following the Waples & Larson code it was adapted from. 08b now does
#       the same (lib/ldne_direct.R, sample_size_correction = TRUE). With it,
#       08b reproduces NeEstimator's point estimates on simulated data
#       (n = 45, 94, 166) and on the goby files run under Wine.
#     * strataG's own Ne remains lower than 08b and NeEstimator. Part of the
#       gap is its MAF filtering, which does not behave as documented in
#       version 2.5.5; the rest we could not trace. We therefore report
#       NeEstimator-equivalent values from 08b and keep this script as a
#       record of the cross-check, not as the source of any reported number.
#
# Order of operations:
#   1. For each check group, read the Genepop file written by 07.
#   2. Draw n_replicates random subsets of n_loci_subset loci.
#   3. On each subset, compute 08b's all-pairs estimate (uncorrected r^2 for
#      the r^2 comparison; corrected r^2 for Ne) and strataG::ldNe.
#   4. Write one row per subset and log the Ne ratio.
#
# Install strataG (not on CRAN):
#   options(repos = c(zkamvar = "https://zkamvar.r-universe.dev",
#                     CRAN = "https://cloud.r-project.org"))
#   install.packages("strataG")
#
# Inputs:   results/07_ldne/genepop/<group>.gen (from 07)
# Outputs:  results/07_ldne/stratag/ldne_stratag_check.csv
#
# References:
#   Archer FI, Adams PE, Schneiders BB (2017) strataG: an R package for
#     manipulating, summarizing and analysing population genetic data.
#     Molecular Ecology Resources 17:5-11.
#   Waples RS, Do C (2008) LDNE: a program for estimating effective population
#     size from data on linkage disequilibrium. Molecular Ecology Resources
#     8:753-756.
# =============================================================================

source("scripts/00_config.R")
need_pkgs(c("strataG", "dplyr", "readr"))
suppressPackageStartupMessages({ library(dplyr); library(readr) })

log_start("08c_ldne_stratag")

genepop_dir <- project_path("results/07_ldne/genepop")
out_dir     <- project_path("results/07_ldne/stratag")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pcrit       <- ldne_settings$pcrit
block_size  <- 1500
source("scripts/lib/ldne_direct.R", local = TRUE)

check_groups  <- c("Cohort_2023", "Cohort_2024", "NearComplete_All")
n_loci_subset <- 1000
n_replicates  <- 3
n_cores       <- max(1, parallel::detectCores() - 1)

msg(sprintf("  strataG %s | %d replicate subsets of %d loci | groups: %s",
            as.character(packageVersion("strataG")), n_replicates, n_loci_subset,
            paste(check_groups, collapse = ", ")))

to_gtypes <- function(dosage, ids) {
  # two allele columns per locus: dosage 0 -> 1/1, 1 -> 1/2, 2 -> 2/2, NA -> NA/NA
  a1 <- ifelse(is.na(dosage), NA, ifelse(dosage == 2, 2L, 1L))
  a2 <- ifelse(is.na(dosage), NA, ifelse(dosage == 0, 1L, 2L))
  alleles <- matrix(NA_integer_, nrow(dosage), 2 * ncol(dosage))
  alleles[, seq(1, ncol(alleles), 2)] <- a1
  alleles[, seq(2, ncol(alleles), 2)] <- a2
  colnames(alleles) <- paste0(rep(paste0("L", seq_len(ncol(dosage))), each = 2), c(".1", ".2"))
  df <- data.frame(id = ids, stratum = "all", alleles, check.names = FALSE)
  strataG::df2gtypes(df, ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3)
}

rows <- list()
for (g in check_groups) {
  gen_file <- file.path(genepop_dir, paste0(g, ".gen"))
  if (!file.exists(gen_file)) { msg("  (", g, ".gen not found; skipped)"); next }
  dat <- read_genepop(gen_file)
  gen_lines <- readLines(gen_file)
  ids <- trimws(sub(",.*$", "", gen_lines[-(1:3)]))
  ids <- ids[nzchar(ids) & toupper(ids) != "POP"]

  for (replicate in seq_len(n_replicates)) {
    set.seed(random_seed + replicate)
    cols <- sort(sample(ncol(dat$dosage), n_loci_subset))
    sub_dosage <- dat$dosage[, cols, drop = FALSE]

    ours     <- ld_direct(sub_dosage, dat$loci[cols])$all
    ours_raw <- ld_direct(sub_dosage, dat$loci[cols], sample_size_correction = FALSE)$all
    t0 <- Sys.time()
    theirs <- strataG::ldNe(to_gtypes(sub_dosage, ids), maf.threshold = pcrit,
                            num.cores = n_cores)
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    theirs <- as.data.frame(theirs)[1, ]

    rows[[length(rows) + 1]] <- data.frame(
      group = g, replicate = replicate, n_fish = nrow(sub_dosage), loci = n_loci_subset,
      r2_08b_uncorrected = ours_raw[["r2"]], r2_08b_corrected = ours[["r2"]],
      r2_strataG = theirs$mean.rsq,
      E_r2_08b = ours[["E"]], E_r2_strataG = theirs$mean.E.rsq,
      Ne_08b = ours[["Ne"]], Ne_strataG = theirs$Ne,
      strataG_L95 = theirs$param.lci, strataG_U95 = theirs$param.uci,
      strataG_minutes = minutes)
    msg(sprintf("  %-17s rep %d: Ne 08b %7.1f | strataG %7.1f (%.1f-%.1f) | r2 %.6f vs %.6f  [%.1f min]",
                g, replicate, ours[["Ne"]], theirs$Ne, theirs$param.lci, theirs$param.uci,
                ours_raw[["r2"]], theirs$mean.rsq, minutes))
  }
}

out <- bind_rows(rows) %>% mutate(ratio = Ne_08b / Ne_strataG)
write_csv(out, file.path(out_dir, "ldne_stratag_check.csv"))
msg(sprintf("\n  08b / strataG Ne ratio: median %.3f, range %.3f-%.3f over %d comparisons",
            median(out$ratio), min(out$ratio), max(out$ratio), nrow(out)))
msg("  wrote ", out_dir)
msg("\n=== 08c_ldne_stratag.R complete ===")
log_end()
