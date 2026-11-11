# =============================================================================
# 08c_ldne_stratag.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Cross-check of 08b against an independent R implementation: strataG::ldNe
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-11
#
# Why this exists:
#   The Table 3 LD estimates come from our own code (08b, scripts/lib/
#   ldne_direct.R). A reader will ask how we know it is right. Three answers:
#     1. It recovers known Ne in simulations (true Ne = 100; 94-102 across
#        1,200-10,000 loci, random or clustered missing data; 07b_ldne_check.R).
#     2. It matches NeEstimator 2.1 run on the same files (all pairs; e.g.
#        near-complete set 503 vs 476.9, 2023 cohort 587 vs 516.3).
#     3. It matches strataG::ldNe (Archer et al. 2017), an independent,
#        published R implementation of the same Waples & Do (2008) estimator
#        "adapted from code by R. Waples and W. Larson". This script.
#
#   strataG computes r^2 for every locus pair in R, which is slow, so the
#   comparison uses replicate random locus subsets. Both methods see exactly
#   the same fish and loci. The comparison is on all locus pairs, which is
#   what strataG computes; excluding same-scaffold pairs (our primary
#   estimate) is a filter applied on top of the same r^2 values.
#
# Install strataG (not on CRAN):
#   options(repos = c(zkamvar = "https://zkamvar.r-universe.dev",
#                     CRAN = "https://cloud.r-project.org"))
#   install.packages("strataG")
#
# Inputs:   results/07_ldne/genepop/<group>.gen (from 07)
# Outputs:  results/07_ldne/stratag/ldne_stratag_check.csv
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
source("scripts/lib/ldne_direct.R")

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

    ours <- ld_direct(sub_dosage, dat$loci[cols])$all
    t0 <- Sys.time()
    theirs <- strataG::ldNe(to_gtypes(sub_dosage, ids), maf.threshold = pcrit,
                            num.cores = n_cores)
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    theirs <- as.data.frame(theirs)[1, ]

    rows[[length(rows) + 1]] <- data.frame(
      group = g, replicate = replicate, n_fish = nrow(sub_dosage), loci = n_loci_subset,
      r2_08b = ours[["r2"]], r2_strataG = theirs$mean.rsq,
      E_r2_08b = ours[["E"]], E_r2_strataG = theirs$mean.E.rsq,
      Ne_08b = ours[["Ne"]], Ne_strataG = theirs$Ne,
      strataG_L95 = theirs$param.lci, strataG_U95 = theirs$param.uci,
      strataG_minutes = minutes)
    msg(sprintf("  %-17s rep %d: Ne 08b %7.1f | strataG %7.1f (%.1f-%.1f) | r2 %.6f vs %.6f  [%.1f min]",
                g, replicate, ours[["Ne"]], theirs$Ne, theirs$param.lci, theirs$param.uci,
                ours[["r2"]], theirs$mean.rsq, minutes))
  }
}

out <- bind_rows(rows) %>% mutate(ratio = Ne_08b / Ne_strataG)
write_csv(out, file.path(out_dir, "ldne_stratag_check.csv"))
msg(sprintf("\n  08b / strataG Ne ratio: median %.3f, range %.3f-%.3f over %d comparisons",
            median(out$ratio), min(out$ratio), max(out$ratio), nrow(out)))
msg("  wrote ", out_dir)
msg("\n=== 08c_ldne_stratag.R complete ===")
log_end()
