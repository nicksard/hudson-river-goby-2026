# =============================================================================
# lib/ldne_direct.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Linkage-disequilibrium Ne computed directly from GENEPOP files
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-11 (moved out of 08b_ldne_direct.R)
# Updated:  2026-11-11 (sample-size correction of r^2; readable names)
#
# Used by:  08b_ldne_direct.R (Table 3), 08c_ldne_stratag.R (cross-check).
#           Both set `pcrit` and `block_size` before sourcing this file with
#           local = TRUE.
#
# Method (Waples & Do 2008; Waples et al. 2016):
#   1. Keep loci whose minor allele frequency in the group is >= pcrit.
#   2. For every locus pair, r^2 = squared Pearson correlation of allele
#      dosages over the fish genotyped at both loci (S of them). For
#      biallelic loci this is the composite r^2 used by LDNe/NeEstimator.
#   3. Multiply r^2 by (S / (S - 1))^2. This sample-size correction comes
#      from the Waples & Larson code adapted in strataG::ldNe (Archer et al.
#      2017); with it this estimator matches NeEstimator 2.1 run on the same
#      files (near-complete set 476.3 vs 476.9; 2023 cohort 503.6 vs 516.3).
#   4. Expected r^2 from sampling alone, per pair (Waples & Do 2008):
#      S >= 30: 1/S + 3.19/S^2;  S < 30: 0.0018 + 0.907/S + 4.44/S^2.
#   5. Drift r^2 (r^2') = S^2-weighted mean of (r^2 - expected r^2).
#   6. Ne from r^2' with the Waples & Do (2008) formula for S >= 30 or < 30.
#   7. 95% CI (primary): delete-a-group jackknife over FISH. Fish are split
#      at random into 10 groups; r^2' is recomputed with each group left out
#      and its jackknife standard error gives r^2' +/- 1.96 SE, converted to
#      Ne (a lower r^2' limit <= 0 gives an infinite upper Ne limit). Most of
#      the uncertainty in LD Ne comes from which fish were sampled; a
#      jackknife over loci or contigs alone was overconfident in simulations
#      (true Ne = 500: contig-jackknife CIs excluded the truth in 2 of 3
#      samples, fish-group CIs contained it in 15 of 16, with widths close
#      to NeEstimator's). The contig jackknife is kept as a secondary CI.
#   Two versions are returned: all locus pairs, and pairs on different
#   contigs (scaffolds) only, which removes physically linked pairs.
#   Work is done in blocks of loci, accumulating contig x contig sums, so
#   memory stays bounded for the full 10,943-locus panel.
#
# References:
#   Archer FI, Adams PE, Schneiders BB (2017) strataG: an R package for
#     manipulating, summarizing and analysing population genetic data.
#     Molecular Ecology Resources 17:5-11.
#   Waples RS, Do C (2008) LDNE: a program for estimating effective population
#     size from data on linkage disequilibrium. Molecular Ecology Resources
#     8:753-756.
#   Waples RK, Larson WA, Waples RS (2016) Estimating contemporary effective
#     population size in non-model species using linkage disequilibrium
#     across thousands of loci. Heredity 117:233-240.
# =============================================================================

# GENEPOP (6-digit, one population) -> fish x loci matrix of allele-2 dosages
read_genepop <- function(path) {
  gen_lines  <- readLines(path)
  locus_names <- trimws(strsplit(gen_lines[2], ",")[[1]])
  fish_lines <- gen_lines[-(1:3)]
  fish_lines <- fish_lines[nzchar(trimws(fish_lines)) & toupper(trimws(fish_lines)) != "POP"]
  genotype_codes <- do.call(rbind, strsplit(trimws(sub("^[^,]*,", "", fish_lines)), "\\s+"))
  need(ncol(genotype_codes) == length(locus_names),
       "Genotype rows do not match the locus line in ", basename(path))
  dosage <- matrix(c("001001" = 0, "001002" = 1, "002001" = 1, "002002" = 2,
                     "000000" = NA)[genotype_codes], nrow = nrow(genotype_codes))
  list(dosage = dosage, loci = locus_names)
}

# Waples & Do (2008) Ne from drift r^2
ne_from_r2_drift <- function(r2_drift, small_sample) {
  if (!is.finite(r2_drift) || r2_drift <= 0) return(Inf)
  if (small_sample) (0.308 + sqrt(0.308^2 - 2.08 * r2_drift)) / (2 * r2_drift)
  else              (1/3 + sqrt(1/9 - 2.76 * r2_drift)) / (2 * r2_drift)
}

ld_direct <- function(dosage, loci, sample_size_correction = TRUE) {
  allele_freq   <- colMeans(dosage, na.rm = TRUE) / 2
  passes_maf    <- !is.na(allele_freq) & pmin(allele_freq, 1 - allele_freq) >= pcrit
  genotypes     <- dosage[, passes_maf, drop = FALSE]
  locus_contig  <- factor(sub("_[0-9]+$", "", loci[passes_maf]))
  genotyped     <- !is.na(genotypes) * 1
  n_fish        <- nrow(genotypes)
  small_sample  <- n_fish < 30                   # Waples & Do (2008) S < 30 branch
  min_pair_n    <- if (small_sample) 2 else 30
  n_contigs     <- nlevels(locus_contig)

  # contig x contig running sums over locus pairs
  empty <- function() matrix(0, n_contigs, n_contigs)
  sums <- list(weight = empty(), weighted_r2 = empty(), weighted_expected = empty(),
               n_pairs = empty(), weight_over_n = empty())

  locus_blocks <- split(seq_len(ncol(genotypes)), ceiling(seq_len(ncol(genotypes)) / block_size))
  for (block_a in seq_along(locus_blocks)) for (block_b in block_a:length(locus_blocks)) {
    loci_a <- locus_blocks[[block_a]]; loci_b <- locus_blocks[[block_b]]
    correlation <- suppressWarnings(cor(genotypes[, loci_a, drop = FALSE],
                                        genotypes[, loci_b, drop = FALSE],
                                        use = "pairwise.complete.obs"))
    pair_n <- crossprod(genotyped[, loci_a, drop = FALSE], genotyped[, loci_b, drop = FALSE])
    usable <- is.finite(correlation) & pair_n >= min_pair_n
    if (block_a == block_b) usable <- usable & upper.tri(usable)      # each pair once
    r2 <- correlation^2 *                                             # step 3
      (if (sample_size_correction) (pair_n / (pair_n - 1))^2 else 1)
    expected_r2 <- if (small_sample) 0.0018 + 0.907 / pair_n + 4.44 / pair_n^2
                   else 1 / pair_n + 3.19 / pair_n^2
    weight <- ifelse(usable, pair_n^2, 0)

    sum_by_contig <- function(values) {          # rows: contigs of block a; cols: block b
      values[!usable] <- 0
      t(rowsum(t(rowsum(values, locus_contig[loci_a], reorder = FALSE)),
               locus_contig[loci_b], reorder = FALSE))
    }
    contigs_a <- unique(as.integer(locus_contig[loci_a]))
    contigs_b <- unique(as.integer(locus_contig[loci_b]))
    contributions <- list(weight = weight, weighted_r2 = weight * r2,
                          weighted_expected = weight * expected_r2,
                          n_pairs = usable * 1, weight_over_n = weight / pair_n)
    for (quantity in names(sums)) {
      sums[[quantity]][contigs_a, contigs_b] <- sums[[quantity]][contigs_a, contigs_b] +
        sum_by_contig(contributions[[quantity]])
    }
  }
  # symmetric totals: off-diagonal cells hold pairs between two contigs once
  symmetric_sums <- lapply(sums, function(m) { s <- m + t(m); diag(s) <- diag(m); s })

  summarise_pairs <- function(drop_contig = NULL, different_contigs_only) {
    total <- function(m) {
      all_pairs <- (sum(m) + sum(diag(m))) / 2
      kept <- if (different_contigs_only) all_pairs - sum(diag(m)) else all_pairs
      if (!is.null(drop_contig))
        kept <- kept - (sum(m[drop_contig, ]) - if (different_contigs_only) m[drop_contig, drop_contig] else 0)
      kept
    }
    total_weight <- total(symmetric_sums$weight)
    c(r2 = total(symmetric_sums$weighted_r2) / total_weight,
      E  = total(symmetric_sums$weighted_expected) / total_weight,
      r2p = (total(symmetric_sums$weighted_r2) - total(symmetric_sums$weighted_expected)) / total_weight,
      pairs = total(symmetric_sums$n_pairs),
      S_mean = total_weight / total(symmetric_sums$weight_over_n))
  }
  estimate_with_ci <- function(different_contigs_only) {
    full <- summarise_pairs(NULL, different_contigs_only)
    jackknife_r2p <- vapply(seq_len(n_contigs), function(contig)
      summarise_pairs(contig, different_contigs_only)[["r2p"]], numeric(1))
    jackknife_se <- sqrt((n_contigs - 1) / n_contigs *
                           sum((jackknife_r2p - mean(jackknife_r2p))^2))
    c(full,
      Ne     = ne_from_r2_drift(full[["r2p"]], small_sample),
      Ne_L95 = ne_from_r2_drift(full[["r2p"]] + 1.96 * jackknife_se, small_sample),
      Ne_U95 = ne_from_r2_drift(full[["r2p"]] - 1.96 * jackknife_se, small_sample))
  }
  list(all = estimate_with_ci(FALSE), diff = estimate_with_ci(TRUE),
       n_fish = n_fish, loci_used = ncol(genotypes), contigs = n_contigs,
       small_sample = small_sample)
}

# Step 7: ld_direct() plus delete-a-group jackknife CIs over fish. Adds
# Ne_L95_fish / Ne_U95_fish to both the all-pairs and different-contig results.
# Costs (n_groups + 1) runs of ld_direct().
ld_direct_fish_ci <- function(dosage, loci, n_groups = 10, seed = 1, ...) {
  full <- ld_direct(dosage, loci, ...)
  set.seed(seed)
  fish_group <- sample(rep_len(seq_len(n_groups), nrow(dosage)))
  leave_one_out <- lapply(seq_len(n_groups), function(k)
    ld_direct(dosage[fish_group != k, , drop = FALSE], loci, ...))
  for (version in c("all", "diff")) {
    r2p_without <- vapply(leave_one_out, function(res) res[[version]][["r2p"]], numeric(1))
    jackknife_se <- sqrt((n_groups - 1) / n_groups * sum((r2p_without - mean(r2p_without))^2))
    r2p <- full[[version]][["r2p"]]
    full[[version]] <- c(full[[version]],
      Ne_L95_fish = ne_from_r2_drift(r2p + 1.96 * jackknife_se, full$small_sample),
      Ne_U95_fish = ne_from_r2_drift(r2p - 1.96 * jackknife_se, full$small_sample))
  }
  full
}
