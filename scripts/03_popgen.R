# =============================================================================
# 03_popgen.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Genetic diversity and structure by birth cohort and by river location
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-09-08
# Updated:  2026-09-10
#
# Description:
#   Runs per-group diversity and pairwise Weir & Cockerham Fst under every
#   grouping listed in popgen_groupings (00_config.R). One code path serves all
#   of them, so a cohort result and a location result differ only in how
#   individuals were labelled, and the two can be compared directly.
#
#   The PCA is computed ONCE, on the full analysis set, and each grouping
#   recolours it. PCA is unsupervised -- the labels play no part in the
#   calculation -- so there is one ordination of these individuals and several
#   ways of colouring it, not several analyses.
#
#   Two analytical choices worth stating plainly:
#
#   1. Grouping is by BIRTH COHORT, not collection year. A collection year
#      contains both Age-0 and Age-1 fish, so consecutive collection years
#      share spawners. That shared drift history deflates pairwise Fst between
#      adjacent years, and anything derived from that Fst inherits the bias.
#
#   2. Fis is reported as a ratio of sums, 1 - sum(Ho)/sum(Hs), across loci --
#      not as the median of the per-locus Fis values. Per-locus Fis is badly
#      behaved at low minor allele frequency (the denominator approaches zero),
#      so the distribution of per-locus values has heavy tails and its median
#      is pulled around by sample size rather than by inbreeding. The ratio of
#      sums weights each locus by how much information it actually carries.
#      The per-locus values are still plotted, so the spread is visible.
#
#   A group with fewer than cohort_settings$min_n_estimator individuals is
#   analysed and reported, but flagged. Fst between small groups is dominated by
#   sampling noise, and the noise is not centred on zero in a way that is
#   obvious from the point estimate alone -- comparing the magnitude of the
#   negative Fst estimates in the matrix to the positive ones is the quickest
#   check of whether any of the positives mean anything.
#
# Inputs:
#   paths$filtered_vcf   output of 01_filter_vcf.R
#   paths$sample_meta    output of 02_sample_metadata.R
#
# Outputs, per grouping, under results/03_popgen/<grouping name>/:
#   diversity.csv        N, mean Ho, mean Hs, Fis with bootstrap CI
#   pairwise_fst.csv     lower triangle of the WC Fst matrix
#   pairwise_fst.png
#   pca_scores.csv       shared ordination scores, labelled by this grouping
#   pca_pc1_pc2.png
#   popgen_report_<grouping>.pdf
#
# Outputs shared across groupings, under results/03_popgen/:
#   pca_variance.csv     percent variance on the first twenty axes
#   pca_scores.csv       the shared ordination itself
#   pca_outliers.csv     individuals beyond 4 SD on a leading axis
#   grouping_overview.csv
#
# Usage:
#   Rscript scripts/03_popgen.R
# =============================================================================

source("scripts/00_config.R")
need_pkgs(c("vcfR", "adegenet", "hierfstat", "dplyr", "tidyr", "tibble",
            "readr", "purrr", "ggplot2", "gridExtra"))

suppressPackageStartupMessages({
  library(vcfR); library(adegenet); library(hierfstat); library(dplyr)
  library(tidyr); library(tibble); library(readr); library(purrr)
  library(ggplot2); library(gridExtra); library(grid)
})

log_start("03_popgen")

# Bootstrapping pairwise Fst resamples loci independently. Over the two to three
# generations spanned here recombination has barely begun to break up the 23
# chromosomes, so ~11,000 loci are nowhere near independent and the resulting
# intervals are far too narrow. Off by default; if turned on, report the
# intervals as a lower bound on uncertainty and say so in the text.
run_fst_bootstrap <- FALSE
n_bootstrap_reps  <- 200

# Permutation test for pairwise Fst: shuffle group labels among individuals and
# recompute, giving the distribution of Fst when there is genuinely nothing to
# find. This is the defensible version of the negative-value noise floor, and
# it is what a reviewer will ask for if the paper leans on any Fst value.
#
# Run on a random locus subset because permuting the full ~11,000-locus panel a
# few hundred times is hours of work. Subsetting loci widens the null slightly
# but does not invalidate it, and the observed statistic is computed on the same
# subset so the comparison stays like-for-like. Off by default; turn it on for
# the run that produces the final tables.
run_fst_permutation <- TRUE
n_permutations      <- 200
n_permutation_loci  <- 2000

need(file.exists(paths$filtered_vcf), "Run scripts/01_filter_vcf.R first.")
need(file.exists(paths$sample_meta),  "Run scripts/02_sample_metadata.R first.")


# =============================================================================
# 1. Load genotypes and metadata, and align them
# =============================================================================
goby_vcf <- read.vcfR(paths$filtered_vcf, verbose = FALSE)

# Sample IDs look like 2021001 and would otherwise be read as numbers, which
# then fail to match the character sample names carried by the genind object.
sample_info <- read_csv(paths$sample_meta, col_types = sample_id_col_types) %>%
  filter(analysis_keep) %>%
  mutate(sample = as.character(sample),
         cohort = as.character(cohort))

msg(sprintf("  VCF %d loci x %d samples | metadata analysis set: %d",
            nrow(goby_vcf), ncol(goby_vcf@gt) - 1L, nrow(sample_info)))

genind_all <- vcfR2genind(goby_vcf)

# Intersect in both directions, then reorder the metadata to match the genind.
# Order matters because pop() is assigned positionally: if the two ever drift
# apart, every individual is silently given its neighbour's group label.
genind_all  <- genind_all[indNames(genind_all) %in% sample_info$sample, ]
sample_info <- sample_info[match(indNames(genind_all), sample_info$sample), ]

need(identical(indNames(genind_all), sample_info$sample),
     "The genind object and the metadata table are not in the same order.")
msg("  analysis set after intersecting with the VCF: ", nInd(genind_all), " individuals")


# =============================================================================
# 2. Helpers
# =============================================================================

# Ratio of sums across loci. Equals Weir & Cockerham's f to about four decimal
# places for a single population, and unlike the median of per-locus values it
# is not pulled around by low-MAF loci.
fis_ratio_of_sums <- function(observed_het, expected_het) {
  1 - sum(observed_het, na.rm = TRUE) / sum(expected_het, na.rm = TRUE)
}

# Turns a grouping specification from popgen_groupings into one label per
# individual. NA means "not assignable under this scheme", which happens when a
# fish has no region and the scheme is location-based.
build_group_labels <- function(sample_info, grouping_spec) {
  location_labels <- if (is.null(grouping_spec$location_scheme)) {
    NULL
  } else {
    sample_info[[paste0("loc_", grouping_spec$location_scheme)]]
  }

  switch(grouping_spec$by,
    cohort   = as.character(sample_info$cohort),
    location = location_labels,
    both     = ifelse(is.na(location_labels) | is.na(sample_info$cohort),
                      NA_character_,
                      paste(sample_info$cohort, location_labels, sep = "_")),
    stop("Unknown `by` value in popgen_groupings: ", grouping_spec$by))
}


# =============================================================================
# 3. The analysis, run once per grouping
# =============================================================================

analyse_grouping <- function(grouping_spec) {

  msg("\n", strrep("=", 72))
  msg("  GROUPING: ", grouping_spec$name, "   (by = ", grouping_spec$by,
      if (!is.null(grouping_spec$location_scheme))
        paste0(", location = ", grouping_spec$location_scheme) else "", ")")
  msg(strrep("=", 72))

  output_dir <- file.path(paths$dir_popgen, grouping_spec$name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  group_labels    <- build_group_labels(sample_info, grouping_spec)
  keep_individual <- !is.na(group_labels)

  n_unassigned <- sum(is.na(group_labels))
  if (n_unassigned > 0) {
    msg("  ", n_unassigned, " individual(s) unassigned under this scheme and excluded")
  }

  label_counts <- table(group_labels[keep_individual])
  tiny_groups  <- names(label_counts)[label_counts < cohort_settings$min_n_analysis]
  if (length(tiny_groups)) {
    msg("  dropping group(s) with n < ", cohort_settings$min_n_analysis, ": ",
        paste(tiny_groups, collapse = ", "))
    keep_individual <- keep_individual & !(group_labels %in% tiny_groups)
  }

  genind_group <- genind_all[keep_individual, ]
  group_factor <- factor(group_labels[keep_individual])
  pop(genind_group) <- group_factor

  group_sizes <- as.data.frame(table(group_factor), responseName = "N")
  names(group_sizes)[1] <- "group"

  msg("  ", nlevels(group_factor), " groups, ", nInd(genind_group), " individuals:")
  msg_table(group_sizes)

  well_powered <- group_sizes$group[group_sizes$N >= cohort_settings$min_n_estimator]
  underpowered <- group_sizes$group[group_sizes$N <  cohort_settings$min_n_estimator]
  if (length(underpowered)) {
    msg("  !! below n = ", cohort_settings$min_n_estimator,
        ", so Fst and anything derived from it are low-precision here: ",
        paste(underpowered, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Diversity
  # ---------------------------------------------------------------------------
  hierfstat_data <- genind2hierfstat(genind_group)
  basic_stats    <- basic.stats(hierfstat_data, diploid = TRUE)

  observed_het_long <- as.data.frame(basic_stats$Ho) %>%
    rownames_to_column("locus") %>%
    pivot_longer(-locus, names_to = "group", values_to = "Ho")
  expected_het_long <- as.data.frame(basic_stats$Hs) %>%
    rownames_to_column("locus") %>%
    pivot_longer(-locus, names_to = "group", values_to = "Hs")
  heterozygosity_long <- left_join(observed_het_long, expected_het_long,
                                   by = c("locus", "group"))

  diversity_summary <- heterozygosity_long %>%
    group_by(group) %>%
    summarise(mean_Ho = mean(Ho, na.rm = TRUE),
              mean_Hs = mean(Hs, na.rm = TRUE),
              Fis     = fis_ratio_of_sums(Ho, Hs),
              .groups = "drop") %>%
    left_join(group_sizes %>% mutate(group = as.character(group)), by = "group") %>%
    relocate(N, .after = group)

  # Per-group Fis with bootstrap confidence intervals. This one is cheap because
  # it does not loop over pairs, so it runs by default unlike the pairwise Fst
  # bootstrap. Same caveat about locus non-independence applies to its width.
  diversity_summary <- tryCatch({
    ### CHANGED (patch 14): bootstrap one group at a time. boot.ppfis returns
    ### its CI rows in its own population order, which is not
    ### levels(group_factor); labelling them positionally attached CIs to the
    ### wrong groups (e.g. reach: Above_dam 0.0046 with CI 0.0077-0.0105).
    ### Rows of hierfstat_data follow the individuals of genind_group, which
    ### are the individuals of group_factor in the same order.
    need(nrow(hierfstat_data) == length(group_factor),
         "hierfstat_data and group_factor have different lengths.")
    fis_intervals <- bind_rows(lapply(levels(group_factor), function(level) {
      rows <- which(group_factor == level)
      fis_bootstrap <- boot.ppfis(dat = hierfstat_data[rows, , drop = FALSE],
                                  nboot = n_bootstrap_reps, quant = c(0.025, 0.975))
      tibble(group = level, Fis_L95 = fis_bootstrap$fis.ci[1, 1],
             Fis_U95 = fis_bootstrap$fis.ci[1, 2])
    }))
    out <- left_join(diversity_summary, fis_intervals, by = "group")
    outside <- out$group[!is.na(out$Fis_L95) & (out$Fis < out$Fis_L95 | out$Fis > out$Fis_U95)]
    if (length(outside) > 0)
      msg("  !! Fis outside its own bootstrap CI (check): ", paste(outside, collapse = ", "))
    out
  }, error = function(e) {
    msg("  (per-group Fis bootstrap skipped: ", conditionMessage(e), ")")
    diversity_summary
  })

  write_csv(diversity_summary, file.path(output_dir, "diversity.csv"))
  msg("\n  diversity:")
  msg_table(as.data.frame(diversity_summary %>%
                        mutate(across(where(is.numeric), ~round(.x, 4)))))

  plot_heterozygosity <- heterozygosity_long %>%
    pivot_longer(c(Ho, Hs), names_to = "statistic", values_to = "value") %>%
    filter(!is.na(value)) %>%
    ggplot(aes(group, value, fill = statistic)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_manual(name = "",
                      values = c(Ho = "cornflowerblue", Hs = "lightgoldenrodyellow"),
                      labels = c("Observed (Ho)", "Expected (Hs)")) +
    labs(title = paste0("Heterozygosity (", grouping_spec$name, ")"),
         x = NULL, y = "Heterozygosity") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Loci that are monomorphic within a group have an undefined Fis. They are
  # dropped here rather than left to raise a warning from the plotting layer,
  # and the count is reported so the reader knows how many were set aside.
  per_locus_fis <- as.data.frame(basic_stats$Fis) %>%
    rownames_to_column("locus") %>%
    pivot_longer(-locus, names_to = "group", values_to = "Fis")
  n_undefined_fis <- sum(!is.finite(per_locus_fis$Fis))
  if (n_undefined_fis > 0) {
    msg(sprintf("  %d of %d per-locus Fis values are undefined (locus monomorphic within group)",
                n_undefined_fis, nrow(per_locus_fis)))
  }

  plot_per_locus_fis <- per_locus_fis %>%
    filter(is.finite(Fis)) %>%
    ggplot(aes(group, Fis, fill = group)) +
    geom_boxplot(outlier.alpha = 0.3) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
    labs(title = paste0("Per-locus Fis (", grouping_spec$name, ")"),
         subtitle = "point estimates in diversity.csv are the ratio of sums, not the median of these",
         x = NULL, y = "Fis") +
    theme_bw() + theme(legend.position = "none",
                       axis.text.x = element_text(angle = 45, hjust = 1))

  # ---------------------------------------------------------------------------
  # Pairwise Fst
  # ---------------------------------------------------------------------------
  plot_fst <- NULL
  fst_pairs <- NULL

  if (nlevels(group_factor) >= 2) {

    fst_matrix <- as.matrix(pairwise.WCfst(hierfstat_data))
    fst_matrix[upper.tri(fst_matrix)] <- NA     # keep one triangle only
    fst_pairs <- as.data.frame(fst_matrix) %>%
      rownames_to_column("group1") %>%
      pivot_longer(-group1, names_to = "group2", values_to = "Fst") %>%
      filter(!is.na(Fst))

    if (run_fst_bootstrap) {
      msg("  bootstrapping pairwise Fst (", n_bootstrap_reps, " reps; slow) ...")
      fst_bootstrap <- boot.ppfst(dat = hierfstat_data, nboot = n_bootstrap_reps,
                                  quant = c(0.025, 0.975))
      # The bootstrap fills one triangle, which is not necessarily the same
      # triangle kept above, so both orientations are joined and coalesced.
      fst_intervals <- expand.grid(group1 = rownames(fst_bootstrap$ll),
                                   group2 = colnames(fst_bootstrap$ll),
                                   stringsAsFactors = FALSE) %>%
        mutate(Fst_L95 = as.vector(fst_bootstrap$ll),
               Fst_U95 = as.vector(fst_bootstrap$ul)) %>%
        filter(!is.na(Fst_L95))
      fst_pairs <- fst_pairs %>%
        left_join(fst_intervals, by = c("group1", "group2")) %>%
        left_join(fst_intervals %>% rename(group1 = group2, group2 = group1,
                                           L95_flipped = Fst_L95,
                                           U95_flipped = Fst_U95),
                  by = c("group1", "group2")) %>%
        mutate(Fst_L95 = coalesce(Fst_L95, L95_flipped),
               Fst_U95 = coalesce(Fst_U95, U95_flipped)) %>%
        select(-L95_flipped, -U95_flipped)
    }

    write_csv(fst_pairs, file.path(output_dir, "pairwise_fst.csv"))
    msg("\n  pairwise Fst:")
    msg_table(as.data.frame(fst_pairs %>% mutate(across(where(is.numeric), ~signif(.x, 4)))))

    # A negative Weir & Cockerham Fst estimate means the between-group variance
    # component came out below zero, which cannot be true and simply says the
    # groups are indistinguishable. The largest negative value in the matrix is
    # therefore a rough empirical noise floor: any positive estimate smaller
    # than its magnitude should not be read as differentiation.
    #
    # Rough is the operative word. The floor is set by whichever pair is
    # noisiest, and that is usually a pair involving the smallest group -- so a
    # single tiny group can inflate the floor and make everything look like
    # noise. The floor is therefore computed twice: over all pairs, and over
    # pairs where both groups clear min_n_estimator. The second is the one to
    # quote. Neither is a substitute for the permutation test below.
    report_noise_floor <- function(pairs, label) {
      if (nrow(pairs) == 0) return(invisible(NULL))
      floor_value <- min(pairs$Fst)
      if (floor_value >= 0) {
        msg(sprintf("  %-22s no negative estimates, so no empirical floor", label))
        return(invisible(NULL))
      }
      n_above <- sum(pairs$Fst > abs(floor_value))
      msg(sprintf("  %-22s floor |%.5f| -> %d of %d pairs exceed it",
                  label, floor_value, n_above, nrow(pairs)))
    }

    msg("\n  noise-floor check (magnitude of the most negative estimate):")
    report_noise_floor(fst_pairs, "all pairs")
    report_noise_floor(fst_pairs %>% filter(group1 %in% well_powered,
                                            group2 %in% well_powered),
                       "well-powered pairs")

    # Permutation test. Shuffling group labels among individuals and recomputing
    # gives the distribution of Fst under no differentiation, which is what the
    # negative-value heuristic is only gesturing at. Run on a random subset of
    # loci because the full panel makes this far too slow, and the subset costs
    # precision in the null rather than validity.
    if (run_fst_permutation) {
      msg(sprintf("\n  permutation test: %d reps on %d random loci ...",
                  n_permutations, n_permutation_loci))
      locus_subset <- sample(nLoc(genind_group),
                             min(n_permutation_loci, nLoc(genind_group)))
      genind_subset <- genind_group[loc = locus_subset]

      observed_fst <- as.matrix(pairwise.WCfst(genind2hierfstat(genind_subset)))

      permuted_fst <- vapply(seq_len(n_permutations), function(rep_i) {
        shuffled <- genind_subset
        pop(shuffled) <- sample(pop(genind_subset))
        as.vector(as.matrix(pairwise.WCfst(genind2hierfstat(shuffled))))
      }, numeric(nlevels(group_factor)^2))

      # One-sided: only an excess of differentiation is evidence of structure.
      p_values <- vapply(seq_len(nlevels(group_factor)^2), function(cell_i) {
        observed <- as.vector(observed_fst)[cell_i]
        if (is.na(observed)) return(NA_real_)
        (1 + sum(permuted_fst[cell_i, ] >= observed, na.rm = TRUE)) / (n_permutations + 1)
      }, numeric(1))

      p_matrix <- matrix(p_values, nrow = nlevels(group_factor),
                         dimnames = dimnames(observed_fst))
      p_long <- as.data.frame(p_matrix) %>%
        rownames_to_column("group1") %>%
        pivot_longer(-group1, names_to = "group2", values_to = "p_perm") %>%
        filter(!is.na(p_perm))
      fst_pairs <- fst_pairs %>% left_join(p_long, by = c("group1", "group2"))

      msg("  permutation p-values (one-sided; p > 0.05 means indistinguishable):")
      msg_table(as.data.frame(fst_pairs %>%
                            mutate(Fst = signif(Fst, 3), p_perm = round(p_perm, 3))))
      write_csv(fst_pairs, file.path(output_dir, "pairwise_fst.csv"))
    }

    # Pull out the contrasts where only one factor varies. Everything else in
    # the matrix confounds cohort with location and should not be read alone.
    if (exists("clean_contrasts") && grouping_spec$name == "cohort_by_site") {
      msg("\n  --- unconfounded contrasts ---")
      for (contrast in clean_contrasts) {
        this_pair <- fst_pairs %>%
          filter((group1 == contrast$a & group2 == contrast$b) |
                 (group1 == contrast$b & group2 == contrast$a))
        if (nrow(this_pair) == 1) {
          msg(sprintf("    %-24s Fst = %9.6f   %s%s", contrast$label,
                      this_pair$Fst[1], contrast$type,
                      if ("p_perm" %in% names(this_pair))
                        sprintf("   p = %.3f", this_pair$p_perm[1]) else ""))
        } else {
          msg(sprintf("    %-24s not available (a group fell below the size cutoff)",
                      contrast$label))
        }
      }
    }

    plot_fst <- ggplot(fst_pairs,
                       aes(group1, group2, fill = Fst, label = sprintf("%.4f", Fst))) +
      geom_tile(colour = "black") +
      geom_text(size = 4) +
      scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick",
                           midpoint = 0, name = "Fst") +
      labs(title = paste0("Pairwise Fst (", grouping_spec$name, ")"),
           subtitle = "diverging scale centred on zero: blue cells are negative estimates, i.e. no differentiation",
           x = NULL, y = NULL) +
      theme_bw() +
      theme(axis.text = element_text(size = 12, colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(output_dir, "pairwise_fst.png"), plot_fst,
           width = 11, height = 8, dpi = 300)
  } else {
    msg("  only one group under this scheme; Fst skipped.")
  }

  # ---------------------------------------------------------------------------
  # PCA
  # ---------------------------------------------------------------------------
  # The ordination itself is computed ONCE, on the full analysis set, before any
  # grouping is applied (see section 3 below). PCA is unsupervised: group labels
  # do not enter the calculation, so recomputing it per grouping produces either
  # the identical ordination with different colours, or -- worse -- a slightly
  # different one because a few individuals were unassignable under that scheme,
  # which then reads as though several separate analyses had been run.
  #
  # So each grouping reuses the shared scores and only recolours them.
  # Individuals with no label under this scheme are drawn in grey rather than
  # dropped, so the reader can see how much of the sample the scheme covers.
  pca_scores <- shared_pca$scores %>%
    left_join(tibble(sample = sample_info$sample, group = group_labels),
              by = "sample") %>%
    mutate(group = ifelse(is.na(group), "unassigned", group))
  variance_explained <- shared_pca$variance

  write_csv(pca_scores, file.path(output_dir, "pca_scores.csv"))

  pca_panel <- function(x_axis, y_axis) {
    axis_number <- function(axis_name) as.integer(substr(axis_name, 5, 5))
    labelled <- pca_scores %>% filter(group != "unassigned")
    group_centroids <- labelled %>%
      group_by(group) %>%
      summarise(x = mean(.data[[x_axis]]), y = mean(.data[[y_axis]]), .groups = "drop")
    ggplot(labelled, aes(.data[[x_axis]], .data[[y_axis]], colour = group)) +
      geom_point(data = pca_scores %>% filter(group == "unassigned"),
                 colour = "grey80", size = 2, alpha = 0.5, show.legend = FALSE) +
      geom_point(size = 2.5, alpha = 0.35) +
      stat_ellipse(type = "t", level = 0.95, linetype = "dashed",
                   linewidth = 0.7, show.legend = FALSE) +
      geom_point(data = group_centroids, aes(x, y), size = 5, shape = 18,
                 show.legend = FALSE) +
      labs(title = paste0("PCA coloured by ", grouping_spec$name),
           subtitle = "one shared ordination of the full analysis set; grey = unassigned under this scheme",
           x = sprintf("PC%d (%.1f%%)", axis_number(x_axis),
                       variance_explained[axis_number(x_axis)]),
           y = sprintf("PC%d (%.1f%%)", axis_number(y_axis),
                       variance_explained[axis_number(y_axis)])) +
      theme_bw()
  }
  plot_pc1_pc2 <- pca_panel("Axis1", "Axis2")
  plot_pc1_pc3 <- pca_panel("Axis1", "Axis3")
  plot_pc2_pc3 <- pca_panel("Axis2", "Axis3")
  ggsave(file.path(output_dir, "pca_pc1_pc2.png"), plot_pc1_pc2,
         width = 8, height = 6, dpi = 300)

  # ---------------------------------------------------------------------------
  # Report
  # ---------------------------------------------------------------------------
  pdf(file.path(output_dir, paste0("popgen_report_", grouping_spec$name, ".pdf")),
      width = 14, height = 11)
  grid.arrange(
    plot_heterozygosity, plot_per_locus_fis,
    if (!is.null(plot_fst)) plot_fst else nullGrob(), plot_pc1_pc2,
    ncol = 2,
    top = textGrob(sprintf("%s -- %d loci, %d individuals, %d groups",
                           grouping_spec$name, nLoc(genind_group),
                           nInd(genind_group), nlevels(group_factor)),
                   gp = gpar(fontsize = 14))
  )
  grid.arrange(plot_pc1_pc3, plot_pc2_pc3, ncol = 2,
               top = textGrob("Supplementary PCA axes", gp = gpar(fontsize = 14)))
  grid.arrange(
    tableGrob(as.data.frame(diversity_summary %>%
                              mutate(across(where(is.numeric), ~round(.x, 4)))),
              rows = NULL),
    top = textGrob("Diversity summary", gp = gpar(fontsize = 14))
  )
  dev.off()

  msg("  wrote ", output_dir)

  # PCA columns are deliberately absent here: the ordination is shared across
  # groupings, so its variance decomposition belongs to the run, not to any one
  # grouping. It is reported once, above, and in pca_variance.csv.
  well_powered_pairs <- if (!is.null(fst_pairs)) {
    fst_pairs %>% filter(group1 %in% well_powered, group2 %in% well_powered)
  } else NULL

  tibble(grouping        = grouping_spec$name,
         n_groups        = nlevels(group_factor),
         n_indiv         = nInd(genind_group),
         n_unassigned    = n_unassigned,
         n_loci          = nLoc(genind_group),
         min_fst         = if (!is.null(fst_pairs)) min(fst_pairs$Fst) else NA_real_,
         max_fst         = if (!is.null(fst_pairs)) max(fst_pairs$Fst) else NA_real_,
         n_pairs_powered = if (!is.null(well_powered_pairs)) nrow(well_powered_pairs) else NA_integer_,
         min_fst_powered = if (!is.null(well_powered_pairs) && nrow(well_powered_pairs) > 0)
           min(well_powered_pairs$Fst) else NA_real_,
         max_fst_powered = if (!is.null(well_powered_pairs) && nrow(well_powered_pairs) > 0)
           max(well_powered_pairs$Fst) else NA_real_)
}


# =============================================================================
# 4. One shared ordination
# =============================================================================
# Computed once, on the full analysis set, and reused by every grouping. PCA is
# unsupervised, so the labels play no part in it -- there is one ordination of
# these individuals and several ways of colouring it, not several analyses.
#
# Centred but not scaled, which is standard for allele-frequency data: scaling
# would upweight rare variants, which are the noisiest. Missing genotypes are
# replaced with the locus mean, which pulls individuals with a lot of missing
# data toward the origin, so PCA outliers are genotype outliers rather than
# missing-data artefacts.
msg("\n--- shared PCA on the full analysis set ---")

scaled_genotypes <- scaleGen(genind_all, center = TRUE, scale = FALSE,
                             NA.method = "mean")
pca_result <- dudi.pca(scaled_genotypes, cent = FALSE, scale = FALSE,
                       nf = 3, scannf = FALSE)

shared_pca <- list(
  variance = 100 * pca_result$eig / sum(pca_result$eig),
  scores   = as_tibble(pca_result$li) %>% mutate(sample = indNames(genind_all))
)

msg(sprintf("  %d individuals, %d loci | PC1 %.2f%% | PC2 %.2f%% | PC3 %.2f%%",
            nInd(genind_all), nLoc(genind_all), shared_pca$variance[1],
            shared_pca$variance[2], shared_pca$variance[3]))

# With this many loci and no real structure the leading axes carry barely more
# variance than the ones behind them, so PC1 on its own is uninterpretable. The
# ratio to a mid-pack noise axis is the number to quote: real structure puts PC1
# many times above the noise floor, drift-free noise puts it just above.
msg(sprintf("  axis 20 carries %.2f%%, so PC1 is %.2fx the noise floor.",
            shared_pca$variance[20],
            shared_pca$variance[1] / shared_pca$variance[20]))
msg("  A ratio near 1 means the leading axes are noise, not structure.")

write_csv(tibble(axis = seq_along(shared_pca$variance),
                 pct_variance = shared_pca$variance) %>% slice_head(n = 20),
          file.path(paths$dir_popgen, "pca_variance.csv"))
write_csv(shared_pca$scores, file.path(paths$dir_popgen, "pca_scores.csv"))

# Flag individuals sitting far off the cloud. Because missing data pulls toward
# the origin, a large score is a genotype signal, not a coverage artefact --
# worth passing through the same allele-balance check that identified the
# contaminated sample before trusting it.
outlier_threshold <- 4
pca_outliers <- shared_pca$scores %>%
  mutate(across(starts_with("Axis"), ~abs(.x - mean(.x)) / sd(.x), .names = "z_{.col}")) %>%
  filter(if_any(starts_with("z_Axis"), ~.x > outlier_threshold))
if (nrow(pca_outliers) > 0) {
  msg(sprintf("\n  !! %d individual(s) beyond %d SD on a leading axis:",
              nrow(pca_outliers), outlier_threshold))
  msg_table(as.data.frame(pca_outliers %>%
                        mutate(across(where(is.numeric), ~round(.x, 2)))))
  msg("  Run the allele-balance diagnostic on these before interpreting them.")
  write_csv(pca_outliers, file.path(paths$dir_popgen, "pca_outliers.csv"))
}


# =============================================================================
# 5. Run every grouping and summarise across them
# =============================================================================
grouping_overview <- map_dfr(popgen_groupings, analyse_grouping)

msg("\n", strrep("=", 72))
msg("  OVERVIEW")
msg(strrep("=", 72))
msg_table(as.data.frame(grouping_overview %>%
                      mutate(across(where(is.numeric), ~signif(.x, 4)))))
write_csv(grouping_overview, file.path(paths$dir_popgen, "grouping_overview.csv"))

msg("\n=== 03_popgen.R complete ===")
log_end()
