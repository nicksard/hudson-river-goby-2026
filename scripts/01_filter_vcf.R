# =============================================================================
# 01_filter_vcf.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Raw Stacks VCF -> final filtered SNP panel
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-09-08
# Updated:  2026-09-10
#
# Description:
#   Takes the VCF written by the Stacks populations module and applies a
#   loci-first filtering chain, ending in the SNP panel every downstream
#   analysis uses. Loci-first means the cheap, decisive locus filters (MAF,
#   thinning to one SNP per RAD locus) run before the expensive per-genotype
#   work, which cuts the dataset by roughly 90% before anything slow happens.
#
#   Order of operations, and why it is this order:
#     1. Diagnose index hopping in the negative controls, then remove them.
#        The diagnosis has to happen while the controls are still in the file.
#     2. MAF, then thin to one SNP per RAD locus. Both are locus-level and
#        neither depends on genotype quality, so they go first.
#     3. HWE and paralog screens, on the thinned set. These need per-locus
#        genotype counts, so they are slower and benefit from the smaller panel.
#     4. Mask low-quality genotypes, then cull individuals, then cull loci.
#        Individuals before loci: a locus should be judged on the samples that
#        are actually going to be analysed, not on samples about to be dropped.
#     5. Remove loci that became monomorphic once individuals were culled.
#
#   Every step appends a row to a filtering ledger, written out as
#   filtering_summary.csv. That table is the Methods filtering paragraph and
#   Table 1; it is not retyped by hand anywhere.
#
# Inputs:
#   paths$raw_vcf   raw Stacks VCF (or paths$test_vcf with --test)
#
# Outputs:
#   results/01_filtering/populations.filtered.final.vcf.gz
#   results/01_filtering/filtering_summary.csv
#   results/01_filtering/hwe_per_locus.csv
#   results/01_filtering/paralog_per_locus.csv
#   results/01_filtering/missingness_per_individual.csv
#   results/01_filtering/missingness_summary_by_year.csv
#   results/01_filtering/ntc_index_hopping_stats.csv
#   results/01_filtering/VCF_filtering_report.pdf
#
# Usage:
#   Rscript scripts/01_filter_vcf.R            full run
#   Rscript scripts/01_filter_vcf.R --test     smoke test on the 2,000-locus subset
# =============================================================================

source("scripts/00_config.R")
need_pkgs(c("vcfR", "adegenet", "hierfstat", "pegas", "dplyr", "tidyr",
            "tibble", "readr", "stringr", "purrr", "ggplot2", "gridExtra"))

suppressPackageStartupMessages({
  library(vcfR); library(dplyr); library(tidyr); library(tibble)
  library(readr); library(stringr); library(purrr); library(ggplot2)
  library(adegenet); library(hierfstat); library(gridExtra)
})

log_start("01_filter_vcf")

command_args <- commandArgs(trailingOnly = TRUE)
run_on_test_subset <- "--test" %in% command_args
input_vcf_path <- if (run_on_test_subset) paths$test_vcf else paths$raw_vcf

if (run_on_test_subset) {
  msg("!! --test: running on ", basename(paths$test_vcf), " (2,000 loci).")
  msg("!! This exercises the code path only. Do not report these numbers.")
}
need(file.exists(input_vcf_path),
     "Input VCF not found: ", input_vcf_path, "\n",
     "  Set the RAW_VCF environment variable, or put the file in data/raw/.")

diagnostic_plots <- list()   # collected as we go, written to one PDF at the end


# -----------------------------------------------------------------------------
# Filtering ledger
# -----------------------------------------------------------------------------
# One row per filtering step. Kept as a running object rather than reconstructed
# at the end so that a run that dies halfway still leaves a record of how far it
# got and what it had removed.

filtering_ledger <- tibble(
  step = character(), filter = character(),
  loci_removed = integer(),  loci_remaining = integer(),
  indiv_removed = integer(), indiv_remaining = integer()
)

count_loci        <- function(vcf) nrow(vcf)
count_individuals <- function(vcf) ncol(vcf@gt) - 1L   # column 1 is FORMAT

record_step <- function(step, filter, vcf, loci_removed = 0L, indiv_removed = 0L) {
  filtering_ledger <<- add_row(
    filtering_ledger,
    step = step, filter = filter,
    loci_removed   = as.integer(loci_removed),  loci_remaining  = count_loci(vcf),
    indiv_removed  = as.integer(indiv_removed), indiv_remaining = count_individuals(vcf)
  )
  msg(sprintf("  %-22s %-38s loci %7d (-%6d) | indiv %4d (-%3d)",
              step, filter, count_loci(vcf), loci_removed,
              count_individuals(vcf), indiv_removed))
}


# =============================================================================
# 1. Load, rename the mislabelled control well, diagnose index hopping
# =============================================================================
msg("\n--- 1. Load and pre-filter ---")

goby_vcf <- read.vcfR(input_vcf_path, verbose = FALSE)
msg(sprintf("  raw: %d variants x %d samples",
            nrow(goby_vcf), ncol(goby_vcf@gt) - 1L))
record_step("Initial state", "none", goby_vcf)

for (old_name in names(control_well_renames)) {
  if (old_name %in% colnames(goby_vcf@gt)) {
    colnames(goby_vcf@gt)[colnames(goby_vcf@gt) == old_name] <- control_well_renames[[old_name]]
    msg("  renamed '", old_name, "' -> '", control_well_renames[[old_name]], "'")
  }
}

ntc_sample_names <- grep("NTC", colnames(goby_vcf@gt), value = TRUE)
msg("  negative controls found: ", length(ntc_sample_names))

# Genotypes appearing in a well that contained only water are, by construction,
# not real. The question is whether they came from physical contamination
# between wells or from index hopping on the flow cell. The two leave different
# signatures in the allele depths: contamination transfers whole genotypes and
# so produces balanced heterozygotes, while index hopping misassigns individual
# reads and so produces very lopsided allele balances at trivial depth. This
# block measures which one it was.
if (length(ntc_sample_names) > 0) {

  ntc_allele_depths <- vcfR2tidy(goby_vcf, format_fields = c("GT", "AD"),
                                 verbose = FALSE)$gt %>%
    filter(Indiv %in% ntc_sample_names,
           !is.na(gt_GT), gt_GT != "./.", !is.na(gt_AD)) %>%
    separate(gt_AD, into = c("depth_ref", "depth_alt"),
             sep = ",", extra = "drop", fill = "right") %>%
    mutate(across(c(depth_ref, depth_alt), as.numeric),
           total_depth    = depth_ref + depth_alt,
           allele_balance = depth_alt / total_depth) %>%
    filter(!is.na(allele_balance))

  index_hopping_stats <- tibble(
    n_spurious_genotypes = nrow(ntc_allele_depths),
    max_depth            = max(ntc_allele_depths$total_depth, na.rm = TRUE),
    pct_balance_below_25 = round(100 * mean(ntc_allele_depths$allele_balance < 0.25), 1)
  )
  write_csv(index_hopping_stats,
            file.path(paths$dir_filtering, "ntc_index_hopping_stats.csv"))
  msg(sprintf("  %d spurious NTC genotypes | max depth %sx | %s%% below 0.25 allele balance",
              index_hopping_stats$n_spurious_genotypes,
              index_hopping_stats$max_depth,
              index_hopping_stats$pct_balance_below_25))

  diagnostic_plots$index_hopping <-
    ggplot(ntc_allele_depths, aes(allele_balance)) +
    geom_histogram(binwidth = 0.05, fill = "tomato", colour = "black", boundary = 0) +
    labs(title = "Allele balance in negative controls",
         subtitle = sprintf("max depth %sx; %s%% below 0.25 -- index hopping, not contamination",
                            index_hopping_stats$max_depth,
                            index_hopping_stats$pct_balance_below_25),
         x = "Alternate allele balance (alt reads / total depth)",
         y = "Spurious genotypes") +
    theme_bw()

  goby_vcf <- goby_vcf[, c("FORMAT", setdiff(colnames(goby_vcf@gt)[-1], ntc_sample_names))]
  record_step("NTC removal", "drop negative-control wells", goby_vcf,
              indiv_removed = length(ntc_sample_names))
}

# Read depth by collection year. This is the evidence behind ranking SNPs on
# 2022 completeness in the thinning step below.
raw_depth_matrix <- extract.gt(goby_vcf, element = "DP", as.numeric = TRUE)
diagnostic_plots$read_depth <- as.data.frame(raw_depth_matrix) %>%
  rownames_to_column("locus") %>%
  pivot_longer(-locus, names_to = "sample", values_to = "depth") %>%
  mutate(collection_year = substr(sample, 1, 4)) %>%
  ggplot(aes(collection_year, depth, fill = collection_year)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Raw read depth by collection year",
       subtitle = "2022 is the low-coverage year; the thinning step protects it",
       x = "Collection year", y = "Read depth (DP)") +
  theme_bw() + theme(legend.position = "none")
rm(raw_depth_matrix)


# =============================================================================
# 2. Minor allele frequency
# =============================================================================
msg("\n--- 2. Locus culling ---")

# Allele frequency is read from the INFO AF= field rather than recomputed. It
# is what Stacks called across the full sample set, and at this stage no
# individuals have been removed, so recomputing would give the same answer far
# more slowly.
loci_before <- count_loci(goby_vcf)
allele_freq_info <- as_tibble(goby_vcf@fix) %>%
  mutate(alt_freq = as.numeric(str_match(INFO, "AF=([0-9.]+)")[, 2]),
         minor_allele_freq = pmin(alt_freq, 1 - alt_freq))

keep_locus <- !is.na(allele_freq_info$minor_allele_freq) &
  allele_freq_info$minor_allele_freq >= filter_settings$min_maf
goby_vcf <- goby_vcf[keep_locus, ]
record_step("MAF filter", sprintf("MAF >= %.2f", filter_settings$min_maf),
            goby_vcf, loci_removed = loci_before - count_loci(goby_vcf))


# =============================================================================
# 3. Thin to one SNP per RAD locus
# =============================================================================
# SNPs on the same ~500 bp RAD locus are completely linked, so keeping more than
# one inflates every LD-based statistic and adds no independent information.
# Which SNP is kept is chosen on genotyping completeness in the low-coverage
# priority year first and overall second, so that the thinned panel is the one
# that best preserves the weakest year.

loci_before <- count_loci(goby_vcf)

genotypes_long <- vcfR2tidy(goby_vcf, format_fields = "GT", verbose = FALSE)$gt %>%
  mutate(is_called = !is.na(gt_GT) & gt_GT != "./.")

completeness_priority_year <- genotypes_long %>%
  filter(substr(Indiv, 1, 4) == filter_settings$thin_priority_year) %>%
  group_by(ChromKey) %>%
  summarise(n_called_priority = sum(is_called), .groups = "drop")

completeness_overall <- genotypes_long %>%
  group_by(ChromKey) %>%
  summarise(n_called_overall = sum(is_called), .groups = "drop")

# The Stacks ID field looks like "5:174:+" -- RAD locus, position in the read,
# strand. Splitting on the first colon gives the RAD locus we group by.
best_snp_per_rad_locus <- as_tibble(goby_vcf@fix) %>%
  mutate(row_index = row_number()) %>%
  separate(ID, into = "rad_locus_id", sep = ":", remove = FALSE, extra = "drop") %>%
  left_join(completeness_overall,       by = c("row_index" = "ChromKey")) %>%
  left_join(completeness_priority_year, by = c("row_index" = "ChromKey")) %>%
  mutate(n_called_priority = coalesce(n_called_priority, 0L),
         n_called_overall  = coalesce(n_called_overall,  0L)) %>%
  group_by(rad_locus_id) %>%
  arrange(desc(n_called_priority), desc(n_called_overall), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  arrange(row_index) %>%
  pull(row_index)

goby_vcf <- goby_vcf[best_snp_per_rad_locus, ]
record_step("Thinning",
            sprintf("1 SNP/RAD locus (%s priority)", filter_settings$thin_priority_year),
            goby_vcf, loci_removed = loci_before - count_loci(goby_vcf))

need(!anyDuplicated(str_extract(goby_vcf@fix[, "ID"], "^[^:]+")),
     "Thinning failed: duplicate RAD locus IDs remain in the panel.")
rm(genotypes_long, completeness_priority_year, completeness_overall,
   best_snp_per_rad_locus)


# =============================================================================
# 4. Hardy-Weinberg screen, via per-locus Fis
# =============================================================================
msg("\n--- 3. Locus quality control ---")

# Fis is used as the filter rather than the HWE p-value because with ~370
# samples the chi-square test has enough power to reject trivially small
# departures, and rejecting on significance would remove hundreds of loci whose
# actual deviation is biologically irrelevant. Fis measures the size of the
# departure, which is what we care about. The p-values are computed and written
# out anyway so the choice can be checked.

loci_before <- count_loci(goby_vcf)

goby_genind <- vcfR2genind(goby_vcf)
pop(goby_genind) <- rep("Overall", nInd(goby_genind))   # one pool: HWE is global here

per_locus_fis <- basic.stats(genind2hierfstat(goby_genind), diploid = TRUE)$perloc %>%
  as_tibble() %>%
  mutate(locus = locNames(goby_genind), .before = 1)

per_locus_hwe_p <- as_tibble(pegas::hw.test(goby_genind, B = 0), rownames = "locus") %>%
  select(locus, p_value = `Pr(chi^2 >)`)

hwe_table <- left_join(per_locus_fis, per_locus_hwe_p, by = "locus")
write_csv(hwe_table, file.path(paths$dir_filtering, "hwe_per_locus.csv"))

bonferroni_alpha <- 0.05 / nrow(hwe_table)
msg(sprintf("  |Fis| > %.2f: %d loci | p < Bonferroni (%.2e): %d loci",
            filter_settings$max_abs_fis,
            sum(abs(hwe_table$Fis) > filter_settings$max_abs_fis, na.rm = TRUE),
            bonferroni_alpha,
            sum(hwe_table$p_value <= bonferroni_alpha, na.rm = TRUE)))

diagnostic_plots$hwe <- ggplot(hwe_table, aes(Fis, -log10(p_value))) +
  geom_point(alpha = 0.4, colour = "midnightblue") +
  geom_vline(xintercept = c(-filter_settings$max_abs_fis, filter_settings$max_abs_fis),
             colour = "red", linetype = "dashed") +
  labs(title = "HWE deviation: effect size against significance",
       subtitle = "filtered on Fis (dashed), not on p, because n is large enough to reject trivial departures",
       x = "Inbreeding coefficient (Fis)", y = expression(-log[10](p))) +
  theme_bw()

loci_failing_hwe <- hwe_table$locus[
  abs(hwe_table$Fis) > filter_settings$max_abs_fis & !is.na(hwe_table$Fis)
]
goby_vcf <- goby_vcf[!(locNames(goby_genind) %in% loci_failing_hwe), ]
record_step("HWE filter", sprintf("|Fis| <= %.2f", filter_settings$max_abs_fis),
            goby_vcf, loci_removed = loci_before - count_loci(goby_vcf))
rm(goby_genind, per_locus_fis, per_locus_hwe_p, hwe_table)


# =============================================================================
# 5. Paralog screen (McKinney et al. 2017)
# =============================================================================
# A locus assembled from two duplicated genomic regions behaves like a permanent
# heterozygote: every individual carries one copy of each "allele", so observed
# heterozygosity is high, and because the two paralogs rarely amplify equally the
# read ratio at those apparent heterozygotes is skewed away from 50:50. Screening
# on both together separates real high-heterozygosity loci from collapsed
# duplicates.

loci_before <- count_loci(goby_vcf)

paralog_input <- vcfR2tidy(goby_vcf, format_fields = c("GT", "AD"), verbose = FALSE)$gt %>%
  rename(genotype = gt_GT, locus_index = ChromKey) %>%
  separate(gt_AD, into = c("depth_ref", "depth_alt"),
           sep = ",", extra = "drop", fill = "right")

# Genotype classes are counted explicitly rather than reshaped with pivot_wider.
# A locus that happens to have no homozygous-alternate individuals would leave
# pivot_wider without that column, and the arithmetic that follows would fail.
observed_heterozygosity <- paralog_input %>%
  filter(!is.na(genotype), genotype != "./.") %>%
  group_by(locus_index) %>%
  summarise(n_hom_ref = sum(genotype == "0/0"),
            n_het     = sum(genotype %in% c("0/1", "1/0")),
            n_hom_alt = sum(genotype == "1/1"),
            .groups = "drop") %>%
  mutate(n_total = n_hom_ref + n_het + n_hom_alt,
         hobs    = ifelse(n_total > 0, n_het / n_total, NA_real_))

# Read-ratio deviation, pooled across all heterozygotes at the locus. Under a
# fair 50:50 expectation the reference read count is binomial, so this is a
# standard normal deviate: (observed - expected) / sd.
read_ratio_deviation <- paralog_input %>%
  filter(genotype %in% c("0/1", "1/0"), !is.na(depth_ref), !is.na(depth_alt)) %>%
  mutate(across(c(depth_ref, depth_alt), as.numeric)) %>%
  group_by(locus_index) %>%
  summarise(total_reads = sum(depth_ref) + sum(depth_alt),
            z_score = (sum(depth_ref) - total_reads / 2) / sqrt(total_reads * 0.25),
            .groups = "drop") %>%
  select(locus_index, z_score)

paralog_table <- left_join(observed_heterozygosity, read_ratio_deviation,
                           by = "locus_index") %>%
  mutate(is_singleton = hobs <= filter_settings$max_paralog_hobs &
           abs(z_score) <= filter_settings$max_paralog_zscore & !is.na(z_score))
write_csv(paralog_table, file.path(paths$dir_filtering, "paralog_per_locus.csv"))

diagnostic_plots$paralog <- ggplot(paralog_table, aes(hobs, z_score, colour = is_singleton)) +
  geom_point(alpha = 0.6) +
  labs(title = "Paralog screen (McKinney et al. 2017)",
       subtitle = "duplicated loci sit high on both axes: excess heterozygotes with skewed read ratios",
       x = "Observed heterozygosity", y = "Read-ratio deviation (z)",
       colour = "Retained as\nsingle-copy") +
  theme_bw()

paralog_row_indices <- paralog_table$locus_index[!paralog_table$is_singleton]
goby_vcf <- goby_vcf[!(seq_len(nrow(goby_vcf)) %in% paralog_row_indices), ]
record_step("Paralog filter",
            sprintf("Hobs > %.2f or |z| > %g",
                    filter_settings$max_paralog_hobs, filter_settings$max_paralog_zscore),
            goby_vcf, loci_removed = loci_before - count_loci(goby_vcf))
rm(paralog_input, observed_heterozygosity, read_ratio_deviation, paralog_table)


# =============================================================================
# 6. Mask low-quality genotypes
# =============================================================================
msg("\n--- 4. Genotype and individual filtering ---")

depth_matrix   <- extract.gt(goby_vcf, element = "DP", as.numeric = TRUE)
quality_matrix <- extract.gt(goby_vcf, element = "GQ", as.numeric = TRUE)

low_quality_mask <- (depth_matrix   < filter_settings$min_depth |
                     quality_matrix < filter_settings$min_gq)
low_quality_mask[is.na(low_quality_mask)] <- FALSE

# Threshold-justification curves: how much data does each candidate cutoff cost?
# Denominator is called genotypes, so the y axis reads as "percent of the data
# you would be throwing away".
genotype_block <- goby_vcf@gt[, -1, drop = FALSE]
n_called_genotypes <- sum(!is.na(genotype_block) & genotype_block != "./.")

threshold_curve <- function(value_matrix, thresholds, axis_label, chosen_value) {
  map_dfr(thresholds, ~tibble(
    threshold  = .x,
    pct_masked = 100 * sum(value_matrix < .x, na.rm = TRUE) / n_called_genotypes
  )) %>%
    ggplot(aes(threshold, pct_masked)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    geom_vline(xintercept = chosen_value, colour = "red", linetype = "dashed") +
    labs(title = paste("Cost of the minimum", axis_label, "filter"),
         subtitle = paste("dashed line is the chosen cutoff:", chosen_value),
         x = paste("Minimum", axis_label),
         y = "Cumulative % of called genotypes masked") +
    theme_bw()
}
diagnostic_plots$depth_cost <- threshold_curve(
  depth_matrix, seq(0, 50, 5), "read depth (DP)", filter_settings$min_depth)
diagnostic_plots$quality_cost <- threshold_curve(
  quality_matrix, seq(0, 40, 5), "genotype quality (GQ)", filter_settings$min_gq)

# The mask covers the genotype columns only. vcf@gt carries an extra leading
# FORMAT column, so the two objects differ by one column and the mask must be
# applied to the genotype block, never to vcf@gt directly. The assertions make
# that requirement explicit rather than trusting it.
genotype_block  <- goby_vcf@gt[, -1, drop = FALSE]
format_column_before <- goby_vcf@gt[, 1]

need(identical(dim(genotype_block), dim(low_quality_mask)),
     "Genotype block is ", paste(dim(genotype_block), collapse = " x "),
     " but the mask is ", paste(dim(low_quality_mask), collapse = " x "), ".")

genotype_block[low_quality_mask] <- "./."
goby_vcf@gt[, -1] <- genotype_block

need(identical(goby_vcf@gt[, 1], format_column_before),
     "The FORMAT column was altered by genotype masking.")

msg(sprintf("  masked %s genotypes with DP < %d or GQ < %d (%.1f%% of called)",
            format(sum(low_quality_mask), big.mark = ","),
            filter_settings$min_depth, filter_settings$min_gq,
            100 * sum(low_quality_mask) / n_called_genotypes))
record_step("Genotype masking",
            sprintf("DP < %d or GQ < %d -> ./.",
                    filter_settings$min_depth, filter_settings$min_gq),
            goby_vcf)
rm(depth_matrix, quality_matrix, low_quality_mask, genotype_block)


# =============================================================================
# 7. Individual call rate
# =============================================================================
individuals_before <- count_individuals(goby_vcf)

genotype_block <- goby_vcf@gt[, -1, drop = FALSE]
missing_per_individual <- colSums(is.na(genotype_block) | genotype_block == "./.") /
  nrow(genotype_block)

individual_missingness <- tibble(
  sample          = names(missing_per_individual),
  prop_missing    = as.numeric(missing_per_individual),
  collection_year = substr(names(missing_per_individual), 1, 4)
)
write_csv(individual_missingness,
          file.path(paths$dir_filtering, "missingness_per_individual.csv"))

diagnostic_plots$individual_missingness <-
  ggplot(individual_missingness, aes(collection_year, prop_missing, fill = collection_year)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_hline(yintercept = filter_settings$max_missing_per_individual,
             colour = "red", linetype = "dashed") +
  labs(title = "Individual missing data after genotype masking",
       subtitle = sprintf("cutoff %.0f%%", 100 * filter_settings$max_missing_per_individual),
       x = "Collection year", y = "Proportion missing") +
  theme_bw() + theme(legend.position = "none")

samples_to_drop <- individual_missingness$sample[
  individual_missingness$prop_missing > filter_settings$max_missing_per_individual
]
msg("  dropping ", length(samples_to_drop), " individuals over ",
    100 * filter_settings$max_missing_per_individual, "% missing")

goby_vcf <- goby_vcf[, c("FORMAT", setdiff(colnames(goby_vcf@gt)[-1], samples_to_drop))]
record_step("Individual call rate",
            sprintf("> %.0f%% missing", 100 * filter_settings$max_missing_per_individual),
            goby_vcf, indiv_removed = individuals_before - count_individuals(goby_vcf))
rm(genotype_block, missing_per_individual)


# =============================================================================
# 8. Locus call rate
# =============================================================================
# Both NA and the literal "./." count as missing here. The masking step above
# writes the string, not NA, so a check that tested only is.na() would see
# almost none of the data it just removed.

loci_before <- count_loci(goby_vcf)
genotype_block <- goby_vcf@gt[, -1, drop = FALSE]
missing_per_locus <- rowSums(is.na(genotype_block) | genotype_block == "./.") /
  ncol(genotype_block)

goby_vcf <- goby_vcf[missing_per_locus <= filter_settings$max_missing_per_locus, ]
record_step("Locus call rate",
            sprintf("> %.0f%% missing", 100 * filter_settings$max_missing_per_locus),
            goby_vcf, loci_removed = loci_before - count_loci(goby_vcf))
rm(genotype_block, missing_per_locus)


# =============================================================================
# 9. Monomorphic cleanup
# =============================================================================
# The MAF filter in step 2 was applied to the full sample set. After ~50
# individuals have been culled some loci retain only one allele, and a
# monomorphic locus contributes nothing but does inflate the locus count and
# breaks per-locus Fis with a divide-by-zero.

loci_before <- count_loci(goby_vcf)
final_maf <- vcfR::maf(goby_vcf)[, "Frequency"]
goby_vcf <- goby_vcf[!is.na(final_maf) & final_maf > 0, ]
record_step("Monomorphic cleanup", "MAF == 0 or NA after individual cull",
            goby_vcf, loci_removed = loci_before - count_loci(goby_vcf))


# =============================================================================
# 10. Verify, export, report
# =============================================================================
msg("\n--- 5. Finalisation ---")

final_genotypes <- goby_vcf@gt[, -1, drop = FALSE]
final_missing_by_individual <- colSums(is.na(final_genotypes) | final_genotypes == "./.") /
  nrow(final_genotypes)
final_missing_by_locus <- rowSums(is.na(final_genotypes) | final_genotypes == "./.") /
  ncol(final_genotypes)

# These four checks restate the filter's contract. If any fails, something
# upstream did not do what its ledger row claims it did.
need(all(substr(goby_vcf@gt[, 1], 1, 2) == "GT"),
     "FORMAT column is malformed in the final VCF.")
need(max(final_missing_by_individual) <= filter_settings$max_missing_per_individual + 1e-9,
     "An individual exceeds the missingness cutoff in the final VCF.")
need(max(final_missing_by_locus) <= filter_settings$max_missing_per_locus + 1e-9,
     "A locus exceeds the missingness cutoff in the final VCF.")
need(min(vcfR::maf(goby_vcf)[, "Frequency"], na.rm = TRUE) > 0,
     "A monomorphic locus survived the cleanup.")

msg(sprintf("  max missing per individual %.3f | per locus %.3f",
            max(final_missing_by_individual), max(final_missing_by_locus)))
msg(sprintf("\n  FINAL PANEL: %d loci x %d individuals",
            count_loci(goby_vcf), count_individuals(goby_vcf)))

# Missingness broken out by year, for the supplement. Locus-level missingness is
# computed within each year separately: a locus that failed only in 2022 looks
# fine when averaged over all four years, and that is exactly the failure mode
# worth reporting given the coverage difference.
collection_years <- substr(colnames(final_genotypes), 1, 4)
bind_rows(
  tibble(sample = names(final_missing_by_individual),
         prop_missing = as.numeric(final_missing_by_individual),
         year = substr(names(final_missing_by_individual), 1, 4)) %>%
    group_by(year) %>%
    summarise(type = "Individual", mean = mean(prop_missing),
              median = median(prop_missing), sd = sd(prop_missing),
              min = min(prop_missing), max = max(prop_missing), .groups = "drop"),
  map_dfr(unique(collection_years), function(this_year) {
    year_block <- final_genotypes[, collection_years == this_year, drop = FALSE]
    prop_missing <- rowSums(is.na(year_block) | year_block == "./.") / ncol(year_block)
    tibble(year = this_year, type = "Locus", mean = mean(prop_missing),
           median = median(prop_missing), sd = sd(prop_missing),
           min = min(prop_missing), max = max(prop_missing))
  })
) %>%
  write_csv(file.path(paths$dir_filtering, "missingness_summary_by_year.csv"))

write.vcf(goby_vcf, file = paths$filtered_vcf)
msg("  wrote ", paths$filtered_vcf)

write_csv(filtering_ledger, paths$filter_log)
msg("  wrote ", paths$filter_log)
print(as.data.frame(filtering_ledger))

pdf(paths$filter_report, width = 14, height = 11)
if (!is.null(diagnostic_plots$index_hopping)) print(diagnostic_plots$index_hopping)
print(diagnostic_plots$read_depth)
grid.arrange(diagnostic_plots$hwe, diagnostic_plots$paralog, ncol = 2)
grid.arrange(diagnostic_plots$depth_cost, diagnostic_plots$quality_cost,
             diagnostic_plots$individual_missingness, ncol = 2)
grid.arrange(tableGrob(as.data.frame(filtering_ledger), rows = NULL))
dev.off()
msg("  wrote ", paths$filter_report)

msg("\n=== 01_filter_vcf.R complete ===")
log_end()
