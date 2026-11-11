# =============================================================================
# 04_temporal_ne.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Temporal Ne from pairwise Fst between birth cohorts, plus the
# generation-length sensitivity table
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
#
# Description:
#   Reads the pairwise Fst tables that 03_popgen.R writes. It does not recompute
#   any genetics, so it runs in seconds and needs only dplyr and readr.
#
#   Method (Waples 1989; same as the September cohort script): for two temporal
#   samples F ~= 2 * Fst and Ne = t / (2F), so
#
#       Ne = t / (4 * Fst),     t = years_apart / generation_length
#
#   Weir & Cockerham Fst is already corrected for sample size, so no separate
#   1/(2S) sampling term is subtracted.
#
#   Only comparisons at least temporal_ne_settings$min_years_apart years apart
#   are kept (Waples 1990; Waples & Yokota 2007). Fst <= 0 means the two cohorts
#   cannot be told apart, so Ne is reported as Inf (undefined), never negative.
#
#   Two sources of comparisons:
#     cohort          every cohort pair (the September framing)
#     cohort_by_site  same-site pairs only, so site is held fixed (the
#                     clean_contrasts in 00_config.R)
#
#   Intervals: 03 writes Fst_L95/Fst_U95 only if run_fst_bootstrap = TRUE, and
#   p_perm only if run_fst_permutation = TRUE. This script carries whichever is
#   present and says in the log which it used. Locus-bootstrap intervals are a
#   lower bound on uncertainty (loci are linked; see 03).
#
# Inputs:   results/03_popgen/{cohort,cohort_by_site}/pairwise_fst.csv
#           results/03_popgen/{cohort,cohort_by_site}/diversity.csv   (for N)
# Outputs:  results/04_temporal_ne/temporal_ne.csv
#           results/04_temporal_ne/generation_length_sensitivity.csv
#
# Usage:    Rscript scripts/04_temporal_ne.R      # from the project root
# =============================================================================

source("scripts/00_config.R")

need_pkgs(c("dplyr", "readr", "tidyr"))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})

log_start("04_temporal_ne")

output_dir <- project_path("results", "04_temporal_ne")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

settings <- temporal_ne_settings
msg("  min years apart: ", settings$min_years_apart,
    " | generation lengths: ", paste(settings$generation_lengths, collapse = ", "),
    " | primary: ", settings$primary_generation_length)


# -----------------------------------------------------------------------------
# Read one grouping's Fst table and attach N
# -----------------------------------------------------------------------------
read_grouping <- function(grouping_name) {
  grouping_dir <- file.path(paths$dir_popgen, grouping_name)
  fst_file <- file.path(grouping_dir, "pairwise_fst.csv")
  div_file <- file.path(grouping_dir, "diversity.csv")
  need(file.exists(fst_file), "Missing ", fst_file, " -- run scripts/03_popgen.R first.")
  need(file.exists(div_file), "Missing ", div_file, " -- run scripts/03_popgen.R first.")

  group_n <- read_csv(div_file, col_types = cols(group = col_character(),
                                                 .default = col_guess())) %>%
    select(group, N)

  read_csv(fst_file, col_types = cols(group1 = col_character(),
                                      group2 = col_character(),
                                      .default = col_guess())) %>%
    left_join(group_n %>% rename(group1 = group, N1 = N), by = "group1") %>%
    left_join(group_n %>% rename(group2 = group, N2 = N), by = "group2") %>%
    mutate(source = grouping_name)
}

# Group labels are "2021" (cohort) or "2021_Albany" (cohort_by_site).
split_label <- function(fst_pairs) {
  fst_pairs %>%
    mutate(cohort1 = as.integer(sub("_.*$", "", group1)),
           cohort2 = as.integer(sub("_.*$", "", group2)),
           site1   = ifelse(grepl("_", group1), sub("^[0-9]+_", "", group1), NA_character_),
           site2   = ifelse(grepl("_", group2), sub("^[0-9]+_", "", group2), NA_character_))
}

comparisons <- bind_rows(
  read_grouping("cohort") %>% split_label(),
  read_grouping("cohort_by_site") %>% split_label() %>%
    filter(site1 == site2),                    # temporal only: site held fixed
  ### ADDED (patch 18): age-ambiguity sensitivity, if 03 produced it
  if (dir.exists(file.path(paths$dir_popgen, "cohort_strict")))
    read_grouping("cohort_strict") %>% split_label()
) %>%
  mutate(site        = coalesce(site1, "all sites"),
         years_apart = abs(cohort2 - cohort1),
         # order each pair earlier -> later so the table reads in time order
         cohort_a    = pmin(cohort1, cohort2),
         cohort_b    = pmax(cohort1, cohort2),
         N_a         = ifelse(cohort1 < cohort2, N1, N2),
         N_b         = ifelse(cohort1 < cohort2, N2, N1),
         min_N       = pmin(N1, N2)) %>%
  select(source, site, cohort_a, cohort_b, N_a, N_b, min_N, years_apart,
         Fst, any_of(c("Fst_L95", "Fst_U95", "p_perm", "Fst_null95")))

n_dropped <- sum(comparisons$years_apart < settings$min_years_apart)
comparisons <- comparisons %>% filter(years_apart >= settings$min_years_apart)
msg("  ", nrow(comparisons), " comparisons kept; ", n_dropped,
    " dropped for spanning < ", settings$min_years_apart, " years")

has_ci   <- all(c("Fst_L95", "Fst_U95") %in% names(comparisons))
has_perm <- "p_perm" %in% names(comparisons)
msg("  Fst intervals: ", ifelse(has_ci, "locus bootstrap (lower bound on uncertainty)",
                                "none -- set run_fst_bootstrap = TRUE in 03 to get them"))
msg("  permutation p: ", ifelse(has_perm, "present", "none -- set run_fst_permutation = TRUE in 03"))


# -----------------------------------------------------------------------------
# Ne at a given generation length
# -----------------------------------------------------------------------------
# Low Fst means high Ne, so the Fst upper limit gives the Ne lower limit and
# vice versa. Any Fst limit at or below zero puts that Ne limit at Inf.
fst_to_ne <- function(t_gen, fst) ifelse(fst > 0, t_gen / (4 * fst), Inf)
alpha_temporal <- if (is.null(settings$alpha)) 0.05 else settings$alpha

ne_at <- function(comparison_table, generation_length) {
  out <- comparison_table %>%
    mutate(generation_length = generation_length,
           t_gen  = years_apart / generation_length,
           Ne     = fst_to_ne(t_gen, Fst),
           # Waples & Yokota (2007) want the samples at least one generation
           # apart; flag rather than drop so the sensitivity table stays square
           t_ge_g = t_gen >= 1)
  if (has_ci) {
    out <- out %>% mutate(Ne_L95 = fst_to_ne(t_gen, Fst_U95),
                          Ne_U95 = fst_to_ne(t_gen, Fst_L95))
  }
  ### ADDED (patch 16): detectable-Ne floor. Fst_null95 is the largest Fst that
  ### sampling alone produces for this pair (95th percentile of the permutation
  ### null). Drift from an Ne below t / (4 * Fst_null95) would be expected to
  ### exceed it, so that value is an approximate lower bound on Ne when the
  ### observed Fst is not significant. Approximate: at the bound itself the
  ### power to detect is only about 50%.
  if ("Fst_null95" %in% names(out)) {
    ### CHANGED (patch 20): a floor only means something where drift was not
    ### detected; where Fst is significant the estimate replaces it.
    out <- out %>% mutate(Ne_floor = fst_to_ne(t_gen, Fst_null95))
    if ("p_perm" %in% names(out))
      out <- out %>% mutate(Ne_floor = ifelse(p_perm < alpha_temporal, NA_real_, Ne_floor))
  }
  out
}

temporal_ne <- ne_at(comparisons, settings$primary_generation_length) %>%
  ### CHANGED (patch 16): Ne is only estimated when Fst is significantly > 0.
  ### Otherwise the drift signal is indistinguishable from sampling noise and
  ### only the floor (Ne_floor) is reported.
  mutate(significant = if (has_perm) p_perm < alpha_temporal else NA,
         Ne = if (has_perm) ifelse(significant %in% TRUE, Ne, NA_real_) else Ne,
         flag = case_when(has_perm & !(significant %in% TRUE) ~
                            paste0("Fst not > 0 (p = ", if (has_perm) round(p_perm, 3) else NA,
                                   "): Ne not estimable; floor only"),
                          Fst <= 0 ~ "Fst <= 0: Ne undefined",
                          min_N < cohort_settings$min_n_estimator ~ "low precision (min N < 30)",
                          TRUE ~ "")) %>%
  arrange(source, site, cohort_a, cohort_b)

write_csv(temporal_ne, file.path(output_dir, "temporal_ne.csv"))
msg("\n  temporal Ne at generation length = ", settings$primary_generation_length, " yr:")
msg_table(as.data.frame(temporal_ne %>%
  mutate(Fst = signif(Fst, 3),
         across(any_of(c("Ne", "Ne_L95", "Ne_U95", "Ne_floor")), ~round(.x, 0))) %>%
  select(-any_of(c("Fst_L95", "Fst_U95")))), row.names = FALSE)

informative <- temporal_ne %>% filter(is.finite(Ne), flag == "")
if (nrow(informative)) {
  msg(sprintf("\n  informative comparisons (Fst significant, finite Ne, min N >= %d): Ne %.0f - %.0f",
              cohort_settings$min_n_estimator, min(informative$Ne), max(informative$Ne)))
} else {
  msg("\n  no comparison has Fst significantly > 0: temporal Ne is not estimable.")
}
if ("Ne_floor" %in% names(temporal_ne)) {
  floors <- temporal_ne %>% filter(min_N >= cohort_settings$min_n_estimator, is.finite(Ne_floor))
  if (!nrow(floors)) msg("  detectable-Ne floor: no well-powered pair lacks significant drift")
  if (nrow(floors))
    msg(sprintf("  detectable-Ne floor (well-powered pairs): Ne below ~%.0f would have produced detectable drift (range %.0f - %.0f)",
                max(floors$Ne_floor), min(floors$Ne_floor), max(floors$Ne_floor)))
}


# -----------------------------------------------------------------------------
# Generation-length sensitivity: one row per comparison, one column per T
# -----------------------------------------------------------------------------
sensitivity_long <- bind_rows(lapply(settings$generation_lengths,
                                     function(g) ne_at(comparisons, g)))

### CHANGED (patch 16): where Fst is not significant, the sensitivity table
### carries the detectable-Ne floor instead of a point estimate (basis column).
use_floor <- has_perm && "Ne_floor" %in% names(sensitivity_long)
sensitivity_wide <- sensitivity_long %>%
  mutate(basis  = if (use_floor) ifelse(p_perm < alpha_temporal, "estimate", "floor") else "estimate",
         Ne     = if (use_floor) ifelse(basis == "estimate", Ne, Ne_floor) else Ne,
         column = paste0("Ne_T", generation_length),
         Ne     = ifelse(t_ge_g, round(Ne, 0), NA)) %>%     # NA where t < g
  select(source, site, cohort_a, cohort_b, years_apart, Fst, basis, column, Ne) %>%
  pivot_wider(names_from = column, values_from = Ne) %>%
  arrange(source, site, cohort_a, cohort_b)

write_csv(sensitivity_wide, file.path(output_dir, "generation_length_sensitivity.csv"))
msg("\n  generation-length sensitivity (basis = estimate or detectable floor; NA = samples less than one generation apart):")
msg_table(as.data.frame(sensitivity_wide %>% mutate(Fst = signif(Fst, 3))), row.names = FALSE)

msg("\n  wrote ", output_dir)
msg("\n=== 04_temporal_ne.R complete ===")
log_end()
