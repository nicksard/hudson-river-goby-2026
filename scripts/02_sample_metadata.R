# =============================================================================
# 02_sample_metadata.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Biological data workbook -> one authoritative sample table
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-09-08
# Updated:  2026-09-10
#
# Description:
#   Builds sample_metadata.csv: one row per fish, carrying birth cohort,
#   location under each grouping scheme, QC status, and a single analysis_keep
#   flag. Every downstream script reads this file and nothing else for its
#   sample set, so that cohort membership, region labels and QC exclusions are
#   defined once and cannot drift between analyses.
#
#   Three things happen here that are worth understanding before trusting the
#   output:
#
#   1. Cohort is derived, not recorded. birth_year = collection_year - age,
#      where collection year is the first four digits of the plate ID and age
#      is the otolith reading. If age is missing the fish gets no cohort and
#      drops out of every cohort-based analysis.
#
#   2. The workbook has two different identifier columns. ID is the plate/well
#      identifier and is what collection year is parsed from. Fish_ID is what
#      the VCF sample columns are named. They are not interchangeable.
#
#   3. Region is free text with seventeen distinct values, including four
#      spellings of the same Lock C1 site. It is mapped through
#      region_patterns in 00_config.R, and any unmapped level is a hard error
#      rather than a silent NA -- see the note at that step.
#
# Inputs:
#   paths$bio_data       GobyRad_tabularPlateWell with Bio Data.xlsx
#   paths$filtered_vcf   output of 01_filter_vcf.R
#
# Outputs:
#   results/02_metadata/sample_metadata.csv      <- the sample set
#   results/02_metadata/region_level_mapping.csv
#   results/02_metadata/length_qc_flags.csv      (only if flags are raised)
#
# Usage:
#   Rscript scripts/02_sample_metadata.R
# =============================================================================

source("scripts/00_config.R")
need_pkgs(c("vcfR", "readxl", "dplyr", "tidyr", "readr", "stringr", "tibble"))

suppressPackageStartupMessages({
  library(vcfR); library(readxl); library(dplyr); library(tidyr)
  library(readr); library(stringr); library(tibble)
})

log_start("02_sample_metadata")

need(file.exists(paths$bio_data),
     "Biological data workbook not found: ", paths$bio_data, "\n",
     "  Set the BIO_DATA environment variable, or put the file in data/raw/.\n",
     "  Cohort and location exist nowhere else, so nothing downstream can run without it.")
need(file.exists(paths$filtered_vcf),
     "Filtered VCF not found: ", paths$filtered_vcf,
     "\n  Run scripts/01_filter_vcf.R first.")


# =============================================================================
# 1. Which samples survived filtering
# =============================================================================
genotyped_samples <- colnames(read.vcfR(paths$filtered_vcf, verbose = FALSE)@gt)[-1]
msg("  filtered VCF carries ", length(genotyped_samples), " samples")


# =============================================================================
# 2. Read the workbook and normalise its column names
# =============================================================================
bio_data <- read_xlsx(paths$bio_data, sheet = 1)

# Column headers in the workbook carry spaces, parentheses and stray brackets
# that make them awkward to reference. Normalising once here means the rest of
# the script can use plain names.
names(bio_data) <- names(bio_data) %>%
  trimws() %>%
  str_replace_all(" ", "_") %>%
  str_remove_all("[()#`']") %>%
  str_remove_all("_\\[|\\]_") %>%
  trimws()
msg("  workbook: ", nrow(bio_data), " rows, ", ncol(bio_data), " columns")

for (required_column in c("ID", "Fish_ID", "Age", "Region")) {
  need(required_column %in% names(bio_data),
       "Workbook is missing the required column '", required_column, "'.\n",
       "  Columns present: ", paste(names(bio_data), collapse = ", "))
}


# =============================================================================
# 3. Remove control wells
# =============================================================================
# Two kinds of control row: wells labelled NTC, and the 2024999-style IDs used
# for the well that was relabelled during filtering. The control-well pattern is
# anchored (four digits then exactly 999) rather than matched loosely anywhere
# in the string, so that a real fish whose sequence number happened to contain
# 999 is not swept up with them. Anything caught by the loose pattern but not
# the anchored one is reported, so the difference is visible rather than assumed.

rows_before_control_removal <- nrow(bio_data)
is_ntc_well     <- str_detect(bio_data$ID, regex("NTC", ignore_case = TRUE))
is_control_well <- str_detect(bio_data$ID, "^\\d{4}999$")

if (any(is_ntc_well | is_control_well)) {
  msg("  dropping ", sum(is_ntc_well | is_control_well), " control rows: ",
      paste(head(bio_data$ID[is_ntc_well | is_control_well], 20), collapse = ", "))
}

contains_999_but_not_control <- str_detect(bio_data$ID, "999") &
  !is_control_well & !is_ntc_well
if (any(contains_999_but_not_control)) {
  msg("  NOTE: ", sum(contains_999_but_not_control),
      " row(s) contain '999' but are not control wells, and are KEPT: ",
      paste(bio_data$ID[contains_999_but_not_control], collapse = ", "))
}

bio_data <- bio_data[!(is_ntc_well | is_control_well), ]
msg("  bio rows: ", rows_before_control_removal, " -> ", nrow(bio_data))


# =============================================================================
# 4. Derive birth cohort
# =============================================================================
bio_data <- bio_data %>%
  mutate(
    plate_id        = as.character(ID),
    sample          = as.character(Fish_ID),
    collection_year = suppressWarnings(as.numeric(substr(plate_id, 1, 4))),
    age             = suppressWarnings(as.numeric(Age)),
    cohort          = collection_year - age
  )

msg(sprintf("  age data: %d of %d rows usable (%d missing)",
            sum(!is.na(bio_data$age)), nrow(bio_data), sum(is.na(bio_data$age))))
if (any(is.na(bio_data$age))) {
  msg("  !! fish with no age get no cohort and drop out of cohort analyses")
}

duplicated_fish_ids <- bio_data$sample[duplicated(bio_data$sample) & !is.na(bio_data$sample)]
if (length(duplicated_fish_ids)) {
  msg("  !! Fish_ID appears more than once in the workbook: ",
      paste(unique(duplicated_fish_ids), collapse = ", "))
}


# =============================================================================
# 5. Map Region to analysis units
# =============================================================================
bio_data$region_raw <- trimws(as.character(bio_data$Region))

# First matching pattern in region_patterns wins, which is why the Lock C1
# pattern is listed before anything that might also appear in those long
# free-text strings.
map_region_to_unit <- function(region_strings) {
  mapped <- rep(NA_character_, length(region_strings))
  for (pattern in names(region_patterns)) {
    matches_this_pattern <- is.na(mapped) & !is.na(region_strings) &
      str_detect(region_strings, regex(pattern, ignore_case = TRUE))
    mapped[matches_this_pattern] <- region_patterns[[pattern]]
  }
  mapped
}
bio_data$region <- map_region_to_unit(bio_data$region_raw)

msg("\n  Region levels found -> analysis unit:")
region_mapping_summary <- bio_data %>%
  count(region_raw, region, name = "n") %>%
  arrange(desc(n))
for (row_i in seq_len(nrow(region_mapping_summary))) {
  msg(sprintf("    %-60s -> %-14s n = %3d%s",
              ifelse(is.na(region_mapping_summary$region_raw[row_i]), "<NA>",
                     region_mapping_summary$region_raw[row_i]),
              ifelse(is.na(region_mapping_summary$region[row_i]), "UNMAPPED",
                     region_mapping_summary$region[row_i]),
              region_mapping_summary$n[row_i],
              ifelse(is.na(region_mapping_summary$region[row_i]) &&
                       !is.na(region_mapping_summary$region_raw[row_i]),
                     "   <-- ADD TO region_patterns", "")))
}
write_csv(region_mapping_summary,
          file.path(paths$dir_metadata, "region_level_mapping.csv"))

# An unmapped level is a hard stop rather than a silent NA. An NA region does
# not remove the fish from the dataset, it removes it from every location
# analysis while leaving it in the cohort analyses -- so the two would be run on
# different sample sets without anything in the output saying so.
unmapped_levels <- region_mapping_summary %>%
  filter(is.na(region), !is.na(region_raw))
need(nrow(unmapped_levels) == 0,
     "Unmapped Region level(s): ", paste(unmapped_levels$region_raw, collapse = ", "), "\n",
     "  Add a pattern for each to `region_patterns` in scripts/00_config.R and rerun.\n",
     "  Full listing: results/02_metadata/region_level_mapping.csv")

# Order sites north to south so that figures read like the river rather than
# like a dictionary.
bio_data$region <- factor(bio_data$region,
                          levels = intersect(region_downstream_order,
                                             unique(bio_data$region)))

for (scheme_name in names(location_schemes)) {
  bio_data[[paste0("loc_", scheme_name)]] <-
    location_schemes[[scheme_name]](bio_data)
}


# =============================================================================
# 6. Length plausibility check
# =============================================================================
# Diagnostic only. Length is never used to assign a cohort here -- otolith age
# is -- but an impossible length is a flag that the whole metadata row may be
# unreliable, which matters when that row is also supplying a region.

length_column <- intersect(
  c("Length_mm", "Total_Length_mm", "TL_mm", "Length", "TL", "Length_MM"),
  names(bio_data)
)
bio_data$length_mm <- if (length(length_column)) {
  suppressWarnings(as.numeric(bio_data[[length_column[1]]]))
} else {
  NA_real_
}
msg(if (length(length_column)) paste0("\n  length column: ", length_column[1])
    else "\n  no length column found")

# Round goby hatch at roughly 4-5 mm and top out near 250 mm, and neither
# extreme is catchable with the gears used here.
implausible_lengths <- bio_data %>%
  filter(!is.na(length_mm), length_mm < 12 | length_mm > 250)
if (nrow(implausible_lengths)) {
  msg("  !! ", nrow(implausible_lengths), " implausible length value(s):")
  for (row_i in seq_len(nrow(implausible_lengths))) {
    msg(sprintf("     %s  region %s  length %s mm",
                implausible_lengths$sample[row_i],
                implausible_lengths$region_raw[row_i],
                implausible_lengths$length_mm[row_i]))
  }
  msg("     A length below hatch size cannot be a real capture; treat the whole")
  msg("     metadata row for these fish as suspect, region included.")
  write_csv(implausible_lengths,
            file.path(paths$dir_metadata, "length_qc_flags.csv"))
}


# =============================================================================
# 6b. Age-ambiguity flag for the cohort sensitivity analysis (patch 18)
# =============================================================================
# Otolith ages in 12-126 mm fish are interpretations, and fish grow through a
# sampling season that runs April to October. So a fish is compared only with
# fish of the other age class caught within cohort_settings$age_window_days of
# it in the same year. Flagged fish stay in the main analyses; cohort_strict
# drops them for the sensitivity runs.
parse_collection_date <- function(x) {
  number <- suppressWarnings(as.numeric(as.character(x)))
  out <- as.Date(number, origin = "1899-12-30")            # Excel serial days
  as_text <- is.na(out) & !is.na(x)
  out[as_text] <- as.Date(as.character(x[as_text]), optional = TRUE,
                          tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d"))
  out
}

flag_age_ambiguous <- function(collection_year, collection_date, age, length_mm,
                               window_days) {
  day <- as.numeric(format(collection_date, "%j"))
  ambiguous <- rep(FALSE, length(age))
  for (i in seq_along(age)) {
    if (is.na(age[i]) || is.na(length_mm[i]) || is.na(day[i])) { ambiguous[i] <- NA; next }
    nearby <- which(collection_year == collection_year[i] & !is.na(day) &
                    abs(day - day[i]) <= window_days & !is.na(length_mm) &
                    !is.na(age) & age != age[i])
    if (!length(nearby)) next
    older   <- nearby[age[nearby] > age[i]]
    younger <- nearby[age[nearby] < age[i]]
    ambiguous[i] <- (length(older)   && length_mm[i] >= min(length_mm[older])) ||
                    (length(younger) && length_mm[i] <= max(length_mm[younger]))
  }
  ambiguous
}

date_column <- names(bio_data)[grepl("^date", trimws(names(bio_data)), ignore.case = TRUE)][1]
bio_data$collection_date <- if (!is.na(date_column)) {
  parse_collection_date(bio_data[[date_column]])
} else as.Date(NA)
bio_data$age_ambiguous <- flag_age_ambiguous(bio_data$collection_year, bio_data$collection_date,
                                             bio_data$age, bio_data$length_mm,
                                             cohort_settings$age_window_days)
bio_data$cohort_strict <- ifelse(bio_data$age_ambiguous %in% FALSE, bio_data$cohort, NA)
msg(sprintf("\n  age ambiguity (+/- %d days): %d fish flagged; date column: %s",
            cohort_settings$age_window_days, sum(bio_data$age_ambiguous %in% TRUE),
            ifelse(is.na(date_column), "none found", date_column)))


# =============================================================================
# 7. Assemble the table and decide analysis_keep
# =============================================================================
sample_metadata <- bio_data %>%
  transmute(sample, plate_id, collection_year, collection_date, age, cohort,
            age_ambiguous, cohort_strict,
            region_raw, region = as.character(region),
            across(starts_with("loc_")),
            length_mm) %>%
  mutate(in_vcf     = sample %in% genotyped_samples,
         has_cohort = !is.na(cohort)) %>%
  left_join(qc_exclusions %>% rename(qc_reason = reason), by = "sample") %>%
  mutate(qc_exclude = !is.na(qc_reason))

# A genotyped sample with no workbook row is carried through as a row with
# has_cohort = FALSE rather than dropped, so that it shows up in the accounting
# below instead of quietly reducing n.
samples_without_metadata <- setdiff(genotyped_samples, sample_metadata$sample)
if (length(samples_without_metadata)) {
  msg("\n  !! ", length(samples_without_metadata),
      " genotyped sample(s) have no workbook row: ",
      paste(samples_without_metadata, collapse = ", "))
  sample_metadata <- bind_rows(
    sample_metadata,
    tibble(sample = samples_without_metadata, in_vcf = TRUE,
           has_cohort = FALSE, qc_exclude = FALSE)
  )
}

cohort_sizes <- sample_metadata %>%
  filter(in_vcf, has_cohort, !qc_exclude) %>%
  count(cohort, name = "n")
undersized_cohorts <- cohort_sizes$cohort[cohort_sizes$n < cohort_settings$min_n_analysis]
if (length(undersized_cohorts)) {
  msg("\n  dropping cohort(s) with n < ", cohort_settings$min_n_analysis, ": ",
      paste(undersized_cohorts, collapse = ", "))
}

sample_metadata <- sample_metadata %>%
  mutate(analysis_keep = in_vcf & has_cohort & !qc_exclude &
           !(cohort %in% undersized_cohorts))


# =============================================================================
# 8. Report the accounting
# =============================================================================
# Printed in full because the final n is a number that appears in the abstract,
# and it should be derivable from this block rather than asserted.

msg("\n  --- sample accounting ---")
msg(sprintf("    in filtered VCF                          %4d", sum(sample_metadata$in_vcf)))
msg(sprintf("      - no cohort (missing age / no record)  %4d",
            sum(sample_metadata$in_vcf & !sample_metadata$has_cohort)))
msg(sprintf("      - QC excluded                          %4d",
            sum(sample_metadata$in_vcf & sample_metadata$qc_exclude)))
for (row_i in seq_len(nrow(qc_exclusions))) {
  this_sample <- qc_exclusions$sample[row_i]
  msg(sprintf("          %-9s %-12s %s", this_sample,
              ifelse(this_sample %in% genotyped_samples, "present", "already gone"),
              substr(qc_exclusions$reason[row_i], 1, 58)))
}
msg(sprintf("      - cohort below n = %-2d                   %4d",
            cohort_settings$min_n_analysis,
            sum(sample_metadata$in_vcf & sample_metadata$has_cohort &
                  !sample_metadata$qc_exclude &
                  sample_metadata$cohort %in% undersized_cohorts)))
msg(sprintf("    ANALYSIS SET                             %4d", sum(sample_metadata$analysis_keep)))

analysis_set <- sample_metadata %>% filter(analysis_keep)

msg("\n  cohort sizes:")
msg_table(as.data.frame(analysis_set %>% count(cohort, name = "n")))

# This cross-tabulation is the age structure of the sample, and it is also the
# check that grouping by cohort is doing something different from grouping by
# collection year. Off-diagonal cells are Age-1 fish.
msg("\n  collection year x birth cohort (off-diagonal cells are Age-1 fish):")
msg_table(table(analysis_set$collection_year, analysis_set$cohort,
            dnn = c("collected", "cohort")))

# The cohort-by-location tables below are worth reading carefully before
# interpreting any cohort result. If a cohort is drawn overwhelmingly from one
# location, then "cohort" and "location" are not separable in these data and a
# temporal difference cannot be told apart from a spatial one.
for (scheme_name in names(location_schemes)) {
  msg("\n  cohort x location [", scheme_name, "]  (cells >= ",
      cohort_settings$min_n_estimator, " are usable for Fst):")
  msg_table(table(analysis_set[[paste0("loc_", scheme_name)]], analysis_set$cohort,
              dnn = c(scheme_name, "cohort")))
}

write_csv(sample_metadata, paths$sample_meta)
msg("\n  wrote ", paths$sample_meta)
msg("\n=== 02_sample_metadata.R complete ===")
log_end()
