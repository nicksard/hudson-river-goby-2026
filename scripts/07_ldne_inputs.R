# =============================================================================
# 07_ldne_inputs.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# GENEPOP inputs and a NeEstimator batch for LDNe (Table 3)
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
#
# Description:
#   Replaces scripts/archive/vcf2genepop_cohorts_v2.R and
#   scripts/archive/07_neestimator_batch.R for the phase 1 panel.
#
#   * Samples: the analysis set in sample_metadata.csv (analysis_keep == TRUE),
#     so QC exclusions come from 00_config.R, not a second hard-coded list.
#   * Groups: ldne_settings$groupings. Every cohort is written regardless of N;
#     cohorts under cohort_settings$min_n_estimator are flagged in the manifest
#     rather than dropped, so Table 3 can show them with a caveat.
#   * Genotype coding, same as v2:  0/0 -> 001001, 0/1 -> 001002,
#     1/1 -> 002002, ./. -> 000000
#   * NeEstimator: one info file per .gen (LD method only, Pcrit and mating
#     from config) and run_all.bat. Same info-file layout as the old 07.
#   * A locus -> contig map is written for the Thursday between-contig LD check.
#
# Inputs:   paths$filtered_vcf, paths$sample_meta
# Outputs:  results/07_ldne/genepop/<group>.gen, ldne_manifest.csv,
#           locus_contig_map.csv, neestimator/<group>_info.txt, run_all.bat
#
# Usage:    Rscript scripts/07_ldne_inputs.R
#           then double-click results/07_ldne/neestimator/run_all.bat
# =============================================================================

source("scripts/00_config.R")

need_pkgs(c("dplyr", "readr"))
suppressPackageStartupMessages({ library(dplyr); library(readr) })

log_start("07_ldne_inputs")

settings    <- ldne_settings
output_dir  <- project_path("results", "07_ldne")
genepop_dir <- file.path(output_dir, "genepop")
ne_dir      <- file.path(output_dir, "neestimator")
for (d in c(genepop_dir, ne_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

need(file.exists(paths$filtered_vcf), "Run scripts/01_filter_vcf.R first.")
need(file.exists(paths$sample_meta),  "Run scripts/02_sample_metadata.R first.")


# -----------------------------------------------------------------------------
# 1. Analysis set and genotypes
# -----------------------------------------------------------------------------
sample_metadata <- read_csv(paths$sample_meta, col_types = sample_id_col_types)
analysis_set    <- sample_metadata %>% filter(analysis_keep)

vcf_table <- read_tsv(paths$filtered_vcf, comment = "##", show_col_types = FALSE,
                      col_types = cols(.default = col_character()))
names(vcf_table)[1] <- "CHROM"
need(startsWith(vcf_table$FORMAT[1], "GT"), "Expected GT as the first FORMAT field.")
need(all(analysis_set$sample %in% names(vcf_table)),
     "Some analysis-set samples are missing from the VCF.")

gt_calls <- sub(":.*$", "", as.matrix(vcf_table[, analysis_set$sample]))
genepop_code <- c("0/0" = "001001", "0/1" = "001002", "1/0" = "001002",
                  "1/1" = "002002", "./." = "000000")
unknown_calls <- setdiff(unique(as.vector(gt_calls)), names(genepop_code))
need(length(unknown_calls) == 0, "Unexpected genotype calls: ", paste(unknown_calls, collapse = ", "))
genepop_matrix <- matrix(genepop_code[gt_calls], nrow = nrow(gt_calls),
                         dimnames = dimnames(gt_calls))

locus_names <- paste0(vcf_table$CHROM, "_", vcf_table$POS)
msg("  ", nrow(vcf_table), " loci x ", nrow(analysis_set), " fish (analysis set)")

write_csv(data.frame(locus = locus_names, contig = vcf_table$CHROM, pos = vcf_table$POS),
          file.path(output_dir, "locus_contig_map.csv"))
msg("  ", length(unique(vcf_table$CHROM)), " contigs in locus_contig_map.csv")


# -----------------------------------------------------------------------------
# 2. Groups
# -----------------------------------------------------------------------------
group_members <- list()
if ("all" %in% settings$groupings) {
  group_members[["All"]] <- analysis_set$sample
}
if ("cohort" %in% settings$groupings) {
  for (cohort_year in sort(unique(analysis_set$cohort))) {
    group_members[[paste0("Cohort_", cohort_year)]] <-
      analysis_set$sample[analysis_set$cohort == cohort_year]
  }
}
### ADDED (patch 18): cohorts without age-ambiguous fish
if ("cohort_strict" %in% settings$groupings && "cohort_strict" %in% names(analysis_set)) {
  strict <- analysis_set %>% filter(!is.na(cohort_strict))
  for (cohort_year in sort(unique(strict$cohort_strict))) {
    group_members[[paste0("CohortStrict_", cohort_year)]] <-
      strict$sample[strict$cohort_strict == cohort_year]
  }
}
if ("albany_troy_cohort" %in% settings$groupings) {
  pooled <- analysis_set %>% filter(loc_albany_troy == "Albany_Troy")
  group_members[["AlbanyTroy_All"]] <- pooled$sample
  for (cohort_year in sort(unique(pooled$cohort))) {
    group_members[[paste0("AlbanyTroy_Cohort_", cohort_year)]] <-
      pooled$sample[pooled$cohort == cohort_year]
  }
}


### NEW -- near-complete subset (robustness check on the missing-data revision)
group_loci <- list()
if ("near_complete" %in% settings$groupings) {
  is_missing    <- gt_calls == "./."
  fish_missing  <- colMeans(is_missing)
  kept_fish     <- names(fish_missing)[fish_missing <= settings$near_complete_max_fish_missing]
  locus_missing <- rowMeans(is_missing[, kept_fish, drop = FALSE])
  kept_loci     <- which(locus_missing <= settings$near_complete_max_locus_missing)
  residual      <- mean(is_missing[kept_loci, kept_fish])
  group_members[["NearComplete_All"]] <- kept_fish
  group_loci[["NearComplete_All"]]    <- kept_loci
  msg(sprintf("  near-complete subset: %d fish (<= %.0f%% missing), %d loci (<= %.0f%%), residual missing %.2f%%",
              length(kept_fish), 100 * settings$near_complete_max_fish_missing,
              length(kept_loci), 100 * settings$near_complete_max_locus_missing, 100 * residual))
}


# -----------------------------------------------------------------------------
# 3. Write .gen files and NeEstimator info files
# -----------------------------------------------------------------------------
write_genepop <- function(group_name, sample_ids, loci = seq_len(nrow(genepop_matrix))) {   ### CHANGED: optional locus subset
  gen_file <- file.path(genepop_dir, paste0(group_name, ".gen"))
  con <- file(gen_file, "w")
  on.exit(close(con))
  writeLines(sprintf("Hudson River round goby | %s | n = %d | %d loci | phase 1 panel",
                     group_name, length(sample_ids), length(loci)), con)
  writeLines(paste(locus_names[loci], collapse = ", "), con)   # same layout as v2, which NeEstimator read
  writeLines("Pop", con)
  for (id in sample_ids) {
    writeLines(paste0(id, " , ", paste(genepop_matrix[loci, id], collapse = " ")), con)
  }
  basename(gen_file)
}

# NeEstimator wants Windows paths with a trailing backslash
win_dir <- function(p) paste0(gsub("/", "\\\\", normalizePath(p, winslash = "/", mustWork = FALSE)), "\\")

manifest <- data.frame()
info_files <- character()
for (group_name in names(group_members)) {
  ids <- group_members[[group_name]]
  loci <- if (is.null(group_loci[[group_name]])) seq_len(nrow(genepop_matrix)) else group_loci[[group_name]]
  gen_name <- write_genepop(group_name, ids, loci)

  info_file <- file.path(ne_dir, paste0(group_name, "_info.txt"))
  writeLines(c(
    "1",                                   # method: 1 = LD only
    win_dir(genepop_dir),                  # input directory
    gen_name,                              # input file
    "2",                                   # format: 2 = GENEPOP
    win_dir(ne_dir),                       # output directory
    paste0(group_name, "_Ne.txt"),         # output file
    "1",                                   # number of critical values
    as.character(settings$pcrit),          # critical value(s)
    as.character(settings$mating)          # mating model
  ), info_file)
  info_files <- c(info_files, basename(info_file))
  ### NEW -- one .bat per group, so a single run can be repeated or added alone
  writeLines(c("@echo off", "cd /d \"%~dp0\"",
               sprintf("\"%s\" i:%s", gsub("/", "\\\\", settings$ne_exe), basename(info_file)),
               "pause"),
             file.path(ne_dir, paste0("run_", group_name, ".bat")))

  manifest <- rbind(manifest, data.frame(
    group = group_name, n = length(ids), n_loci = length(loci), gen_file = gen_name,
    ne_output = paste0(group_name, "_Ne.txt"),
    flag = ifelse(length(ids) < cohort_settings$min_n_estimator,
                  paste0("n < ", cohort_settings$min_n_estimator, ": low precision"), "")))
  msg(sprintf("  %-22s n = %3d  %s", group_name, length(ids), tail(manifest$flag, 1)))
}

write_csv(manifest, file.path(output_dir, "ldne_manifest.csv"))
msg_table(manifest, row.names = FALSE)

writeLines(c("@echo off",
             "cd /d \"%~dp0\"",
             sprintf("\"%s\" i:%s", gsub("/", "\\\\", settings$ne_exe), info_files),
             "echo ALL NEESTIMATOR RUNS COMPLETE", "pause"),
           file.path(ne_dir, "run_all.bat"))

msg("\n  wrote ", length(info_files), " runs; start with ", file.path(ne_dir, "run_all.bat"))
msg("\n=== 07_ldne_inputs.R complete ===")
log_end()
