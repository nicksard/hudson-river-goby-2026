# =============================================================================
# 06_colony_inputs.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Write COLONY sibship input files from the phase 1 panel
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
#
# Description:
#   Replaces the September driver (scripts/archive/Pedigree Pipeline for
#   loop.V6.R). Same COLONY writer (now scripts/lib/colony_dat_writer.R) and
#   the same run settings, read from colony_settings in 00_config.R. What
#   changed:
#
#   * Input is the phase 1 filtered VCF and the analysis set from
#     sample_metadata.csv (analysis_keep == TRUE), so the contaminated sample
#     and the 2024221/2023016 duplicate are out. KING keeps the duplicate as
#     its positive control; COLONY must not, or it is scored as a full-sib pair.
#   * The VCF is parsed directly (readr), not through vcfR. GT 0/1 becomes
#     COLONY alleles 1/2 and missing becomes 0 -- equivalent to the old
#     A/C/G/T -> integer map, since COLONY only needs distinct codes per locus.
#   * "maf" subsets rank loci by minor-allele frequency computed in the
#     analysis set. The September driver ranked by maf(vcf, element = 2)[, 2];
#     in vcfR that column is the count of missing genotypes ("NA"), not the
#     frequency, so those subsets were likely ordered by missingness. Check
#     with head(vcfR::maf(vcf, element = 2)) before saying so in print.
#   * One run per locus subset at the Figure 4 settings, not the 448-run grid.
#     The grid showed run length, precision and ALF updating barely matter.
#
# Inputs:   paths$filtered_vcf, paths$sample_meta
# Outputs:  results/06_colony/<run>.dat, <run>_LOCUS_MAP.csv,
#           run_plan.csv, run_all.bat (Windows), job_list.txt (cluster array)
#
# Usage:    Rscript scripts/06_colony_inputs.R      # from the project root
#           then run_all.bat on Windows, or job_list.txt with an array job
# =============================================================================

source("scripts/00_config.R")
source(project_path("scripts", "lib", "colony_dat_writer.R"))   # colony.dat.create_FLEX.R()

need_pkgs(c("dplyr", "readr"))
suppressPackageStartupMessages({ library(dplyr); library(readr) })

log_start("06_colony_inputs")

# COLONY is Fortran and reads "1e+06" as garbage, so no scientific notation
# anywhere in the .dat (the writer prints the monitor interval as 1e+06 otherwise).
options(scipen = 999)

output_dir <- project_path("results", "06_colony")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
settings <- colony_settings

need(file.exists(paths$filtered_vcf), "Run scripts/01_filter_vcf.R first.")
need(file.exists(paths$sample_meta),  "Run scripts/02_sample_metadata.R first.")


# -----------------------------------------------------------------------------
# 1. Analysis set and genotypes
# -----------------------------------------------------------------------------
sample_metadata <- read_csv(paths$sample_meta, col_types = sample_id_col_types)
analysis_ids    <- sample_metadata$sample[sample_metadata$analysis_keep]

vcf_table <- read_tsv(paths$filtered_vcf, comment = "##", show_col_types = FALSE,
                      col_types = cols(.default = col_character()))
names(vcf_table)[1] <- "CHROM"
need(all(vcf_table$FORMAT == vcf_table$FORMAT[1]) && startsWith(vcf_table$FORMAT[1], "GT"),
     "Expected GT as the first FORMAT field on every line.")

missing_ids <- setdiff(analysis_ids, names(vcf_table))
need(length(missing_ids) == 0,
     "Analysis-set samples missing from the VCF: ", paste(missing_ids, collapse = ", "))
excluded <- setdiff(names(vcf_table)[-(1:9)], analysis_ids)
msg("  VCF ", nrow(vcf_table), " loci x ", ncol(vcf_table) - 9, " samples | analysis set ",
    length(analysis_ids), " | left out: ", paste(excluded, collapse = ", "))

# Genotype calls only ("0/1" from "0/1:37:37,0:40:..."), loci in rows
gt_calls <- as.matrix(vcf_table[, analysis_ids])
gt_calls <- sub(":.*$", "", gt_calls)
allele_1 <- substr(gt_calls, 1, 1)
allele_2 <- substr(gt_calls, 3, 3)
to_colony <- function(a) ifelse(a == "0", "1", ifelse(a == "1", "2", "0"))   # missing -> 0
allele_1 <- to_colony(allele_1)
allele_2 <- to_colony(allele_2)
need(all(nchar(gt_calls) == 3), "Unexpected genotype string (multi-allelic or haploid call?).")

locus_ids <- paste0(vcf_table$CHROM, ":", vcf_table$POS)

# Minor-allele frequency in the analysis set, for the "maf" subsets
n_alt    <- rowSums(allele_1 == "2") + rowSums(allele_2 == "2")
n_called <- rowSums(allele_1 != "0") + rowSums(allele_2 != "0")
alt_freq <- n_alt / n_called
locus_maf <- pmin(alt_freq, 1 - alt_freq)
msg(sprintf("  locus MAF in analysis set: min %.3f, median %.3f", min(locus_maf), median(locus_maf)))


# -----------------------------------------------------------------------------
# 2. Run plan
# -----------------------------------------------------------------------------
run_plan <- data.frame(subset_type = "all", num_loci = nrow(vcf_table), stringsAsFactors = FALSE)
if (isTRUE(settings$run_locus_sweep)) {
  run_plan <- bind_rows(run_plan,
                        expand.grid(subset_type = settings$sweep_subsets,
                                    num_loci = settings$sweep_loci,
                                    stringsAsFactors = FALSE))
}
run_plan <- run_plan %>%
  mutate(run_name = sprintf("%s_%s_%05d", settings$project_name, subset_type, num_loci))

set.seed(random_seed)   # random subsets reproducible from the config seed


# -----------------------------------------------------------------------------
# 3. Write one .dat per run
# -----------------------------------------------------------------------------
write_colony_run <- function(locus_index, run_name) {
  locus_names <- paste0("LOCUS", seq_along(locus_index))

  genotypes <- data.frame(Sample_ID = analysis_ids, stringsAsFactors = FALSE)
  allele_cols <- vector("list", 2 * length(locus_index))
  allele_cols[c(TRUE, FALSE)] <- lapply(locus_index, function(i) allele_1[i, ])
  allele_cols[c(FALSE, TRUE)] <- lapply(locus_index, function(i) allele_2[i, ])
  names(allele_cols) <- as.vector(rbind(paste0(locus_names, "_A"), paste0(locus_names, "_B")))
  genotypes <- cbind(genotypes, as.data.frame(allele_cols, stringsAsFactors = FALSE))

  markers <- as.data.frame(matrix(c(0, settings$error_rate_geno, settings$error_rate_other),
                                  nrow = 3, ncol = length(locus_index)))
  colnames(markers) <- locus_names

  dat_file <- file.path(output_dir, paste0(run_name, ".dat"))
  suppressMessages(colony.dat.create_FLEX.R(
    kids = genotypes, markers = markers, moms = NA, dads = NA,
    analysis_type = "sibship",
    output_filename = dat_file, project_name = run_name,
    update.alfs = settings$update_alfs, run.length = settings$run_length,
    likelihood.precision = settings$likelihood_precision,
    spp.type = settings$spp_type, inbreeding = settings$inbreeding, ploidy = 0,
    fem.gamy = settings$female_gamy, mal.gamy = settings$male_gamy,
    clone = 0, sib.scale = settings$sib_scale, sib.prior = settings$sib_prior,
    known.alfs = 0, run.number = 1, monitor = 0, windows.version = 0,
    full.likelihood = settings$full_likelihood, prob.mom = 0, prob.dad = 0
  ))

  write_csv(data.frame(generic_name = locus_names, original_name = locus_ids[locus_index]),
            file.path(output_dir, paste0(run_name, "_LOCUS_MAP.csv")))
  basename(dat_file)
}

dat_files <- character(nrow(run_plan))
for (run_i in seq_len(nrow(run_plan))) {
  run <- run_plan[run_i, ]
  locus_index <- switch(run$subset_type,
    all    = seq_len(nrow(vcf_table)),
    maf    = order(locus_maf, decreasing = TRUE)[seq_len(run$num_loci)],
    random = sort(sample(nrow(vcf_table), run$num_loci))
  )
  dat_files[run_i] <- write_colony_run(locus_index, run$run_name)
  msg(sprintf("  wrote %-34s %5d loci x %d fish", dat_files[run_i], length(locus_index),
              length(analysis_ids)))
}

run_plan$dat_file <- dat_files
write_csv(run_plan, file.path(output_dir, "run_plan.csv"))
writeLines(dat_files, file.path(output_dir, "job_list.txt"))

# Windows batch: runs the primary (full panel) first, then the sweep
writeLines(c("@echo off",
             "cd /d \"%~dp0\"",
             sprintf("\"%s\" IFN:%s", gsub("/", "\\\\", settings$colony_exe), dat_files),
             "echo ALL COLONY RUNS COMPLETE", "pause"),
           file.path(output_dir, "run_all.bat"))

msg("\n  ", length(dat_files), " runs written to ", output_dir)
msg("  primary run: ", dat_files[1])
msg("\n=== 06_colony_inputs.R complete ===")
log_end()
