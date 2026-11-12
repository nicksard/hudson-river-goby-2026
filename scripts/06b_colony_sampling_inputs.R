# =============================================================================
# 06b_colony_sampling_inputs.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# COLONY inputs for the sample-size sensitivity analysis (2023 cohort)
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-11
# Updated:  2026-11-11
#
# Description:
#   Sibship Nb (COLONY, ~41,000) and LD Nb (~500-1,000) disagree by about
#   fortyfold, whereas Wang (2016) found the two similar when their
#   assumptions hold, with precision driven mainly by sample size. This
#   script sets up the sibship half of a sample-size sensitivity test: fish
#   from the 2023 cohort (the largest, n = 166) are subsampled at several
#   sizes, with replicate draws, and each subsample gets its own COLONY run.
#   The LD half reuses the same fish lists (fish_subsets.csv), so both methods
#   see exactly the same fish.
#
#   Fixed choices, so that only the fish change between runs:
#     * Loci: one random 8,000-locus subset, the same draw as the replicate-1
#       8,000-locus run of 06's locus sweep (results had settled by 8,000).
#     * Settings: the primary-run settings from colony_settings (medium
#       length, medium precision, allele-frequency updating on).
#     * Fish: independent random draws without replacement for each
#       replicate at each size. At n = 166 every replicate has the same fish,
#       so those 5 runs differ only in COLONY's random seed and show its
#       run-to-run variation.
#
#   This script is not part of run_all.R: like 06, its outputs go to the
#   cluster, and it writes to its own folder so it cannot disturb 06/09.
#
# Order of operations:
#   1. Read the analysis set and keep the 2023 cohort.
#   2. Read genotypes from the filtered VCF and draw the 8,000-locus subset.
#   3. Draw fish subsets (sizes x replicates); write fish_subsets.csv.
#   4. Write one COLONY .dat per subset, plus run_plan.csv and job_list.txt.
#
# Inputs:   paths$filtered_vcf, paths$sample_meta (from 01 and 02)
# Outputs:  results/06b_colony_sampling/<run>.dat, run_plan.csv,
#           fish_subsets.csv, locus_subset.csv, job_list.txt
#
# Cluster:  upload the .dat files to PROJ_DIR/colony_dat_files/ and
#           job_list.txt to PROJ_DIR/ in a NEW folder (e.g.
#           rg_sampling_sensitivity/), then from that folder:
#             mkdir -p logs
#             sbatch --array=1-$(wc -l < job_list.txt) colony_array.sh
#           (cluster/colony_array.sh unchanged). Download
#           colony_final_outputs/ into results/06b_colony_sampling/outputs/.
#
# References:
#   Jones OR, Wang J (2010) COLONY: a program for parentage and sibship
#     inference from multilocus genotype data. Molecular Ecology Resources
#     10:551-555.
#   Wang J (2016) A comparison of single-sample estimators of effective
#     population sizes from genetic marker data. Molecular Ecology
#     25:4692-4711.
# =============================================================================

source("scripts/00_config.R")
source(project_path("scripts", "lib", "colony_dat_writer.R"), local = TRUE)   # colony.dat.create_FLEX.R()

need_pkgs(c("dplyr", "readr"))
suppressPackageStartupMessages({ library(dplyr); library(readr) })

log_start("06b_colony_sampling_inputs")
options(scipen = 999)   # COLONY (Fortran) cannot read 1e+06

# ---- Design -----------------------------------------------------------------
focal_cohort    <- 2023
sample_sizes    <- c(25, 50, 75, 100, 125)   # plus the whole cohort, added below
n_replicates    <- 5
n_loci          <- 8000
settings        <- colony_settings
output_dir      <- project_path("results", "06b_colony_sampling")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

need(file.exists(paths$filtered_vcf), "Run scripts/01_filter_vcf.R first.")
need(file.exists(paths$sample_meta),  "Run scripts/02_sample_metadata.R first.")

# ---- 1. Analysis set, 2023 cohort -------------------------------------------
sample_metadata <- read_csv(paths$sample_meta, col_types = sample_id_col_types)
analysis_ids    <- sample_metadata$sample[sample_metadata$analysis_keep]
cohort_ids      <- sample_metadata$sample[sample_metadata$analysis_keep &
                                          !is.na(sample_metadata$cohort) &
                                          sample_metadata$cohort == focal_cohort]
cohort_size     <- length(cohort_ids)
msg(sprintf("  %d cohort: %d fish in the analysis set", focal_cohort, cohort_size))
need(cohort_size > max(sample_sizes), "Cohort smaller than the largest subsample.")
sample_sizes <- c(sample_sizes, cohort_size)

# ---- 2. Genotypes and the locus subset --------------------------------------
vcf_table <- read_tsv(paths$filtered_vcf, comment = "##", show_col_types = FALSE,
                      col_types = cols(.default = col_character()))
names(vcf_table)[1] <- "CHROM"
need(all(cohort_ids %in% names(vcf_table)), "Cohort fish missing from the VCF.")
n_all <- nrow(vcf_table)

# Same draw as 06's locus sweep (replicate 1, 8,000 loci); see 06 section 3
set.seed(random_seed + 1000L * 1L + n_loci)
locus_index <- sort(sample(n_all, n_loci))
locus_ids   <- paste0(vcf_table$CHROM, ":", vcf_table$POS)[locus_index]
write_csv(data.frame(generic_name = paste0("LOCUS", seq_along(locus_index)),
                     original_name = locus_ids),
          file.path(output_dir, "locus_subset.csv"))

gt_calls  <- sub(":.*$", "", as.matrix(vcf_table[locus_index, cohort_ids]))
need(all(nchar(gt_calls) == 3), "Unexpected genotype string (multi-allelic or haploid call?).")
to_colony <- function(a) ifelse(a == "0", "1", ifelse(a == "1", "2", "0"))   # missing -> 0
allele_1  <- to_colony(substr(gt_calls, 1, 1))     # loci in rows, fish in columns
allele_2  <- to_colony(substr(gt_calls, 3, 3))

# ---- 3. Fish subsets --------------------------------------------------------
run_plan <- expand.grid(replicate = seq_len(n_replicates), n_fish = sample_sizes) %>%
  arrange(n_fish, replicate) %>%
  mutate(order    = row_number(),
         run_name = sprintf("s%02d_n%03d_r%d", order, n_fish, replicate),
         num_loci = n_loci, cohort = focal_cohort,
         run_length = settings$run_length, precision = settings$likelihood_precision,
         update_alfs = settings$update_alfs)

fish_subsets <- bind_rows(lapply(seq_len(nrow(run_plan)), function(i) {
  run <- run_plan[i, ]
  set.seed(random_seed + 100L * run$n_fish + run$replicate)
  chosen <- if (run$n_fish == cohort_size) cohort_ids else sort(sample(cohort_ids, run$n_fish))
  data.frame(run_name = run$run_name, n_fish = run$n_fish, replicate = run$replicate,
             sample = chosen)
}))
write_csv(fish_subsets, file.path(output_dir, "fish_subsets.csv"))

# ---- 4. One .dat per subset -------------------------------------------------
old_files <- list.files(output_dir, pattern = "\\.dat$", full.names = TRUE)
if (length(old_files)) { file.remove(old_files); msg("  removed ", length(old_files), " old .dat files") }

locus_names <- paste0("LOCUS", seq_len(n_loci))
markers <- as.data.frame(matrix(c(0, settings$error_rate_geno, settings$error_rate_other),
                                nrow = 3, ncol = n_loci))
colnames(markers) <- locus_names

write_colony_run <- function(run, fish) {
  columns <- match(fish, cohort_ids)
  genotypes <- data.frame(Sample_ID = fish, stringsAsFactors = FALSE)
  allele_cols <- vector("list", 2 * n_loci)
  allele_cols[c(TRUE, FALSE)] <- lapply(seq_len(n_loci), function(l) allele_1[l, columns])
  allele_cols[c(FALSE, TRUE)] <- lapply(seq_len(n_loci), function(l) allele_2[l, columns])
  names(allele_cols) <- as.vector(rbind(paste0(locus_names, "_A"), paste0(locus_names, "_B")))
  genotypes <- cbind(genotypes, as.data.frame(allele_cols, stringsAsFactors = FALSE))

  dat_file <- file.path(output_dir, paste0(run$run_name, ".dat"))
  suppressMessages(colony.dat.create_FLEX.R(
    kids = genotypes, markers = markers, moms = NA, dads = NA,
    analysis_type = "sibship",
    output_filename = dat_file, project_name = run$run_name,
    update.alfs = run$update_alfs, run.length = run$run_length,
    likelihood.precision = run$precision,
    spp.type = settings$spp_type, inbreeding = settings$inbreeding, ploidy = 0,
    fem.gamy = settings$female_gamy, mal.gamy = settings$male_gamy,
    clone = 0, sib.scale = settings$sib_scale, sib.prior = settings$sib_prior,
    known.alfs = 0, run.number = 1, monitor = 0, windows.version = 0,
    full.likelihood = settings$full_likelihood, prob.mom = 0, prob.dad = 0
  ))
  basename(dat_file)
}

run_plan$dat_file <- vapply(seq_len(nrow(run_plan)), function(i) {
  run <- run_plan[i, ]
  set.seed(random_seed + 5000L + run$order)    # COLONY's own seed, drawn inside the writer
  write_colony_run(run, fish_subsets$sample[fish_subsets$run_name == run$run_name])
}, character(1))

# Largest runs first, so the long ones start at once on the cluster
run_plan <- run_plan %>% arrange(desc(n_fish), replicate)
write_csv(run_plan, file.path(output_dir, "run_plan.csv"))
writeLines(run_plan$dat_file, file.path(output_dir, "job_list.txt"))

msg(sprintf("  wrote %d runs: sizes %s x %d replicates, %d loci, RL%d P%d ALF%d",
            nrow(run_plan), paste(sample_sizes, collapse = "/"), n_replicates, n_loci,
            settings$run_length, settings$likelihood_precision, settings$update_alfs))
msg("  ", output_dir)
msg("\n=== 06b_colony_sampling_inputs.R complete ===")
log_end()
