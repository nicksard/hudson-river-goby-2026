# ==============================================================================
#
# R SCRIPT FOR COHORT & LOCATION GENEPOP (.gen) EXPORT (NeEstimator)  -- v2
#
# WHAT CHANGED FROM v1
#   1. QC exclusions applied BEFORE grouping (Part 0). Two samples out.
#   2. Region collapse fixed: only Albany|Troy are relabelled, and it happens
#      BEFORE the n >= 30 filter. v1 overwrote every region unconditionally
#      after filtering, which only worked because nothing else reached 30.
#   3. Writes BOTH grouping schemes -- "All" (all regions) and "Albany_Troy" --
#      so Table 3 comes from one script.
#   4. Optional exclusion of above-confluence sites (Mohawk_River, Lock_C1).
#   5. N in / N out printed at every step, and a manifest CSV written at the
#      end listing every group with its N before and after exclusions, so it
#      is obvious which LDNe runs need repeating.
#   6. Writes a locus -> scaffold map in NeEstimator chromosome-file format,
#      for the "pairs on different chromosomes only" LD option.
#   7. Writes a one-SNP-per-scaffold genepop for each group -- strictly
#      between-scaffold by construction -- as the second between-scaffold set.
#
# Genepop Mapping:
#   VCF 0/0  -> 001001 | 0/1 -> 001002 | 1/1 -> 002002 | ./. -> 000000
#
# ==============================================================================

library(vcfR); library(dplyr); library(readr); library(tidyr)
library(tibble); library(readxl); library(beepr)

setwd("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/pedigree reconstruction")

vcf_file_path <- "../post_RAD analysis/Output/populations.filtered.final.vcf.gz"
bio_data_path <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/plate_barcodes/Input/GobyRad_tabularPlateWell with Bio Data.xlsx"
output_dir    <- "Output/genepop_cohorts_v2"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ==============================================================================
# PART 0: QC DECISIONS
# ==============================================================================

# --- Samples excluded on QC grounds (see 04_vcf_sample_diagnostics.R) ---
#   2023018  contaminated: trimodal allele balance at het sites, highest het in
#            dataset (z = 2.52), 19% missing. Two-individual mixture.
#   2024221  duplicate of 2023016: 99.58% concordance, 0 opposite homozygotes.
#            No extraction record in master metadata; Lock C1 metadata block
#            unreliable. 2023016 (Troy, has extraction record) is retained.
EXCLUDE_SAMPLES <- c("2023018", "2024221")

# --- Above-confluence sites ---
# Mohawk_River (n = 15, rkm 250) and Lock_C1 (n = 8, rkm 253) are upstream of
# the Federal Dam. Set TRUE to drop them from the "All" grouping. Does not
# affect Albany_Troy. Decision pending with Rich.
EXCLUDE_ABOVE_CONFLUENCE <- FALSE
ABOVE_CONFLUENCE_REGIONS <- c("Mohawk_River", "Mohawk River", "Lock_C1", "Lock C1")

MIN_N <- 30


# ==============================================================================
# SECTION 1: GENEPOP WRITER
# ==============================================================================
write_ne_genepop_fast <- function(mat, filename, pop_labels, description,
                                  loci_names, sample_names) {
  con <- file(filename, "w")
  writeLines(description, con)
  writeLines(paste(loci_names, collapse = ", "), con)
  for (p in unique(pop_labels)) {
    writeLines("Pop", con)
    for (i in which(pop_labels == p)) {
      writeLines(paste0(sample_names[i], " , ", paste(mat[, i], collapse = " ")), con)
    }
  }
  close(con)
  message(paste("  wrote:", filename))
}


# ==============================================================================
# SECTION 2: BIO DATA, COHORTS, REGIONS
# ==============================================================================
message("\n--- Loading Biological Data & Inferring Cohorts/Regions ---")
bio_data <- read_xlsx(path = bio_data_path, sheet = 1)
names(bio_data) <- trimws(names(bio_data))
names(bio_data) <- gsub(" ", "_", names(bio_data))
names(bio_data) <- gsub("\\(|\\)|#|`|'", "", names(bio_data))

n0 <- nrow(bio_data)
bio_data <- bio_data[!grepl("NTC", bio_data$ID), ]
bio_data <- bio_data[!grepl("999", bio_data$ID), ]
message(sprintf("  NTC/999 filter: %d -> %d", n0, nrow(bio_data)))

bio_data$Region      <- trimws(as.character(bio_data$Region))
bio_data$Sample_Year <- as.numeric(substr(bio_data$ID, 1, 4))
bio_data$Age_num     <- suppressWarnings(as.numeric(bio_data$Age))
bio_data$Birth_Year  <- bio_data$Sample_Year - bio_data$Age_num
bio_data$Fish_ID     <- as.character(bio_data$Fish_ID)

bio_data_clean <- bio_data[!is.na(bio_data$Birth_Year) & !is.na(bio_data$Fish_ID) &
                             !is.na(bio_data$Region) & bio_data$Region != "", ]
message(sprintf("  usable cohort+region: %d", nrow(bio_data_clean)))

# --- FIX: targeted collapse, BEFORE any filtering ---
bio_data_clean$Region_AT <- ifelse(grepl("Albany|Troy", bio_data_clean$Region),
                                   "Albany_Troy", NA_character_)

# --- QC exclusions ---
n_pre <- nrow(bio_data_clean)
bio_before <- bio_data_clean                       # keep for the manifest
bio_data_clean <- bio_data_clean[!bio_data_clean$Fish_ID %in% EXCLUDE_SAMPLES, ]
message(sprintf("  QC exclusions (%s): %d -> %d",
                paste(EXCLUDE_SAMPLES, collapse = ", "), n_pre, nrow(bio_data_clean)))

if (EXCLUDE_ABOVE_CONFLUENCE) {
  n_pre <- nrow(bio_data_clean)
  bio_data_clean <- bio_data_clean[!bio_data_clean$Region %in% ABOVE_CONFLUENCE_REGIONS, ]
  message(sprintf("  above-confluence exclusion: %d -> %d", n_pre, nrow(bio_data_clean)))
}


# ==============================================================================
# SECTION 3: DEFINE GROUPS  (both schemes)
# ==============================================================================
# Scheme "All":         one group per birth cohort, all regions
# Scheme "Albany_Troy": one group per birth cohort, Albany + Troy only
groups <- bind_rows(
  bio_data_clean %>% mutate(Scheme = "All", Group = as.character(Birth_Year)),
  bio_data_clean %>% filter(!is.na(Region_AT)) %>%
    mutate(Scheme = "Albany_Troy", Group = paste0(Birth_Year, "_AlbanyTroy"))
) %>% select(Scheme, Group, Birth_Year, Fish_ID, Region)

groups_before <- bind_rows(
  bio_before %>% mutate(Scheme = "All", Group = as.character(Birth_Year)),
  bio_before %>% filter(!is.na(Region_AT)) %>%
    mutate(Scheme = "Albany_Troy", Group = paste0(Birth_Year, "_AlbanyTroy"))
) %>% count(Scheme, Group, name = "N_before_QC")


# ==============================================================================
# SECTION 4: VCF, SCAFFOLD MAP, EXPORT
# ==============================================================================
message("\n--- Loading VCF ---")
vcf <- read.vcfR(vcf_file_path, verbose = FALSE)
vcf_samples <- colnames(vcf@gt)[-1]
message(sprintf("  VCF: %d loci x %d samples", nrow(vcf), length(vcf_samples)))

loci_all  <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_")   # unique locus name
scaf_all  <- vcf@fix[, "CHROM"]
rownames(vcf@gt) <- loci_all   # not strictly needed; keeps names consistent

# --- NeEstimator chromosome file: "locus chromosome", one per line ---
# Used with the LD option that restricts to pairs on different chromosomes.
# Check the exact option name in your NeEstimator version's manual.
chrom_map <- file.path(output_dir, "neestimator_chromosome_map.txt")
writeLines(paste(loci_all, scaf_all), chrom_map)
message(sprintf("  wrote chromosome map: %s (%d loci on %d scaffolds)",
                chrom_map, length(loci_all), n_distinct(scaf_all)))

# --- one-SNP-per-scaffold index (highest global MAF within each scaffold) ---
gt_all   <- extract.gt(vcf, element = "GT")
alt_ct   <- (gt_all %in% c("0/1","1/0")) * 1 + (gt_all == "1/1") * 2
alt_ct   <- matrix(alt_ct, nrow = nrow(gt_all))
called   <- !(is.na(gt_all) | gt_all == "./.")
af_alt   <- rowSums(alt_ct * called, na.rm = TRUE) / (2 * rowSums(called))
maf_all  <- pmin(af_alt, 1 - af_alt)
one_per_scaf <- tibble(idx = seq_along(loci_all), scaf = scaf_all, maf = maf_all) %>%
  group_by(scaf) %>% slice_max(maf, n = 1, with_ties = FALSE) %>% ungroup() %>%
  arrange(idx) %>% pull(idx)
message(sprintf("  one-SNP-per-scaffold set: %d loci", length(one_per_scaf)))

# --- export loop ---
manifest <- list()
message(sprintf("\n--- Generating Genepop files (Min N = %d) ---", MIN_N))

for (g in unique(groups$Group)) {
  meta      <- groups %>% filter(Group == g)
  scheme    <- meta$Scheme[1]
  fish      <- meta$Fish_ID
  valid     <- intersect(vcf_samples, fish)
  n_before  <- groups_before$N_before_QC[groups_before$Group == g]
  n_before  <- if (length(n_before)) n_before else NA

  changed   <- !is.na(n_before) && n_before != length(valid)
  message(sprintf("\n%-22s  N before QC = %3s | N now = %3d  %s",
                  g, as.character(n_before), length(valid),
                  ifelse(changed, "<-- CHANGED, rerun LDNe", "unchanged")))

  manifest[[g]] <- tibble(Scheme = scheme, Group = g, N_before_QC = n_before,
                          N_after_QC = length(valid), Changed = changed,
                          Exported = length(valid) >= MIN_N)

  if (length(valid) < MIN_N) { message("  -> below MIN_N, skipped"); next }

  vcf_g <- vcf[, c("FORMAT", valid)]
  gt    <- extract.gt(vcf_g, element = "GT")
  keep  <- rowMeans(is.na(gt) | gt == "./.") < 1.0
  gt    <- gt[keep, , drop = FALSE]
  loci  <- loci_all[keep]

  rec <- gt
  rec[rec == "0/0"] <- "001001"
  rec[rec == "0/1" | rec == "1/0"] <- "001002"
  rec[rec == "1/1"] <- "002002"
  rec[is.na(rec) | rec == "./."] <- "000000"

  # (a) full locus set -- the standard run
  write_ne_genepop_fast(rec, file.path(output_dir, paste0("round_goby_", g, ".gen")),
                        rep(g, ncol(rec)), paste("Round_Goby", g), loci, colnames(rec))

  # (b) one SNP per scaffold -- strictly between-scaffold pairs
  ops <- which(keep)[which(keep) %in% one_per_scaf]
  ops_rows <- match(loci_all[ops], loci)
  write_ne_genepop_fast(rec[ops_rows, , drop = FALSE],
                        file.path(output_dir, paste0("round_goby_", g, "_onePerScaffold.gen")),
                        rep(g, ncol(rec)), paste("Round_Goby", g, "onePerScaffold"),
                        loci[ops_rows], colnames(rec))
}

manifest <- bind_rows(manifest)
write_csv(manifest, file.path(output_dir, "ldne_rerun_manifest.csv"))

message("\n\n==================== RERUN MANIFEST ====================")
print(as.data.frame(manifest))
message("\nGroups marked Changed = TRUE need LDNe rerun. Others are unaffected.")
message("Files per group:")
message("   round_goby_<group>.gen                 -> standard run (all loci)")
message("   round_goby_<group>_onePerScaffold.gen  -> between-scaffold only")
message("   neestimator_chromosome_map.txt         -> for the different-chromosomes LD option")
beep(sound = 8)
