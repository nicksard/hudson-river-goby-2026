################################################################################
### KING-ROBUST KINSHIP DIAGNOSTIC -- COLONY FALSE-NEGATIVE TEST
###
### Author: Nicholas Sard, Ph.D.
### Assistant: Claude
### Date: September 2026
###
### Description: Model-free test of whether COLONY is missing sib relationships.
### Computes KING-robust pairwise kinship straight from the filtered VCF -- no
### pedigree reconstruction, no error model, no sibship prior -- and asks whether
### there is an excess of dyads near phi = 0.125 (half-sib).
###
### ---------------------------------------------------------------------------
### WHY
### ---------------------------------------------------------------------------
### The parental detection spectrum (f1 ~ 634, f2 ~ 6, f3 = f4 = 0) is jointly
### consistent with EITHER:
###     (a) Np ~ 34,000 and no false negatives, OR
###     (b) Np ~ 2,500 and a ~90% per-dyad false-negative rate in COLONY
### and these are non-identifiable from the pedigree alone. KING sidesteps
### COLONY entirely.
###
### PREDICTION under (b): a modest but detectable excess of dyads near
### phi = 0.125 that COLONY did not assign. A clean unimodal distribution
### centred on zero kills the false-negative hypothesis outright.
###
### ---------------------------------------------------------------------------
### CONVENTIONS INHERITED FROM 02_vcf_basic_pop_gen_analysis_cohorts.R
### ---------------------------------------------------------------------------
###   - ID is the plate/well identifier; Fish_ID is the VCF sample name.
###     They are DIFFERENT columns. Sample_Year comes from substr(ID, 1, 4).
###   - Age is stored as character and needs suppressWarnings(as.numeric()).
###   - Birth_Year <- Sample_Year - Age_num
###   - Region, not Location. Albany|Troy are collapsed to "Albany_Troy".
###   - NTC and 999 rows are dropped from bio_data by grepl on ID.
###   - Cohorts with N < DROP_COHORT_BELOW are dropped (removes the lone 2020 fish).
################################################################################


# ==============================================================================
# Part 0: Options
# ==============================================================================

setwd("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis")

vcf_file_path <- "Output/populations.filtered.final.vcf.gz"
bio_data_path <- paste0("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/",
                        "plate_barcodes/Input/",
                        "GobyRad_tabularPlateWell with Bio Data.xlsx")

# Optional: COLONY dyad export, so the per-dyad false-negative rate can be
# measured rather than inferred. Needs two ID columns holding Fish_IDs.
# Set to NULL if you don't have it handy yet.
colony_pairs_path <- NULL
# colony_pairs_path <- "../pedigree reconstruction/Output/colony_dyads.csv"
colony_id_cols    <- c("OffspringID1", "OffspringID2")

output_dir        <- "Output/king_kinship"

DROP_COHORT_BELOW <- 5             # mirrors script 02
MAF_THRESHOLDS    <- c(0.05, 0.10) # KING degrades at low MAF; run both
MISSING_RATE      <- 0.10
DO_LD_PRUNE       <- TRUE
LD_THRESH         <- 0.2
NUM_THREADS       <- 4

# Np grid spanning the Np / false-negative trade-off curve
NP_GRID <- c(400, 642, 1000, 2000, 2500, 4000, 8000, 17000, 34000)

EXPECTED_COHORT_N <- c("2021" = 44, "2022" = 16, "2023" = 166, "2024" = 94)


# ==============================================================================
# Part 1: Setup & Data Loading
# ==============================================================================

### Step 1.1: Load Libraries ----
# ------------------------------------------------------------------------------
message("--- Part 1: Setup & Data Loading ---")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("SNPRelate", "gdsfmt")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
}

library(SNPRelate)
library(gdsfmt)
library(tidyverse)
library(readxl)
library(beepr)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260906)


### Step 1.2: Convert VCF to GDS ----
# ------------------------------------------------------------------------------
# SNPRelate needs its own on-disk format. Conversion is done once and cached.
message("\n--- Step 1.2: Converting VCF to GDS ---")

gds_file <- file.path(output_dir, "goby_king.gds")
if (!file.exists(gds_file)) {
  snpgdsVCF2GDS(vcf_file_path, gds_file, method = "biallelic.only", verbose = TRUE)
} else {
  message("  Existing GDS found; delete it to force a rebuild.")
}
genofile <- snpgdsOpen(gds_file)

vcf_samples <- read.gdsn(index.gdsn(genofile, "sample.id"))
message(sprintf("  GDS holds %d samples and %d SNPs.",
                length(vcf_samples),
                length(read.gdsn(index.gdsn(genofile, "snp.id")))))


### Step 1.3: Load Bio Data & Infer Birth Cohorts ----
# ------------------------------------------------------------------------------
# Name cleaning and cohort derivation copied verbatim from script 02 so all
# scripts agree on what a cohort is.
message("\n--- Step 1.3: Loading Biological Data & Inferring Cohorts ---")

bio_data <- read_xlsx(path = bio_data_path, sheet = 1)
names(bio_data) <- trimws(names(bio_data))
names(bio_data) <- trimws(gsub(pattern = " ", replacement = "_", x = names(bio_data)))
names(bio_data) <- trimws(gsub(pattern = "\\(|\\)", replacement = "", x = names(bio_data)))
names(bio_data) <- trimws(gsub(pattern = "_\\[|\\]_", replacement = "", x = names(bio_data)))

# Clean bad IDs. NOTE: grepl("999") is a wide net -- it removes ANY ID containing
# 999, not just 2024999. Kept identical to 02 for consistency, but the count is
# printed so a silent over-drop is visible.
n_before <- nrow(bio_data)
bio_data <- bio_data[!grepl(pattern = "NTC", x = bio_data$ID), ]
bio_data <- bio_data[!grepl(pattern = "999", x = bio_data$ID), ]
message(sprintf("  Dropped %d NTC/999 bio records (%d -> %d).",
                n_before - nrow(bio_data), n_before, nrow(bio_data)))

# Infer cohorts
bio_data$Sample_Year <- as.numeric(substr(x = bio_data$ID, start = 1, stop = 4))
bio_data$Age_num     <- suppressWarnings(as.numeric(bio_data$Age))
bio_data$Birth_Year  <- bio_data$Sample_Year - bio_data$Age_num

message(sprintf("  Bio records: %d | with usable Age: %d | missing Age: %d",
                nrow(bio_data), sum(!is.na(bio_data$Age_num)),
                sum(is.na(bio_data$Age_num))))

# Region cleaning, same as 02
bio_data$Region <- trimws(bio_data$Region)
bio_data$Region[grepl(pattern = "Albany|Troy", x = bio_data$Region)] <- "Albany_Troy"

# Length, for annotating candidate dyads only -- not used to assign cohorts
len_col <- intersect(c("Length_mm","Total_Length_mm","TL_mm","Length","TL","Length_MM"),
                     names(bio_data))
bio_data$Length_QC <- if (length(len_col) > 0) {
  suppressWarnings(as.numeric(bio_data[[len_col[1]]]))
} else NA_real_


### Step 1.4: Build the sample lookup and subset ----
# ------------------------------------------------------------------------------
message("\n--- Step 1.4: Matching VCF samples to cohorts ---")

samp_lu <- tibble(Sample = vcf_samples) %>%
  left_join(
    bio_data %>%
      transmute(Sample       = as.character(Fish_ID),
                Collect_Year = as.character(Sample_Year),
                Age          = Age_num,
                Cohort       = as.character(Birth_Year),
                Region       = Region,
                Length_mm    = Length_QC) %>%
      distinct(Sample, .keep_all = TRUE),
    by = "Sample"
  )

n_unmatched <- sum(is.na(samp_lu$Cohort))
if (n_unmatched > 0) {
  message(sprintf("  !! %d genotyped fish have no cohort and will be dropped:", n_unmatched))
  print(samp_lu$Sample[is.na(samp_lu$Cohort)])
}

cohort_n <- samp_lu %>% filter(!is.na(Cohort)) %>% count(Cohort, name = "N")
message("\n  Cohort sample sizes:")
print(as.data.frame(cohort_n))
for (ch in names(EXPECTED_COHORT_N)) {
  got <- cohort_n$N[cohort_n$Cohort == ch]; got <- if (length(got) == 0) 0 else got
  message(sprintf("  Cohort %s: observed %3d | Table 3 expects %3d  [%s]",
                  ch, got, EXPECTED_COHORT_N[[ch]],
                  ifelse(got == EXPECTED_COHORT_N[[ch]], "OK", "MISMATCH")))
}

small <- cohort_n$Cohort[cohort_n$N < DROP_COHORT_BELOW]
if (length(small) > 0) {
  message(sprintf("  Dropping cohort(s) with N < %d: %s",
                  DROP_COHORT_BELOW, paste(small, collapse = ", ")))
}

samp_lu  <- samp_lu %>% filter(!is.na(Cohort), !(Cohort %in% small))
keep_ids <- samp_lu$Sample
message(sprintf("  Retained %d individuals across %d cohorts.",
                length(keep_ids), n_distinct(samp_lu$Cohort)))

write_csv(samp_lu, file.path(output_dir, "king_sample_lookup.csv"))


# ==============================================================================
# Part 2: KING-Robust Kinship
# ==============================================================================
message("\n--- Part 2: KING-robust kinship ---")

### Step 2.1: Run KING across MAF x pruning settings ----
# ------------------------------------------------------------------------------
# Running the grid guards against the obvious reviewer objection that the result
# is an artefact of marker choice. KING-robust is used rather than KING-homo
# because it does not assume homogeneous ancestry.

run_king <- function(maf, prune) {
  
  snps <- NULL
  if (prune) {
    set.seed(1)
    pr <- snpgdsLDpruning(genofile, sample.id = keep_ids, snp.id = NULL,
                          autosome.only = FALSE,          ### FIX
                          maf = maf, missing.rate = MISSING_RATE,
                          ld.threshold = LD_THRESH,
                          num.thread = NUM_THREADS, verbose = FALSE)
    snps <- unlist(pr, use.names = FALSE)
  }
  
  king <- snpgdsIBDKING(genofile, sample.id = keep_ids, snp.id = snps,
                        autosome.only = FALSE,            ### FIX
                        maf = maf, missing.rate = MISSING_RATE,
                        type = "KING-robust",
                        num.thread = NUM_THREADS, verbose = FALSE)
  
  n_snp <- length(king$snp.id)
  
  snpgdsIBDSelection(king, kinship.cutoff = NaN) %>%
    as_tibble() %>%
    mutate(maf = maf, pruned = prune, n_snp = n_snp)
}

grid <- expand.grid(maf   = MAF_THRESHOLDS,
                    prune = if (DO_LD_PRUNE) c(FALSE, TRUE) else FALSE,
                    stringsAsFactors = FALSE)


k <- snpgdsIBDKING(genofile, sample.id = keep_ids, autosome.only = FALSE,
                   maf = 0.05, missing.rate = 0.10,
                   type = "KING-robust", num.thread = 4, verbose = TRUE)
length(k$snp.id)

king_long <- purrr::pmap_dfr(list(grid$maf, grid$prune), function(m, p) {
  message(sprintf("  MAF >= %.2f | pruned = %s", m, p))
  run_king(m, p)
})
beep(sound = 2)

# Reference setting used for the headline figures and tables
REF_MAF <- min(MAF_THRESHOLDS)


### Step 2.2: Annotate dyads with cohort and region ----
# ------------------------------------------------------------------------------
king_long <- king_long %>%
  left_join(samp_lu %>% select(Sample, Cohort1 = Cohort, Region1 = Region,
                               Len1 = Length_mm),
            by = c("ID1" = "Sample")) %>%
  left_join(samp_lu %>% select(Sample, Cohort2 = Cohort, Region2 = Region,
                               Len2 = Length_mm),
            by = c("ID2" = "Sample")) %>%
  mutate(same_cohort = Cohort1 == Cohort2,
         cohort      = ifelse(same_cohort, Cohort1, NA_character_),
         by_gap      = abs(as.numeric(Cohort1) - as.numeric(Cohort2)),
         same_region = Region1 == Region2,
         pair_key    = paste(pmin(ID1, ID2), pmax(ID1, ID2), sep = "_"))

# Optional: flag dyads COLONY already assigned
if (!is.null(colony_pairs_path) && file.exists(colony_pairs_path)) {
  cp <- read_csv(colony_pairs_path, show_col_types = FALSE)
  a <- as.character(cp[[colony_id_cols[1]]])
  b <- as.character(cp[[colony_id_cols[2]]])
  colony_keys <- paste(pmin(a, b), pmax(a, b), sep = "_")
  king_long$colony_assigned <- king_long$pair_key %in% colony_keys
  message(sprintf("  Loaded %d COLONY dyads.", length(colony_keys)))
} else {
  king_long$colony_assigned <- NA
  message("  No COLONY dyad file supplied; false-negative rate will be inferred, not measured.")
}

write_csv(king_long, file.path(output_dir, "king_pairwise_all.csv"))


# --- effective locus count and null SD ---
af    <- snpgdsSNPRateFreq(genofile, sample.id = keep_ids, snp.id = k$snp.id)$AlleleFreq
h     <- 2 * af * (1 - af)
L_eff <- sum(h)^2 / sum(h^2)
sd0   <- sqrt(3) / (2 * sqrt(L_eff))
cat(sprintf("L = %d | L_eff = %.0f | null SD = %.4f | half-sib at %.1f SD\n",
            length(h), L_eff, sd0, 0.125 / sd0))

# --- pairwise, annotated with cohort ---
kl <- snpgdsIBDSelection(k, kinship.cutoff = NaN) %>%
  as_tibble() %>%
  left_join(samp_lu %>% select(Sample, C1 = Cohort), by = c("ID1" = "Sample")) %>%
  left_join(samp_lu %>% select(Sample, C2 = Cohort), by = c("ID2" = "Sample")) %>%
  mutate(same_cohort = C1 == C2)

# --- excess via negative-tail mirroring ---
for (cut in c(0.0442, 0.0625, 0.0884, 0.125)) {
  a <- sum(kl$kinship > cut); b <- sum(kl$kinship < -cut)
  cat(sprintf("all pairs  cutoff %.4f: pos %5d | mirror %5d | excess %5d\n", cut, a, b, a - b))
}
cat("\n")
for (ch in c("2021","2023","2024")) {
  s <- kl %>% filter(same_cohort, C1 == ch)
  for (cut in c(0.0625, 0.0884, 0.125)) {
    a <- sum(s$kinship > cut); b <- sum(s$kinship < -cut)
    cat(sprintf("cohort %s (%d pairs) cutoff %.4f: excess %4d\n",
                ch, nrow(s), cut, a - b))
  }
}

cat(sprintf("\nmax kinship = %.4f | min = %.4f\n", max(kl$kinship), min(kl$kinship)))

# ==============================================================================
# Part 3: Excess Relatedness via Negative-Tail Mirroring
# ==============================================================================
# For truly unrelated pairs KING kinship is centred on zero and approximately
# symmetric. Real relatives can only push mass to the RIGHT, so
#
#       excess(x) = count(phi > x) - count(phi < -x)
#
# is an assumption-free lower bound on the number of genuinely related dyads
# above x. No distributional assumption, no simulation, and it is immune to the
# "your markers are uninformative" objection -- weak markers widen BOTH tails.
message("\n--- Part 3: Excess relatedness ---")

excess_at <- function(df, cut) {
  pos <- sum(df$kinship >  cut, na.rm = TRUE)
  neg <- sum(df$kinship < -cut, na.rm = TRUE)
  tibble(cutoff = cut, n_pos = pos, n_mirror = neg,
         excess = pos - neg, excess_se = sqrt(pos + neg))
}

# Standard KING class boundaries (midpoints on the log2 scale)
cuts <- c(0.0442, 0.0625, 0.0884, 0.125, 0.177, 0.25)

excess_overall <- king_long %>%
  group_by(maf, pruned, n_snp) %>%
  group_modify(~ map_dfr(cuts, function(cc) excess_at(.x, cc))) %>%
  ungroup()

excess_by_cohort <- king_long %>%
  filter(same_cohort) %>%
  group_by(maf, pruned, cohort) %>%
  group_modify(~ map_dfr(cuts, function(cc) excess_at(.x, cc))) %>%
  ungroup()

write_csv(excess_overall,   file.path(output_dir, "king_excess_overall.csv"))
write_csv(excess_by_cohort, file.path(output_dir, "king_excess_by_cohort.csv"))

message("\n--- Excess dyads above threshold (all pairs) ---")
print(as.data.frame(excess_overall))

message("\n--- Excess dyads within cohort ---")
print(as.data.frame(excess_by_cohort %>% filter(maf == REF_MAF, !pruned)))


# ==============================================================================
# Part 4: Expected Half-Sib Dyads Across the Np Trade-Off Curve
# ==============================================================================
# Two offspring drawn from Np parents (Np/2 per sex) share at least one parent
# with probability ~ 4/Np, and share both with probability ~ 4/Np^2.
# Comparing these to the observed excess near phi = 0.125 reads off which point
# on the Np / false-negative trade-off curve the data actually support.
#
# CAVEAT: assumes random sampling of offspring. Sampling here is spatially
# clustered, which inflates sib co-detection, so these are conservative.
message("\n--- Part 4: Expected dyads under competing Np ---")

expected_dyads <- expand_grid(Cohort = cohort_n$Cohort[cohort_n$N >= DROP_COHORT_BELOW],
                              Np = NP_GRID) %>%
  left_join(cohort_n, by = "Cohort") %>%
  mutate(n_pairs   = choose(N, 2),
         E_halfsib = n_pairs * 4 / Np,
         E_fullsib = n_pairs * 4 / Np^2)

write_csv(expected_dyads, file.path(output_dir, "expected_dyads_by_Np.csv"))

print(expected_dyads %>%
        select(Cohort, N, Np, E_halfsib) %>%
        pivot_wider(names_from = Np, values_from = E_halfsib, names_prefix = "Np=") %>%
        mutate(across(where(is.numeric), ~round(.x, 1))) %>%
        as.data.frame())


# ==============================================================================
# Part 5: Plots
# ==============================================================================
message("\n--- Part 5: Generating plots ---")

class_lines <- tibble(
  phi   = c(0.0625, 0.125, 0.25),
  label = c("1st cousin", "half-sib / GP-GO", "full-sib / parent-offspring")
)

ref <- king_long %>% filter(maf == REF_MAF, !pruned)

# 5.1 Main histogram, log1p y so the right tail is legible
p_king_all <- ggplot(ref, aes(x = kinship)) +
  geom_histogram(binwidth = 0.004, fill = "grey35", colour = NA) +
  geom_vline(data = class_lines, aes(xintercept = phi),
             linetype = "dashed", colour = "firebrick") +
  geom_text(data = class_lines, aes(x = phi, y = Inf, label = label),
            angle = 90, vjust = -0.4, hjust = 1.05, size = 2.8, colour = "firebrick") +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 10, 100, 1000, 10000, 50000)) +
  labs(x = "KING-robust kinship coefficient", y = "Number of dyads (log1p)",
       title = "Pairwise kinship, all dyads",
       subtitle = sprintf("MAF >= %.2f, unpruned. Right-tail mass above ~0.09 is the COLONY false-negative signal.",
                          REF_MAF)) +
  theme_bw()
print(p_king_all)

# 5.2 Within-cohort dyads, faceted -- this is the relationship class COLONY was
#     asked to find, so it is the primary panel
p_king_cohort <- ggplot(ref %>% filter(same_cohort), aes(x = kinship)) +
  geom_histogram(binwidth = 0.004, fill = "steelblue4", colour = NA) +
  geom_vline(xintercept = c(0.0625, 0.125, 0.25),
             linetype = "dashed", colour = "firebrick") +
  coord_cartesian(xlim = c(-0.15, 0.35)) +
  facet_wrap(~ cohort, scales = "free_y", ncol = 2) +
  labs(x = "KING-robust kinship coefficient", y = "Number of dyads",
       title = "Within-cohort dyads by birth cohort",
       subtitle = "Any asymmetry above zero, judged against the mirrored negative tail, is real relatedness") +
  theme_bw()
print(p_king_cohort)

# 5.3 Sensitivity across MAF and LD pruning
p_king_sens <- ggplot(king_long, aes(x = kinship,
                                     colour = interaction(maf, pruned, sep = " / prune="))) +
  stat_bin(binwidth = 0.004, geom = "step", position = "identity") +
  scale_y_continuous(trans = "log1p") +
  coord_cartesian(xlim = c(-0.15, 0.35)) +
  geom_vline(xintercept = c(0.0625, 0.125, 0.25), linetype = "dashed", colour = "grey40") +
  labs(x = "KING-robust kinship", y = "Dyads (log1p)", colour = "MAF / pruned",
       title = "Kinship distribution across MAF and LD-pruning settings") +
  theme_bw()
print(p_king_sens)

# 5.4 IBS0 vs kinship. KING cannot separate full-sib from parent-offspring on
#     kinship alone (both ~0.25), but PO pairs have IBS0 near zero and FS pairs
#     do not. Relevant because Age-1 fish could in principle be parents of Age-0
#     fish, which bears directly on the grandparent / half-sib ambiguity.
p_king_ibs0 <- ggplot(ref %>% filter(kinship > 0.05),
                      aes(x = IBS0, y = kinship, colour = same_cohort)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_hline(yintercept = c(0.0884, 0.177), linetype = "dashed", colour = "grey40") +
  labs(x = "Proportion of loci sharing zero alleles (IBS0)",
       y = "KING-robust kinship", colour = "Same cohort",
       title = "IBS0 vs kinship for candidate relatives",
       subtitle = "Parent-offspring: high kinship, IBS0 near zero. Full-sib: high kinship, IBS0 > 0.") +
  theme_bw()
print(p_king_ibs0)

ggsave(file.path(output_dir, "Fig_KING_histogram_all.png"),    p_king_all,
       width = 9, height = 5.5, dpi = 300)
ggsave(file.path(output_dir, "Fig_KING_histogram_cohort.png"), p_king_cohort,
       width = 9, height = 6.5, dpi = 300)
ggsave(file.path(output_dir, "Fig_KING_sensitivity.png"),      p_king_sens,
       width = 9, height = 5.5, dpi = 300)
ggsave(file.path(output_dir, "Fig_KING_IBS0_vs_kinship.png"),  p_king_ibs0,
       width = 8, height = 5.5, dpi = 300)


# ==============================================================================
# Part 6: Candidate Relatives COLONY Missed
# ==============================================================================
message("\n--- Part 6: Candidate relatives ---")

candidates <- ref %>%
  filter(kinship > 0.0884) %>%
  arrange(desc(kinship)) %>%
  select(ID1, ID2, kinship, IBS0, Cohort1, Cohort2, by_gap,
         Region1, Region2, same_region, Len1, Len2, colony_assigned, n_snp)

write_csv(candidates, file.path(output_dir, "king_candidate_relatives.csv"))

message(sprintf("  %d dyads above phi = 0.0884 (half-sib lower bound).", nrow(candidates)))
if (any(!is.na(candidates$colony_assigned)) && nrow(candidates) > 0) {
  message(sprintf("  COLONY assigned %d, missed %d -> measured per-dyad FN rate = %.2f",
                  sum(candidates$colony_assigned, na.rm = TRUE),
                  sum(!candidates$colony_assigned, na.rm = TRUE),
                  mean(!candidates$colony_assigned, na.rm = TRUE)))
}

# Are the candidates clustered by sampling event? Clustered sampling inflates
# sib co-detection relative to the random-sampling model in Part 4.
if (nrow(candidates) > 0) {
  message("\n  Candidate dyads by region pairing:")
  print(candidates %>% count(Region1, Region2, sort = TRUE) %>% as.data.frame())
  message("\n  Candidate dyads by cohort gap (0 = within cohort):")
  print(candidates %>% count(by_gap) %>% as.data.frame())
}


# ==============================================================================
# Part 7: PDF Report
# ==============================================================================
pdf(file.path(output_dir, "KING_Kinship_Diagnostic_Report.pdf"), width = 14, height = 11)
print(p_king_all)
print(p_king_cohort)
print(p_king_sens)
print(p_king_ibs0)
dev.off()

message("\nReport saved to: ", file.path(output_dir, "KING_Kinship_Diagnostic_Report.pdf"))

snpgdsClose(genofile)


################################################################################
### HOW TO READ THE OUTPUT
###
### Compare the excess at cutoff = 0.0884 in king_excess_by_cohort.csv against
### expected_dyads_by_Np.csv. For the 2023 cohort (N = 166, 13,695 pairs) the
### expected half-sib counts are roughly:
###     Np =    400 -> 137
###     Np =  2,500 ->  22
###     Np = 34,000 ->   1.6
###
### OUTCOME A -- excess of roughly 20-40 in the 2023 cohort
###   Consistent with Np ~ 2,000-2,700 and a ~90% COLONY false-negative rate.
###   This is the working hypothesis. Section 5 becomes the non-identifiability
###   argument and Section 7's connectivity story is corroborated, not
###   contradicted.
###
### OUTCOME B -- excess near zero, distribution symmetric about zero
###   No hidden relatives. COLONY was right, Np really is very large, and the
###   sibship and iChao estimates should be taken at face value. Section 5 has
###   to be rewritten around WHERE the parents are rather than around estimator
###   failure, and the census reductio carries the whole argument.
###
### OUTCOME C -- excess well above 100
###   Np is genuinely in the low hundreds and COLONY failed almost completely.
###   Go straight to the error-rate sensitivity runs; the assumed 1-2% genotyping
###   error is the prime suspect.
###
### In every case, check the region and cohort-gap tables in Part 6. If the
### excess dyads pile up within single sampling events, the random-sampling
### model in Part 4 understates expected co-detection and the Np read-off shifts.
################################################################################
beep(sound = 8)
