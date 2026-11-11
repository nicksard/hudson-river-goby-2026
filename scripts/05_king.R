# =============================================================================
# 05_king.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# KING-robust kinship: a Mendelian, allele-frequency-free check on COLONY
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-11
#
# Why this exists:
#   COLONY's sibship likelihoods depend on allele frequencies, genotyping-error
#   rates and its mating-system model. KING-robust kinship (Manichaikul et al.
#   2010) does not use allele frequencies at all: it counts, across loci typed
#   in both fish, shared heterozygous genotypes and opposite homozygotes. So it
#   is an independent test of how many close relatives the sample contains.
#
#     phi_ij = (N_het,het - 2 * N_opposite_hom) / (N_het_i + N_het_j)
#
#   counted over loci genotyped in both fish. Expected values: duplicate /
#   identical 0.5, first degree (full sib, parent-offspring) 0.25, second
#   degree (half sib, grandparent, avuncular) 0.125, third degree 0.0625,
#   unrelated 0. Bins follow Manichaikul et al. (2010): > 0.354, 0.177-0.354,
#   0.0884-0.177, 0.0442-0.0884.
#
#   First-degree pairs are split by IBS0 (share of loci with opposite
#   homozygotes): parent-offspring pairs have IBS0 ~ 0 apart from genotyping
#   error; full sibs have about a quarter of the unrelated rate.
#
#   KING-robust is inflated for individuals with excess heterozygosity (here
#   mean kinship correlates r = 0.71 with Ho and only 0.11 with missingness).
#   Each fish's mean kinship to everyone else is reported, and fish more than
#   2 SD above the mean are flagged; third-degree pairs involving them should
#   not be read as relatives. Their elevated Ho may also deserve the
#   allele-balance check used for 2023018.
#
# Inputs:   results/07_ldne/genepop/All.gen   (the 321-fish analysis set, from 07)
#           results/02_metadata/sample_metadata.csv
#           results/09_colony/primary_sibling_pairs.csv  (optional, for comparison)
# Outputs:  results/05_king/king_pairs_related.csv   (pairs above the 3rd-degree cut)
#           results/05_king/king_individuals.csv     (per-fish missingness, Ho, mean kinship)
#           results/05_king/king_summary.csv         (counts per bin)
#           results/05_king/FigureS_king_kinship.png
# =============================================================================

source("scripts/00_config.R")

need_pkgs(c("dplyr", "readr", "ggplot2"))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2) })

log_start("05_king")

gen_file   <- project_path("results/07_ldne/genepop/All.gen")
output_dir <- project_path("results/05_king")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
need(file.exists(gen_file), "Run scripts/07_ldne_inputs.R first (All.gen missing).")

kin_breaks <- c(-Inf, 0.0442, 0.0884, 0.177, 0.354, Inf)
kin_labels <- c("unrelated", "third degree", "second degree (half-sib)",
                "first degree (full-sib / parent-offspring)", "duplicate")

# --- genotypes ----------------------------------------------------------------
gen_lines <- readLines(gen_file)
body <- gen_lines[-(1:3)]
body <- body[nzchar(trimws(body)) & toupper(trimws(body)) != "POP"]
ids  <- trimws(sub(",.*$", "", body))
geno <- do.call(rbind, strsplit(trimws(sub("^[^,]*,", "", body)), "\\s+"))
typed    <- (geno != "000000") * 1
het      <- (geno == "001002" | geno == "002001") * 1
hom_1    <- (geno == "001001") * 1
hom_2    <- (geno == "002002") * 1
msg(sprintf("  %d fish x %d loci; missing %.2f%%", nrow(geno), ncol(geno),
            100 * (1 - mean(typed))))

# --- kinship (all pairs at once) ----------------------------------------------
het_het    <- tcrossprod(het)
opposite   <- tcrossprod(hom_1, hom_2) + tcrossprod(hom_2, hom_1)
het_shared <- tcrossprod(het, typed)            # het count of i at loci typed in both
both_typed <- tcrossprod(typed)
kinship    <- (het_het - 2 * opposite) / (het_shared + t(het_shared))
ibs0       <- opposite / both_typed
dimnames(kinship) <- dimnames(ibs0) <- list(ids, ids)

upper <- which(upper.tri(kinship), arr.ind = TRUE)
pairs <- tibble(sample1 = ids[upper[, 1]], sample2 = ids[upper[, 2]],
                kinship = kinship[upper], ibs0 = ibs0[upper],
                loci_both = both_typed[upper]) %>%
  mutate(bin = cut(kinship, kin_breaks, labels = kin_labels))

# --- per-fish quality: does a fish inflate kinship with everyone? ---------------
diag(kinship) <- NA
individuals <- tibble(sample = ids,
                      missing = 1 - rowMeans(typed),
                      Ho = rowSums(het) / rowSums(typed),
                      mean_kinship = rowMeans(kinship, na.rm = TRUE)) %>%
  mutate(z_mean_kinship = (mean_kinship - mean(mean_kinship)) / sd(mean_kinship),
         inflated = z_mean_kinship > 2)
inflated_ids <- individuals$sample[individuals$inflated]
msg(sprintf("  fish with inflated mean kinship (z > 2): %d %s", length(inflated_ids),
            if (length(inflated_ids)) paste0("(", paste(inflated_ids, collapse = ", "), ")") else ""))
msg(sprintf("  cor(mean kinship, Ho) = %.2f; cor(mean kinship, missingness) = %.2f",
            cor(individuals$mean_kinship, individuals$Ho),
            cor(individuals$mean_kinship, individuals$missing)))
msg(sprintf("  their median missingness %.1f%% and Ho %.3f vs %.1f%% and %.3f overall",
            100 * median(individuals$missing[individuals$inflated]),
            median(individuals$Ho[individuals$inflated]),
            100 * median(individuals$missing), median(individuals$Ho)))

# Noise in unrelated pairs, from the bulk of the distribution
noise_sd <- sd(pairs$kinship[pairs$kinship < 0.0442])
unrelated_unrelated <- median(pairs$ibs0[pairs$kinship < 0.0442])

summary_table <- pairs %>%
  mutate(involves_inflated = sample1 %in% inflated_ids | sample2 %in% inflated_ids) %>%
  count(bin, involves_inflated, name = "pairs", .drop = FALSE)
write_csv(summary_table, file.path(output_dir, "king_summary.csv"))
msg(sprintf("\n  %d pairs | kinship median %.4f, SD of unrelated pairs %.4f",
            nrow(pairs), median(pairs$kinship), noise_sd))
msg("  pairs per bin (split by whether a flagged fish is involved):")
msg_table(as.data.frame(summary_table), row.names = FALSE)

# --- related pairs, with metadata -------------------------------------------------
related <- pairs %>% filter(kinship >= 0.0442) %>%
  mutate(sd_above_noise = kinship / noise_sd,
         involves_inflated = sample1 %in% inflated_ids | sample2 %in% inflated_ids,
         first_degree_type = ifelse(kinship >= 0.177,
                                    ifelse(ibs0 < 0.25 * unrelated_unrelated / 2,
                                           "parent-offspring?", "full sib?"), NA)) %>%
  arrange(desc(kinship))

if (file.exists(paths$sample_meta)) {
  meta <- read_csv(paths$sample_meta, col_types = sample_id_col_types) %>%
    select(sample, any_of(c("region", "cohort", "collection_year", "length_mm", "age")))
  related <- related %>%
    left_join(meta %>% rename_with(~paste0(.x, "1"), -sample), by = c("sample1" = "sample")) %>%
    left_join(meta %>% rename_with(~paste0(.x, "2"), -sample), by = c("sample2" = "sample"))
}

colony_file <- project_path("results/09_colony/primary_sibling_pairs.csv")
if (file.exists(colony_file)) {
  colony <- read_csv(colony_file, col_types = cols(.default = col_character())) %>%
    transmute(key = paste(pmin(sample1, sample2), pmax(sample1, sample2)), colony = relationship)
  related <- related %>%
    mutate(key = paste(pmin(sample1, sample2), pmax(sample1, sample2))) %>%
    left_join(colony, by = "key") %>% select(-key)
  colony_missing <- setdiff(colony$key, paste(pmin(related$sample1, related$sample2),
                                              pmax(related$sample1, related$sample2)))
  msg(sprintf("\n  COLONY primary pairs: %d; with KING kinship < 0.0442 (not supported): %d",
              nrow(colony), length(colony_missing)))
} else {
  msg("\n  (COLONY primary_sibling_pairs.csv not found yet; rerun 05 after 09)")
}

write_csv(related, file.path(output_dir, "king_pairs_related.csv"))
write_csv(individuals, file.path(output_dir, "king_individuals.csv"))

msg("\n  first- and second-degree pairs:")
close_pairs <- related %>% filter(kinship >= 0.0884)
if (nrow(close_pairs)) {
  msg_table(as.data.frame(close_pairs %>% mutate(across(c(kinship, ibs0, sd_above_noise), ~round(.x, 4))) %>%
                            select(-bin)), row.names = FALSE)
} else msg("  none")

# --- figure -------------------------------------------------------------------------
figure <- ggplot(pairs, aes(kinship)) +
  geom_histogram(binwidth = 0.005, fill = "grey60", colour = "white") +
  geom_vline(xintercept = c(0.0442, 0.0884, 0.177, 0.354), linetype = 2, colour = "grey30") +
  annotate("text", x = c(0.0663, 0.1327, 0.2655), y = Inf, vjust = 1.5, size = 3,
           label = c("3rd", "2nd", "1st")) +
  scale_y_log10() +
  labs(x = "KING-robust kinship", y = "Pairs (log scale)",
       title = sprintf("KING-robust kinship, %d fish (%d pairs)", length(ids), nrow(pairs))) +
  theme_bw()
ggsave(file.path(output_dir, "FigureS_king_kinship.png"), figure, width = 7, height = 4.5, dpi = 300)

msg("\n  wrote ", output_dir)
msg("\n=== 05_king.R complete ===")
log_end()
