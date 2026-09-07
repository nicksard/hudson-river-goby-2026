################################################################################
### 06_ldne_scaffold_partition.R
###
### Within-scaffold vs between-scaffold vs all-pairs LD-Ne, computed in R.
###
### WHY
###   NeEstimator can restrict LD to pairs on DIFFERENT chromosomes (with the
###   chromosome map), but has no "same chromosome only" option. This script
###   computes all three partitions directly so they can be compared on the
###   same footing. It uses the Waples (2006) / Waples & Do (2010) bias-corrected
###   formula for S >= 30. It is an APPROXIMATION to NeEstimator -- genotypic
###   correlation rather than the Burrows composite, unweighted mean r2 -- so
###   use it for the within/between CONTRAST, and report NeEstimator's numbers
###   as the headline values.
###
### READING THE RESULT
###   If between-scaffold Ne >> within-scaffold Ne, physical linkage is
###   inflating LD and the corrected LDNe is still too low.
###   If they agree, the small Ne is not a linkage artefact.
###   Nick's prior: the Waples 2016 correction moved things ~25%, so expect
###   agreement. This is the check that lets you say so with a number.
################################################################################

library(vcfR); library(tidyverse); library(readxl)
setwd("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/post_RAD analysis")

vcf_file_path <- "Output/populations.filtered.final.vcf.gz"
bio_data_path <- paste0("G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/",
                        "plate_barcodes/Input/GobyRad_tabularPlateWell with Bio Data.xlsx")
output_dir    <- "Output/ldne_scaffold"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

EXCLUDE_SAMPLES <- c("2023018", "2024221")
PCRIT           <- 0.05      # match NeEstimator run
MIN_N           <- 30
COHORTS         <- c("2023", "2024")          # the informative ones
SCHEMES         <- c("All", "Albany_Troy")


## ---------------------------------------------------------------------------
## 1. Cohort assignment (same logic as 02 / genepop v2)
## ---------------------------------------------------------------------------
bio <- read_xlsx(bio_data_path, sheet = 1)
names(bio) <- gsub("\\(|\\)|#|`|'", "", gsub(" ", "_", trimws(names(bio))))
bio <- bio[!grepl("NTC|999", bio$ID), ]
bio$Region     <- trimws(as.character(bio$Region))
bio$Birth_Year <- as.numeric(substr(bio$ID, 1, 4)) - suppressWarnings(as.numeric(bio$Age))
bio$Fish_ID    <- as.character(bio$Fish_ID)
bio <- bio %>% filter(!is.na(Birth_Year), !Fish_ID %in% EXCLUDE_SAMPLES) %>%
  mutate(is_AT = grepl("Albany|Troy", Region))


## ---------------------------------------------------------------------------
## 2. VCF -> numeric matrix, scaffold vector
## ---------------------------------------------------------------------------
vcf  <- read.vcfR(vcf_file_path, verbose = FALSE)
scaf <- vcf@fix[, "CHROM"]
gt   <- extract.gt(vcf, element = "GT")
num  <- matrix(NA_real_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
num[gt %in% c("0/0","0|0")] <- 0
num[gt %in% c("0/1","1/0","0|1","1|0")] <- 1
num[gt %in% c("1/1","1|1")] <- 2


## ---------------------------------------------------------------------------
## 3. LD-Ne from mean r2  (Waples 2006; Waples & Do 2010, S >= 30)
## ---------------------------------------------------------------------------
ldne_from_r2 <- function(mean_r2, S) {
  r2_drift <- mean_r2 - (1/S + 3.19/S^2)
  if (r2_drift <= 0) return(Inf)
  disc <- 1/9 - 2.76 * r2_drift
  if (disc < 0) return(NA_real_)
  (1/3 + sqrt(disc)) / (2 * r2_drift)
}

partition_ldne <- function(sample_ids, label) {
  X <- num[, sample_ids, drop = FALSE]
  S <- ncol(X)

  # Pcrit screen on this sample, as NeEstimator does
  af  <- rowMeans(X, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  ok  <- !is.na(maf) & maf >= PCRIT & rowMeans(is.na(X)) < 0.5
  X   <- X[ok, ]; sc <- scaf[ok]
  L   <- nrow(X)

  # mean-impute missing for the correlation (small effect at <5% missing)
  Xi <- t(apply(X, 1, function(v) { v[is.na(v)] <- mean(v, na.rm = TRUE); v }))

  # genotypic correlation matrix; loci x loci
  R  <- cor(t(Xi))
  r2 <- R^2
  rm(R)

  same <- outer(sc, sc, "==")
  ut   <- upper.tri(r2)

  within  <- r2[ut & same]
  between <- r2[ut & !same]
  allp    <- r2[ut]

  tibble(
    Group     = label, S = S, L = L,
    n_within  = length(within), n_between = length(between),
    r2_within = mean(within), r2_between = mean(between), r2_all = mean(allp),
    Ne_within  = ldne_from_r2(mean(within),  S),
    Ne_between = ldne_from_r2(mean(between), S),
    Ne_all     = ldne_from_r2(mean(allp),    S)
  )
}


## ---------------------------------------------------------------------------
## 4. Run
## ---------------------------------------------------------------------------
res <- list()
for (sch in SCHEMES) for (ch in COHORTS) {
  ids <- bio %>% filter(Birth_Year == as.numeric(ch))
  if (sch == "Albany_Troy") ids <- ids %>% filter(is_AT)
  ids <- intersect(ids$Fish_ID, colnames(num))
  lab <- paste(ch, sch, sep = "_")
  if (length(ids) < MIN_N) { message(lab, ": N=", length(ids), " < ", MIN_N, ", skipped"); next }
  message("running ", lab, " (N=", length(ids), ") ...")
  res[[lab]] <- partition_ldne(ids, lab)
}
res <- bind_rows(res)

write_csv(res, file.path(output_dir, "ldne_scaffold_partition.csv"))
print(res %>% mutate(across(starts_with("r2"), ~signif(.x, 4)),
                     across(starts_with("Ne"), ~round(.x))) %>% as.data.frame())

cat("\nNe_between / Ne_within ratio (>> 1 means linkage is inflating LD):\n")
print(res %>% transmute(Group, ratio = round(Ne_between / Ne_within, 2)) %>% as.data.frame())


## ---------------------------------------------------------------------------
## 5. Figure
## ---------------------------------------------------------------------------
p <- res %>%
  pivot_longer(c(Ne_within, Ne_between, Ne_all), names_to = "partition", values_to = "Ne") %>%
  mutate(partition = factor(sub("Ne_", "", partition), c("within","between","all"))) %>%
  ggplot(aes(x = Group, y = Ne, fill = partition)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = c(within = "firebrick", between = "steelblue4", all = "grey50")) +
  labs(y = "LD-Ne (approximate)", x = NULL, fill = "locus pairs",
       title = "LD-Ne by scaffold partition",
       subtitle = "within-scaffold pairs carry physical linkage; between-scaffold pairs are (approximately) unlinked") +
  theme_bw()
ggsave(file.path(output_dir, "Fig_LDNe_scaffold_partition.png"), p, width = 8, height = 5, dpi = 300)
