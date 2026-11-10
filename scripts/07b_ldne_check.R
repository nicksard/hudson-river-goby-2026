# =============================================================================
# 07b_ldne_check.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Independent check on NeEstimator's LD Ne: which of its numbers is right?
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
#
# Why this exists:
#   NeEstimator 2.1 prints a point estimate that its own CIs don't contain,
#   and the near-complete run (256 fish x 8,775 loci, 0.71% missing) still
#   shows the gap (306.1 vs CI ~111), so missing data is NOT the explanation.
#   The printed CI is exactly Waples & Do (2008) applied to the printed r^2
#   minus the printed "Expected r^2 Sample", and that printed expectation is
#   quantised (1/64 or 1/256) rather than 1/S + 3.19/S^2.
#
# What it does:
#   1. Reads results/07_ldne/genepop/NearComplete_All.gen.
#   2. Draws random locus subsets (sizes below; fixed seed), keeps loci with
#      MAF >= pcrit in the subset fish.
#   3. Computes LD Ne directly: r^2 per locus pair = squared Pearson
#      correlation of allele counts (a close approximation to the Burrows
#      composite r^2 NeEstimator uses), pair-specific S, E(r^2) = 1/S + 3.19/S^2,
#      r^2' = S^2-weighted mean of (r^2 - E), then the S >= 30 formula.
#   4. Writes each subset as a GENEPOP file + NeEstimator info file + .bat,
#      so NeEstimator can be run on exactly the same loci.
#
# Reading the result:
#   Our value ~ NeEstimator's point  -> the point is right; its CIs are wrong.
#   Our value ~ NeEstimator's CI     -> the point is the odd one out.
#   If the NeEstimator gap grows with locus count, it is numerical, not biology.
#
# Outputs (results/07_ldne/check/):
#   ld_check_direct.csv; Check_<n>loci.gen, _info.txt, run_check_<n>.bat
# =============================================================================

source("scripts/00_config.R")
suppressPackageStartupMessages({ library(readr); library(dplyr) })

check_dir  <- project_path("results/07_ldne/check")
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
gen_file   <- project_path("results/07_ldne/genepop/NearComplete_All.gen")
subset_sizes <- c(1000, 4000)
pcrit      <- ldne_settings$pcrit

msg("07b: independent LD Ne check on ", basename(gen_file))

# --- read GENEPOP (layout written by 07: title, comma locus line, Pop, rows) --
gen_lines  <- readLines(gen_file)
locus_names <- trimws(strsplit(gen_lines[2], ",")[[1]])
body       <- gen_lines[-(1:3)]
body       <- body[nzchar(trimws(body))]
ids        <- trimws(sub(",.*$", "", body))
geno_txt   <- strsplit(trimws(sub("^[^,]*,", "", body)), "\\s+")
need(all(lengths(geno_txt) == length(locus_names)), "Genotype rows do not match the locus line.")
geno_code  <- do.call(rbind, geno_txt)                       # fish x loci, "001002"
dosage     <- matrix(c("001001" = 0, "001002" = 1, "002001" = 1, "002002" = 2,
                       "000000" = NA)[geno_code], nrow = nrow(geno_code))
msg(sprintf("  %d fish x %d loci; missing %.2f%%", nrow(dosage), ncol(dosage),
            100 * mean(is.na(dosage))))

ld_ne_direct <- function(x) {
  p    <- colMeans(x, na.rm = TRUE) / 2
  keep <- pmin(p, 1 - p) >= pcrit
  x    <- x[, keep, drop = FALSE]
  obs  <- !is.na(x) * 1
  S    <- crossprod(obs)                                      # pair-specific S
  r    <- suppressWarnings(cor(x, use = "pairwise.complete.obs"))
  up   <- upper.tri(r) & S >= 30 & is.finite(r)
  r2   <- r[up]^2; s <- S[up]
  E    <- 1 / s + 3.19 / s^2
  w    <- s^2
  r2p  <- sum(w * (r2 - E)) / sum(w)
  ne   <- if (r2p > 0) (1/3 + sqrt(1/9 - 2.76 * r2p)) / (2 * r2p) else Inf
  data.frame(loci_used = sum(keep), pairs = length(r2),
             S_harmonic = 1 / mean(1 / s),
             r2_weighted = sum(w * r2) / sum(w),
             E_r2_weighted = sum(w * E) / sum(w),
             r2_prime = r2p, Ne_direct = ne)
}

results <- list()
for (n_loci in subset_sizes) {
  set.seed(random_seed + n_loci)
  cols <- sort(sample(ncol(dosage), min(n_loci, ncol(dosage))))
  res  <- ld_ne_direct(dosage[, cols])
  results[[length(results) + 1]] <- cbind(subset = paste0("Check_", n_loci, "loci"),
                                          loci_drawn = length(cols), res)
  msg(sprintf("  %5d loci drawn, %5d used: r2 = %.6f, E = %.6f, r2' = %.6f -> Ne = %.1f",
              length(cols), res$loci_used, res$r2_weighted, res$E_r2_weighted,
              res$r2_prime, res$Ne_direct))

  # same loci for NeEstimator
  stem <- paste0("Check_", n_loci, "loci")
  con  <- file(file.path(check_dir, paste0(stem, ".gen")), "wb")
  writeLines(sprintf("Hudson River round goby | NearComplete random %d loci | 07b check", length(cols)), con)
  writeLines(paste(locus_names[cols], collapse = ", "), con)
  writeLines("Pop", con)
  for (i in seq_along(ids))
    writeLines(paste0(ids[i], " , ", paste(geno_code[i, cols], collapse = " ")), con)
  close(con)
  win <- function(p) paste0(gsub("/", "\\\\", normalizePath(p, winslash = "/", mustWork = FALSE)), "\\")
  writeLines(c("1", win(check_dir), paste0(stem, ".gen"), "2", win(check_dir),
               paste0(stem, "_Ne.txt"), "1", as.character(pcrit),
               as.character(ldne_settings$mating)),
             file.path(check_dir, paste0(stem, "_info.txt")))
  writeLines(c("@echo off", "cd /d \"%~dp0\"",
               sprintf("\"%s\" i:%s", gsub("/", "\\\\", ldne_settings$ne_exe),
                       paste0(stem, "_info.txt"))),
             file.path(check_dir, paste0("run_check_", n_loci, ".bat")))
}

out <- bind_rows(results)
write_csv(out, file.path(check_dir, "ld_check_direct.csv"))
msg_table(out)
msg("  Now run run_check_1000.bat and run_check_4000.bat in ", check_dir,
    " and compare NeEstimator's point and CI with Ne_direct.")
