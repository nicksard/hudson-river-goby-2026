# =============================================================================
# 08b_ldne_direct.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# LD Ne computed directly, excluding locus pairs on the same contig
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
#
# Why this exists (findings of 2026-11-10):
#   * NeEstimator 2.1 itself is sound. On simulated data (true Ne = 100,
#     1,200-10,000 loci, random or clustered missingness) it returns 94-100
#     with CIs that contain the point. On the goby near-complete file the same
#     exe returns 476.9 (461.2-493.7) when run under Wine on Linux. The Windows
#     runs of the morning of Nov 10 (306.1 with CI ~111, etc.) were corrupted
#     in the run itself, not in the data.
#   * Locus pairs on the same contig are physically linked: r^2 = 0.061 vs
#     0.0044 for pairs on different contigs (near-complete set). Including
#     them roughly halves Ne (503 vs 974). NeEstimator has no option to drop
#     them, so this script does the estimate itself.
#
# Method:
#   r^2 per locus pair = squared Pearson correlation of allele dosages (the
#   composite r^2 for biallelic loci), loci with MAF >= pcrit in the group,
#   pairs with S >= 30 (S < 30 for groups with fewer fish: see below).
#   E(r^2) per pair from Waples & Do (2008): S >= 30: 1/S + 3.19/S^2;
#   S < 30: 0.0018 + 0.907/S + 4.44/S^2. r^2' = S^2-weighted mean of
#   (r^2 - E), then Ne from the matching Waples & Do (2008) formula.
#   CI: delete-one-contig jackknife on r^2' (pairs are not independent, and
#   contigs are the natural blocks), converted to Ne; a lower r^2' limit <= 0
#   gives an infinite upper Ne limit.
#   Computed in blocks of loci, accumulating contig x contig sums, so memory
#   stays bounded for the full 10,943-locus panel.
#
#   Two estimates per group:
#     Ne_all_pairs       comparable to NeEstimator (all pairs)
#     Ne_diff_contig     primary: pairs on different contigs only
#   The Waples et al. (2016) chromosome correction (divisor 0.7847) is reported
#   for the all-pairs estimate only. It corrects the bias from linked pairs on
#   the same chromosome, which the different-contig estimate already removes
#   in part (contigs on one chromosome are still loosely linked). Which to
#   report is a Methods decision: see the log.
#
# Inputs:  results/07_ldne/genepop/<group>.gen (written by 07)
# Outputs: results/07_ldne/direct/ld_ne_direct.csv
# =============================================================================

source("scripts/00_config.R")
suppressPackageStartupMessages({ library(readr); library(dplyr) })

genepop_dir <- project_path("results/07_ldne/genepop")
out_dir     <- project_path("results/07_ldne/direct")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pcrit       <- ldne_settings$pcrit
block_size  <- 1500
waples_divisor <- 0.098 + 0.219 * log(ldne_settings$n_chromosomes)

groups <- sub("\\.gen$", "", list.files(genepop_dir, pattern = "\\.gen$"))
if (exists("ld_direct_groups")) groups <- intersect(groups, ld_direct_groups)
msg("08b: direct LD Ne for ", length(groups), " group(s): ", paste(groups, collapse = ", "))

read_genepop <- function(path) {
  gen_lines <- readLines(path)
  loci <- trimws(strsplit(gen_lines[2], ",")[[1]])
  body <- gen_lines[-(1:3)]
  body <- body[nzchar(trimws(body)) & toupper(trimws(body)) != "POP"]
  geno <- do.call(rbind, strsplit(trimws(sub("^[^,]*,", "", body)), "\\s+"))
  need(ncol(geno) == length(loci), "Genotype rows do not match the locus line in ", basename(path))
  dosage <- matrix(c("001001" = 0, "001002" = 1, "002001" = 1, "002002" = 2,
                     "000000" = NA)[geno], nrow = nrow(geno))
  list(dosage = dosage, loci = loci)
}

ne_from_r2p <- function(r2p, small_sample) {
  if (!is.finite(r2p) || r2p <= 0) return(Inf)
  if (small_sample) (0.308 + sqrt(0.308^2 - 2.08 * r2p)) / (2 * r2p)
  else              (1/3 + sqrt(1/9 - 2.76 * r2p)) / (2 * r2p)
}

ld_direct <- function(dosage, loci) {
  p    <- colMeans(dosage, na.rm = TRUE) / 2
  keep <- !is.na(p) & pmin(p, 1 - p) >= pcrit
  x    <- dosage[, keep, drop = FALSE]
  contig <- factor(sub("_[0-9]+$", "", loci[keep]))
  obs  <- !is.na(x) * 1
  n_fish <- nrow(x)
  small_sample <- n_fish < 30                     # Waples & Do (2008) S < 30 branch
  min_S <- if (small_sample) 2 else 30
  K <- nlevels(contig)
  sums <- list(w = matrix(0, K, K), wr2 = matrix(0, K, K), wE = matrix(0, K, K),
               n = matrix(0, K, K), winvS = matrix(0, K, K))
  blocks <- split(seq_len(ncol(x)), ceiling(seq_len(ncol(x)) / block_size))
  for (a in seq_along(blocks)) for (b in a:length(blocks)) {
    ia <- blocks[[a]]; ib <- blocks[[b]]
    r <- suppressWarnings(cor(x[, ia, drop = FALSE], x[, ib, drop = FALSE],
                              use = "pairwise.complete.obs"))
    S <- crossprod(obs[, ia, drop = FALSE], obs[, ib, drop = FALSE])
    ok <- is.finite(r) & S >= min_S
    if (a == b) ok <- ok & upper.tri(ok)
    E  <- if (small_sample) 0.0018 + 0.907 / S + 4.44 / S^2 else 1 / S + 3.19 / S^2
    w  <- ifelse(ok, S^2, 0)
    add <- function(M) {                            # contig(ia) x contig(ib) sums
      M[!ok] <- 0
      t(rowsum(t(rowsum(M, contig[ia], reorder = FALSE)), contig[ib], reorder = FALSE))
    }
    ca <- unique(as.integer(contig[ia])); cb <- unique(as.integer(contig[ib]))
    for (nm in names(sums)) {
      M <- switch(nm, w = w, wr2 = w * r^2, wE = w * E, n = ok * 1, winvS = w / S)
      sums[[nm]][ca, cb] <- sums[[nm]][ca, cb] + add(M)
    }
  }
  # symmetric totals over unordered pairs: off-diagonal cells counted once
  sym <- lapply(sums, function(C) { T <- C + t(C); diag(T) <- diag(C); T })

  estimate <- function(drop_k = NULL, diff_only) {
    tot <- function(T) {
      s_all <- (sum(T) + sum(diag(T))) / 2           # unordered pairs incl. same-contig
      s <- if (diff_only) s_all - sum(diag(T)) else s_all
      if (!is.null(drop_k)) s <- s - (sum(T[drop_k, ]) - if (diff_only) T[drop_k, drop_k] else 0)
      s
    }
    w <- tot(sym$w)
    c(r2 = tot(sym$wr2) / w, E = tot(sym$wE) / w, r2p = (tot(sym$wr2) - tot(sym$wE)) / w,
      pairs = tot(sym$n), S_mean = w / tot(sym$winvS))
  }
  one <- function(diff_only) {
    full <- estimate(NULL, diff_only)
    jk   <- vapply(seq_len(K), function(k) estimate(k, diff_only)["r2p"], numeric(1))
    se   <- sqrt((K - 1) / K * sum((jk - mean(jk))^2))
    lo_r2p <- full["r2p"] - 1.96 * se; hi_r2p <- full["r2p"] + 1.96 * se
    c(full, Ne = ne_from_r2p(unname(full["r2p"]), small_sample),
      Ne_L95 = ne_from_r2p(unname(hi_r2p), small_sample),
      Ne_U95 = ne_from_r2p(unname(lo_r2p), small_sample))
  }
  list(all = one(FALSE), diff = one(TRUE), n_fish = n_fish, loci_used = ncol(x),
       contigs = K, small_sample = small_sample)
}

rows <- list()
for (g in groups) {
  t0 <- Sys.time()
  dat <- read_genepop(file.path(genepop_dir, paste0(g, ".gen")))
  res <- ld_direct(dat$dosage, dat$loci)
  rows[[g]] <- data.frame(
    group = g, n_fish = res$n_fish, loci_used = res$loci_used, contigs = res$contigs,
    formula = if (res$small_sample) "S<30" else "S>=30",
    pairs_all = res$all[["pairs"]], r2_all = res$all[["r2"]], E_r2 = res$all[["E"]],
    Ne_all_pairs = res$all[["Ne"]], L95_all = res$all[["Ne_L95"]], U95_all = res$all[["Ne_U95"]],
    Ne_all_waples2016 = res$all[["Ne"]] / waples_divisor,
    L95_all_waples2016 = res$all[["Ne_L95"]] / waples_divisor,
    U95_all_waples2016 = res$all[["Ne_U95"]] / waples_divisor,
    pairs_diff = res$diff[["pairs"]], r2_diff = res$diff[["r2"]],
    Ne_diff_contig = res$diff[["Ne"]], L95_diff = res$diff[["Ne_L95"]], U95_diff = res$diff[["Ne_U95"]])
  msg(sprintf("  %-22s n=%3d loci=%5d  all pairs: %7.1f (%s-%s)  different contigs: %7.1f (%s-%s)  [%.1f min]",
              g, res$n_fish, res$loci_used,
              res$all[["Ne"]], format(round(res$all[["Ne_L95"]], 1)), format(round(res$all[["Ne_U95"]], 1)),
              res$diff[["Ne"]], format(round(res$diff[["Ne_L95"]], 1)), format(round(res$diff[["Ne_U95"]], 1)),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
out <- bind_rows(rows)
write_csv(out, file.path(out_dir, "ld_ne_direct.csv"))
msg_table(out %>% select(group, n_fish, loci_used, Ne_all_pairs, L95_all, U95_all,
                         Ne_diff_contig, L95_diff, U95_diff))
msg("  Methods decision: report Ne_diff_contig (linked pairs removed) or",
    " Ne_all_waples2016 (all pairs, Waples et al. 2016 correction), not both corrections.")
