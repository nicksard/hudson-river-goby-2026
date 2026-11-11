# =============================================================================
# lib/ldne_direct.R
#
# LD Ne computed directly from GENEPOP files (Waples & Do 2008 bias correction,
# pair-specific S, S^2 weighting, delete-one-contig jackknife). Shared by
# 08b_ldne_direct.R and 08c_ldne_stratag.R. Expects pcrit and block_size in
# the calling environment (08b sets them from ldne_settings).
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-11 (moved out of 08b unchanged)
# =============================================================================

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
