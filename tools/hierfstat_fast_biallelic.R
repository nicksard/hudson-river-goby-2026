# =============================================================================
# tools/hierfstat_fast_biallelic.R
#
# Proposal for hierfstat: matrix-algebra fast paths for biallelic (SNP) data
#
# Author:   Nick Sard (with Claude, Anthropic)
# Created:  2026-11-11
#
# For a GitHub issue on jgx65/hierfstat. Self-contained: simulates a SNP data
# set in hierfstat format (first column = population, genotypes coded 11 /
# 12 / 22, NA missing), then compares hierfstat's basic.stats and
# pairwise.WCfst with fast versions built from per-population counts.
#
# The idea: for biallelic loci, everything basic.stats and pairwise.WCfst need
# comes from three populations x loci count matrices (genotyped individuals,
# copies of allele 1, heterozygotes). Each is one crossproduct of a one-hot
# population matrix with an individuals x loci indicator matrix, so the cost
# is a few matrix multiplications instead of loops over loci and pairs.
# Permutation tests then cost the same per replicate.
#
# Run:  Rscript tools/hierfstat_fast_biallelic.R
# =============================================================================

library(hierfstat)
set.seed(1)

# --- simulate: 6 populations x 50 individuals, 10,000 SNPs, 5% missing ---------
n_pops <- 6; n_per_pop <- 50; n_loci <- 10000
p_anc  <- runif(n_loci, 0.05, 0.95)
geno <- do.call(rbind, lapply(seq_len(n_pops), function(k) {
  p_k <- pmin(pmax(p_anc + rnorm(n_loci, 0, 0.03), 0.01), 0.99)
  t(replicate(n_per_pop, rbinom(n_loci, 2, p_k)))        # copies of allele 2
}))
geno[sample(length(geno), 0.05 * length(geno))] <- NA
coded <- ifelse(is.na(geno), NA, c(11L, 12L, 22L)[geno + 1])
dat <- data.frame(pop = rep(seq_len(n_pops), each = n_per_pop), coded)

# --- fast versions ------------------------------------------------------------
counts_from_hierfstat <- function(dat) {
  g <- as.matrix(dat[, -1])
  copies1 <- (g %/% 10 == 1) + (g %% 10 == 1)               # copies of allele 1
  typed <- !is.na(g) * 1
  allele <- ifelse(is.na(copies1), 0, copies1)
  het <- (allele == 1) * typed
  pops <- factor(dat[, 1])
  z <- outer(as.integer(pops), seq_len(nlevels(pops)), "==") * 1
  colnames(z) <- levels(pops)
  list(n = crossprod(z, typed), allele = crossprod(z, allele), het = crossprod(z, het))
}

fast_Ho_Hs <- function(counts) {                    # as basic.stats()$Ho and $Hs
  n <- counts$n; p <- counts$allele / (2 * n); Ho <- counts$het / n
  Hs <- n / (n - 1) * (1 - p^2 - (1 - p)^2 - Ho / (2 * n))
  Ho[n < 1] <- NA; Hs[n < 2] <- NA
  list(Ho = t(Ho), Hs = t(Hs))
}

fast_pairwise_WCfst <- function(counts) {          # as pairwise.WCfst()
  n_all <- counts$n; groups <- rownames(n_all)
  out <- matrix(NA_real_, length(groups), length(groups), dimnames = list(groups, groups))
  for (i in seq_along(groups)[-1]) for (j in seq_len(i - 1)) {
    n1 <- n_all[i, ]; n2 <- n_all[j, ]
    ok <- n1 > 0 & n2 > 0 & (n1 + n2) > 2
    n1 <- n1[ok]; n2 <- n2[ok]
    p1 <- counts$allele[i, ok] / (2 * n1); p2 <- counts$allele[j, ok] / (2 * n2)
    h1 <- counts$het[i, ok] / n1;          h2 <- counts$het[j, ok] / n2
    nbar <- (n1 + n2) / 2
    nc   <- 2 * nbar - (n1^2 + n2^2) / (2 * nbar)
    pbar <- (n1 * p1 + n2 * p2) / (2 * nbar)
    s2   <- (n1 * (p1 - pbar)^2 + n2 * (p2 - pbar)^2) / nbar
    hbar <- (n1 * h1 + n2 * h2) / (2 * nbar)
    pq   <- pbar * (1 - pbar)
    a <- nbar / nc * (s2 - (pq - s2 / 2 - hbar / 4) / (nbar - 1))
    b <- nbar / (nbar - 1) * (pq - s2 / 2 - (2 * nbar - 1) / (4 * nbar) * hbar)
    c <- hbar / 2
    out[i, j] <- out[j, i] <- sum(a) / sum(a + b + c)
  }
  out
}

# --- compare --------------------------------------------------------------------
timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; value <- expr
  list(value = value, seconds = proc.time()[["elapsed"]] - t0) }

hier_bs  <- timed(basic.stats(dat, diploid = TRUE))
hier_fst <- timed(as.matrix(pairwise.WCfst(dat)))
fast     <- timed({ cnt <- counts_from_hierfstat(dat)
                    list(hs = fast_Ho_Hs(cnt), fst = fast_pairwise_WCfst(cnt)) })

max_diff <- function(a, b) max(abs(as.matrix(a) - as.matrix(b)), na.rm = TRUE)
cat(sprintf("\n%d populations x %d individuals x %d SNPs, %.0f%% missing\n",
            n_pops, n_per_pop, n_loci, 100 * mean(is.na(geno))))
cat(sprintf("basic.stats:      %7.1f s | fast Ho/Hs: max |difference| Ho %.1e, Hs %.1e\n",
            hier_bs$seconds, max_diff(hier_bs$value$Ho, fast$value$hs$Ho),
            max_diff(hier_bs$value$Hs, fast$value$hs$Hs)))
cat(sprintf("pairwise.WCfst:   %7.1f s | fast WC Fst: max |difference| %.1e\n",
            hier_fst$seconds, max_diff(hier_fst$value, fast$value$fst)))
cat(sprintf("fast (both):      %7.1f s\n", fast$seconds))

# one permutation replicate of all pairwise Fst, fast path
perm <- timed({ d2 <- dat; d2[, 1] <- sample(d2[, 1]); fast_pairwise_WCfst(counts_from_hierfstat(d2)) })
cat(sprintf("one permutation replicate (all %d pairs, fast): %.2f s\n",
            choose(n_pops, 2), perm$seconds))
cat("\nsessionInfo:\n"); print(sessionInfo()$otherPkgs$hierfstat$Version)
