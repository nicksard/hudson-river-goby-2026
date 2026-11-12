# =============================================================================
# scripts/lib/colony_parse.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Readers for COLONY output and the sibship/Np arithmetic built on them
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10 (in 09_colony_summary.R)
# Updated:  2026-11-11 (moved here so 09 and 10 share one copy)
#
# Description:
#   find_output()      first file in a run folder with a given extension
#   colony_number()    COLONY numbers; 2147483647 (its "infinite") -> Inf
#   parse_ne()         random-mating Nb and 95% CI from the .Ne file
#   read_best_config() .BestConfig (offspring, father, mother, cluster)
#   np_estimates()     S_obs, f1-f4, bias-corrected Chao1, iChao1
#   sibling_pairs()    every full- and half-sib pair implied by a BestConfig
#
#   Parsing follows the September 02.2 script (compiled_colony_results.csv),
#   which was checked against COLONY's own output.
#
# Used by:  09_colony_summary.R, 10_sampling_sensitivity.R
#
# References:
#   Chao A (2005) Species estimation and applications. In: Encyclopedia of
#     Statistical Sciences, 2nd edn. Wiley, New York, pp 7907-7916.
#   Chiu C-H, Wang Y-T, Walther BA, Chao A (2014) An improved nonparametric
#     lower bound of species richness via a modified Good-Turing frequency
#     formula. Biometrics 70:671-682.
#   Jones OR, Wang J (2010) COLONY: a program for parentage and sibship
#     inference from multilocus genotype data. Molecular Ecology Resources
#     10:551-555.
# =============================================================================

find_output <- function(run_dir, extension) {
  hit <- list.files(run_dir, pattern = paste0("\\.", extension, "$"), full.names = TRUE)
  if (length(hit)) hit[1] else NA_character_
}

colony_number <- function(x) {
  value <- suppressWarnings(as.numeric(trimws(x)))
  ifelse(!is.na(value) & value >= 2147483647, Inf, value)
}

# First Ne / CI95(L) / CI95(U) lines = the random-mating block (02.2 logic)
parse_ne <- function(ne_file) {
  empty <- c(Nb = NA_real_, Nb_L95 = NA_real_, Nb_U95 = NA_real_)
  if (is.na(ne_file)) return(empty)
  lines <- readLines(ne_file, warn = FALSE)
  grab <- function(pattern) {
    hit <- grep(pattern, lines)
    if (!length(hit)) return(NA_real_)
    colony_number(sub("^[^=]*=", "", lines[hit[1]]))
  }
  c(Nb = grab("^\\s*Ne\\s*="), Nb_L95 = grab("^\\s*CI95\\(L\\)\\s*="),
    Nb_U95 = grab("^\\s*CI95\\(U\\)\\s*="))
}

read_best_config <- function(config_file) {
  config <- read.table(config_file, header = TRUE, stringsAsFactors = FALSE,
                       comment.char = "", colClasses = "character")
  names(config)[1:4] <- c("sample", "father", "mother", "cluster")
  config
}

# Chao (2005) bias-corrected Chao1, and iChao1 (Chiu et al. 2014; f4 -> 1 if 0)
np_estimates <- function(parent_counts) {
  n  <- sum(parent_counts)
  f  <- function(k) sum(parent_counts == k)
  f1 <- f(1); f2 <- f(2); f3 <- f(3); f4 <- f(4)
  s_obs <- length(parent_counts)
  chao1 <- s_obs + ((n - 1) / n) * f1 * (f1 - 1) / (2 * (f2 + 1))
  f4_used <- max(f4, 1)
  ichao1 <- chao1 + ((n - 3) / (4 * n)) * (f3 / f4_used) *
    max(f1 - ((n - 3) / (2 * (n - 1))) * f2 * f3 / f4_used, 0)
  c(S_obs = s_obs, f1 = f1, f2 = f2, f3 = f3, f4 = f4, Np_Chao1 = chao1, Np_iChao1 = ichao1)
}

sibling_pairs <- function(config) {
  pairs_for <- function(groups, type) {
    out <- lapply(split(config$sample, groups), function(ids) {
      if (length(ids) < 2) return(NULL)
      combos <- t(utils::combn(sort(ids), 2))
      data.frame(sample1 = combos[, 1], sample2 = combos[, 2], stringsAsFactors = FALSE)
    })
    res <- bind_rows(out)
    if (nrow(res)) res$share <- type
    res
  }
  fs  <- pairs_for(paste(config$father, config$mother), "both")
  dad <- pairs_for(config$father, "father")
  mom <- pairs_for(config$mother, "mother")
  all_pairs <- bind_rows(fs, dad, mom)
  if (!nrow(all_pairs)) {
    return(data.frame(sample1 = character(), sample2 = character(), relationship = character()))
  }
  all_pairs %>%
    distinct(sample1, sample2, share) %>%
    group_by(sample1, sample2) %>%
    summarise(relationship = ifelse(any(share == "both"), "full-sib", "half-sib"), .groups = "drop")
}
