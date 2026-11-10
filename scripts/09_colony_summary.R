# =============================================================================
# 09_colony_summary.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Parse COLONY outputs from the 06 run design: sibship Nb, sibling dyads,
# parent-count Np estimators, Figure 3 (locus + settings sweeps), Figure 4
# (spatial pedigree of the primary run)
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
#
# Description:
#   Adapted from the September cluster scripts (scripts/archive/
#   cluster_rgsibship/: 02.2_colony_data_production.R for the parsing and dyad
#   arithmetic, 03.5_unified.figures_png_pub.R for Figure 3,
#   06.2_all_loci_consensus_spatial_family_plot_png_pub.R for Figure 4).
#   Changes from those scripts:
#
#   * Run identity comes from results/06_colony/run_plan.csv (run name ->
#     loci, replicate, settings), not analysis_plan_key.csv.
#   * Nb is the random-mating block of the .Ne file (as in 02.2). Note
#     02.1/v8 took the non-random block and its CI regex returned "95" from
#     "CI95(L)" for every run; colony_final_metrics_v8.csv is not usable.
#   * COLONY's 2147483647 (its "infinite") is read as Inf.
#   * Np: Chao1 (bias-corrected, Chao 2005) AND iChao1 (Chiu et al. 2014).
#     The September Figure 3 panel B was labelled "iChao" but plotted Chao1.
#     iChao1 equals Chao1 whenever f3 = 0. Cross-check one run against
#     SpadeR::ChaoSpecies before quoting either in the text.
#   * Sibling pairs keep site, cohort and length from sample_metadata.csv,
#     and sites are ordered by river km from 00_config.R (site_rkm).
#   * Full-sib check: September found exactly one FS pair in every run with
#     >= 5,000 loci, almost certainly the 2024221/2023016 duplicate, which is
#     no longer in the analysis set. The log reports the FS pairs found.
#
# Inputs:   results/06_colony/run_plan.csv
#           results/06_colony/outputs/<run_name>/   one folder per run, as the
#               cluster script leaves them in colony_final_outputs/ (download
#               the whole folders; any run not yet finished is skipped)
#           results/02_metadata/sample_metadata.csv
# Outputs:  results/09_colony/colony_run_metrics.csv
#           results/09_colony/locus_sweep_summary.csv
#           results/09_colony/settings_sweep.csv
#           results/09_colony/primary_sibling_pairs.csv
#           results/09_colony/Figure3_locus_sweep.png
#           results/09_colony/FigureS_settings_sweep.png
#           results/09_colony/Figure4_spatial_pedigree.png
#
# Usage:    Rscript scripts/09_colony_summary.R   (rerun as more runs finish)
# =============================================================================

source("scripts/00_config.R")

need_pkgs(c("dplyr", "readr", "tidyr", "ggplot2", "gridExtra"))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})

log_start("09_colony_summary")

colony_dir  <- project_path("results", "06_colony")
outputs_dir <- file.path(colony_dir, "outputs")
output_dir  <- project_path("results", "09_colony")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plan_file <- file.path(colony_dir, "run_plan.csv")
need(file.exists(plan_file), "Run scripts/06_colony_inputs.R first (run_plan.csv missing).")
need(dir.exists(outputs_dir),
     "Download the finished run folders from the cluster's colony_final_outputs/ into ", outputs_dir)

run_plan <- read_csv(plan_file, show_col_types = FALSE)
sample_metadata <- read_csv(paths$sample_meta, col_types = sample_id_col_types)


# -----------------------------------------------------------------------------
# Parsers
# -----------------------------------------------------------------------------
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


# -----------------------------------------------------------------------------
# Parse every finished run
# -----------------------------------------------------------------------------
metrics <- list()
primary_pairs <- NULL

for (i in seq_len(nrow(run_plan))) {
  run <- run_plan[i, ]
  run_dir <- file.path(outputs_dir, run$run_name)
  config_file <- find_output(run_dir, "BestConfig")
  if (is.na(config_file)) next

  config  <- read_best_config(config_file)
  pairs   <- sibling_pairs(config)
  parents <- table(c(config$father, config$mother))
  ne      <- parse_ne(find_output(run_dir, "Ne"))

  metrics[[length(metrics) + 1]] <- bind_cols(
    run,
    tibble(n_offspring = nrow(config),
           n_fs = sum(pairs$relationship == "full-sib"),
           n_hs = sum(pairs$relationship == "half-sib"),
           n_fs_families = length(unique(paste(config$father, config$mother))),
           fs_pairs = paste(apply(pairs[pairs$relationship == "full-sib", 1:2, drop = FALSE], 1,
                                  paste, collapse = "-"), collapse = "; ")),
    as_tibble(as.list(np_estimates(as.numeric(parents)))),
    as_tibble(as.list(ne))
  )

  if (isTRUE(run$is_primary)) primary_pairs <- pairs
}

need(length(metrics) > 0, "No finished runs found under ", outputs_dir)
run_metrics <- bind_rows(metrics)
write_csv(run_metrics, file.path(output_dir, "colony_run_metrics.csv"))
msg(sprintf("  parsed %d of %d planned runs", nrow(run_metrics), nrow(run_plan)))

# Degenerate runs (everyone in one family), as seen with MAF-ranked subsets
collapsed <- run_metrics %>% filter(n_fs_families == 1)
if (nrow(collapsed)) {
  msg("  !! single-family (collapsed) runs: ", paste(collapsed$run_name, collapse = ", "))
}

msg("\n  per-run summary:")
msg_table(as.data.frame(run_metrics %>%
  transmute(run = sub("^(c[0-9]+)_.*", "\\1", run_name), loci = num_loci, rep = replicate,
            RL = run_length, P = precision, ALF = update_alfs, n_fs, n_hs, S_obs, f1, f2,
            Chao1 = round(Np_Chao1), iChao1 = round(Np_iChao1),
            Nb = Nb, Nb_L95, Nb_U95)), row.names = FALSE)

fs_any <- run_metrics %>% filter(n_fs > 0, num_loci >= 5000)
msg("\n  full-sib pairs in runs with >= 5,000 loci: ",
    if (nrow(fs_any)) paste(unique(unlist(strsplit(fs_any$fs_pairs, "; "))), collapse = ", ")
    else "none (the September single FS pair was the 2024221/2023016 duplicate)")


# -----------------------------------------------------------------------------
# Locus sweep (Figure 3) and settings sweep (supplement)
# -----------------------------------------------------------------------------
locus_runs <- run_metrics %>% filter(part == "locus_sweep" | is_primary)
locus_summary <- locus_runs %>%
  group_by(num_loci) %>%
  summarise(n_runs = n(),
            across(c(Nb, Np_Chao1, Np_iChao1, n_fs, n_hs),
                   list(mean = ~mean(.x[is.finite(.x)]), sd = ~sd(.x[is.finite(.x)])),
                   .names = "{.col}_{.fn}"),
            .groups = "drop")
write_csv(locus_summary, file.path(output_dir, "locus_sweep_summary.csv"))
msg("\n  locus sweep (mean across replicates):")
msg_table(as.data.frame(locus_summary %>% mutate(across(-c(num_loci, n_runs), ~signif(.x, 4)))),
          row.names = FALSE)

settings_runs <- run_metrics %>% filter(part == "settings_sweep" |
                                          (num_loci %in% c(colony_settings$design$settings_loci_small,
                                                           max(run_plan$num_loci)) & replicate == 1))
write_csv(settings_runs, file.path(output_dir, "settings_sweep.csv"))

panel <- function(y, label, tag) {
  ggplot(locus_runs, aes(num_loci, .data[[y]])) +
    geom_point(size = 2, alpha = 0.7, colour = "#377EB8") +
    stat_summary(fun = mean, geom = "line", colour = "black", linewidth = 0.8) +
    scale_x_log10(breaks = sort(unique(locus_runs$num_loci))) +
    labs(x = "Number of loci (log10)", y = label, tag = tag) +
    theme_bw() +
    theme(plot.tag = element_text(face = "bold", size = 16),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 11, colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
}
figure3 <- gridExtra::arrangeGrob(
  panel("Nb", "Sibship Nb (random mating)", "A"),
  panel("Np_iChao1", "Np (iChao1)", "B"),
  panel("n_fs", "Full-sib dyads", "C"),
  panel("n_hs", "Half-sib dyads", "D"),
  ncol = 2)
ggsave(file.path(output_dir, "Figure3_locus_sweep.png"), figure3, width = 11, height = 9, dpi = 300)

if (nrow(settings_runs)) {
  figure_settings <- ggplot(settings_runs,
                            aes(factor(run_length), Nb, colour = factor(update_alfs),
                                shape = factor(precision))) +
    geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
    facet_wrap(~ paste(num_loci, "loci"), scales = "free_y") +
    scale_colour_manual(values = c(`0` = "#E41A1C", `1` = "#377EB8"),
                        name = "Update allele freqs", labels = c(`0` = "No", `1` = "Yes")) +
    scale_shape_discrete(name = "Precision (0-3)") +
    labs(x = "Run length (1 short - 4 very long)", y = "Sibship Nb") +
    theme_bw()
  ggsave(file.path(output_dir, "FigureS_settings_sweep.png"), figure_settings,
         width = 10, height = 5, dpi = 300)
}


# -----------------------------------------------------------------------------
# Primary run: sibling pairs with site / cohort, and Figure 4
# -----------------------------------------------------------------------------
if (is.null(primary_pairs)) {
  msg("\n  primary run not finished yet -- sibling-pair table and Figure 4 skipped")
} else {
  site_info <- sample_metadata %>% select(sample, region, cohort, collection_year, length_mm)
  primary_pairs <- primary_pairs %>%
    left_join(site_info %>% rename_with(~paste0(.x, "1"), -sample), by = c("sample1" = "sample")) %>%
    left_join(site_info %>% rename_with(~paste0(.x, "2"), -sample), by = c("sample2" = "sample")) %>%
    mutate(rkm1 = site_rkm[region1], rkm2 = site_rkm[region2],
           dist_rkm = abs(rkm1 - rkm2),
           spatial = ifelse(region1 == region2, "same site", "different sites"),
           cohorts = ifelse(cohort1 == cohort2, "same cohort", "different cohorts"))
  write_csv(primary_pairs, file.path(output_dir, "primary_sibling_pairs.csv"))

  msg("\n  primary run sibling pairs:")
  msg_table(as.data.frame(primary_pairs %>% count(relationship, spatial, cohorts)), row.names = FALSE)
  dispersed <- primary_pairs %>% filter(spatial == "different sites")
  if (nrow(dispersed)) {
    msg("\n  dispersed pairs (different sites):")
    msg_table(as.data.frame(dispersed %>%
      select(sample1, sample2, relationship, region1, region2, dist_rkm, cohort1, cohort2,
             length_mm1, length_mm2)), row.names = FALSE)
  }

  # Figure 4 (06.2 layout): one row per multi-offspring family, sites N -> S
  primary_config <- read_best_config(find_output(
    file.path(outputs_dir, run_plan$run_name[run_plan$is_primary][1]), "BestConfig"))
  related <- unique(c(primary_pairs$sample1, primary_pairs$sample2))
  family_edges <- primary_pairs %>% select(sample1, sample2)
  # connected components of the sibling graph = families for plotting
  family_id <- setNames(seq_along(related), related)
  repeat {
    changed <- FALSE
    for (k in seq_len(nrow(family_edges))) {
      a <- family_edges$sample1[k]; b <- family_edges$sample2[k]
      m <- min(family_id[a], family_id[b])
      if (family_id[a] != m || family_id[b] != m) {
        family_id[family_id == family_id[a] | family_id == family_id[b]] <- m
        changed <- TRUE
      }
    }
    if (!changed) break
  }
  families <- tibble(sample = names(family_id), family = as.integer(factor(family_id))) %>%
    left_join(site_info, by = "sample") %>%
    group_by(family) %>%
    mutate(family_size = n(),
           spatial_status = ifelse(n_distinct(region) > 1, "Dispersed (multiple sites)",
                                   "Sympatric (same site)")) %>%
    ungroup() %>%
    mutate(site = factor(region, levels = region_downstream_order[region_downstream_order %in% region]))
  family_order <- families %>% distinct(family, spatial_status, family_size) %>%
    arrange(spatial_status, desc(family_size)) %>% pull(family)
  families$row <- match(families$family, family_order)
  write_csv(families, file.path(output_dir, "primary_families.csv"))

  status_colours <- c("Sympatric (same site)" = "#2B83BA", "Dispersed (multiple sites)" = "#D7191C")
  figure4 <- ggplot(families, aes(site, row, group = row)) +
    geom_line(aes(colour = spatial_status, linewidth = family_size), alpha = 0.75) +
    geom_point(aes(fill = spatial_status, shape = factor(cohort)), size = 3.8,
               colour = "black", stroke = 0.8) +
    scale_shape_manual(name = "Birth cohort", values = c(`2021` = 21, `2022` = 22, `2023` = 23, `2024` = 24)) +
    scale_colour_manual(values = status_colours, guide = "none") +
    scale_fill_manual(values = status_colours, guide = "none") +
    scale_linewidth_continuous(range = c(0.9, 2.8), guide = "none") +
    scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
    guides(shape = guide_legend(override.aes = list(size = 4.5, fill = "grey40"))) +
    labs(x = "Sampling location (north to south)", y = "Reconstructed families") +
    theme_classic() +
    theme(axis.title = element_text(face = "bold", size = 16),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13, colour = "black"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          legend.position = "bottom")
  ggsave(file.path(output_dir, "Figure4_spatial_pedigree.png"), figure4,
         width = 10.5, height = 8.5, dpi = 300)
  msg(sprintf("\n  Figure 4: %d multi-offspring families (%d dispersed)",
              length(family_order), n_distinct(families$family[families$spatial_status != "Sympatric (same site)"])))
}

msg("\n  wrote ", output_dir)
msg("\n=== 09_colony_summary.R complete ===")
log_end()
