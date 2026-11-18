# =============================================================================
# 11_figures.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Manuscript Figures 1-5, numbered as the text cites them, in one folder
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-18
#
# Description:
#   The analysis scripts each draw the figures for their own step, under their
#   own names. This script collects the five figures the manuscript cites into
#   results/figures/ as Figure1.png ... Figure5.png, so there is one folder to
#   share with coauthors and, later, to upload to the journal. Three are drawn
#   here and two are copied:
#
#     Figure 1  Collection sites by river km and year     drawn here (02)
#     Figure 2  Pairwise Fst among cohorts + PCA          drawn here (03)
#     Figure 3  Locus sweep                                copied from 09
#     Figure 4  LD and sibship Nb against sample size      copied from 10
#     Figure 5  Length by date of capture and age class    drawn here (02)
#
#   Figures 3 and 4 are copied rather than redrawn so each figure has one
#   piece of code behind it. If 09 or 10 has not produced its figure yet, that
#   figure is reported as missing and the others are still written.
#
#   Figure 1 is a schematic: sites placed in river order, labelled with river
#   km, with the number of fish collected in each year. It needs nothing beyond
#   site_rkm in 00_config.R. A geographic map (with a regional inset) needs
#   site coordinates, which are not in the metadata yet; it would replace this
#   panel or sit beside it.
#
#   Figures 1 and 5 use every fish collected (the 371 in Table 1), not only
#   the 321 analysed, because they describe the sampling. The negative-control
#   well renamed in 00_config.R (2024999) is not a fish and is left out.
#
# Inputs:   results/02_metadata/sample_metadata.csv                   (02)
#           results/03_popgen/cohort/pairwise_fst.csv                  (03)
#           results/03_popgen/cohort/pca_scores.csv                    (03)
#           results/03_popgen/pca_variance.csv                         (03)
#           results/09_colony/Figure3_locus_sweep.png                  (09, optional)
#           results/10_sampling_sensitivity/Figure4_sampling_sensitivity.png
#                                                                      (10, optional)
# Outputs:  results/figures/Figure1.png ... Figure5.png
#           results/figures/figure_list.csv   file, source and caption of each
#
# Usage:    Rscript scripts/11_figures.R
# =============================================================================

source("scripts/00_config.R")

need_pkgs(c("dplyr", "readr", "tidyr", "ggplot2", "gridExtra"))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})

log_start("11_figures")

output_dir <- project_path("results/figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fst_file      <- file.path(paths$dir_popgen, "cohort", "pairwise_fst.csv")
pca_file      <- file.path(paths$dir_popgen, "cohort", "pca_scores.csv")
variance_file <- file.path(paths$dir_popgen, "pca_variance.csv")
figure3_file  <- project_path("results/09_colony/Figure3_locus_sweep.png")
figure4_file  <- project_path("results/10_sampling_sensitivity/Figure4_sampling_sensitivity.png")

for (f in c(paths$sample_meta, fst_file, pca_file, variance_file)) {
  need(file.exists(f), "Input not found: ", f, "\n  Run 02 and 03 first (scripts/run_all.R).")
}

# One palette for birth cohort in every figure drawn here (Okabe-Ito, which
# stays distinguishable in greyscale and for the common colour-vision types).
cohort_colours <- c(`2021` = "#E69F00", `2022` = "#56B4E9",
                    `2023` = "#009E73", `2024` = "#CC79A7")
age_colours    <- c(`0` = "#0072B2", `1` = "#D55E00")

theme_figure <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black"),
        plot.tag = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "grey92"))

site_label <- function(site) gsub("_", " ", site)


# =============================================================================
# Fish collected (Figures 1 and 5)
# =============================================================================
sample_metadata <- read_csv(paths$sample_meta, col_types = sample_id_col_types,
                            show_col_types = FALSE)
not_fish <- qc_exclusions$sample[grepl("not a fish", qc_exclusions$reason)]
collected <- sample_metadata %>%
  filter(!(sample %in% not_fish), !is.na(collection_year)) %>%
  mutate(collection_year = as.integer(collection_year))
msg(sprintf("  fish collected: %d (Table 1 has 371)", nrow(collected)))
if (nrow(collected) != 371)
  msg("  !! count differs from Table 1; check the workbook before sharing Figures 1 and 5")


# =============================================================================
# Figure 1. Collection sites by river km and year
# =============================================================================
site_counts <- collected %>%
  filter(!is.na(region)) %>%
  count(region, collection_year, name = "n")
unplaced <- setdiff(unique(site_counts$region), names(site_rkm))
need(!length(unplaced), "Sites without a river km in site_rkm (00_config.R): ",
     paste(unplaced, collapse = ", "))

# Sites in river order, top to bottom. Sites share river km (Lloyd and
# Poughkeepsie are both 121), so position is rank, and the km is in the label.
sites_present <- region_downstream_order[region_downstream_order %in% site_counts$region]
site_position <- setNames(seq_along(sites_present), sites_present)
site_counts <- site_counts %>%
  mutate(y = site_position[region],
         year = collection_year)
# The Federal Dam at Troy sits between the sites above it and Troy
above_dam <- intersect(c("Lock_C1", "Mohawk_River"), sites_present)
dam_y <- if (length(above_dam)) max(site_position[above_dam]) + 0.5 else NA

figure1 <- ggplot(site_counts, aes(year, y)) +
  { if (!is.na(dam_y)) geom_hline(yintercept = dam_y, linetype = "dashed", colour = "grey40") } +
  { if (!is.na(dam_y)) annotate("text", x = min(site_counts$year) - 0.45, y = dam_y - 0.15, label = "Federal Dam, Troy",
                                hjust = 0, vjust = 0, size = 3, colour = "grey30") } +
  geom_point(aes(size = n), shape = 21, fill = "#56B4E9", colour = "black", alpha = 0.8) +
  geom_text(aes(label = n), size = 2.8, nudge_x = 0.28, hjust = 0) +
  scale_y_reverse(breaks = site_position,
                  labels = sprintf("%s (%d)", site_label(names(site_position)),
                                   site_rkm[names(site_position)]),
                  expand = expansion(add = 0.6)) +
  scale_x_continuous(breaks = sort(unique(site_counts$year)),
                     expand = expansion(add = c(0.5, 0.6))) +
  scale_size_area(max_size = 12, name = "Fish collected", breaks = c(1, 10, 50, 150)) +
  labs(x = "Year of collection", y = "Site (river km), upstream to downstream") +
  theme_figure + theme(panel.grid.major.y = element_line(colour = "grey92"))
ggsave(file.path(output_dir, "Figure1.png"), figure1, width = 7, height = 6, dpi = 300)
msg(sprintf("  Figure 1: %d sites, %d fish", length(sites_present), sum(site_counts$n)))


# =============================================================================
# Figure 2. Pairwise Fst among cohorts (A) and PCA of all fish (B)
# =============================================================================
fst_pairs <- read_csv(fst_file, col_types = cols(group1 = "c", group2 = "c"),
                      show_col_types = FALSE)
cohort_levels <- sort(unique(c(fst_pairs$group1, fst_pairs$group2)))
n_pairs <- nrow(fst_pairs)
alpha_corrected <- 0.05 / n_pairs
fst_pairs <- fst_pairs %>%
  mutate(group1 = factor(group1, levels = cohort_levels),
         group2 = factor(group2, levels = rev(cohort_levels)),
         significant = if ("p_perm" %in% names(.)) p_perm < alpha_corrected else FALSE,
         label = paste0(sprintf("%.4f", Fst), ifelse(significant, "*", "")))
if (!"p_perm" %in% names(fst_pairs))
  msg("  !! no permutation p-values in ", fst_file, "; Figure 2A has no significance marks")

fst_limit <- max(abs(fst_pairs$Fst))
panel_fst <- ggplot(fst_pairs, aes(group1, group2, fill = Fst)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = label), size = 3.6) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0,
                       limits = c(-fst_limit, fst_limit), name = expression(F[ST])) +
  labs(x = "Birth cohort", y = "Birth cohort", tag = "A") +
  coord_fixed() +
  theme_figure + theme(panel.grid = element_blank())

pca_scores <- read_csv(pca_file, col_types = sample_id_col_types, show_col_types = FALSE) %>%
  filter(group != "unassigned") %>%
  mutate(group = factor(group, levels = cohort_levels))
pca_variance <- read_csv(variance_file, show_col_types = FALSE)$pct_variance
panel_pca <- ggplot(pca_scores, aes(Axis1, Axis2, colour = group)) +
  geom_point(size = 1.8, alpha = 0.7) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.6, show.legend = FALSE) +
  scale_colour_manual(values = cohort_colours, name = "Birth cohort") +
  labs(x = sprintf("PC1 (%.1f%%)", pca_variance[1]),
       y = sprintf("PC2 (%.1f%%)", pca_variance[2]), tag = "B") +
  theme_figure + theme(legend.position = "bottom")

figure2 <- gridExtra::arrangeGrob(panel_fst, panel_pca, ncol = 2, widths = c(1, 1.1))
ggsave(file.path(output_dir, "Figure2.png"), figure2, width = 11, height = 5, dpi = 300)
msg(sprintf("  Figure 2: %d cohort pairs (* = permutation p < %.4f), %d fish in the PCA",
            n_pairs, alpha_corrected, nrow(pca_scores)))


# =============================================================================
# Figures 3 and 4. Copied from 09 and 10
# =============================================================================
copy_figure <- function(from, to, step) {
  if (!file.exists(from)) {
    msg("  !! ", basename(to), " missing: ", step, " has not written ", basename(from), " yet")
    if (file.exists(to)) file.remove(to)   # never leave a stale copy behind
    return(FALSE)
  }
  file.copy(from, to, overwrite = TRUE)
  msg("  ", basename(to), ": copied from ", sub(project_path(""), "", from, fixed = TRUE))
  TRUE
}
invisible(copy_figure(figure3_file, file.path(output_dir, "Figure3.png"), "09_colony_summary.R"))
invisible(copy_figure(figure4_file, file.path(output_dir, "Figure4.png"), "10_sampling_sensitivity.R"))


# =============================================================================
# Figure 5. Length by date of capture and age class
# =============================================================================
# Plotted against day of the year so the four years share one axis. Open
# symbols are the age-ambiguous fish removed from the strict cohorts.
length_data <- collected %>%
  filter(!is.na(length_mm), !is.na(age), !is.na(collection_date)) %>%
  mutate(collection_date = as.Date(collection_date),
         day = as.Date(format(collection_date, "2000-%m-%d")),   # common year
         age = factor(age),
         ambiguous = ifelse(age_ambiguous %in% TRUE, "Age ambiguous", "Age clear"))
n_dropped <- nrow(collected) - nrow(length_data)
if (n_dropped) msg("  Figure 5: ", n_dropped, " fish lack length, age or date and are not plotted")

figure5 <- ggplot(length_data, aes(day, length_mm, colour = age, shape = ambiguous)) +
  geom_point(size = 1.8, alpha = 0.8, stroke = 0.7) +
  facet_wrap(~ collection_year, nrow = 1) +
  scale_colour_manual(values = age_colours, name = "Otolith age") +
  scale_shape_manual(values = c(`Age clear` = 16, `Age ambiguous` = 1), name = NULL) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  labs(x = "Date of capture", y = "Total length (mm)") +
  theme_figure + theme(legend.position = "bottom")
ggsave(file.path(output_dir, "Figure5.png"), figure5, width = 10, height = 4, dpi = 300)
msg(sprintf("  Figure 5: %d fish, %d-%d mm, %d age-ambiguous",
            nrow(length_data), min(length_data$length_mm), max(length_data$length_mm),
            sum(length_data$ambiguous == "Age ambiguous")))


# =============================================================================
# Figure list
# =============================================================================
figure_list <- tibble(
  figure = paste("Figure", 1:5),
  file   = paste0("Figure", 1:5, ".png"),
  source = c("11_figures.R (02 metadata)", "11_figures.R (03 popgen)",
             "09_colony_summary.R", "10_sampling_sensitivity.R",
             "11_figures.R (02 metadata)"),
  caption = c(
    "Round goby collection sites on the Hudson River and lower Mohawk River, 2021-2024, in river order with river kilometre, and the number of fish collected at each site in each year. The dashed line marks the Federal Dam at Troy.",
    "(A) Pairwise FST between birth cohorts; * marks pairs significant by permutation after Bonferroni correction. (B) Principal components analysis of all analysed fish, coloured by birth cohort, with 95% ellipses.",
    "Sibship results against the number of loci: (A) sibship Nb, (B) number of parents (Chao1), (C) full-sibling pairs and (D) half-sibling pairs. Points are replicate runs; lines join the means.",
    "LD and sibship estimates of Nb against the number of fish sampled from the 2023 cohort. Lines join medians across five random draws; open triangles are draws with no finite sibship estimate.",
    "Total length against date of capture for every fish collected, by year, coloured by otolith age. Open symbols are fish whose length overlapped the other age class within 30 days (removed from the strict cohorts)."),
  present = file.exists(file.path(output_dir, paste0("Figure", 1:5, ".png")))
)
write_csv(figure_list, file.path(output_dir, "figure_list.csv"))
msg("\n  figures written to ", output_dir, ": ",
    paste(figure_list$figure[figure_list$present], collapse = ", "))

msg("\n=== 11_figures.R complete ===")
log_end()
