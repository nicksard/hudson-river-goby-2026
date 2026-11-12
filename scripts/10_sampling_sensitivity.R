# =============================================================================
# 10_sampling_sensitivity.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# How LD Nb and sibship Nb respond to sample size (2023 cohort)
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-11
# Updated:  2026-11-11
#
# Description:
#   For the full 2023 cohort, LD Nb is about 500-1,000 and sibship Nb about
#   41,000. Wang (2016) found the two estimators similar when their
#   assumptions hold, with precision governed mostly by sample size. This
#   script asks how each estimate behaves as the sample grows, using the fish
#   subsets written by 06b (25, 50, 75, 100, 125 and all 166 fish; 5 draws
#   each). Both methods see exactly the same fish in each subset.
#
#   What to look for:
#     * An estimate that keeps rising with n instead of levelling off means
#       the sample is still small relative to Nb.
#     * Half-sib pairs should grow roughly as n^2 / Nb; at small n most draws
#       should find none, and sibship Nb is then unbounded.
#     * The 5 runs at n = 166 have identical fish, so their spread is COLONY's
#       run-to-run (seed) variation; the LD values there are identical.
#
#   LD uses all 10,943 loci (the Table 3 estimator: scripts/lib/ldne_direct.R,
#   (S/(S-1))^2 r^2 correction), point estimates only, all pairs (raw and
#   with the Waples et al. 2016 linkage correction) and different-contig
#   pairs. COLONY uses 06b's fixed 8,000-locus subset. COLONY runs that have
#   not been downloaded yet are skipped, so the LD half can run first.
#
# Order of operations:
#   1. Read 06b's run plan and fish lists, and the 2023 cohort Genepop file.
#   2. LD Nb for every subset (cached in ld_by_subset.csv).
#   3. Parse every finished COLONY run: sibling pairs, Nb, Np (Chao1, iChao1).
#   4. Summarise by sample size; write the table and the figure.
#
# Inputs:   results/06b_colony_sampling/run_plan.csv, fish_subsets.csv (06b)
#           results/07_ldne/genepop/Cohort_2023.gen (07)
#           results/06b_colony_sampling/outputs/<run>/ (cluster, optional)
# Outputs:  results/10_sampling_sensitivity/ld_by_subset.csv (cache; redone
#           automatically if lib/ldne_direct.R changes),
#           sampling_by_run.csv, sampling_summary.csv,
#           Figure_S_sampling_sensitivity.png
#
# References:
#   Waples RS, Do C (2008) LDNE: a program for estimating effective
#     population size from data on linkage disequilibrium. Molecular Ecology
#     Resources 8:753-756.
#   Waples RS, Larson WA, Waples RK (2016) Estimating contemporary effective
#     population size in non-model species using linkage disequilibrium across
#     thousands of loci. Heredity 117:233-240.
#   Wang J (2016) A comparison of single-sample estimators of effective
#     population sizes from genetic marker data. Molecular Ecology
#     25:4692-4711.
# =============================================================================

source("scripts/00_config.R")
need_pkgs(c("dplyr", "readr", "tidyr", "ggplot2"))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2) })

log_start("10_sampling_sensitivity")

sampling_dir   <- project_path("results", "06b_colony_sampling")
colony_out_dir <- file.path(sampling_dir, "outputs")
genepop_file   <- project_path("results", "07_ldne", "genepop", "Cohort_2023.gen")
output_dir     <- project_path("results", "10_sampling_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pcrit          <- ldne_settings$pcrit        # used inside ld_direct()
block_size     <- 1500                       # used inside ld_direct()
waples_divisor <- 0.098 + 0.219 * log(ldne_settings$n_chromosomes)
source("scripts/lib/ldne_direct.R",  local = TRUE)   # read_genepop(), ld_direct()
source("scripts/lib/colony_parse.R", local = TRUE)   # parse_ne(), sibling_pairs(), ...

need(file.exists(file.path(sampling_dir, "fish_subsets.csv")),
     "Run scripts/06b_colony_sampling_inputs.R first.")
need(file.exists(genepop_file), "Run scripts/07_ldne_inputs.R first (Cohort_2023.gen missing).")

# ---- 1. Plan, fish lists, genotypes -----------------------------------------
run_plan     <- read_csv(file.path(sampling_dir, "run_plan.csv"), show_col_types = FALSE)
fish_subsets <- read_csv(file.path(sampling_dir, "fish_subsets.csv"),
                         col_types = cols(sample = col_character(), .default = col_guess()))

genepop   <- read_genepop(genepop_file)
gen_lines <- readLines(genepop_file)
fish_ids  <- trimws(sub(",.*$", "", gen_lines[-(1:3)]))
fish_ids  <- fish_ids[nzchar(fish_ids) & toupper(fish_ids) != "POP"]
need(length(fish_ids) == nrow(genepop$dosage), "Fish IDs and genotype rows do not line up.")
missing_fish <- setdiff(fish_subsets$sample, fish_ids)
need(length(missing_fish) == 0,
     "Subset fish not in Cohort_2023.gen: ", paste(head(missing_fish), collapse = ", "))
msg(sprintf("  %d runs | Cohort_2023.gen: %d fish x %d loci",
            nrow(run_plan), length(fish_ids), ncol(genepop$dosage)))

# ---- 2. LD Nb per subset (cached) -------------------------------------------
ld_cache <- file.path(output_dir, "ld_by_subset.csv")
ld_code  <- unname(tools::md5sum(project_path("scripts", "lib", "ldne_direct.R")))
ld_done  <- if (file.exists(ld_cache)) read_csv(ld_cache, show_col_types = FALSE) else NULL
if (!is.null(ld_done) && !identical(unique(ld_done$ld_code), ld_code)) {
  msg("  LD estimator code changed since the cache was written; recomputing every subset")
  ld_done <- NULL
}
ld_todo  <- setdiff(run_plan$run_name, ld_done$run_name)
msg(sprintf("  LD: %d subsets cached, %d to compute", if (is.null(ld_done)) 0L else nrow(ld_done), length(ld_todo)))

for (run_name in ld_todo) {
  started <- Sys.time()
  rows    <- match(fish_subsets$sample[fish_subsets$run_name == run_name], fish_ids)
  result  <- ld_direct(genepop$dosage[rows, , drop = FALSE], genepop$loci)
  ld_done <- bind_rows(ld_done, tibble(
    run_name = run_name, ld_code = ld_code, ld_n_fish = length(rows),
    LD_Nb_all = result$all[["Ne"]], LD_Nb_all_waples = result$all[["Ne"]] / waples_divisor,
    LD_Nb_diff = result$diff[["Ne"]], r2_all = result$all[["r2"]], E_r2_all = result$all[["E"]]))
  write_csv(ld_done, ld_cache)        # saved after each subset, so a stop loses little
  msg(sprintf("    %-14s n = %3d  all %8.1f  diff %8.1f  [%.1f min]", run_name, length(rows),
              result$all[["Ne"]], result$diff[["Ne"]],
              as.numeric(difftime(Sys.time(), started, units = "mins"))))
}

# ---- 3. COLONY runs that have finished --------------------------------------
colony_rows <- lapply(seq_len(nrow(run_plan)), function(i) {
  run_dir     <- file.path(colony_out_dir, run_plan$run_name[i])
  config_file <- find_output(run_dir, "BestConfig")
  if (is.na(config_file)) return(NULL)
  config  <- read_best_config(config_file)
  pairs   <- sibling_pairs(config)
  parents <- table(c(config$father, config$mother))
  bind_cols(tibble(run_name = run_plan$run_name[i],
                   n_fs = sum(pairs$relationship == "full-sib"),
                   n_hs = sum(pairs$relationship == "half-sib")),
            as_tibble(as.list(np_estimates(as.numeric(parents)))),
            as_tibble(as.list(parse_ne(find_output(run_dir, "Ne")))))
})
colony_done <- bind_rows(colony_rows)
msg(sprintf("  COLONY: %d of %d runs finished and downloaded", nrow(colony_done), nrow(run_plan)))

by_run <- run_plan %>%
  select(run_name, n_fish, replicate) %>%
  left_join(ld_done, by = "run_name") %>%
  left_join(colony_done, by = "run_name") %>%
  arrange(n_fish, replicate)
write_csv(by_run, file.path(output_dir, "sampling_by_run.csv"))

# ---- 4. Summary and figure --------------------------------------------------
# Medians and ranges use finite runs only; n_infinite counts unbounded runs.
long <- by_run %>%
  select(n_fish, replicate, LD_Nb_all, LD_Nb_all_waples, LD_Nb_diff,
         any_of(c("Nb", "Np_Chao1", "n_hs"))) %>%
  rename(any_of(c(COLONY_Nb = "Nb"))) %>%
  pivot_longer(-c(n_fish, replicate), names_to = "estimate", values_to = "value")

finite_stat <- function(x, f) {
  x <- x[!is.na(x) & is.finite(x) & x >= 0]
  if (length(x)) f(x) else NA_real_
}
summary_table <- long %>%
  group_by(estimate, n_fish) %>%
  summarise(runs = sum(!is.na(value)),
            median = finite_stat(value, median),
            min = finite_stat(value, min),
            max = finite_stat(value, max),
            n_infinite = sum(!is.na(value) & (!is.finite(value) | value < 0)),
            .groups = "drop") %>%
  filter(runs > 0)
write_csv(summary_table, file.path(output_dir, "sampling_summary.csv"))
msg("\n  by sample size (median, range over finite runs; n_infinite = unbounded runs):")
msg_table(as.data.frame(summary_table %>% mutate(across(c(median, min, max), ~ round(.x)))),
          row.names = FALSE)

labels <- c(LD_Nb_all = "LD, all pairs", LD_Nb_all_waples = "LD, all pairs (Waples 2016)",
            LD_Nb_diff = "LD, different contigs", COLONY_Nb = "Sibship (COLONY)")
plot_data <- long %>%
  filter(estimate %in% names(labels), !is.na(value)) %>%
  mutate(unbounded = !is.finite(value) | value < 0,
         estimate  = factor(labels[estimate], levels = labels))
ceiling_value <- 10^ceiling(log10(max(plot_data$value[!plot_data$unbounded], 1e4)))
plot_data$shown <- ifelse(plot_data$unbounded, ceiling_value, plot_data$value)
# Small horizontal offsets so the four estimates do not sit on top of each other
plot_data$x_shown <- plot_data$n_fish + c(-4.5, -1.5, 1.5, 4.5)[as.integer(plot_data$estimate)]
medians <- plot_data %>% filter(!unbounded) %>%
  group_by(estimate, x_shown) %>% summarise(value = median(value), .groups = "drop")

figure <- ggplot(plot_data, aes(x_shown, shown, colour = estimate)) +
  geom_line(data = medians, aes(y = value)) +
  geom_point(aes(shape = unbounded), size = 2, alpha = 0.8) +
  scale_y_log10() +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2),
                     labels = c(`FALSE` = "finite", `TRUE` = "unbounded (plotted at top)")) +
  labs(x = "Fish sampled from the 2023 cohort", y = "Nb (log scale; lines = medians)",
       colour = NULL, shape = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", legend.box = "vertical")
ggsave(file.path(output_dir, "Figure_S_sampling_sensitivity.png"), figure,
       width = 7, height = 5, dpi = 300)

msg("\n  wrote ", output_dir)
msg("\n=== 10_sampling_sensitivity.R complete ===")
log_end()
