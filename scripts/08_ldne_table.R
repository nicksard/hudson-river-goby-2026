# =============================================================================
# 08_ldne_table.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Parse NeEstimator LD output into Table 3, with the Waples et al. (2016)
# finite-chromosome correction applied to the point estimate AND both CI limits
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-10
#
# Description:
#   Reads every <group>_Ne.txt that run_all.bat wrote (listed in
#   ldne_manifest.csv), takes the column for ldne_settings$pcrit, and writes
#   raw and corrected values side by side.
#
#   Correction (Waples et al. 2016): divisor = 0.098 + 0.219 * ln(n_chr);
#   n_chr = 23 gives 0.7847. Nb, L95 and U95 are each divided by it; an
#   infinite limit stays infinite. This replaces the old Table 3, where some
#   rows paired corrected CIs with uncorrected points.
#
#   Consistency check. Each row also gets Nb_from_r2: the Waples & Do (2008)
#   S >= 30 formula applied to the r^2 and E(r^2) values NeEstimator printed,
#
#       r2' = r2 - E(r2_sample),  Nb = (1/3 + sqrt(1/9 - 2.76 r2')) / (2 r2')
#
#   and two flags: whether that reproduces the printed point estimate, and
#   whether the printed CI contains the printed point estimate. In the
#   September 2023-cohort output the CI (60.7-60.8) did NOT contain the point
#   estimate (214.7), and the CI matched this formula while the point estimate
#   did not. Do not lock Table 3 while either flag is FALSE.
#
# Inputs:   results/07_ldne/ldne_manifest.csv, results/07_ldne/neestimator/*_Ne.txt
# Outputs:  results/08_ldne_table/table3_ldne.csv
#
# Usage:    Rscript scripts/08_ldne_table.R
# =============================================================================

source("scripts/00_config.R")

need_pkgs(c("dplyr", "readr"))
suppressPackageStartupMessages({ library(dplyr); library(readr) })

log_start("08_ldne_table")

settings   <- ldne_settings
ldne_dir   <- project_path("results", "07_ldne")
ne_dir     <- file.path(ldne_dir, "neestimator")
output_dir <- project_path("results", "08_ldne_table")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest_file <- file.path(ldne_dir, "ldne_manifest.csv")
need(file.exists(manifest_file), "Run scripts/07_ldne_inputs.R first.")
manifest <- read_csv(manifest_file, show_col_types = FALSE)

chr_divisor <- 0.098 + 0.219 * log(settings$n_chromosomes)
msg(sprintf("  Waples et al. 2016 divisor (n_chr = %d): %.4f", settings$n_chromosomes, chr_divisor))


# -----------------------------------------------------------------------------
# Parse one NeEstimator LD output file
# -----------------------------------------------------------------------------
# Values are laid out in columns, one per "Lowest Allele Frequency Used".
# "Infinite" is NeEstimator's word for an unbounded estimate or limit.
as_value <- function(x) ifelse(grepl("^Inf", x, ignore.case = TRUE), Inf, suppressWarnings(as.numeric(x)))

parse_ne_file <- function(ne_file, pcrit) {
  lines <- readLines(ne_file, warn = FALSE)
  values_after <- function(pattern, offset = 0) {
    hit <- grep(pattern, lines)
    need(length(hit) >= 1, "Could not find '", pattern, "' in ", basename(ne_file))
    line <- lines[hit[1] + offset]
    fields <- strsplit(trimws(sub("^.*=", "", sub("^\\*[^0-9I]*", "", line))), "\\s+")[[1]]
    as_value(fields)
  }
  crit_line <- lines[grep("^Lowest Allele Frequency Used", lines)[1]]
  crit <- strsplit(trimws(sub("Lowest Allele Frequency Used", "", crit_line)), "\\s+")[[1]]
  col <- match(sprintf("%.3f", pcrit), crit)
  need(!is.na(col), "Pcrit ", pcrit, " not among the columns of ", basename(ne_file))

  param_l <- values_after("^\\* Parametric")
  param_u <- as_value(strsplit(trimws(lines[grep("^\\* Parametric", lines)[1] + 1]), "\\s+")[[1]])
  jack_l  <- values_after("^\\* JackKnife")
  jack_u  <- as_value(strsplit(trimws(lines[grep("^\\* JackKnife", lines)[1] + 1]), "\\s+")[[1]])

  data.frame(
    harmonic_S  = values_after("^Harmonic Mean Sample Size")[col],
    comparisons = values_after("^Independent Comparisons")[col],
    r2          = values_after("^OverAll r\\^2")[col],
    E_r2        = values_after("^Expected r\\^2 Sample")[col],
    Nb          = values_after("^Estimated Ne\\^")[col],
    Nb_L95_param = param_l[col], Nb_U95_param = param_u[col],
    Nb_L95_jack  = jack_l[col],  Nb_U95_jack  = jack_u[col]
  )
}

# Waples & Do (2008), S >= 30
nb_from_r2 <- function(r2_prime) {
  ifelse(r2_prime <= 0, Inf, (1/3 + sqrt(pmax(1/9 - 2.76 * r2_prime, 0))) / (2 * r2_prime))
}


# -----------------------------------------------------------------------------
# Build the table
# -----------------------------------------------------------------------------
rows <- list()
for (i in seq_len(nrow(manifest))) {
  ne_file <- file.path(ne_dir, manifest$ne_output[i])
  if (!file.exists(ne_file)) {
    msg("  !! missing ", manifest$ne_output[i], " (run not finished?) -- skipped")
    next
  }
  rows[[i]] <- cbind(group = manifest$group[i], n = manifest$n[i],
                     parse_ne_file(ne_file, settings$pcrit))
}
need(length(rows) > 0, "No NeEstimator outputs found in ", ne_dir)

table3 <- bind_rows(rows) %>%
  mutate(
    Nb_from_r2        = nb_from_r2(r2 - E_r2),
    point_reproduces  = abs(Nb_from_r2 - Nb) / Nb < 0.02,
    ci_contains_point = Nb >= Nb_L95_jack & Nb <= Nb_U95_jack,
    # Waples et al. 2016 correction: point and both limits
    Nb_adj           = Nb / chr_divisor,
    Nb_L95_jack_adj  = Nb_L95_jack / chr_divisor,
    Nb_U95_jack_adj  = ifelse(is.infinite(Nb_U95_jack), Inf, Nb_U95_jack / chr_divisor),
    flag = ifelse(n < cohort_settings$min_n_estimator, "low precision (n < 30)", "")
  )

write_csv(table3, file.path(output_dir, "table3_ldne.csv"))

msg("\n  Table 3 (Pcrit = ", settings$pcrit, "; jackknife CIs; adj = / ", round(chr_divisor, 4), "):")
msg_table(as.data.frame(table3 %>%
  transmute(group, n, Nb = round(Nb, 1),
            CI = sprintf("%.1f-%.1f", Nb_L95_jack, Nb_U95_jack),
            Nb_adj = round(Nb_adj, 1),
            CI_adj = sprintf("%.1f-%.1f", Nb_L95_jack_adj, Nb_U95_jack_adj), flag)),
  row.names = FALSE)

msg("\n  consistency check (Waples & Do 2008 formula on the printed r^2):")
msg_table(as.data.frame(table3 %>%
  transmute(group, r2 = signif(r2, 5), E_r2 = signif(E_r2, 5), Nb = round(Nb, 1),
            Nb_from_r2 = round(Nb_from_r2, 1), point_reproduces, ci_contains_point)),
  row.names = FALSE)

if (!all(table3$point_reproduces & table3$ci_contains_point)) {
  msg("\n  !! NeEstimator's point estimate and CI disagree for at least one group.")
  msg("  !! Do not lock Table 3 until this is resolved.")
}

msg("\n  wrote ", file.path(output_dir, "table3_ldne.csv"))
msg("\n=== 08_ldne_table.R complete ===")
log_end()
