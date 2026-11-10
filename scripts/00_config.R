# =============================================================================
# 00_config.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Shared configuration: paths, filtering thresholds, grouping definitions
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-09-08
# Updated:  2026-09-10
#
# Description:
#   Every other script in scripts/ opens with source("scripts/00_config.R").
#   Nothing else in the pipeline contains a file path or a filtering threshold.
#   Two reasons for that. First, the Methods section of the manuscript quotes
#   these numbers, and if they live in five places they will eventually
#   disagree with each other and with the text. Second, a sensitivity run
#   ("what if MAF were 0.01?") should be a one-line edit here, not a hunt
#   through several hundred lines of analysis code.
#
#   There is no setwd() anywhere in the pipeline. Scripts find the project root
#   by walking up from the working directory looking for the .Rproj file, which
#   means the same code runs on the Windows laptop, on a shared drive, and on
#   the cluster without edits.
#
# Usage:
#   Rscript scripts/run_phase1.R          # from the project root
#   source("scripts/run_phase1.R")        # or in RStudio with the .Rproj open
# =============================================================================


# -----------------------------------------------------------------------------
# Locate the project root
# -----------------------------------------------------------------------------
# Walks up at most six levels from the working directory looking for the .Rproj
# file. Failing loudly here is deliberate: the alternative is a script that
# quietly writes results into a home directory and is never found again.

find_project_root <- function(max_levels = 6) {
  candidate_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (level in seq_len(max_levels)) {
    if (file.exists(file.path(candidate_dir, "hudson-river-goby-2026.Rproj"))) {
      return(candidate_dir)
    }
    parent_dir <- dirname(candidate_dir)
    if (identical(parent_dir, candidate_dir)) break   # hit the filesystem root
    candidate_dir <- parent_dir
  }
  stop("Cannot find the project root (hudson-river-goby-2026.Rproj).\n",
       "  Run from the project root, or open the .Rproj in RStudio.\n",
       "  Working directory is currently: ", getwd(), call. = FALSE)
}

project_root <- find_project_root()
project_path <- function(...) file.path(project_root, ...)

project_path()
# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# The two large inputs are not tracked in git and are not redistributed with
# the project. Point at your own copies either by dropping them in data/raw/ or
# by setting the environment variable named in the comment.
#
# On the shared drive these live at:
#   raw VCF   .../SOAC Lab/Round Goby/RADseq/R/post_RAD analysis/Input/
#             populations.snps.vcf
#   bio data  .../SOAC Lab/Round Goby/RADseq/R/plate_barcodes/Input/
#             GobyRad_tabularPlateWell with Bio Data.xlsx

paths <- list(

  # --- inputs (external; gitignored) ---
  raw_vcf  = Sys.getenv("RAW_VCF",  project_path("data/raw/populations.snps.vcf")),
  bio_data = Sys.getenv("BIO_DATA", project_path("data/raw/GobyRad_tabularPlateWell with Bio Data.xlsx")),

  # 2,000-locus x 384-sample slice of the raw VCF, tracked in git. Exists only
  # so the filtering script can be exercised without the full 1.23 GB file.
  test_vcf = project_path("data/test_subset.vcf.gz"),

  # --- outputs (produced here; safe to delete and regenerate) ---
  filtered_vcf  = project_path("results/01_filtering/populations.filtered.final.vcf.gz"),
  filter_log    = project_path("results/01_filtering/filtering_summary.csv"),
  filter_report = project_path("results/01_filtering/VCF_filtering_report.pdf"),
  sample_meta   = project_path("results/02_metadata/sample_metadata.csv"),

  dir_filtering = project_path("results/01_filtering"),
  dir_metadata  = project_path("results/02_metadata"),
  dir_popgen    = project_path("results/03_popgen")
)

for (output_dir in c(paths$dir_filtering, paths$dir_metadata, paths$dir_popgen)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}


# -----------------------------------------------------------------------------
# Filtering thresholds
# -----------------------------------------------------------------------------
# These are the values quoted in the Methods. Each is annotated with why it was
# chosen, because "MAF >= 0.05" on its own tells a reader (or a future us)
# nothing about whether 0.05 was reasoned or inherited.

filter_settings <- list(

  # Rare variants at this depth are dominated by sequencing error rather than
  # real polymorphism, and they carry almost no information for Fst or LD.
  min_maf = 0.05,

  # SNPs within one RAD locus are completely linked, so only one is kept. Which
  # one matters: the 2022 collection has markedly lower coverage than the other
  # years, so SNPs are ranked by how completely they genotyped in 2022 first and
  # overall second. Picking by MAF instead would systematically favour SNPs that
  # happen to be well covered in the three good years and would push the 2022
  # samples toward the missing-data cutoff.
  thin_priority_year = "2022",

  # Extreme heterozygote excess or deficit at a single locus is far more often a
  # genotyping artefact (a collapsed repeat, a null allele) than real biology in
  # a panmictic sample. Symmetric so both directions are caught.
  max_abs_fis = 0.25,

  # McKinney et al. (2017) paralog screen. Loci assembled from two duplicated
  # genomic regions look like permanent heterozygotes with skewed read ratios.
  # A true singleton locus should sit below both thresholds.
  max_paralog_hobs   = 0.60,
  max_paralog_zscore = 7.00,

  # Genotype-level quality. Below ~10 reads a heterozygote is frequently called
  # as a homozygote simply because one allele was not sampled, which biases Ho
  # downward and inflates Fis. GQ 20 is the matching Phred-scaled confidence.
  min_depth = 10,
  min_gq    = 20,

  # Individual call rate. Stricter than a typical population-genetics cutoff
  # because these genotypes also feed pedigree reconstruction, where missing
  # data at a locus removes that locus from every pairwise comparison involving
  # the individual.
  max_missing_per_individual = 0.40,

  # Locus call rate, applied after the individual cull so that loci are judged
  # on the samples actually being analysed.
  max_missing_per_locus = 0.20
)

# One plate well was labelled as a negative control but in fact held a fish.
# It is renamed so that it survives the NTC strip and its behaviour stays
# visible in the filtering log, then excluded by name in 02 (see below).
control_well_renames <- c("NTC4a" = "2024999")


# -----------------------------------------------------------------------------
# Quality-control exclusions
# -----------------------------------------------------------------------------
# Individuals dropped for reasons that have nothing to do with call rate.
#
# These are applied downstream of filtering, through sample_metadata.csv, and
# never inside 01_filter_vcf.R. Two reasons. First, these individuals were
# identified by diagnostics run on the already-filtered VCF, so folding the
# result back into the filter would make the pipeline circular. Second, the
# KING kinship screen needs the duplicate pair present in its input, because
# that pair is the positive control that shows the screen works.
#
# Downstream scripts subset on analysis_keep == TRUE. A script that wants the
# excluded individuals back (KING) reads the metadata and opts out explicitly.

qc_exclusions <- data.frame(
  sample = c("2024999", "2023018", "2024221"),
  reason = c(
    "negative-control well (renamed from NTC4a); not a fish",
    "contaminated: trimodal allele balance at heterozygous sites, highest observed heterozygosity in the dataset, elevated missingness; consistent with a two-individual DNA mixture",
    "duplicate genotype of 2023016 (99.6% concordance, zero opposite homozygotes); no extraction record and unreliable Lock C1 metadata, so 2023016 is the one retained"
  ),
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# Birth cohorts
# -----------------------------------------------------------------------------
# Cohort is derived, not recorded:
#
#     birth_year = collection_year - otolith_age
#
# Grouping by birth cohort rather than collection year matters here. A single
# collection year contains both Age-0 and Age-1 fish, so consecutive collection
# years share spawners. That shared drift history deflates pairwise Fst between
# adjacent years, and anything derived from that Fst inherits the bias.

cohort_settings <- list(

  # Below this a group cannot support an allele-frequency estimate at all; it
  # exists mainly to drop the single 2020-born fish.
  min_n_analysis = 5,

  # Below this a group can be described but its Fst and any Ne derived from it
  # are too imprecise to interpret. Groups under this are reported with a
  # warning rather than silently dropped.
  min_n_estimator = 30
)


### NEW -- temporal Ne (04_temporal_ne.R)
# Comparisons spanning fewer years than this are dropped (Waples 1990; Waples &
# Yokota 2007). The Methods quote this value.
temporal_ne_settings <- list(
  min_years_apart           = 2,
  primary_generation_length = 1,                    # years; the value in the text
  generation_lengths        = c(1.0, 1.5, 2.0, 2.5) # sensitivity table
)


### NEW -- COLONY sibship inputs (06_colony_inputs.R)
# Values carried over unchanged from the September runs (Pedigree Pipeline for
# loop.V6.R and the Figure 4 run, "Medium_Medium_ALF_Yes"), so the rerun differs
# from the stand-in results only in the data. The Methods quote these.
colony_settings <- list(
  project_name         = "rgSibship",
  error_rate_geno      = 0.02,    # allelic dropout
  error_rate_other     = 0.001,   # other genotyping error
  run_length           = 2,       # 1 short, 2 medium, 3 long, 4 very long
  likelihood_precision = 1,       # 0 low, 1 medium, 2 high, 3 very high
  update_alfs          = 1,       # 1 = update allele frequencies
  spp_type             = 2,       # dioecious
  inbreeding           = 0,
  female_gamy          = 0,       # 0 = polygamy
  male_gamy            = 0,
  full_likelihood      = 1,
  sib_scale            = 1,
  sib_prior            = 0,

  # Optional locus sweep (rebuilds the locus-number axis of Figure 3) at the
  # same settings. The full panel is always written as the primary run.
  run_locus_sweep = TRUE,
  sweep_loci      = c(250, 500, 1000, 2000, 5000, 8000),
  sweep_subsets   = c("maf", "random"),

  # Path to the COLONY executable, used only to write run_all.bat. Edit to
  # match your install.
  colony_exe = "C:/Program Files (x86)/COLONY2/colony2s.exe"
)


### NEW -- LDNe inputs (07_ldne_inputs.R) and Table 3 corrections (08)
ldne_settings <- list(
  pcrit  = 0.05,       # the critical value reported in Table 3
  mating = 0,          # 0 = random mating, 1 = monogamy

  # Groups written as one .gen each. "all" = the whole analysis set; "cohort" =
  # one file per birth cohort. Add "albany_troy_cohort" to also write the
  # Albany+Troy-pooled cohort files the earlier LDNe tables used.
  groupings = c("all", "cohort", "near_complete"),

  ### NEW -- "near_complete": a robustness subset with almost no missing data,
  # where NeEstimator's missing-data revision (Peel et al. 2013) should change
  # nothing. Fish are filtered first, then loci are judged on the kept fish.
  # 0.05 / 0.02 keeps ~256 fish x ~8,775 loci at ~0.7% residual missingness.
  near_complete_max_fish_missing  = 0.05,
  near_complete_max_locus_missing = 0.02,

  # Waples et al. (2016) finite-chromosome correction: divide Nb and both CI
  # limits by 0.098 + 0.219 * ln(n_chromosomes). Round goby n = 23 (cite karyotype).
  n_chromosomes = 23,

  ne_exe = "G:/My Drive/Programs/NeEstimator/Ne2-1.exe"
)


# -----------------------------------------------------------------------------
# Location
# -----------------------------------------------------------------------------
# The Region field in the workbook is free text, with 17 distinct values
# including four different spellings of the Lock C1 site. `region_patterns`
# collapses them into analysis units.
#
# Keys are case-insensitive regular expressions matched against the trimmed
# Region string; FIRST MATCH WINS, so order matters. 02_sample_metadata.R
# prints every raw level next to the label it received and stops if any level
# is unmatched, so a new site appearing in a future version of the workbook
# cannot be quietly folded into an existing group.
#
# NOTE: the sites below Albany were added on 2026-09-10 and their assignment to
# river reaches is Claude's proposal from place names alone -- Nick to confirm,
# in particular Mill Creek, which is ambiguous (several in the Hudson valley).

region_patterns <- c(
  # above the Federal Dam at Troy
  "Lock[ _]?C1"  = "Lock_C1",
  "Mohawk"       = "Mohawk_River",

  # the invasion front itself
  "Poesten"      = "Poesten_Kill",   # tributary entering at Troy
  "Troy"         = "Troy",
  "Albany"       = "Albany",

  # below Albany, within the established range
  "Coxsackie"    = "Coxsackie",
  "Hannacroix"   = "Hannacroix",
  "Catskill"     = "Catskill",
  "Black Creek"  = "Black_Creek",
  "Twaalfkill"   = "Twaalfkill",
  "Mill Creek"   = "Mill_Creek",     
  "Lloyd"        = "Lloyd",
  "Poughkeepsie" = "Poughkeepsie",
  "Roseton"      = "Roseton"
)

# Approximate downstream ordering, used only to order factor levels so that
# figures read north-to-south instead of alphabetically. Not used to compute
# anything.
region_downstream_order <- c(
  "Lock_C1", "Mohawk_River",
  "Poesten_Kill", "Troy", "Albany",
  "Coxsackie", "Hannacroix", "Catskill",
  "Black_Creek", "Twaalfkill", "Mill_Creek",
  "Lloyd", "Poughkeepsie", "Roseton"
)

# Named ways of turning region into a grouping variable. Each is a function of
# the metadata table returning one label per row; NA means "leave this fish out
# of this particular scheme". Adding an entry here makes it available to every
# downstream script without touching their code.
location_schemes <- list(

  # Every site kept separate. Most are far too small for Fst, but this is the
  # view that shows where the sampling effort actually went.
  site = function(sample_info) sample_info$region,

  # Three reaches relative to the Federal Dam at Troy, which is the barrier the
  # invasion front is pressing against. This is the scheme that uses all the
  # data, and it is the one to use as primary unless a reach split proves
  # uninformative.
  reach = function(sample_info) {
    dplyr::case_when(
      sample_info$region %in% c("Lock_C1", "Mohawk_River")            ~ "Above_dam",
      sample_info$region %in% c("Poesten_Kill", "Troy", "Albany")     ~ "Front",
      is.na(sample_info$region)                                       ~ NA_character_,
      TRUE                                                            ~ "Downstream"
    )
  },

  # Albany and Troy pooled, everything else left as-is. This is the grouping the
  # existing LDNe tables use, kept so the new numbers stay comparable to them.
  albany_troy = function(sample_info) {
    ifelse(sample_info$region %in% c("Albany", "Troy"), "Albany_Troy", as.character(sample_info$region))   ### CHANGED: region is a factor, so ifelse() was returning its integer codes
  }
)

# The groupings 03_popgen.R actually runs. `by` is one of "cohort", "location"
# or "both"; `location_scheme` names an entry in location_schemes above.
#
# cohort_by_site is the one that matters analytically. Sampling effort moved
# upstream over the study: the 2021 cohort is almost entirely Albany fish and
# the 2024 cohort is almost entirely Troy fish, so a plain cohort contrast
# cannot separate a temporal signal from a spatial one. Crossing cohort with
# SITE (not reach, and not the Albany+Troy pool, both of which hide the problem
# rather than fixing it) isolates the cells where only one factor varies:
#
#     Albany 2021 (n = 31) vs Albany 2023 (n = 62)   temporal, site held fixed
#     Troy   2023 (n = 75) vs Troy   2024 (n = 92)   temporal, site held fixed
#     Albany 2023 (n = 62) vs Troy   2023 (n = 75)   spatial, cohort held fixed
#
# Those three comparisons are the interpretable ones. Everything else in the
# matrix confounds the two factors.
popgen_groupings <- list(
  list(name = "cohort",            by = "cohort",   location_scheme = NULL),
  list(name = "reach",             by = "location", location_scheme = "reach"),
  list(name = "site",              by = "location", location_scheme = "site"),
  list(name = "cohort_by_reach",   by = "both",     location_scheme = "reach"),
  list(name = "cohort_by_site",    by = "both",     location_scheme = "site")
)

# Pairs to report as the primary, unconfounded contrasts. Named so that
# 03_popgen.R can pull them out of the cohort_by_site matrix and print them
# separately from the rest, which are not interpretable on their own.
clean_contrasts <- list(
  list(label = "Albany: 2021 vs 2023", a = "2021_Albany", b = "2023_Albany",
       type = "temporal (site fixed)"),
  list(label = "Troy: 2023 vs 2024",   a = "2023_Troy",   b = "2024_Troy",
       type = "temporal (site fixed)"),
  list(label = "2023: Albany vs Troy", a = "2023_Albany", b = "2023_Troy",
       type = "spatial (cohort fixed)")
)


# -----------------------------------------------------------------------------
# Session settings and small shared helpers
# -----------------------------------------------------------------------------

random_seed <- 20260908
set.seed(random_seed)

# Sample IDs look like 2021001 and will be silently read as numbers unless told
# otherwise, which then fails to match the character sample names in the VCF.
# Every read and write of the metadata table goes through these.
sample_id_col_types <- readr::cols(sample = readr::col_character(),
                                   .default = readr::col_guess())

# Console messages are also written to results/log_<script>.txt so that a run
# leaves a transcript, rather than the record living only in whatever console
# buffer happens to still be open.
log_connection <- NULL

log_start <- function(script_name) {
  log_file <- project_path("results", paste0("log_", script_name, ".txt"))
  log_connection <<- file(log_file, open = "wt")
  msg("=== ", script_name, " | ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ===")
  msg("R ", getRversion(), " | seed ", random_seed)
  invisible(log_file)
}

msg <- function(...) {
  line <- paste0(...)
  message(line)
  if (!is.null(log_connection)) writeLines(line, log_connection)
  invisible(line)
}

### NEW -- print() output goes to the console only, so tables never reached the
### log. msg_table() prints as before and also writes each line to the log.
msg_table <- function(x, ...) {
  table_lines <- utils::capture.output(print(x, ...))
  for (line in table_lines) msg(line)
  invisible(x)
}

log_end <- function() {
  if (!is.null(log_connection)) {
    close(log_connection)
    log_connection <<- NULL
  }
}

# Stop with a readable explanation instead of letting a bad assumption surface
# as a cryptic error two hundred lines later.
need <- function(condition, ...) {
  if (!isTRUE(condition)) stop(paste0(...), call. = FALSE)
}

need_pkgs <- function(package_names) {
  missing_pkgs <- package_names[
    !vapply(package_names, requireNamespace, logical(1), quietly = TRUE)
  ]
  need(length(missing_pkgs) == 0,
       "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
       '  install.packages(c("', paste(missing_pkgs, collapse = '", "'), '"))')
}

