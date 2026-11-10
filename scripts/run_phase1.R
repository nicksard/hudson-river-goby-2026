# =============================================================================
# run_phase1.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# Phase 1 driver: raw VCF -> filtered panel -> sample table -> population genetics
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-09-08
# Updated:  2026-09-10
#
# Description:
#   Runs the three Phase 1 scripts in order, timing each. Each script is sourced
#   into its own environment so that an object left behind by one cannot be
#   picked up accidentally by the next -- the scripts communicate only through
#   files on disk, which is what makes them individually rerunnable.
#
# Usage:
#   Rscript scripts/run_phase1.R          full run
#   Rscript scripts/run_phase1.R --test   filtering only, on the 2,000-locus subset
#                                         (02 and 03 need the bio workbook)
# =============================================================================

command_args <- commandArgs(trailingOnly = TRUE)
run_on_test_subset <- "--test" %in% command_args

phase1_scripts <- if (run_on_test_subset) {
  "scripts/01_filter_vcf.R"
} else {
  c("scripts/01_filter_vcf.R",
    "scripts/02_sample_metadata.R",
    "scripts/03_popgen.R")
}

for (script_path in phase1_scripts) {

  message("\n", strrep("#", 72))
  message("#  ", script_path)
  message(strrep("#", 72))

  started_at <- Sys.time()

  if (run_on_test_subset) {
    # Shadow commandArgs inside the sourced environment so that 01 sees --test
    # without the flag having to be plumbed through as a global.
    local({
      commandArgs <- function(trailingOnly = FALSE) "--test"
      source(script_path, local = TRUE)
    })
  } else {
    source(script_path, local = new.env())
  }

  elapsed_minutes <- as.numeric(difftime(Sys.time(), started_at, units = "mins"))
  message(sprintf("#  done in %.1f min", elapsed_minutes))
}

message("\nPhase 1 complete. Outputs are under results/.")
