# =============================================================================
# run_phase1.R -- raw VCF to population genetics, in order
#
#   Rscript R/run_phase1.R           full run
#   Rscript R/run_phase1.R --test    smoke test on data/test_subset.vcf.gz
#                                    (filtering only; the metadata and popgen
#                                     steps need the bio workbook)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
test <- "--test" %in% args

steps <- if (test) {
  "R/01_filter_vcf.R"
} else {
  c("R/01_filter_vcf.R", "R/02_sample_metadata.R", "R/03_popgen.R")
}

for (s in steps) {
  message("\n", strrep("#", 72), "\n#  ", s, "\n", strrep("#", 72))
  t0 <- Sys.time()
  if (test) {
    # forward --test to the filtering script
    local({ commandArgs <- function(trailingOnly = FALSE) "--test"; source(s, local = TRUE) })
  } else {
    source(s, local = new.env())
  }
  message(sprintf("#  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

message("\nPhase 1 complete. Outputs under results/.")
