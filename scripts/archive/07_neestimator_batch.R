################################################################################
### 07_neestimator_batch.R
###
### Writes one NeEstimator "info" file per .gen and a .bat that runs them all.
### The GUI does exactly this behind the scenes for one file at a time.
###
### VERIFY THE INFO-FILE LAYOUT ONCE before batching:
###   run one file through the GUI with your usual settings, then look in the
###   output folder for the info file the GUI wrote (usually named after your
###   output). Match its line order here if it differs. The layout below is
###   NeEstimator 2.1's for LD-only runs; the temporal method adds lines.
################################################################################

library(tidyverse)

gen_dir   <- "G:/Shared drives/SOAC Lab/Round Goby/RADseq/R/pedigree reconstruction/Output/genepop_cohorts_v2"
out_dir   <- file.path(gen_dir, "neestimator_out")
ne_exe <- "G:/My Drive/Programs/NeEstimator/Ne2-1.exe"
crit      <- c(0.05)              # Pcrit values; 0.05 is the one in Table 3
mating    <- 0                                # 0 = random mating, 1 = monogamy

dir.create(out_dir, showWarnings = FALSE)

# NeEstimator wants trailing slashes and Windows-style paths
win <- function(p) paste0(gsub("/", "\\\\", p), "\\")

gen_files <- list.files(gen_dir, pattern = "\\.gen$", full.names = FALSE)
if (!length(gen_files)) stop("no .gen files found in ", gen_dir)

info_files <- character()
for (f in gen_files) {
  stem <- sub("\\.gen$", "", f)
  info <- file.path(out_dir, paste0(stem, "_info.txt"))
  writeLines(c(
    "1",                                        # method: 1 = LD only
    win(gen_dir),                               # input directory
    f,                                          # input file
    "2",                                        # format: 2 = GENEPOP
    win(out_dir),                               # output directory
    paste0(stem, "_Ne.txt"),                    # output file
    as.character(length(crit)),                 # number of critical values
    paste(crit, collapse = " "),                # critical values
    as.character(mating)                        # mating model
  ), info)
  info_files <- c(info_files, info)
}

# --- Windows batch file ---
bat <- file.path(out_dir, "run_all.bat")
writeLines(c(
  "@echo off",
  sprintf("cd /d \"%s\"", gsub("/", "\\\\", out_dir)),
  sprintf("\"%s\" i:\"%s\"", gsub("/", "\\\\", ne_exe), basename(info_files)),
  "echo ALL RUNS COMPLETE",
  "pause"
), bat)

cat(sprintf("Wrote %d info files and %s\n", length(info_files), bat))
cat("Double-click run_all.bat, or from a cmd prompt:  ", bat, "\n")
