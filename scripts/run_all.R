# =============================================================================
# run_all.R
#
# Hudson River round goby (Neogobius melanostomus) -- RADseq
# One entry point: runs each numbered script only if it needs to run
#
# Author:   Nick Sard
# Coauthor: Claude (Anthropic), coding assistant
# Created:  2026-11-11
#
# How it decides:
#   Each step lists its inputs and outputs (table below). A step runs when
#     * any of its outputs is missing, or
#     * any input (including the scripts/lib files it uses), its own script
#       or 00_config.R has changed since the step last ran successfully.
#   "Changed" means different file contents (MD5), recorded in
#   results/.pipeline_state.rds after each successful step. Contents rather
#   than timestamps, because a git checkout or a copy touches timestamps
#   without changing anything. The first time (no record yet) timestamps are
#   used: the step runs if an input is newer than its oldest output.
#   A step that reruns rewrites its outputs, so everything downstream whose
#   inputs actually changed runs too.
#
# External steps (cluster COLONY runs, NeEstimator) cannot run from here. The
# runner says what to do and skips the steps that depend on their outputs.
#
# Usage (from the project root, or source() in RStudio):
#   Rscript scripts/run_all.R              run whatever is out of date
#   Rscript scripts/run_all.R --dry-run    only report what would run and why
#   Rscript scripts/run_all.R --force 03   rerun 03 (and whatever that changes)
#   Rscript scripts/run_all.R --only 08b   consider only that step
#   Rscript scripts/run_all.R --record     mark every step whose outputs exist
#                                          as up to date now, without running
#                                          (set a baseline once, e.g. after
#                                          editing config or moving code)
#
# Manual steps (01 filtering, 06 COLONY inputs) never run on their own: if
# they are out of date the runner says so, and --force 01 / --force 06 runs
# them. Re-filtering or regenerating cluster inputs should be deliberate.
# =============================================================================

source("scripts/00_config.R")

args    <- if (exists("run_all_args")) run_all_args else commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
force   <- if ("--force" %in% args) args[which(args == "--force") + 1] else character()
only    <- if ("--only"  %in% args) args[which(args == "--only")  + 1] else character()
record  <- "--record" %in% args

p <- function(...) project_path(...)

# -----------------------------------------------------------------------------
# The pipeline. inputs/outputs are files or directories (a directory stands for
# every file in it). optional inputs count when present and never block.
# -----------------------------------------------------------------------------
steps <- list(
  list(id = "01", script = "01_filter_vcf.R", manual = TRUE,
       inputs = c(paths$raw_vcf),
       outputs = c(paths$filtered_vcf, paths$filter_log)),
  list(id = "02", script = "02_sample_metadata.R",
       inputs = c(paths$filtered_vcf, paths$bio_data),
       outputs = c(paths$sample_meta)),
  list(id = "03", script = "03_popgen.R",
       inputs = c(paths$filtered_vcf, paths$sample_meta),
       outputs = c(p("results/03_popgen/grouping_overview.csv"))),
  list(id = "04", script = "04_temporal_ne.R",
       inputs = c(p("results/03_popgen/grouping_overview.csv")),
       outputs = c(p("results/04_temporal_ne/temporal_ne.csv"))),
  list(id = "07", script = "07_ldne_inputs.R",
       inputs = c(paths$filtered_vcf, paths$sample_meta),
       outputs = c(p("results/07_ldne/ldne_manifest.csv"))),
  list(id = "08b", script = "08b_ldne_direct.R",
       inputs = c(p("results/07_ldne/ldne_manifest.csv"), p("scripts/lib/ldne_direct.R")),
       outputs = c(p("results/07_ldne/direct/ld_ne_direct.csv"))),
  list(id = "08c", script = "08c_ldne_stratag.R", needs_package = "strataG",
       inputs = c(p("results/07_ldne/ldne_manifest.csv"), p("scripts/lib/ldne_direct.R")),
       outputs = c(p("results/07_ldne/stratag/ldne_stratag_check.csv"))),
  list(id = "06", script = "06_colony_inputs.R", manual = TRUE,
       inputs = c(paths$filtered_vcf, paths$sample_meta, p("scripts/lib/colony_dat_writer.R")),
       outputs = c(p("results/06_colony/run_plan.csv")),
       then = paste("COLONY runs on the cluster: upload results/06_colony/*.dat,",
                    "submit the array, download finished run folders into",
                    "results/06_colony/outputs/.")),
  list(id = "09", script = "09_colony_summary.R",
       inputs = c(p("results/06_colony/run_plan.csv"), p("scripts/lib/colony_parse.R")),
       external_inputs = c(p("results/06_colony/outputs")),
       outputs = c(p("results/09_colony/colony_run_metrics.csv"))),
  list(id = "05", script = "05_king.R",
       inputs = c(p("results/07_ldne/ldne_manifest.csv")),
       optional = c(p("results/09_colony/primary_sibling_pairs.csv")),
       outputs = c(p("results/05_king/king_summary.csv"))),
  # Sample-size sensitivity: 06b writes the subsets (run it by hand, like 06);
  # 10 runs the LD half at once and adds COLONY runs as they are downloaded.
  list(id = "10", script = "10_sampling_sensitivity.R",
       inputs = c(p("results/06b_colony_sampling/fish_subsets.csv"),
                  p("results/07_ldne/ldne_manifest.csv"),
                  p("scripts/lib/ldne_direct.R"), p("scripts/lib/colony_parse.R")),
       optional = c(p("results/06b_colony_sampling/outputs")),
       outputs = c(p("results/10_sampling_sensitivity/sampling_summary.csv")))
)
# 00_config.R holds every step's settings, so a change there reruns everything.
# Library files are listed as inputs of the steps that use them.
shared_inputs <- c(p("scripts/00_config.R"))

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
present <- function(x) if (length(x)) x[file.exists(x)] else character()
expand <- function(x) unlist(lapply(x, function(f)
  if (dir.exists(f)) list.files(f, recursive = TRUE, full.names = TRUE) else f))

fingerprint <- function(files) {
  files <- sort(expand(files)); files <- files[file.exists(files)]
  if (!length(files)) return(NA_character_)
  sums <- tools::md5sum(files)
  paste(sub(p(""), "", names(sums), fixed = TRUE), unname(sums), collapse = "|")
}

state_file <- p("results/.pipeline_state.rds")
state <- if (file.exists(state_file)) readRDS(state_file) else list()

decide <- function(step) {
  script <- p("scripts", step$script)
  if (!file.exists(script)) return(list(run = FALSE, why = "script not in repo yet"))
  if (!is.null(step$needs_package) && !requireNamespace(step$needs_package, quietly = TRUE))
    return(list(run = FALSE, why = paste0("package ", step$needs_package, " not installed (optional step)")))
  missing_inputs <- step$inputs[!file.exists(step$inputs)]
  if (length(missing_inputs))
    return(list(run = FALSE, blocked = TRUE, why = paste("waiting for", paste(basename(missing_inputs), collapse = ", "))))
  if (length(step$external_inputs) && !length(expand(step$external_inputs)))
    return(list(run = FALSE, blocked = TRUE, why = paste("no files yet in", step$external_inputs)))
  if (step$id %in% force) return(list(run = TRUE, why = "forced"))
  if (any(!file.exists(step$outputs)))
    return(list(run = TRUE, why = paste("missing", paste(basename(step$outputs[!file.exists(step$outputs)]), collapse = ", "))))
  all_inputs <- c(script, shared_inputs, step$inputs, step$external_inputs,
                  present(step$optional))
  now <- fingerprint(all_inputs)
  before <- state[[step$id]]
  if (is.null(before)) {
    newest_input <- max(file.mtime(expand(all_inputs)), na.rm = TRUE)
    oldest_output <- min(file.mtime(step$outputs))
    if (newest_input > oldest_output) return(list(run = TRUE, why = "an input is newer than the outputs (no record yet)", now = now))
    return(list(run = FALSE, why = "up to date (recorded now)", now = now, record = TRUE))
  }
  if (!identical(before, now)) {
    old <- strsplit(before, "|", fixed = TRUE)[[1]]; new <- strsplit(now, "|", fixed = TRUE)[[1]]
    changed <- unique(sub(" .*$", "", setdiff(new, old)))
    return(list(run = TRUE, why = paste("changed:", paste(basename(changed), collapse = ", ")), now = now))
  }
  list(run = FALSE, why = "up to date")
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
if (length(only)) steps <- Filter(function(s) s$id %in% only, steps)
message(strrep("=", 72), "\n run_all.R",
        if (record) " (recording baseline)" else if (dry_run) " (dry run)" else "",
        "\n", strrep("=", 72))

for (step in steps) {
  decision <- decide(step)
  if (record) {
    if (all(file.exists(step$outputs)) && !isTRUE(decision$blocked) &&
        file.exists(p("scripts", step$script))) {
      state[[step$id]] <- fingerprint(c(p("scripts", step$script), shared_inputs, step$inputs,
                                        step$external_inputs, present(step$optional)))
      saveRDS(state, state_file)
      message(sprintf("[%s] %-24s recorded as up to date", step$id, step$script))
    } else message(sprintf("[%s] %-24s not recorded (outputs missing or inputs not ready)",
                           step$id, step$script))
    next
  }
  if (decision$run && isTRUE(step$manual) && !(step$id %in% force)) {
    message(sprintf("\n[%s] %-24s OUT OF DATE (manual step) -- %s\n      run it deliberately with --force %s",
                    step$id, step$script, decision$why, step$id))
    next
  }
  message(sprintf("\n[%s] %-24s %s", step$id, step$script,
                  if (decision$run) paste("RUN  --", decision$why) else paste("skip --", decision$why)))
  if (isTRUE(decision$record) && !dry_run) {
    state[[step$id]] <- decision$now; saveRDS(state, state_file)
  }
  if (!decision$run || dry_run) {
    if (!is.null(step$then) && file.exists(step$outputs[1]))
      message("      next (outside R): ", step$then)
    next
  }
  started <- Sys.time()
  ok <- tryCatch({ source(p("scripts", step$script), local = new.env()); TRUE },
                 error = function(e) { message("      !! failed: ", conditionMessage(e)); FALSE })
  if (!ok) stop("Stopped at step ", step$id, "; fix it and run again (completed steps will be skipped).")
  state[[step$id]] <- fingerprint(c(p("scripts", step$script), shared_inputs, step$inputs,
                                    step$external_inputs, present(step$optional)))
  saveRDS(state, state_file)
  message(sprintf("      done in %.1f min", as.numeric(difftime(Sys.time(), started, units = "mins"))))
  if (!is.null(step$then)) message("      next (outside R): ", step$then)
}
message("\nAll steps checked.")
