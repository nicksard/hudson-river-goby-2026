#!/bin/bash
#SBATCH --job-name=rg_colony
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/colony_%A_%a.out
# =============================================================================
# cluster/colony_array.sh
#
# Runs one COLONY .dat per array task. Adapted from
# nysdec_round_goby/cohort_rgsibship/02_run_colony_array.sh: same sandbox
# pattern, but the run ID is the .dat name minus ".dat", so it works for any
# file 06_colony_inputs.R writes (rgSibship_all_10943, rgSibship_maf_00250, ...).
#
# Setup on the cluster (once per batch):
#   PROJ_DIR/colony_dat_files/   <- copy results/06_colony/*.dat here
#   PROJ_DIR/job_list.txt        <- copy results/06_colony/job_list.txt here
#
# Submit from PROJ_DIR (line 1 of job_list.txt is the full-panel primary run):
#   mkdir -p logs
#   sbatch --array=1-$(wc -l < job_list.txt) colony_array.sh
#
# Outputs land in PROJ_DIR/colony_final_outputs/<run_id>_{Ne,BestConfig,BestCluster}.txt
# plus <run_id>_PairwiseFullSibDyad.txt / _PairwiseHalfSibDyad.txt for the dyad
# counts and the KING comparison.
# =============================================================================

PROJ_DIR="${PROJ_DIR:-$SLURM_SUBMIT_DIR}"
JOB_LIST="${PROJ_DIR}/job_list.txt"
DAT_DIR="${PROJ_DIR}/colony_dat_files"
FINAL_OUTPUT_DIR="${PROJ_DIR}/colony_final_outputs"
mkdir -p "$FINAL_OUTPUT_DIR"

DAT_FILENAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOB_LIST" | tr -d '\r')
if [ -z "$DAT_FILENAME" ]; then
  echo "Error: no filename on line $SLURM_ARRAY_TASK_ID of $JOB_LIST"; exit 1
fi
RUN_ID="${DAT_FILENAME%.dat}"
echo "Task $SLURM_ARRAY_TASK_ID | $(hostname) | $DAT_FILENAME"

RUN_DIR="${PROJ_DIR}/run_${RUN_ID}"
mkdir -p "$RUN_DIR"
cp "${DAT_DIR}/${DAT_FILENAME}" "${RUN_DIR}/colony2.dat"
cd "$RUN_DIR" || exit 1

/usr/bin/time -v colony

# Keep the files later steps parse. COLONY names outputs after the output
# name inside the .dat, so match on extension.
for ext in Ne BestConfig BestCluster PairwiseFullSibDyad PairwiseHalfSibDyad FullSibDyad HalfSibDyad; do
  for f in *."$ext"; do
    [ -e "$f" ] && mv "$f" "${FINAL_OUTPUT_DIR}/${RUN_ID}_${ext}.txt"
  done
done

cd "$PROJ_DIR" || exit 1
rm -rf "$RUN_DIR"
echo "Finished $RUN_ID"
