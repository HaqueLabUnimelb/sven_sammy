#!/bin/bash
#SBATCH --job-name=ss_compare_reps
#SBATCH --account=uom00060
#SBATCH --partition=sapphire
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=logs/compare_replicates_%j.log

# Runs compare_replicates.py on all Hotspot replicates (after the four parts of
# hotspot_replicates.sbatch). Needs only numpy, pandas and scikit-learn; it uses the Python
# of the module Anaconda3/2022.10.

# ---- paths: edit these ----------------------------------------------------------
STAGING=$(cd "${SLURM_SUBMIT_DIR:-.}/.." && pwd)   # submit this script from STAGING/tests
REPO=$(cd "$STAGING/../.." && pwd)
ARCHIVED_HOTSPOT=/data/projects/punim1236/liver_stomicsv1.3/hotspot
UNTITLED3=/data/projects/punim1236/b_cells/Untitled3.ipynb
KEEP_ROOT=/path/to/kept_outputs
OUT=$KEEP_ROOT/hotspot_replicates
A1=$KEEP_ROOT/hotspot_35423407/test_a        # test a of job 35423407
# ---------------------------------------------------------------------------------

for f in /etc/profile.d/*.sh; do source "$f"; done
module purge
module load Anaconda3/2022.10
unset PYTHONPATH
export PYTHONNOUSERSITE=1
T=$STAGING/tests
ls -la "$OUT" "$OUT/reps"
python "$T/compare_replicates.py" --archived "$ARCHIVED_HOTSPOT" \
    --notebook "$REPO/3_STOmics_liver_v1.3/4_Hotspot/hotspot_liverstomics_NI.ipynb" --untitled3 "$UNTITLED3" \
    --auto30 "$OUT"/reps/n30_rep*.csv "$A1/hotspot_results/hotspot_autocorrelation.csv" \
             "$OUT/a2/hotspot_results/hotspot_autocorrelation.csv" "$OUT/a3/hotspot_results/hotspot_autocorrelation.csv" \
             "$OUT/a_seeded/hotspot_results/hotspot_autocorrelation.csv" \
    --auto15 "$OUT"/reps/n15_rep*.csv "$A1/n15/hotspot_autocorrelation.csv" \
             "$OUT/a2/n15/hotspot_autocorrelation.csv" "$OUT/a3/n15/hotspot_autocorrelation.csv" \
             "$OUT/a_seeded/n15/hotspot_autocorrelation.csv" \
    --full "$A1" "$OUT/a2" "$OUT/a3" "$OUT/a_seeded" \
    --json "$T/logs/compare_replicates_${SLURM_JOB_ID}.json" > "$T/logs/compare_replicates_${SLURM_JOB_ID}.txt" 2>&1
echo done
