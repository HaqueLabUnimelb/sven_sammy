#!/bin/bash
#SBATCH --job-name=ss_compare_rerun
#SBATCH --account=uom00060
#SBATCH --partition=sapphire
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=logs/compare_rerun_%j.log

# Re-runs compare_hotspot.py on the tables of test a of job 35423407 (kept outside /tmp),
# after the comparison of the ranked gene tables in compare_hotspot.py was made robust to
# tables of different length (in job 35423407 it stopped at that step). Needs only numpy
# and pandas; it uses the Python of the module Anaconda3/2022.10 (pandas 1.4.4).
# Job 35425958 ran this command (from a copy of this script with literal paths, which
# also made a trial run of compare_replicates.py on duplicated inputs; that output was
# not kept).

# ---- paths: edit these ----------------------------------------------------------
STAGING=$(cd "${SLURM_SUBMIT_DIR:-.}/.." && pwd)   # submit this script from STAGING/tests
REPO=$(cd "$STAGING/../.." && pwd)
ARCHIVED_HOTSPOT=/data/projects/punim1236/liver_stomicsv1.3/hotspot
UNTITLED3=/data/projects/punim1236/b_cells/Untitled3.ipynb
TEST_A=/path/to/kept_outputs/hotspot_35423407/test_a
# ---------------------------------------------------------------------------------

for f in /etc/profile.d/*.sh; do source "$f"; done
module purge
module load Anaconda3/2022.10
unset PYTHONPATH
export PYTHONNOUSERSITE=1
T=$STAGING/tests
python "$T/compare_hotspot.py" --new "$TEST_A" --archived "$ARCHIVED_HOTSPOT" \
    --repo-tables "$REPO/3_STOmics_liver_v1.3/4_Hotspot/result_tables" \
    --notebook "$REPO/3_STOmics_liver_v1.3/4_Hotspot/hotspot_liverstomics_NI.ipynb" --untitled3 "$UNTITLED3" \
    --json "$T/logs/hotspot_35423407_test_a_vs_archived_rerun.json" > "$T/logs/hotspot_35423407_test_a_vs_archived_rerun.txt" 2>&1
echo done
