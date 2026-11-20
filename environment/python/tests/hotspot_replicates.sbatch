#!/bin/bash
#SBATCH --job-name=ss_hotspot_reps
#SBATCH --account=punim2394
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --tmp=60G
#SBATCH --time=01:00:00
#SBATCH --output=logs/hotspot_replicates_%j.log

# Second set of Hotspot jobs. hotspot_env_test.sbatch showed that two runs of the same
# code in the same environment (tests a and b) never write identical tables. Each part
# builds sven_sammy_hotspot again and writes its tables to OUT:
#   knn_n30   shows the cause (hotspot_knn_check.py: the neighbour graph), then runs the
#             autocorrelation step of the notebook 8 times with n_neighbors=30;
#   n15       runs the autocorrelation step 8 times with n_neighbors=15;
#   seeded    runs the script (test a) and a copy of the notebook (test b, the same three
#             changes as before) with numpy's global random state seeded to 0 and numba
#             limited to one thread, and compares their tables byte for byte;
#   unseeded  runs the script twice more without a seed.
# compare_replicates.sbatch then compares the archived run with all replicates.
#
# Usage, from the tests folder:
#   sbatch --cpus-per-task=12 --mem=48G hotspot_replicates.sbatch knn_n30
#   sbatch --cpus-per-task=10 --mem=40G hotspot_replicates.sbatch n15
#   sbatch hotspot_replicates.sbatch seeded
#   sbatch hotspot_replicates.sbatch unseeded

# ---- paths: edit these ----------------------------------------------------------
STAGING=$(cd "${SLURM_SUBMIT_DIR:-.}/.." && pwd)   # folder with the YAML files; submit this script from STAGING/tests
REPO=$(cd "$STAGING/../.." && pwd)     # root of this repository
IN_H5AD=/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.h5ad
ARCHIVED_HOTSPOT=/data/projects/punim1236/liver_stomicsv1.3/hotspot
UNTITLED3=/data/projects/punim1236/b_cells/Untitled3.ipynb
OUT=/path/to/kept_outputs/hotspot_replicates
NAME=sven_sammy_hotspot
N_REP=8
PART=${1:?give the part: knn_n30, n15, seeded or unseeded}
# ---------------------------------------------------------------------------------

for f in /etc/profile.d/*.sh; do source "$f"; done
module purge
module load Miniforge3/24.7.1-2
set -uo pipefail

TESTS=$STAGING/tests
LOGS=$TESTS/logs
ID=${SLURM_JOB_ID}
P=$LOGS/hotspot_${ID}_${PART}
WORK=/tmp/ss_hotspot_${ID}
NB=$REPO/3_STOmics_liver_v1.3/4_Hotspot/hotspot_liverstomics_NI.ipynb
RT=$REPO/3_STOmics_liver_v1.3/4_Hotspot/result_tables
source "$TESTS/env_common.sh"
isolate "$WORK"
trap 'rm -rf "$WORK"' EXIT
TIMES=${P}_times.txt
: > "$TIMES"
ENV=$WORK/envs/$NAME
CONDA=$EBROOTMINIFORGE3/bin/conda
mkdir -p "$OUT"

echo "job $ID part $PART on $(hostname), $(date), ${SLURM_CPUS_PER_TASK} CPUs"
lscpu | grep -E "^Model name|^CPU\(s\)"

# ---- build (same YAML; records checked against the first build) ------------------
timed build "$TIMES" "$CONDA" env create -q -p "$ENV" -f "$STAGING/$NAME.yml" > "${P}_build.log" 2>&1 \
    || { echo "BUILD FAILED"; tail -50 "${P}_build.log"; exit 1; }
mkdir -p "$WORK/records"
record_env "$ENV" "$NAME" "$WORK/records" "$WORK/records" > /dev/null
if diff "$STAGING/$NAME.pip-freeze.txt" "$WORK/records/$NAME.pip-freeze.txt" > /dev/null && \
   diff <(grep -v '^#' "$STAGING/$NAME.lock.txt") <(grep -v '^#' "$WORK/records/$NAME.lock.txt") > /dev/null; then
    echo "REBUILD: lock file and pip freeze identical to the first build" | tee -a "$TIMES"
else
    echo "REBUILD: records DIFFER from the first build" | tee -a "$TIMES"
    diff "$STAGING/$NAME.pip-freeze.txt" "$WORK/records/$NAME.pip-freeze.txt"
fi
export PATH=$ENV/bin:$PATH

autocorr_reps() {   # autocorr_reps N_NEIGHBORS
    local k=$1 t0=$SECONDS i
    mkdir -p "$OUT/reps"
    for i in $(seq 1 $N_REP); do
        python "$TESTS/hotspot_autocorr_replicate.py" "$IN_H5AD" "$OUT/reps/n${k}_rep${i}_job${ID}.csv" $k 1 \
            > "$WORK/n${k}_rep${i}.log" 2>&1 &
    done
    wait
    echo "autocorrelation replicates n_neighbors=$k: $(ls "$OUT"/reps/n${k}_rep*_job${ID}.csv | wc -l) tables, elapsed=$((SECONDS - t0))s" | tee -a "$TIMES"
    grep -h "_rep" "$WORK"/n${k}_rep*.log | grep -v "it/s" >> "${P}_autocorr_replicates.txt"
}

full_run() {   # full_run LABEL JOBS [env assignments...]
    local label=$1 jobs=$2; shift 2
    mkdir -p "$WORK/$label"
    timed "$label" "$TIMES" env "$@" python "$TESTS/hotspot_test_a.py" "$IN_H5AD" "$WORK/$label" "$jobs" > "${P}_${label}.log" 2>&1
}

keep_and_compare() {   # keep_and_compare LABEL DIR
    local label=$1 dir=$2
    rm -rf "$OUT/$label"; mkdir -p "$OUT/$label"
    cp -r "$dir/hotspot_results" "$dir/hotspot_plots" "$OUT/$label/"
    [[ -d "$dir/n15" ]] && cp -r "$dir/n15" "$OUT/$label/"
    rm -f "$OUT/$label/hotspot_results/adata_with_hotspot.h5ad" "$OUT/$label/hotspot_plots/"*.png "$OUT/$label/hotspot_plots/"*.pdf
    python "$TESTS/compare_hotspot.py" --new "$OUT/$label" --archived "$ARCHIVED_HOTSPOT" --repo-tables "$RT" \
        --notebook "$NB" --untitled3 "$UNTITLED3" --json "${P}_${label}_vs_archived.json" > "${P}_${label}_vs_archived.txt" 2>&1
}

case "$PART" in
knn_n30)
    {
        python "$TESTS/hotspot_knn_check.py" "$IN_H5AD" unseeded
        python "$TESTS/hotspot_knn_check.py" "$IN_H5AD" seeded
        NUMBA_NUM_THREADS=1 python "$TESTS/hotspot_knn_check.py" "$IN_H5AD" seeded
    } > "${P}_knn_check.txt" 2>&1
    grep KNNCHECK "${P}_knn_check.txt" | tee -a "$TIMES"
    autocorr_reps 30
    ;;
n15)
    autocorr_reps 15
    ;;
seeded)
    mkdir -p "$WORK/b_seeded/nb" "$WORK/b_seeded/out"
    python "$TESTS/make_notebook_copy.py" "$NB" "$WORK/b_seeded/nb/hotspot_liverstomics_NI_test.ipynb" "$IN_H5AD" "$WORK/b_seeded/out" \
        > "${P}_notebook_diff.txt" 2>&1
    mkdir -p "$IPYTHONDIR/profile_default/startup"
    cat > "$IPYTHONDIR/profile_default/startup/00-seed.py" <<'PY'
import numpy
numpy.random.seed(0)
print("IPython startup: numpy.random.seed(0)")
PY
    (
        cd "$WORK/b_seeded/nb" && timed b_seeded "$TIMES" env NUMBA_NUM_THREADS=1 jupyter nbconvert --to notebook --execute \
            --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=python3 \
            --output hotspot_liverstomics_NI_test_executed.ipynb hotspot_liverstomics_NI_test.ipynb \
            > "${P}_b_seeded_nbconvert.log" 2>&1
    ) &
    full_run a_seeded 8 HOTSPOT_TEST_SEED=0 NUMBA_NUM_THREADS=1
    wait
    {
        echo "# seeded notebook run (b_seeded) against seeded script run (a_seeded), byte for byte"
        for d in hotspot_results hotspot_plots; do
            for f in $(ls "$WORK/b_seeded/out/$d" | grep -E '\.(csv|txt)$'); do
                if [[ -e "$WORK/a_seeded/$d/$f" ]]; then
                    if cmp -s "$WORK/b_seeded/out/$d/$f" "$WORK/a_seeded/$d/$f"; then r=identical; else r=DIFFERENT; fi
                else
                    r="not written by the script"
                fi
                echo "$d/$f $r"
            done
        done
    } > "${P}_b_seeded_vs_a_seeded.txt" 2>&1
    echo "seeded pair: $(grep -c ' identical$' "${P}_b_seeded_vs_a_seeded.txt") files identical, $(grep -c ' DIFFERENT$' "${P}_b_seeded_vs_a_seeded.txt") different" | tee -a "$TIMES"
    python - "$WORK/b_seeded/nb/hotspot_liverstomics_NI_test_executed.ipynb" > "${P}_b_seeded_executed_summary.txt" 2>&1 <<'PY'
import json, sys
nb = json.load(open(sys.argv[1]))
code = [c for c in nb["cells"] if c["cell_type"] == "code"]
errors = [o for c in code for o in c.get("outputs", []) if o.get("output_type") == "error"]
print("code cells:", len(code), "executed:", sum(1 for c in code if c.get("execution_count")), "error outputs:", len(errors))
for c in code[:1]:
    for o in c.get("outputs", []):
        print("".join(o.get("text", ""))[:300])
PY
    cat "${P}_b_seeded_executed_summary.txt"
    keep_and_compare a_seeded "$WORK/a_seeded"
    keep_and_compare b_seeded "$WORK/b_seeded/out"
    cp "$WORK/b_seeded/nb/hotspot_liverstomics_NI_test_executed.ipynb" "$OUT/b_seeded/"
    ;;
unseeded)
    full_run a2 7 &
    full_run a3 7 &
    wait
    keep_and_compare a2 "$WORK/a2"
    keep_and_compare a3 "$WORK/a3"
    ;;
*)
    echo "unknown part $PART"; exit 1
    ;;
esac
echo "done $(date)"
