#!/bin/bash
#SBATCH --job-name=ss_hotspot_env
#SBATCH --account=uom00060
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --tmp=100G
#SBATCH --time=05:00:00
#SBATCH --output=logs/hotspot_env_test_%j.log

# Builds sven_sammy_hotspot on the node-local disk, records it, and runs
#   test c: the pandas command of 2_CellChat/commands_for_two_extra_tables (copy with the
#           renamed input folder and the output in /tmp), compared byte for byte;
#   test a: hotspot_test_a.py (the notebook's Hotspot code, n_neighbors=30, plus the
#           autocorrelation step with n_neighbors=15), compared with the archived tables;
#   test b: a copy of the notebook executed with jupyter nbconvert (input path, os.chdir
#           target and n_neighbors=30 changed), compared with the archived tables and
#           with the tables of test a;
#   lock check: a second environment built from the lock file and the pip freeze list.
# Tests a and b run at the same time.

# ---- paths: edit these ----------------------------------------------------------
STAGING=$(cd "${SLURM_SUBMIT_DIR:-.}/.." && pwd)   # folder with the YAML files; submit this script from STAGING/tests
REPO=$(cd "$STAGING/../.." && pwd)     # root of this repository
IN_H5AD=/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.h5ad
ARCHIVED_HOTSPOT=/data/projects/punim1236/liver_stomicsv1.3/hotspot
UNTITLED3=/data/projects/punim1236/b_cells/Untitled3.ipynb
CELLCHAT_RENAMED_DIR=/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/MNP2_CD4_20261023
KEEP_ROOT=/path/to/kept_outputs
NAME=sven_sammy_hotspot
JOBS_A=23            # Hotspot worker processes in test a (test b keeps the notebook's jobs=8)
# ---------------------------------------------------------------------------------

for f in /etc/profile.d/*.sh; do source "$f"; done
module purge
module load Miniforge3/24.7.1-2
set -uo pipefail

TESTS=$STAGING/tests
LOGS=$TESTS/logs
ID=${SLURM_JOB_ID}
P=$LOGS/hotspot_${ID}
WORK=/tmp/ss_hotspot_${ID}
KEEP=$KEEP_ROOT/hotspot_${ID}
NB=$REPO/3_STOmics_liver_v1.3/4_Hotspot/hotspot_liverstomics_NI.ipynb
CMD=$REPO/2_CellChat/commands_for_two_extra_tables/make_MNP2_vs_original_combinedMNP_comparison.sh
CMD_EXPECTED=$REPO/2_CellChat/MNP2_CD4_20261023/MNP2_vs_original_combinedMNP_comparison.csv
source "$TESTS/env_common.sh"
isolate "$WORK"
trap 'rm -rf "$WORK"' EXIT
TIMES=${P}_times.txt
: > "$TIMES"
ENV=$WORK/envs/$NAME
CONDA=$EBROOTMINIFORGE3/bin/conda

echo "job $ID on $(hostname), $(date)"
lscpu | grep -E "^Model name|^CPU\(s\)"
echo "conda: $("$CONDA" --version)"

# ---- build ----------------------------------------------------------------------
timed build "$TIMES" "$CONDA" env create -q -p "$ENV" -f "$STAGING/$NAME.yml" > "${P}_build.log" 2>&1 \
    || { echo "BUILD FAILED"; tail -50 "${P}_build.log"; exit 1; }
record_env "$ENV" "$NAME" "$STAGING" "$LOGS"
export PATH=$ENV/bin:$PATH
echo "python: $(command -v python) $(python -c 'import sys; print(sys.version)')"
python - <<'PY'
import sys
print("sys.path:", sys.path)
import numpy, pandas, scipy, sklearn, numba, anndata, scanpy, hotspot, statsmodels, h5py, matplotlib, nbconvert, ipykernel
for m in (numpy, pandas, scipy, sklearn, numba, anndata, scanpy, hotspot, statsmodels, h5py, matplotlib, nbconvert, ipykernel):
    print(m.__name__, getattr(m, "__version__", "?"), m.__file__)
numpy.show_config()
PY
jupyter kernelspec list

# ---- test c: pandas command -----------------------------------------------------
mkdir -p "$WORK/test_c"
python "$TESTS/make_command_copy.py" "$CMD" "$WORK/test_c/make_MNP2_vs_original_combinedMNP_comparison_test.sh" \
    "/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_latest_MNP2_CD4_20261023_matched_original_params/03_differential_LR_CD4_to_MNP2_infected_minus_naive.csv" \
    "$CELLCHAT_RENAMED_DIR/03_differential_LR_CD4_to_MNP2_infected_minus_naive.csv" \
    "out=Path('/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_latest_MNP2_CD4_20261023_matched_original_params/MNP2_vs_original_combinedMNP_comparison.csv')" \
    "out=Path('$WORK/test_c/MNP2_vs_original_combinedMNP_comparison.csv')" \
    > "${P}_test_c_command_diff.txt"
( cd "$WORK/test_c" && timed test_c "$TIMES" bash make_MNP2_vs_original_combinedMNP_comparison_test.sh > "${P}_test_c_stdout.txt" 2>&1 )
if cmp "$WORK/test_c/MNP2_vs_original_combinedMNP_comparison.csv" "$CMD_EXPECTED"; then
    echo "TEST C: output identical byte for byte to $CMD_EXPECTED" | tee -a "$TIMES"
else
    echo "TEST C: output DIFFERS from $CMD_EXPECTED" | tee -a "$TIMES"
fi
md5sum "$WORK/test_c/MNP2_vs_original_combinedMNP_comparison.csv" "$CMD_EXPECTED" | tee -a "${P}_test_c_stdout.txt"

# ---- tests a and b --------------------------------------------------------------
mkdir -p "$WORK/test_a" "$WORK/test_b/nb" "$WORK/test_b/out"
python "$TESTS/make_notebook_copy.py" "$NB" "$WORK/test_b/nb/hotspot_liverstomics_NI_test.ipynb" "$IN_H5AD" "$WORK/test_b/out" \
    > "${P}_test_b_notebook_diff.txt" 2>&1 || { echo "notebook copy failed"; cat "${P}_test_b_notebook_diff.txt"; }
(
    cd "$WORK/test_b/nb" && timed test_b "$TIMES" jupyter nbconvert --to notebook --execute \
        --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=python3 \
        --output hotspot_liverstomics_NI_test_executed.ipynb hotspot_liverstomics_NI_test.ipynb \
        > "${P}_test_b_nbconvert.log" 2>&1
) &
PID_B=$!
timed test_a "$TIMES" python "$TESTS/hotspot_test_a.py" "$IN_H5AD" "$WORK/test_a" "$JOBS_A" > "${P}_test_a.log" 2>&1
wait $PID_B
echo "tests a and b finished $(date)"

# ---- comparisons ----------------------------------------------------------------
RT=$REPO/3_STOmics_liver_v1.3/4_Hotspot/result_tables
python "$TESTS/compare_hotspot.py" --new "$WORK/test_a" --archived "$ARCHIVED_HOTSPOT" --repo-tables "$RT" \
    --notebook "$NB" --untitled3 "$UNTITLED3" --json "${P}_test_a_vs_archived.json" > "${P}_test_a_vs_archived.txt" 2>&1
python "$TESTS/compare_hotspot.py" --new "$WORK/test_b/out" --archived "$ARCHIVED_HOTSPOT" --repo-tables "$RT" \
    --notebook "$NB" --json "${P}_test_b_vs_archived.json" > "${P}_test_b_vs_archived.txt" 2>&1
{
    echo "# test b tables against test a tables (cmp, byte for byte)"
    for d in hotspot_results hotspot_plots; do
        for f in $(ls "$WORK/test_b/out/$d" | grep -E '\.(csv|txt)$'); do
            if [[ -e "$WORK/test_a/$d/$f" ]]; then
                if cmp -s "$WORK/test_b/out/$d/$f" "$WORK/test_a/$d/$f"; then r=identical; else r=DIFFERENT; fi
            else
                r="not written by test a"
            fi
            echo "$d/$f $r"
        done
    done
    echo "# files written by test b"
    (cd "$WORK/test_b" && find . -type f -printf '%P %s\n' | sort)
    echo "# md5 of test b tables"
    (cd "$WORK/test_b/out" && md5sum hotspot_results/*.csv hotspot_plots/*.csv hotspot_plots/*.txt)
    echo "# md5 of test a tables"
    (cd "$WORK/test_a" && md5sum hotspot_results/*.csv hotspot_plots/*.csv hotspot_plots/*.txt n15/*.csv)
} > "${P}_test_b_vs_test_a.txt" 2>&1
python - "$WORK/test_b/nb/hotspot_liverstomics_NI_test_executed.ipynb" > "${P}_test_b_executed_summary.txt" 2>&1 <<'PY'
import json, sys
nb = json.load(open(sys.argv[1]))
code = [c for c in nb["cells"] if c["cell_type"] == "code"]
errors = [o for c in code for o in c.get("outputs", []) if o.get("output_type") == "error"]
print("code cells:", len(code), "executed:", sum(1 for c in code if c.get("execution_count")), "error outputs:", len(errors))
for c in code:
    src = "".join(c["source"])
    if "n_neighbors" in src or "pip install" in src or "print(" in src:
        print("--- cell", c.get("execution_count"), "|", src.splitlines()[0] if src else "")
        for o in c.get("outputs", []):
            t = "".join(o.get("text", "")) or "".join(o.get("data", {}).get("text/plain", ""))
            print(t[:3000])
PY

# ---- keep the small outputs for later inspection --------------------------------
mkdir -p "$KEEP/test_a" "$KEEP/test_b"
cp -r "$WORK/test_a/hotspot_results" "$WORK/test_a/hotspot_plots" "$WORK/test_a/n15" "$KEEP/test_a/"
cp "$WORK/test_b/nb/hotspot_liverstomics_NI_test.ipynb" "$WORK/test_b/nb/hotspot_liverstomics_NI_test_executed.ipynb" "$KEEP/test_b/"
cp "$WORK/test_c/MNP2_vs_original_combinedMNP_comparison.csv" "$KEEP/"

# ---- lock check -----------------------------------------------------------------
t0=$SECONDS
lock_check "$WORK/envs/${NAME}_lockcheck" "$NAME" "$STAGING" "$LOGS"
echo "lockcheck elapsed=$((SECONDS - t0))s" | tee -a "$TIMES"
echo "done $(date)"
