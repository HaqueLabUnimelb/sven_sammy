#!/bin/bash
#SBATCH --job-name=ss_stereopy_env
#SBATCH --account=punim1445
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --tmp=60G
#SBATCH --time=03:00:00
#SBATCH --output=logs/stereopy_env_test_%j.log

# Builds an environment from a YAML file on the node-local disk, records it, runs a copy
# of 01_Y01453J8_bin20_gef_to_h5ad.py in which only the output folder is changed (to
# /tmp), and compares the h5ad it writes with the archived h5ad. Then rebuilds the
# environment from its lock file and pip freeze list (lock check).
#
# Usage: sbatch stereopy_env_test.sbatch [YAML_NAME [LABEL "PKG==VER ..."]]
#   YAML_NAME  environment file without .yml, in STAGING (default: sven_sammy_stereopy)
#              or in STAGING/tests.
#   LABEL and the package list are for side checks of other versions: after the build,
#   the listed packages are installed with pip --no-deps over the environment, the
#   records go to tests/logs/LABEL, and the lock check is skipped.

# ---- paths: edit these ----------------------------------------------------------
STAGING=$(cd "${SLURM_SUBMIT_DIR:-.}/.." && pwd)   # folder with the YAML files; submit this script from STAGING/tests
REPO=$(cd "$STAGING/../.." && pwd)     # root of this repository
GEF=/data/projects/punim1236/output_liver_stomics/Y01453J8/outs/feature_expression/Y01453J8.tissue.gef
ARCHIVED_H5AD=/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.h5ad
NAME=${1:-sven_sammy_stereopy}
LABEL=${2:-}
OVERRIDE=${3:-}
# ---------------------------------------------------------------------------------

for f in /etc/profile.d/*.sh; do source "$f"; done
module purge
module load Miniforge3/24.7.1-2
set -uo pipefail

TESTS=$STAGING/tests
LOGS=$TESTS/logs
ID=${SLURM_JOB_ID}
P=$LOGS/stereopy_${ID}
WORK=/tmp/ss_stereopy_${ID}
SCRIPT=$REPO/3_STOmics_liver_v1.3/2_bin20_matrix/01_Y01453J8_bin20_gef_to_h5ad.py
YML=$STAGING/$NAME.yml
[[ -e "$YML" ]] || YML=$TESTS/$NAME.yml
source "$TESTS/env_common.sh"
isolate "$WORK"
trap 'rm -rf "$WORK"' EXIT
TIMES=${P}_times.txt
: > "$TIMES"
ENV=$WORK/envs/$NAME
CONDA=$EBROOTMINIFORGE3/bin/conda

echo "job $ID on $(hostname), $(date); environment file $YML"
lscpu | grep -E "^Model name|^CPU\(s\)"
echo "conda: $("$CONDA" --version)"
[[ "$(grep -c 'output_dir = Path(' "$SCRIPT")" == 1 ]] || { echo "unexpected script"; exit 1; }

# ---- build ----------------------------------------------------------------------
timed build "$TIMES" "$CONDA" env create -q -p "$ENV" -f "$YML" > "${P}_build.log" 2>&1 \
    || { echo "BUILD FAILED"; tail -80 "${P}_build.log"; exit 1; }
if [[ -n "$OVERRIDE" ]]; then
    echo "side check $LABEL: pip install --no-deps $OVERRIDE"
    "$ENV/bin/python" -m pip install --no-deps $OVERRIDE > "${P}_override.log" 2>&1 || { echo "OVERRIDE FAILED"; cat "${P}_override.log"; exit 1; }
fi
if [[ -z "$LABEL" && "$YML" == "$STAGING/$NAME.yml" ]]; then
    record_env "$ENV" "$NAME" "$STAGING" "$LOGS"
else
    mkdir -p "$LOGS/${LABEL:-$NAME}"
    record_env "$ENV" "$NAME" "$LOGS/${LABEL:-$NAME}" "$LOGS/${LABEL:-$NAME}"
fi
"$ENV/bin/python" - <<'PY'
import sys
print(sys.version)
print("sys.path:", sys.path)
import stereo, gefpy, anndata, numpy, scipy, pandas, h5py
for m in (stereo, gefpy, anndata, numpy, scipy, pandas, h5py):
    print(m.__name__, getattr(m, "__version__", "?"), m.__file__)
PY

# ---- run a copy of the script with only the output folder changed ---------------
mkdir -p "$WORK/run" "$WORK/out"
"$ENV/bin/python" "$TESTS/make_command_copy.py" "$SCRIPT" "$WORK/run/01_Y01453J8_bin20_gef_to_h5ad_test.py" \
    $'output_dir = Path(\n    "/data/gpfs/projects/punim1236/"\n    "liver_curio/2026714RCTDnew/bin20"\n)' \
    $'output_dir = Path(\n    "'"$WORK/out"$'"\n)' \
    > "${P}_script_diff.txt"
"$ENV/bin/python" -c "import ast, sys; t = ast.parse(open(sys.argv[1]).read()); v = [ast.literal_eval(n.value) for n in t.body if isinstance(n, ast.Assign) and getattr(n.targets[0], 'id', '') == 'gef_file'][0]; print('script input:', v); sys.exit(v != sys.argv[2])" \
    "$WORK/run/01_Y01453J8_bin20_gef_to_h5ad_test.py" "$GEF" || echo "WARNING: the script reads another input than $GEF"
( cd "$WORK/run" && timed run_script "$TIMES" "$ENV/bin/python" 01_Y01453J8_bin20_gef_to_h5ad_test.py > "${P}_run.log" 2>&1 )
ls -la "$WORK/out"

# ---- compare --------------------------------------------------------------------
timed compare "$TIMES" "$ENV/bin/python" "$TESTS/compare_h5ad.py" "$WORK/out/Y01453J8_tissue_bin20_raw_counts.h5ad" \
    "$ARCHIVED_H5AD" "${P}_compare.json" > "${P}_compare.txt" 2>&1
cat "${P}_compare.txt"

# ---- lock check (main environment only) -----------------------------------------
if [[ -z "$LABEL" && "$YML" == "$STAGING/$NAME.yml" ]]; then
    t0=$SECONDS
    lock_check "$WORK/envs/${NAME}_lockcheck" "$NAME" "$STAGING" "$LOGS"
    echo "lockcheck elapsed=$((SECONDS - t0))s" | tee -a "$TIMES"
fi
echo "done $(date)"
