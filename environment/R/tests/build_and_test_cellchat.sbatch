#!/bin/bash
# =============================================================================
# build_and_test_cellchat.sbatch
#
# Builds the CellChat environment from scratch on node-local /tmp (conda R 4.5.1 from
# cellchat_r.yml, then the R library from cellchat_packages.tsv with install_cellchat_env.R) and
# runs copies of the two archived CellChat scripts, in which only input_rds and outdir are changed,
# on the real input objects. Tables 01 to 04 are compared with the archived tables.
# Everything is written under ${WORK} (node-local /tmp, removed at the end); small result files
# (logs, tables, comparison reports, sessionInfo) are copied to ${RESULTS}.
#
# Submit from this folder:   cd tests && sbatch build_and_test_cellchat.sbatch
# =============================================================================
#SBATCH --job-name=renv_cellchat
#SBATCH --account=punim2394
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --tmp=60G
#SBATCH --time=05:00:00
#SBATCH --export=NONE
#SBATCH --output=logs/%x_%j.out

# ----------------------------------------------------------------------------- paths (edit here)
STAGING=${STAGING:-$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)}          # folder with install_*.R, *.tsv, *.yml
ARCHIVE=${ARCHIVE:-$(cd "${STAGING}/../.." && pwd)}   # root of this repository
INPUT_RUN1=/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS
INPUT_RUN2=/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS
CONDA_MODULE=Miniforge3/24.7.1-2
WORK=/tmp/renv_cellchat_${SLURM_JOB_ID}
RESULTS=${STAGING}/tests/results/cellchat
WORKERS=8      # CellChat workers per run, as in the archived runs (SLURM_CPUS_PER_TASK=8)
STEP_MEM=30G   # memory of each test step

# ----------------------------------------------------------------------------- environment
source /etc/profile.d/*.sh >/dev/null 2>&1
module purge >/dev/null 2>&1
module load ${CONDA_MODULE}
unset PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
[ -n "${SLURM_JOB_ID:-}" ] || { echo "run with sbatch"; exit 1; }
mkdir -p ${WORK}/tmp ${RESULTS}
trap 'rm -rf "${WORK}"' EXIT
export TMPDIR=${WORK}/tmp
export CONDA_PKGS_DIRS=${WORK}/conda_pkgs
export CONDA_SOLVER=libmamba
export CONDA_CHANNEL_PRIORITY=strict
export LANG=en_AU.UTF-8; unset LC_ALL   # locale of the archived runs (sessionInfo)
ulimit -c 0
T=${STAGING}/tests
LIB=${WORK}/Rlib
PREFIX=${WORK}/r451
echo "job ${SLURM_JOB_ID} on $(hostname), $(date), CPUs ${SLURM_CPUS_PER_TASK}"
echo "STAGING=${STAGING}"; echo "ARCHIVE=${ARCHIVE}"

# ----------------------------------------------------------------------------- 1. conda R 4.5.1
t0=$(date +%s)
conda env create -p ${PREFIX} -f ${STAGING}/cellchat_r.yml -q
rc=$?
t1=$(date +%s)
echo "BUILD conda env: exit ${rc}, $(( (t1 - t0) / 60 )) min $(( (t1 - t0) % 60 )) s"
[ ${rc} -eq 0 ] || exit 1
conda list -p ${PREFIX} --explicit --md5 > ${RESULTS}/cellchat_r_conda_explicit.txt
conda list -p ${PREFIX} --show-channel-urls > ${RESULTS}/cellchat_r_conda_list.txt
eval "$(conda shell.bash hook)"
conda activate ${PREFIX}
which Rscript; Rscript --version

set -u
# ----------------------------------------------------------------------------- 2. R library
mkdir -p ${LIB} ${WORK}/build
cd ${WORK}/build
t0=$(date +%s)
LOCK_WORK_DIR=${WORK}/build/lock_build Rscript ${STAGING}/install_cellchat_env.R ${LIB} ${SLURM_CPUS_PER_TASK} \
  > ${RESULTS}/install_cellchat_env.log 2>&1
rc=$?
t1=$(date +%s)
echo "BUILD R library: exit ${rc}, $(( (t1 - t0) / 60 )) min $(( (t1 - t0) % 60 )) s"
cp cellchat_lock_check.tsv ${RESULTS}/ 2>/dev/null
for f in ${WORK}/build/lock_build/install_logs/*.out; do
  if grep -qE "^ERROR|non-zero exit status" "$f"; then mkdir -p ${RESULTS}/failed_build_logs; cp "$f" ${RESULTS}/failed_build_logs/; fi
done
[ ${rc} -eq 0 ] || { echo "library build failed"; exit 1; }

# ----------------------------------------------------------------------------- 3. CellChat runs
export R_LIBS_USER=${LIB} R_LIBS_SITE=${LIB}
run_test () {   # $1 = run name (infected | cytotoxic), $2 = input object, $3 = archived script
  local run=$1 input=$2 script=$3 d=${WORK}/test_$1 r=${RESULTS}/$1
  mkdir -p ${d}/out ${r}
  Rscript ${T}/test_cellchat.R copy ${run} ${ARCHIVE} ${input} ${d}/out ${d}/script_copy.R > ${r}/copy.log 2>&1
  diff ${ARCHIVE}/${script} ${d}/script_copy.R > ${r}/script_copy.diff
  local s=$(date +%s)
  ( cd ${d} && srun --ntasks=1 --cpus-per-task=${WORKERS} --mem=${STEP_MEM} --export=ALL /bin/bash -c \
      "ulimit -c 0; SLURM_CPUS_PER_TASK=${WORKERS} CELLCHAT_LIB=${LIB} SEURAT_FALLBACK_LIB=${LIB} /usr/bin/time -v Rscript ${d}/script_copy.R" \
      > ${r}/run.out 2> ${r}/run.err )
  local rc=$?
  local e=$(date +%s)
  echo "TEST cellchat ${run}: exit ${rc}, $(( (e - s) / 60 )) min $(( (e - s) % 60 )) s"
  Rscript ${T}/test_cellchat.R compare ${run} ${ARCHIVE} ${d}/out ${r}/compare_tables.tsv > ${r}/compare.log 2>&1
  cat ${r}/compare.log
  cp ${d}/out/0*.csv ${d}/out/00_run_config.txt ${d}/out/99_sessionInfo.txt ${r}/ 2>/dev/null
}
run_test infected  ${INPUT_RUN1} 2_CellChat/run_cellchat_MNP2_infected_vs_naive_CD4_latest.R
run_test cytotoxic ${INPUT_RUN2} 2_CellChat/run_cellchat_MNP2_cytotoxic_vs_naive_CD4_latest.R

sacct -j ${SLURM_JOB_ID} --format=JobID,JobName%12,Elapsed,MaxRSS,State,ExitCode -P > ${RESULTS}/sacct_steps.txt
cat ${RESULTS}/sacct_steps.txt
echo "finished $(date)"
