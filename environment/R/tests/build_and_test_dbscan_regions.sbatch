#!/bin/bash
# =============================================================================
# build_and_test_dbscan_regions.sbatch
#
# Builds the dbscan_regions R library from scratch on node-local /tmp (R 4.5.0 module, packages
# from dbscan_regions_packages.tsv with install_dbscan_regions_env.R), then runs test_dbscan_regions.R
# on the real data: copies of liverstomics_hotspot_dbscan_NI.R and colocalisation_script_NI.R with
# only paths changed; DBSCAN selection, DBSCAN, region membership, colocalisation and enrichment;
# comparison with the archived tables and objects.
# Everything is written under ${WORK} (removed at the end); small result files go to ${RESULTS}.
#
# Submit from this folder:   cd tests && sbatch build_and_test_dbscan_regions.sbatch
# =============================================================================
#SBATCH --job-name=renv_dbscan_regions
#SBATCH --account=uom00060
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --tmp=30G
#SBATCH --time=03:00:00
#SBATCH --export=NONE
#SBATCH --output=logs/%x_%j.out

# ----------------------------------------------------------------------------- paths (edit here)
STAGING=${STAGING:-$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)}          # folder with install_*.R and *.tsv
ARCHIVE=${ARCHIVE:-$(cd "${STAGING}/../.." && pwd)}   # root of this repository
RCTD_RDS=/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/rctd_MNP1_2_Tcell_combined_no_sinusoids.rds
HOTSPOT_CSV=/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_results/hotspot_spot_scores.csv
ARCHIVED_OUT=/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot
WORK=/tmp/renv_dbscan_regions_${SLURM_JOB_ID}
RESULTS=${STAGING}/tests/results/dbscan_regions

# ----------------------------------------------------------------------------- environment
source /etc/profile.d/*.sh >/dev/null 2>&1
module purge >/dev/null 2>&1
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0
unset PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
[ -n "${SLURM_JOB_ID:-}" ] || { echo "run with sbatch"; exit 1; }
mkdir -p ${WORK}/tmp ${RESULTS}
trap 'rm -rf "${WORK}"' EXIT
export TMPDIR=${WORK}/tmp
export LANG=en_AU.UTF-8; unset LC_ALL
ulimit -c 0
set -u
T=${STAGING}/tests
LIB=${WORK}/Rlib
echo "job ${SLURM_JOB_ID} on $(hostname), $(date), CPUs ${SLURM_CPUS_PER_TASK}"
echo "STAGING=${STAGING}"; echo "ARCHIVE=${ARCHIVE}"
which Rscript; Rscript --version

# ----------------------------------------------------------------------------- 1. R library
mkdir -p ${LIB} ${WORK}/build
cd ${WORK}/build
t0=$(date +%s)
LOCK_WORK_DIR=${WORK}/build/lock_build Rscript ${STAGING}/install_dbscan_regions_env.R ${LIB} ${SLURM_CPUS_PER_TASK} \
  > ${RESULTS}/install_dbscan_regions_env.log 2>&1
rc=$?
t1=$(date +%s)
echo "BUILD R library: exit ${rc}, $(( (t1 - t0) / 60 )) min $(( (t1 - t0) % 60 )) s"
cp dbscan_regions_lock_check.tsv ${RESULTS}/ 2>/dev/null
for f in ${WORK}/build/lock_build/install_logs/*.out; do
  if grep -qE "^ERROR|non-zero exit status" "$f"; then mkdir -p ${RESULTS}/failed_build_logs; cp "$f" ${RESULTS}/failed_build_logs/; fi
done
[ ${rc} -eq 0 ] || { echo "library build failed"; exit 1; }

# ----------------------------------------------------------------------------- 2. test
export R_LIBS_USER=${LIB} R_LIBS_SITE=${LIB}
OUT=${WORK}/test_dbscan_regions
mkdir -p ${OUT}
s=$(date +%s)
srun --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} --mem=40G --export=ALL --job-name=C_dbscan_regions /bin/bash -c \
  "ulimit -c 0; cd ${OUT} && /usr/bin/time -v Rscript ${T}/test_dbscan_regions.R ${ARCHIVE} ${RCTD_RDS} ${HOTSPOT_CSV} ${ARCHIVED_OUT} ${OUT}" \
  > ${RESULTS}/C_dbscan_regions.log 2>&1
rc=$?; e=$(date +%s)
echo "TEST C_dbscan_regions: exit ${rc}, $(( (e - s) / 60 )) min $(( (e - s) % 60 )) s"
diff ${ARCHIVE}/3_STOmics_liver_v1.3/5_DBSCAN_regions_enrichment/liverstomics_hotspot_dbscan_NI.R \
  ${OUT}/liverstomics_hotspot_dbscan_NI_copy.R > ${RESULTS}/C_script_copy.diff
cp ${OUT}/compare_tables.tsv ${OUT}/compare_objects.txt ${OUT}/sessionInfo_dbscan_regions.txt ${RESULTS}/ 2>/dev/null
mkdir -p ${RESULTS}/tables
cp ${OUT}/enrichment/*.csv ${OUT}/rho_*.csv ${RESULTS}/tables/ 2>/dev/null
cp ${OUT}/background_regions/rho_*.csv ${RESULTS}/tables/ 2>/dev/null
cat ${RESULTS}/compare_objects.txt

sacct -j ${SLURM_JOB_ID} --format=JobID,JobName%18,Elapsed,MaxRSS,State,ExitCode -P > ${RESULTS}/sacct_steps.txt
cat ${RESULTS}/sacct_steps.txt
echo "finished $(date)"
