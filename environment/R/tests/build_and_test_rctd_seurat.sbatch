#!/bin/bash
# =============================================================================
# build_and_test_rctd_seurat.sbatch
#
# Builds the rctd_seurat R library from scratch on node-local /tmp (R 4.5.0 module, packages from
# rctd_seurat_packages.tsv with install_rctd_seurat_env.R), then runs on the real data:
#   D  part 1: load the packages of the Seurat scripts; FindClusters(resolution = 0.7) on the cells
#      of liver_scRNAseq_processed_NI.R line 65, compared with sven_liver_for_paper.RDS;
#   B1 part 3: copy of 02_Y01453J8_bin20_h5ad_to_rds.R (paths changed), compared with the archived
#      bin20 rds that RCTD read;
#   B2 part 3: copy of rctd_script_NI.R (paths changed), compared with the archived RCTD result.
# The copies of B1 and B2 are run with Rscript and the R profile sessioninfo_at_exit.R, which writes
# sessionInfo() when the session ends.
# Everything is written under ${WORK} (removed at the end); small result files go to ${RESULTS}.
#
# Submit from this folder:   cd tests && sbatch build_and_test_rctd_seurat.sbatch
# =============================================================================
#SBATCH --job-name=renv_rctd_seurat
#SBATCH --account=uom00060
#SBATCH --partition=sapphire
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --tmp=60G
#SBATCH --time=06:00:00
#SBATCH --export=NONE
#SBATCH --output=logs/%x_%j.out

# ----------------------------------------------------------------------------- paths (edit here)
STAGING=${STAGING:-$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)}          # folder with install_*.R and *.tsv
ARCHIVE=${ARCHIVE:-$(cd "${STAGING}/../.." && pwd)}   # root of this repository
GEX_HTO=/data/projects/punim1236/liver_scRNAseq/GEX_HTO_Processed.rds
LIVER_FOR_PAPER=/data/projects/punim1236/liver_scRNAseq/sven_liver_for_paper.RDS
H5AD=/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.h5ad
BIN20_RDS=/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.rds
TABULA_MURIS=/data/gpfs/projects/punim1236/liver_curio/tabula_muris/droplet_Liver_seurat_tiss.Robj
IMMUNE=/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS
RCTD_ARCHIVED=/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/rctd_MNP1_2_Tcell_combined_no_sinusoids.rds
WORK=/tmp/renv_rctd_seurat_${SLURM_JOB_ID}
RESULTS=${STAGING}/tests/results/rctd_seurat

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
LOCK_WORK_DIR=${WORK}/build/lock_build Rscript ${STAGING}/install_rctd_seurat_env.R ${LIB} ${SLURM_CPUS_PER_TASK} \
  > ${RESULTS}/install_rctd_seurat_env.log 2>&1
rc=$?
t1=$(date +%s)
echo "BUILD R library: exit ${rc}, $(( (t1 - t0) / 60 )) min $(( (t1 - t0) % 60 )) s"
cp rctd_seurat_lock_check.tsv ${RESULTS}/ 2>/dev/null
for f in ${WORK}/build/lock_build/install_logs/*.out; do
  if grep -qE "^ERROR|non-zero exit status" "$f"; then mkdir -p ${RESULTS}/failed_build_logs; cp "$f" ${RESULTS}/failed_build_logs/; fi
done
[ ${rc} -eq 0 ] || { echo "library build failed"; exit 1; }

# ----------------------------------------------------------------------------- 2. tests
export R_LIBS_USER=${LIB} R_LIBS_SITE=${LIB}
step () {   # $1 = name, $2 = CPUs, $3 = memory, $4 = working dir, rest = command; output to ${RESULTS}/<name>.log
  local name=$1 cpus=$2 mem=$3 dir=$4; shift 4
  local s=$(date +%s)
  srun --ntasks=1 --cpus-per-task=${cpus} --mem=${mem} --export=ALL --job-name=${name} /bin/bash -c "ulimit -c 0; cd ${dir} && /usr/bin/time -v $*" \
    > ${RESULTS}/${name}.log 2>&1
  local rc=$? e=$(date +%s)
  echo "TEST ${name}: exit ${rc}, $(( (e - s) / 60 )) min $(( (e - s) % 60 )) s"
}

# D: part 1 Seurat scripts
D=${WORK}/test_seurat
mkdir -p ${D}
step D_seurat 4 32G ${D} Rscript ${T}/test_rctd_seurat.R seurat ${ARCHIVE} ${GEX_HTO} ${LIVER_FOR_PAPER} ${D}
cp ${D}/compare_clusters.txt ${D}/sessionInfo_seurat_scripts.txt ${RESULTS}/ 2>/dev/null
diff ${ARCHIVE}/1_scRNAseq_preprocessing_clustering_annotation/2_seurat_clustering_annotation/liver_scRNAseq_processed_NI.R \
  ${D}/liver_scRNAseq_processed_NI_copy.R > ${RESULTS}/D_script_copy.diff
cat ${RESULTS}/compare_clusters.txt

# B1: bin20 h5ad -> rds
B1=${WORK}/test_bin20
mkdir -p ${B1}
Rscript ${T}/test_rctd_seurat.R bin20-copy ${ARCHIVE} ${H5AD} ${B1} ${B1}/bin20_copy.R > ${RESULTS}/B1_copy.log 2>&1
diff ${ARCHIVE}/3_STOmics_liver_v1.3/2_bin20_matrix/02_Y01453J8_bin20_h5ad_to_rds.R ${B1}/bin20_copy.R > ${RESULTS}/B1_script_copy.diff
step B1_bin20_run 4 48G ${B1} env R_PROFILE_USER=${T}/sessioninfo_at_exit.R SESSIONINFO_FILE=${RESULTS}/B1_sessionInfo_bin20_run.txt Rscript ${B1}/bin20_copy.R
step B1_bin20_compare 2 48G ${B1} Rscript ${T}/test_rctd_seurat.R bin20-compare ${B1}/Y01453J8_tissue_bin20_raw_counts.rds ${BIN20_RDS} ${RESULTS}/B1_compare.txt
cat ${RESULTS}/B1_compare.txt
rm -f ${B1}/Y01453J8_tissue_bin20_raw_counts.rds

# B2: RCTD (archived wrapper: 16 CPUs, 100 GB; the script uses max_cores = 4)
B2=${WORK}/test_rctd
mkdir -p ${B2}
Rscript ${T}/test_rctd_seurat.R rctd-copy ${ARCHIVE} ${TABULA_MURIS} ${IMMUNE} ${BIN20_RDS} ${B2} ${B2}/rctd_copy.R > ${RESULTS}/B2_copy.log 2>&1
diff ${ARCHIVE}/3_STOmics_liver_v1.3/3_RCTD/rctd_script_NI.R ${B2}/rctd_copy.R > ${RESULTS}/B2_script_copy.diff
step B2_rctd_run 16 80G ${B2} env R_PROFILE_USER=${T}/sessioninfo_at_exit.R SESSIONINFO_FILE=${RESULTS}/B2_sessionInfo_rctd_run.txt Rscript ${B2}/rctd_copy.R
step B2_rctd_compare 2 32G ${B2} Rscript ${T}/test_rctd_seurat.R rctd-compare ${B2}/rctd_MNP1_2_Tcell_combined_no_sinusoids_TEST.rds ${RCTD_ARCHIVED} ${RESULTS}/B2_compare.txt
cat ${RESULTS}/B2_compare.txt

sacct -j ${SLURM_JOB_ID} --format=JobID,JobName%18,Elapsed,MaxRSS,State,ExitCode -P > ${RESULTS}/sacct_steps.txt
cat ${RESULTS}/sacct_steps.txt
echo "finished $(date)"
