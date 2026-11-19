#!/bin/bash
#SBATCH --job-name=Y01453J8_bin20_RCTD
#SBATCH --account=punim1236
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=31-00:00:00
#SBATCH --output=/data/gpfs/projects/punim1236/liver_curio/2026714RCTDnew/bin20/slurm-%x-%j.out
#SBATCH --error=/data/gpfs/projects/punim1236/liver_curio/2026714RCTDnew/bin20/slurm-%x-%j.err

# =============================================================================
# Y01453J8 Stereo-seq bin20 -> H5AD -> Seurat RDS -> RCTD
#
# Step 1:
#   Y01453J8.tissue.gef
#       -> Y01453J8_tissue_bin20_raw_counts.h5ad
#
# Step 2:
#   Y01453J8_tissue_bin20_raw_counts.h5ad
#       -> Y01453J8_tissue_bin20_raw_counts.rds
#
# Step 3:
#   Y01453J8_tissue_bin20_raw_counts.rds
#       -> Y01453J8_tissue_bin20_rctd_results.rds
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------

WORK_DIR="/data/gpfs/projects/punim1236/liver_curio/2026714RCTDnew/bin20"

PYTHON_BIN="/data/scratch/projects/punim1236/AAA_Chen_is_here/conda_envs/stereopy_latest_py38/bin/python3.8"

PYTHON_SCRIPT="${WORK_DIR}/01_Y01453J8_bin20_gef_to_h5ad.py"
H5AD_TO_RDS_SCRIPT="${WORK_DIR}/02_Y01453J8_bin20_h5ad_to_rds.R"
RCTD_SCRIPT="${WORK_DIR}/03_Y01453J8_bin20_run_RCTD.R"

GEF_FILE="/data/projects/punim1236/output_liver_stomics/Y01453J8/outs/feature_expression/Y01453J8.tissue.gef"
H5AD_FILE="${WORK_DIR}/Y01453J8_tissue_bin20_raw_counts.h5ad"
SEURAT_RDS="${WORK_DIR}/Y01453J8_tissue_bin20_raw_counts.rds"
RCTD_RDS="${WORK_DIR}/Y01453J8_tissue_bin20_rctd_results.rds"

# Use the R library in which Seurat, spacexr and related packages are installed.
export R_LIBS_USER="/data/projects/punim1236/liver_stomicsv1.3/RCTD_CLAHE_old/R_library_4.5"

# Prevent BLAS/OpenMP oversubscription. runRctd will use SLURM_CPUS_PER_TASK
# as the number of parallel workers.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

# -----------------------------------------------------------------------------
# Final job report
# -----------------------------------------------------------------------------

finalize_job() {
    exit_status=$?
    trap - EXIT

    echo
    echo "================================================================"
    if [[ ${exit_status} -eq 0 ]]; then
        echo "Pipeline completed successfully."
    else
        echo "Pipeline failed with exit status ${exit_status}."
    fi
    echo "Job ID: ${SLURM_JOB_ID:-not_available}"
    echo "End time: $(date)"
    echo "================================================================"

    if command -v my-job-stats >/dev/null 2>&1 && [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo
        echo "Resource usage:"
        my-job-stats -a -n -s || true
    fi

    exit "${exit_status}"
}

trap finalize_job EXIT

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    echo "ERROR: Submit this file using sbatch, not bash."
    exit 1
fi

echo "================================================================"
echo "Y01453J8 bin20 RCTD pipeline"
echo "================================================================"
echo "Job ID:             ${SLURM_JOB_ID}"
echo "Node:               ${SLURMD_NODENAME:-unknown}"
echo "CPUs:               ${SLURM_CPUS_PER_TASK}"
echo "Memory requested:   250G"
echo "Start time:         $(date)"
echo "Working directory:  ${WORK_DIR}"
echo "R library:          ${R_LIBS_USER}"
echo "================================================================"
echo

# -----------------------------------------------------------------------------
# Check required files
# -----------------------------------------------------------------------------

for required_file in \
    "${GEF_FILE}" \
    "${PYTHON_BIN}" \
    "${PYTHON_SCRIPT}" \
    "${H5AD_TO_RDS_SCRIPT}" \
    "${RCTD_SCRIPT}"
do
    if [[ ! -e "${required_file}" ]]; then
        echo "ERROR: Required file does not exist:"
        echo "  ${required_file}"
        exit 1
    fi
done

# -----------------------------------------------------------------------------
# Step 1: GEF -> bin20 H5AD
# -----------------------------------------------------------------------------

echo "================================================================"
echo "Step 1/3: Convert tissue.gef to bin20 raw-count H5AD"
echo "Start: $(date)"
echo "================================================================"

if [[ -s "${H5AD_FILE}" ]]; then
    echo "Output already exists; skipping Step 1:"
    echo "  ${H5AD_FILE}"
else
    "${PYTHON_BIN}" "${PYTHON_SCRIPT}"

    if [[ ! -s "${H5AD_FILE}" ]]; then
        echo "ERROR: Step 1 finished without creating a valid H5AD file:"
        echo "  ${H5AD_FILE}"
        exit 1
    fi
fi

echo "Step 1 completed: $(date)"
echo

# -----------------------------------------------------------------------------
# Load R environment once for Steps 2 and 3
# -----------------------------------------------------------------------------

module purge
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0

echo "R executable: $(command -v Rscript)"
Rscript --version
echo

# Check required R packages. Do not install packages during the compute job.
Rscript -e '
required <- c(
  "Matrix",
  "hdf5r",
  "Seurat",
  "spacexr",
  "SummarizedExperiment",
  "SpatialExperiment",
  "S4Vectors"
)
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0) {
  stop(
    paste0(
      "Missing R packages in R_LIBS_USER: ",
      paste(missing, collapse = ", ")
    )
  )
}
cat("All required R packages are available.\n")
'

# -----------------------------------------------------------------------------
# Step 2: H5AD -> Seurat RDS
# -----------------------------------------------------------------------------

echo
echo "================================================================"
echo "Step 2/3: Convert bin20 H5AD to Seurat RDS"
echo "Start: $(date)"
echo "================================================================"

if [[ -s "${SEURAT_RDS}" ]]; then
    echo "Output already exists; skipping Step 2:"
    echo "  ${SEURAT_RDS}"
else
    Rscript "${H5AD_TO_RDS_SCRIPT}"

    if [[ ! -s "${SEURAT_RDS}" ]]; then
        echo "ERROR: Step 2 finished without creating a valid Seurat RDS:"
        echo "  ${SEURAT_RDS}"
        exit 1
    fi
fi

echo "Step 2 completed: $(date)"
echo

# -----------------------------------------------------------------------------
# Step 3: Run RCTD
# -----------------------------------------------------------------------------

echo "================================================================"
echo "Step 3/3: Run RCTD on Y01453J8 bin20"
echo "Start: $(date)"
echo "RCTD worker cores: ${SLURM_CPUS_PER_TASK}"
echo "================================================================"

if [[ -s "${RCTD_RDS}" ]]; then
    echo "Final result already exists; skipping Step 3:"
    echo "  ${RCTD_RDS}"
else
    Rscript "${RCTD_SCRIPT}"

    if [[ ! -s "${RCTD_RDS}" ]]; then
        echo "ERROR: Step 3 finished without creating the expected RCTD result:"
        echo "  ${RCTD_RDS}"
        exit 1
    fi
fi

echo "Step 3 completed: $(date)"
echo
echo "Final RCTD result:"
echo "  ${RCTD_RDS}"
