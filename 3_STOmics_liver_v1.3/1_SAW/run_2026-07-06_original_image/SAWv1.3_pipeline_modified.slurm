#!/bin/bash
#SBATCH --job-name=SAW_liver
#SBATCH --partition=sapphire
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/data/projects/punim1236/AAA_Chen_is_here/Ollie_SAW/logs/saw_liver_%j.log
#SBATCH --error=/data/projects/punim1236/AAA_Chen_is_here/Ollie_SAW/logs/saw_liver_%j.err

# Set working directories
WORKDIR=/data/projects/punim1236/AAA_Chen_is_here/Ollie_SAW
SAWDIR=/data/scratch/projects/punim1236/spleen_stomicsv1.3_serial/saw-8.2.2
REFDIR=/data/projects/punim1236/AAA_Chen_is_here/Mus_musculus_index
OUTDIR=${WORKDIR}/output_liver_stomics

# Input directories (from punim1445)
FASTQDIR=/data/scratch/projects/punim1445/liver_stomicsv1.3/raw/fastq
MASKDIR=/data/scratch/projects/punim1445/liver_stomicsv1.3/raw/mask
IMAGEDIR=/data/scratch/projects/punim1445/liver_stomicsv1.3/raw/images

# Sample information
declare -a SAMPLE_IDS=("Y01453J8")

# Image mapping
declare -A IMAGES
IMAGES[Y01453J8]="Y01453J8_DAPI.tif"


# Set PATH to include SAW
export PATH=${SAWDIR}/bin:${PATH}

# Redirect Snakemake cache to work directory (avoid home directory quota issues)
export SNAKEMAKE_OUTPUT_CACHE=${WORKDIR}/.snakemake_cache
export XDG_CACHE_HOME=${WORKDIR}/.cache
export TMPDIR=${WORKDIR}/.tmp
mkdir -p ${SNAKEMAKE_OUTPUT_CACHE} ${XDG_CACHE_HOME} ${TMPDIR}

# Print environment info
echo "================================================"
echo "SAW Stereo-seq Analysis Workflow"
echo "================================================"
echo "Job started at: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Running on node: $(hostname)"
echo "Working directory: ${WORKDIR}"
echo "SAW directory: ${SAWDIR}"
echo "Reference directory: ${REFDIR}"
echo "Output directory: ${OUTDIR}"
echo "FASTQ directory: ${FASTQDIR}"
echo "Snakemake cache: ${SNAKEMAKE_OUTPUT_CACHE}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: 200GB"
echo ""

# Create output and log directories
mkdir -p ${OUTDIR}
mkdir -p ${WORKDIR}/logs
cd ${WORKDIR}

# Verify SAW installation
echo "Verifying SAW installation..."
saw --version
echo ""

# Process each sample
for SN in "${SAMPLE_IDS[@]}"; do
    IMAGE_FILE=${IMAGES[$SN]}
    MASK=${MASKDIR}/${SN}.barcodeToPos.h5
    IMAGE=${IMAGEDIR}/${IMAGE_FILE}

    echo "================================================"
    echo "Processing sample: ${SN}"
    echo "================================================"
    echo "Chip mask: ${MASK}"
    echo "Image: ${IMAGE}"
    echo "FASTQ directory: ${FASTQDIR}"
    echo "Expected FASTQ files: ${SN}_1.fq.gz, ${SN}_2.fq.gz"
    echo ""

    # Check input files
    if [ ! -f "$MASK" ]; then
        echo "ERROR: Mask file not found: ${MASK}"
        continue
    fi

    if [ ! -f "${FASTQDIR}/${SN}_1.fq.gz" ]; then
        echo "ERROR: R1 FASTQ not found: ${FASTQDIR}/${SN}_1.fq.gz"
        continue
    fi

    if [ ! -f "${FASTQDIR}/${SN}_2.fq.gz" ]; then
        echo "ERROR: R2 FASTQ not found: ${FASTQDIR}/${SN}_2.fq.gz"
        continue
    fi

    if [ ! -f "$IMAGE" ]; then
        echo "WARNING: Image file not found: ${IMAGE}"
        echo "Proceeding without image..."
        IMAGE_PARAM=""
    else
        IMAGE_PARAM="--image ${IMAGE}"
    fi

    # Run SAW count
    echo "Starting SAW count for ${SN}..."
    echo ""

    saw count \
        --id ${SN} \
        --sn ${SN} \
        --omics transcriptomics \
        --kit-version "Stereo-seq T FF V1.3" \
        --sequencing-type "PE75_50+100" \
        --chip-mask ${MASK} \
        --organism mouse \
        --tissue liver \
        --reference ${REFDIR} \
        --fastqs ${FASTQDIR} \
        ${IMAGE_PARAM} \
        --output ${OUTDIR} \
        --threads-num ${SLURM_CPUS_PER_TASK} \
        --memory 200 \
        --create-gem \
        --job-mode local

    EXIT_CODE=$?

    if [ $EXIT_CODE -eq 0 ]; then
        echo ""
        echo "SUCCESS: SAW count completed for sample ${SN}"
        echo "Output directory: ${OUTDIR}/${SN}"
    else
        echo ""
        echo "ERROR: SAW count failed for sample ${SN} with exit code ${EXIT_CODE}"
    fi
    echo ""

done

echo "================================================"
echo "All samples processed"
echo "Job completed at: $(date)"
echo "Total runtime: $SECONDS seconds"
echo "================================================"
