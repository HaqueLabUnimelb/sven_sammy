#!/bin/bash
# Standard Cell Ranger 8.0.1 workflow for the scRNA-seq library of part 1: BCL files to FASTQ files split
# by sample index, then counting of the gene expression and hashtag (Antibody Capture) libraries.
#
# Written for this repository (20 November 2026). The commands and settings are the ones Cell Ranger
# recorded for the run of November 2025 (2025-11_mkfastq_and_count/Split_Output_new/_invocation and
# new_cellranger/_cmdline); only the paths are variables here. Reference: refdata-gex-mm10-2020-A
# (mouse GRCm38, GENCODE vM23).
#
# Spartan: module load CellRanger/8.0.1 GCC/11.3.0 bcl2fastq2/2.20.0
# Run in an empty working folder:  bash /path/to/run_cellranger.sh

set -euo pipefail

BCL_RUN=/path/to/250527_VH00915_604_AAH2VJHM5-10042-micheleC   # sequencing run folder (BCL files)
REFERENCE=/path/to/refdata-gex-mm10-2020-A                      # 10x Genomics mouse reference mm10 2020-A
CORES=24
MEM_GB=128

RECORDS=$(cd "$(dirname "$0")" && pwd)/2025-11_mkfastq_and_count

# 1. BCL to FASTQ. Sample indices: SI-TT-C11 = gene expression (GEX_C11), SI-NT-B1 = hashtags (CSP_B1)
cellranger mkfastq --id=Split_Output_new --run="$BCL_RUN" --csv="$RECORDS/samplesheet.csv" \
    --localcores=$CORES --localmem=$MEM_GB

# 2. Count both libraries; the feature reference lists the 8 hashtag oligos
FASTQ_DIR=$PWD/Split_Output_new/outs/fastq_path/AAH2VJHM5
cat > library.csv <<EOF
fastqs,sample,library_type
$FASTQ_DIR,GEX_C11,Gene Expression
$FASTQ_DIR,CSP_B1,Antibody Capture
EOF
cellranger count --id=new_cellranger --libraries=library.csv --transcriptome="$REFERENCE" \
    --feature-ref="$RECORDS/new_cellranger/outs/feature_reference.csv" \
    --create-bam=false --localcores=$CORES --localmem=$MEM_GB
