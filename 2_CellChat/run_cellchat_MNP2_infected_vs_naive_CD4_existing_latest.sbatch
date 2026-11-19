#!/bin/bash
#SBATCH --job-name=cellchat_mnp2_cd4_existing_latest
#SBATCH --account=punim1236
#SBATCH --partition=cascade
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_mnp2_cd4_existing_latest_%j.out
#SBATCH --error=/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_mnp2_cd4_existing_latest_%j.err

set -euo pipefail
export CELLCHAT_LIB=/data/gpfs/projects/punim2394/AAA_Sorry_Its_Chen_Here/Rlib_new
export SEURAT_FALLBACK_LIB=/data/scratch/projects/punim1236/python_envs/rseurat/lib/R/library
export LD_LIBRARY_PATH=/data/scratch/projects/punim1236/python_envs/rseurat/lib:${LD_LIBRARY_PATH:-}
export TMPDIR=/data/scratch/projects/punim1236/AAA_Chen_is_here/tmp
/data/gpfs/projects/punim2394/conda/rctd_env/bin/Rscript /data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/run_cellchat_MNP2_infected_vs_naive_CD4_latest.R
