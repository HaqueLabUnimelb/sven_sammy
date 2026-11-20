#!/bin/bash
# =============================================================================
# make_lock_tables.sh: the commands that wrote the three lock tables (*_packages.tsv) and
# cellchat_r.yml on 20 November 2026, on the Spartan login node (run from this folder). Kept for the record; it reads the
# original R libraries and conda environment, which exist only on Spartan.
#
# make_lock.py follows Depends, Imports and LinkingTo from the root packages through the given
# libraries in order (first hit wins, like .libPaths()), drops base packages, and writes
# package, version, source (CRAN / Bioconductor 3.22 / GitHub user/repo@commit, from the
# DESCRIPTION fields Repository, biocViews, RemoteType, RemoteSha) and origin (the library).
# With --si it stops unless every package of the list resolves to the version given there.
# =============================================================================
set -euo pipefail
OUT=..   # the folder with the lock tables

# A. cellchat: roots = every package of 2_CellChat/MNP2_CD4_20261023/99_sessionInfo.txt (attached and
#    loaded namespaces; the two runs list the same packages) plus BiocVersion. Libraries in the order
#    of the run (.libPaths(c(CELLCHAT_LIB, SEURAT_FALLBACK_LIB))); conda R's own library holds only base packages.
ROOTS=$(mktemp)
(cat cellchat_sessionInfo_packages.txt; echo BiocVersion) > ${ROOTS}
python3 -I make_lock.py --roots ${ROOTS} --si cellchat_sessionInfo_packages.txt \
  --libs /data/gpfs/projects/punim2394/AAA_Sorry_Its_Chen_Here/Rlib_new,/data/scratch/projects/punim1236/python_envs/rseurat/lib/R/library \
  --labels Rlib_new,rseurat > ${OUT}/cellchat_packages.tsv
rm -f ${ROOTS}

# B and D. rctd_seurat: roots = packages loaded by rctd_script_NI.R, 02_Y01453J8_bin20_h5ad_to_rds.R
#    (and checked by its wrapper) and the part 1 Seurat scripts, plus BiocVersion, BiocManager and otel
#    (the one package of R_library_4.5 outside that closure). Libraries in the order of the RCTD job:
#    R_LIBS_USER, then the library of the module R/4.5.0. --si checks the ten versions (tidyverse and its nine core packages) that the
#    start-up message printed in the RCTD log (3_RCTD/slurm-32372823.out).
python3 -I make_lock.py --roots rctd_seurat_roots.txt --si rctd_seurat_logged_versions.txt \
  --libs /data/projects/punim1236/liver_stomicsv1.3/RCTD_CLAHE_old/R_library_4.5,/apps/easybuild-2022/easybuild/software/MPI/GCC/11.3.0/OpenMPI/4.1.4/R/4.5.0/lib/R/library \
  --labels R_library_4.5,R-4.5.0-module > ${OUT}/rctd_seurat_packages.tsv

# C. dbscan_regions (proposal): roots = packages loaded by liverstomics_hotspot_dbscan_NI.R and
#    colocalisation_script_NI.R (dbscan added separately) plus BiocVersion and BiocManager, resolved
#    through the same libraries as B; dbscan 1.2.7 = the CRAN version of 5 to 27 October 2026.
python3 -I make_lock.py --roots dbscan_regions_roots.txt \
  --libs /data/projects/punim1236/liver_stomicsv1.3/RCTD_CLAHE_old/R_library_4.5,/apps/easybuild-2022/easybuild/software/MPI/GCC/11.3.0/OpenMPI/4.1.4/R/4.5.0/lib/R/library \
  --labels rctd_seurat:R_library_4.5,rctd_seurat:R-4.5.0-module \
  --extra "dbscan=1.2.7=CRAN=added:CRAN-current-2026-10-27" > ${OUT}/dbscan_regions_packages.tsv

# cellchat_r.yml: name=version=build of every package in /data/gpfs/projects/punim2394/conda/rctd_env/conda-meta/*.json
# (89 packages), followed by seven unpinned build tools and headers (see the comments in the file).
