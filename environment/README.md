# Environments

Five environments run the analysis scripts of this repository. Each was built from scratch on Spartan
(University of Melbourne; Red Hat Enterprise Linux 9, x86_64) on 20 November 2026 and tested on the
original data.

| Environment | Contents | Runs | Test on the original data |
|---|---|---|---|
| `sven_sammy_stereopy` | conda: Python 3.8.20, stereopy 1.6.2 | part 3, step 2: `01_Y01453J8_bin20_gef_to_h5ad.py` | h5ad file identical byte for byte |
| `sven_sammy_hotspot` | conda: Python 3.9.13, hotspotsc 1.1.3, scanpy 1.10.3, Jupyter | part 3, step 4: Hotspot notebook; part 2: `make_MNP2_vs_original_combinedMNP_comparison.sh` | notebook runs; its tables differ from the archived ones as much as new runs differ from each other (see Notes); the part 2 table is identical byte for byte |
| `cellchat` | R 4.5.1 (conda) and 210 R packages, CellChat 2.2.0.9001 | part 2: both CellChat scripts | tables `01` to `04` of both runs identical byte for byte |
| `rctd_seurat` | R 4.5.0 and 223 R packages, Seurat 5.5.1, spacexr 1.2.0 | part 1: Seurat scripts; part 3, steps 2 and 3: `02_Y01453J8_bin20_h5ad_to_rds.R`, `rctd_script_NI.R` | bin20 Seurat object identical; RCTD weights identical; resolution-0.7 clusters of part 1 reproduced for all 2,934 cells |
| `dbscan_regions` | R 4.5.0 and 76 R packages, dbscan 1.2.7 | part 3, step 5 | enrichment table identical byte for byte; same 12 regions (4,907 bins); 19 of the 24 colocalisation tables that do not depend on random sampling identical, the other 5 differ by at most 1.6e-4 |

Cell Ranger 8.0.1 (part 1; reference `refdata-gex-mm10-2020-A`) and SAW 8.2.2 (part 3, step 1) are installed
from 10x Genomics and STOmics. Not covered by these environments: `make_example_LR_expression_support.sh` of
part 2 (it needs only Seurat; not tested), the RStudio editor buffer in `2_CellChat/earlier_analysis_2025-11/`
and `03_Y01453J8_bin20_run_RCTD.R` (an earlier RCTD run whose output is not used).

Each build block below starts in the root of the repository.

## Python: conda environments

```bash
export PYTHONNOUSERSITE=1      # before building and using the environments; otherwise pip and Python use packages in ~/.local
unset PYTHONPATH
cd environment/python
conda env create -f sven_sammy_stereopy.yml
conda env create -f sven_sammy_hotspot.yml
```

To rebuild the exact tested environments without solving, use the lock files (with the same two settings):

```bash
cd environment/python
conda create -n sven_sammy_hotspot --file sven_sammy_hotspot.lock.txt
conda run -n sven_sammy_hotspot python -m pip install --no-deps -r sven_sammy_hotspot.pip-freeze.txt
# the same for sven_sammy_stereopy
```

The Hotspot notebook can be run without opening it, from the folder of a copy of the notebook whose paths
are edited (see part 3):

```bash
conda activate sven_sammy_hotspot
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 \
    --output hotspot_liverstomics_NI_executed.ipynb hotspot_liverstomics_NI.ipynb
```

## R: package libraries

Each install script downloads every package at the version in its lock table (`*_packages.tsv`) from CRAN,
Bioconductor 3.22 and GitHub, compiles all of them from source into a new library folder and checks the
installed versions at the end (written to `<env>_lock_check.tsv` in the current folder). It needs network
access and takes 6 to 9 minutes with 8 to 16 CPUs. Downloads and build logs go to R's temporary folder, which
R deletes when the script ends; set `LOCK_WORK_DIR=<folder>` to keep them.

```bash
cd environment/R

# cellchat: R 4.5.1 from conda
conda env create -p /path/to/r451 -f cellchat_r.yml    # or, without solving: conda create -p /path/to/r451 --file cellchat_r.explicit.txt
conda activate /path/to/r451
Rscript install_cellchat_env.R /path/to/Rlib_cellchat 16
conda deactivate

# rctd_seurat and dbscan_regions: R 4.5.0
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0    # Spartan; elsewhere R 4.5.0 with the system libraries in R/RESULTS.md
Rscript install_rctd_seurat_env.R /path/to/Rlib_rctd_seurat 16
Rscript install_dbscan_regions_env.R /path/to/Rlib_dbscan_regions 8
```

To use a library:

```bash
# CellChat scripts (lines 6 and 7 read these two variables)
conda activate /path/to/r451
export CELLCHAT_LIB=/path/to/Rlib_cellchat SEURAT_FALLBACK_LIB=/path/to/Rlib_cellchat

# all other R scripts
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0    # Spartan
export R_LIBS_USER=/path/to/Rlib_rctd_seurat    # /path/to/Rlib_dbscan_regions for part 3, step 5
```

## Notes

- **Hotspot gives slightly different tables in every run.** hotspotsc 1.1.3 builds the neighbour graph with
  `pynndescent` without a random seed, and on the bin20 grid many neighbours are equally far away, so each run
  keeps a different set of tied neighbours. The archived tables differ from new runs as much as new runs
  differ from each other; module 5, used in step 5, is the same module in every run (1,854 to 1,891 bins with
  `module_5 > 4` in four new runs, 1,882 in the archived run). With a fixed seed, the notebook and a script
  version of it write identical tables.
- **`dbscan_regions`.** The original session of part 3, step 5 ran R 4.3.1 (recorded in the PDF and `.rds`
  files it wrote). This library uses R 4.5.0 with the versions of the RCTD library plus `dbscan` 1.2.7, the
  CRAN version at the time of the analysis.
- **Part 1.** The original part 1 objects were saved by R 4.4.0. `rctd_seurat` was tested by re-running the
  clustering of the liver cells; `0_make_GEX_HTO_Processed.R` also runs with its package versions.
- The lock tables fix package versions and sources, not compiled binaries.

## Files

| Path | Content |
|---|---|
| [`python/RESULTS.md`](python/RESULTS.md), [`R/RESULTS.md`](R/RESULTS.md) | How each environment was built and tested: job ids, run times, memory, all comparisons, what was not tested |
| `python/*.yml`, `python/*.lock.txt`, `python/*.pip-freeze.txt` | conda environment files, explicit conda lock files, pip lists |
| `R/install_*_env.R`, `R/*_packages.tsv` | install scripts and lock tables |
| `R/cellchat_r.yml`, `R/cellchat_r.explicit.txt` | conda environment with R 4.5.1 for `cellchat`, and its explicit lock file |
| `R/lock_provenance/` | how the versions were read from the original R libraries |
| `python/tests/`, `R/tests/` | SLURM jobs that built each environment from scratch and compared the outputs with the archived results (Spartan paths; submit from the `tests/` folder after editing the paths at the top), their logs and results |
