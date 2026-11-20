# R environments for the archived R scripts: lock tables, install scripts and test results

This folder lets a reader build the R package libraries needed by the R scripts of the archive, with the
package versions of the runs that produced the archived results, and records how these libraries were
built and tested on the real data on Spartan (University of Melbourne) on 20 November 2026.

## Summary

| Environment | Used by | R | Packages (lock table) | Build from scratch | Tests on the real data |
|---|---|---|---|---|---|
| `cellchat` | part 2: `2_CellChat/run_cellchat_MNP2_*_latest.R` | 4.5.1 (conda-forge `r-base`, `cellchat_r.yml`) | 210: 201 CRAN, 8 Bioconductor 3.22, CellChat from GitHub at commit `75253cd` | 8.7 min R library + 33 s conda (16 CPUs) | both runs: tables 01 to 04 identical byte for byte to the archived tables; same 154 package versions in sessionInfo |
| `rctd_seurat` | part 3: `3_RCTD/rctd_script_NI.R`, `2_bin20_matrix/02_Y01453J8_bin20_h5ad_to_rds.R`; part 1: `2_seurat_clustering_annotation/*.R` | 4.5.0 | 223: 205 CRAN, 18 Bioconductor 3.22 | 7.6 min (16 CPUs) | bin20 conversion: identical object. Part 1: the 15 resolution-0.7 clusters reproduced for all 2,934 cells. RCTD: `weights` `identical()` to the archived result (maximum difference 0). |
| `dbscan_regions` (proposal) | part 3: `5_DBSCAN_regions_enrichment/*.R` | 4.5.0 (proposed) | 76: 60 CRAN, 16 Bioconductor 3.22 | 5.7 min (8 CPUs) | enrichment table identical byte for byte; region membership identical (12 regions, 4,907 bins); colocalisation runs (19 of 24 tables that do not depend on the random sampling identical, 5 differ by at most 1.6e-4) |

"Build from scratch" is the time of the install script (download of all source tarballs, compilation of
every package from source, version check) in a SLURM job on the `sapphire` partition, with the library on
node-local `/tmp`. Every installed package matched its locked version (`tests/results/*/*_lock_check.tsv`).
Whole jobs (build and all tests): `cellchat` job 35423699, 28 min 55 s; `rctd_seurat` job 35425723,
28 min 59 s; `dbscan_regions` job 35423372, 7 min 49 s.

## Files

The tests ran in a folder outside the repository, next to a local clone of it; in this repository the
sbatch files take the repository root (`../..` from the folder of the install scripts) as `ARCHIVE`. The
local folder is shown as `<build_dir>` in the logs.

| File | What it is |
|---|---|
| `cellchat_packages.tsv`, `rctd_seurat_packages.tsv`, `dbscan_regions_packages.tsv` | Lock tables: `package`, `version`, `source` (`CRAN`, `Bioconductor 3.22`, or `GitHub user/repo@commit`), `origin` (the original library the version was read from) |
| `install_cellchat_env.R`, `install_rctd_seurat_env.R`, `install_dbscan_regions_env.R` | Install scripts: `Rscript install_<env>_env.R <library_dir> [n_cpus]` |
| `cellchat_r.yml` | Conda environment with R 4.5.1 for `cellchat` (the packages of the original conda environment, plus build tools) |
| `lock_provenance/` | How the lock tables and the yml were made: `make_lock.py`, `make_conda_yml.py`, `make_lock_tables.sh` and the lists of root packages |
| `tests/build_and_test_<env>.sbatch` | SLURM jobs that build one environment from scratch on `/tmp` and run its tests (paths as variables at the top; submit from `tests/`) |
| `tests/test_cellchat.R`, `tests/test_rctd_seurat.R`, `tests/test_dbscan_regions.R` | Test drivers: they write copies of the archived scripts in which only paths are changed, run them or the sections named below, and compare the outputs with the archived results |
| `tests/test_helpers.R`, `tests/sessioninfo_at_exit.R` | Helper functions (path substitution, running line ranges, CSV comparison); an R profile (`R_PROFILE_USER`) that writes `sessionInfo()` at the end of an `Rscript` run of a script that does not write it |
| `tests/logs/` | SLURM logs of the three final jobs (35423699, 35425723, 35423372) |
| `tests/results/<env>/` | Per environment: lock check, test logs, comparison reports, `sessionInfo` of every test run, `diff` of each script copy against the archived script, `sacct` record. Not included in this repository: the install logs, the stderr logs of the CellChat runs (progress bars) and the regenerated tables that are identical to the archived ones. |

## How to build each environment

The install scripts need network access to CRAN (`https://cloud.r-project.org`, or `CRAN_MIRROR`),
Bioconductor (`https://bioconductor.org`, or `BIOC_MIRROR`) and, for CellChat, `codeload.github.com`. They
download the source tarball of every package at exactly the locked version (from the CRAN or Bioconductor
release folder, or from its `Archive/` folder once the version has been superseded), place them in a local
package repository and install them all from it with `install.packages()`, so no other version can be
chosen. Everything is compiled from source. At the end every installed version is checked against the lock
table; the script writes `<env>_lock_check.tsv` into the working directory and exits with status 1 on any
mismatch. Only `<library_dir>` and the base library of R are used (`R_LIBS`, `R_LIBS_USER` and
`R_LIBS_SITE` are set to `<library_dir>` while the script runs), so a personal library is neither read nor
changed. Installing `BiocVersion` 3.22.0 from the lock table has the same effect as
`BiocManager::install(version = "3.22")`: `BiocManager::version()` reports 3.22 in the new library.

### cellchat (R 4.5.1)

```bash
conda env create -p /path/to/r451 -f cellchat_r.yml
conda activate /path/to/r451
Rscript install_cellchat_env.R /path/to/Rlib_cellchat 16
# run the CellChat scripts with
export CELLCHAT_LIB=/path/to/Rlib_cellchat SEURAT_FALLBACK_LIB=/path/to/Rlib_cellchat
Rscript run_cellchat_MNP2_infected_vs_naive_CD4_latest.R     # after editing input_rds and outdir
```

### rctd_seurat (R 4.5.0)

```bash
# Spartan: module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0
# elsewhere: R 4.5.0 with a C/C++/Fortran compiler and the development files of zlib, libcurl, openssl,
# libxml2, freetype, harfbuzz, fribidi, fontconfig, libpng, libtiff, libjpeg, HDF5, ImageMagick, glpk, gmp
Rscript install_rctd_seurat_env.R /path/to/Rlib_rctd_seurat 16
export R_LIBS_USER=/path/to/Rlib_rctd_seurat
```

### dbscan_regions (R 4.5.0, proposal)

```bash
module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0      # Spartan; elsewhere R 4.5.0 with zlib, libcurl, openssl, ImageMagick
Rscript install_dbscan_regions_env.R /path/to/Rlib_dbscan_regions 8
export R_LIBS_USER=/path/to/Rlib_dbscan_regions
```

The `dbscan_regions` packages are a subset of `rctd_seurat` (same versions) plus `dbscan` 1.2.7, so the
`rctd_seurat` library with `dbscan` 1.2.7 added should also run part 3, step 5 (not tested).

## How the versions were determined

All original libraries were read, never changed. `lock_provenance/make_lock_tables.sh` holds the exact
commands. `make_lock.py` starts from a list of root packages and follows `Depends`, `Imports` and
`LinkingTo` through the original libraries in their `.libPaths()` order; `version` is the `Version` field of
the copy that R would load, `source` comes from `Repository`, `biocViews`, `RemoteType` and `RemoteSha`.

- **cellchat.** Roots: all 154 packages (attached and loaded) of `2_CellChat/MNP2_CD4_20261023/99_sessionInfo.txt`
  (the second run lists the same packages and versions), plus `BiocVersion`. Libraries:
  `/data/gpfs/projects/punim2394/AAA_Sorry_Its_Chen_Here/Rlib_new` (188 of the packages; built with R 4.5.0
  between 20 January and 20 March 2026, before the runs of 23 October 2026), then
  `/data/scratch/projects/punim1236/python_envs/rseurat/lib/R/library` (22 packages, conda-forge builds).
  Every sessionInfo version was checked against the resolved library (`--si`). CellChat:
  `RemoteSha 75253cd0c9e68410e6e721a6d3a0419a1d7e358f`, `jinworks/CellChat`. R: the conda environment
  `/data/gpfs/projects/punim2394/conda/rctd_env` (`r-base 4.5.1`, OpenBLAS 0.3.30), whose 89 conda packages
  are pinned in `cellchat_r.yml` by name, version and build.
- **rctd_seurat.** Roots: the packages loaded by `rctd_script_NI.R` (cowplot, tidyverse, spacexr, Seurat,
  SummarizedExperiment, SpatialExperiment), by `02_Y01453J8_bin20_h5ad_to_rds.R` and checked by its wrapper
  (Seurat, Matrix, hdf5r, S4Vectors) and by the part 1 scripts (Seurat, ggplot2, tidyverse, dplyr, ggrepel),
  plus BiocManager, BiocVersion and otel. Libraries in the order of the RCTD job (`R_LIBS_USER`, then the
  module library): `/data/projects/punim1236/liver_stomicsv1.3/RCTD_CLAHE_old/R_library_4.5` (177 packages,
  all installed on 14 July 2026; the closure covers every package of that library) and the library of the
  Spartan module `R/4.5.0` (46 packages that `R_library_4.5` does not hold, e.g. stringr 1.5.0,
  tibble 3.2.1, purrr 1.0.2; also stringi 1.7.12, because the `stringi` folder of `R_library_4.5` is an
  empty leftover of a failed installation, next to `00LOCK-stringi`, so R loads the module's copy). The RCTD log `3_RCTD/slurm-32372823.out` prints ten of the loaded versions
  (tidyverse 2.0.0, dplyr 1.2.1, readr 2.2.0, forcats 1.0.1, stringr 1.5.0, ggplot2 4.0.3, tibble 3.2.1,
  lubridate 1.9.5, tidyr 1.3.2, purrr 1.0.2); all ten agree with the lock table. The same library is
  offered for the part 1 Seurat scripts, whose package versions are not recorded (the objects were saved by
  R 4.4.0 and carry the SeuratObject stamp 5.1.0; the library has Seurat 5.5.1 and SeuratObject 5.4.0).
- **dbscan_regions (proposal).** The original session ran R 4.3.1 (the `Producer` field of all PDF files
  in its output folder and the version stamp of `hs_out.rds` and `sub_df_containsdbscan_info.rds`); its
  package library was not found, and the RCTD library has no `dbscan`. Roots: the packages loaded by the two scripts (data.table, dplyr, tidyr,
  ggplot2, SummarizedExperiment, SpatialExperiment, viridisLite, patchwork, Matrix, MASS, scales, reshape2;
  dbscan added separately), resolved through the same two libraries as `rctd_seurat`, so every version
  equals the `rctd_seurat` one. `dbscan` 1.2.7 is the version on CRAN from 5 October 2026 until after
  27 October 2026, the date of the archived enrichment table (P3M snapshots of 2026-10-06 and 2026-10-27
  both serve 1.2.7; the snapshot of 2026-07-14, the date of the RCTD library, serves 1.2.5).

## cellchat: build and test

- **R.** conda-forge `r-base` 4.5.1 from `cellchat_r.yml`. In the tested environment 87 of the 89 conda
  packages of the original environment have the same version and build; `sysroot_linux-64` 2.12 and
  `kernel-headers_linux-64` 2.6.32 were replaced by 2.17 and 3.10 (see "First attempt" below), and seven
  packages were added for compiling: pkg-config 0.29.2, cmake 4.0.3, zlib 1.3.1, libxml2 2.15.1,
  glpk 5.0, gmp 6.3.0 and nlopt 2.11.0 (nlopt brings Python 3.13 and numpy as conda dependencies; the R
  scripts do not use them). The solved environment, with download URLs and md5 sums of all 105 conda
  packages, is `cellchat_r.explicit.txt`; `conda create -p <prefix> --file cellchat_r.explicit.txt`
  recreates it without solving. OpenBLAS 0.3.30, as in the original run.
- **First attempt** (job 35423193; its files are not included in this repository): with the original
  `sysroot_linux-64` 2.12 (C library headers of CentOS 6) the C++ library of the conda compiler gcc 15.1
  does not compile (`error: 'timespec_get' has not been declared in '::'`): Rcpp, igraph, FNN, RANN,
  rjson and scattermore failed to compile, stringi failed in `configure` ("C++ compiler cannot create
  executables"), 46 packages that depend on them were skipped, and CellChat could not be installed
  (156 of 210 packages installed). Fix: `sysroot_linux-64=2.17` and `kernel-headers_linux-64=3.10` in
  `cellchat_r.yml`. These headers are used only when packages are compiled. (The original CellChat library was not compiled with
  conda compilers: `Rlib_new` was built with the Spartan module R 4.5.0 and GCC 11.3.0, and the 22
  packages from `rseurat` are conda-forge builds for R 4.3.3. The versions are the same; the binaries are
  not.)
- **Build** (job 35423699, account punim2394, 16 CPUs): conda environment 33 s; R library 8.7 min
  (downloads 2.8 min, 209 CRAN and Bioconductor packages 5.4 min, CellChat from GitHub 0.4 min);
  210 of 210 packages at the locked versions; CellChat has `RemoteSha 75253cd0c9e68410e6e721a6d3a0419a1d7e358f`.
- **Tests** (same job). Copies of the two archived scripts in which only `outdir` (line 20) is changed;
  `input_rds` (line 19) is set by the driver and here points to the same input files
  (`tests/results/cellchat/<run>/script_copy.diff`). Run as in the archived wrappers: 8 CPUs and
  `SLURM_CPUS_PER_TASK=8` (8 future workers), locale `en_AU.UTF-8`, `CELLCHAT_LIB` and
  `SEURAT_FALLBACK_LIB` set to the new library.

  | | Run 1: infected vs naive CD4 | Run 2: cytotoxic vs naive CD4 |
  |---|---|---|
  | SLURM step | 35423699.0 | 35423699.1 |
  | Run time (archived run) | 10 min 55 s (12 min 32 s) | 8 min 34 s (10 min 58 s) |
  | Peak memory, sacct MaxRSS (archived) | 13.5 GiB (13.5 GiB) | 13.2 GiB (13.1 GiB) |
  | Exit status (archived) | 0 (139, crash at exit) | 0 (139, crash at exit) |
  | Ligand-receptor pairs kept by `identifyOverExpressedInteractions` | 317 and 347 | 176 and 347 |
  | `01_input_annotation_counts.csv` | identical byte for byte | identical byte for byte |
  | `02_input_annotation_by_HTO.csv` | identical byte for byte | identical byte for byte |
  | `03_differential_LR_CD4_to_MNP2_*.csv` | identical byte for byte (51 pairs: 37 / 14) | identical byte for byte (41 pairs: 23 / 18) |
  | `04_differential_LR_MNP2_to_CD4_*.csv` | identical byte for byte (50 pairs: 39 / 11) | identical byte for byte (42 pairs: 25 / 17) |

  `00_run_config.txt` differs from the archived file only in `output_dir` and `CellChat_path`.
  `99_sessionInfo.txt` of both runs lists the same R version, locale and the same 154 packages at the same
  versions as the archived `99_sessionInfo.txt`; only the BLAS path differs. The rebuilt environment does
  not crash when R exits (the archived jobs ended with a segmentation fault after `DONE`; the cause was not
  investigated). Figures `05` to `13` and the CellChat `.rds` objects were not compared.
- **sessionInfo:** `tests/results/cellchat/infected/99_sessionInfo.txt`, `tests/results/cellchat/cytotoxic/99_sessionInfo.txt`
  (written by the scripts themselves).

## rctd_seurat: build and test

- **R.** The Spartan module R 4.5.0 (`module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0`; EasyBuild, GCC 11.3.0,
  FlexiBLAS with OpenBLAS), the R of the archived RCTD and bin20 runs.
- **Build** (final run, job 35425723, account uom00060, 16 CPUs): 223 of 223 packages at the locked versions
  in 7.6 min (downloads 1.2 min, compilation 6.2 min). No failures. The two earlier jobs below built the same library in 8.5 and 9.0 min.
- **Test D, part 1** (`tests/test_rctd_seurat.R seurat`, step 35425723.0, 11 s, peak memory 2.2 GiB).
  `library()` of Seurat 5.5.1, ggplot2 4.0.3, tidyverse 2.0.0, dplyr 1.2.1 and ggrepel 0.9.8, the packages of
  the three part 1 scripts: all load. Then lines 1-3, 6, 20, 65, 98 and 101 of a copy of
  `liver_scRNAseq_processed_NI.R` (paths changed; `tests/results/rctd_seurat/D_script_copy.diff`): read
  `GEX_HTO_Processed.rds`, select the liver cells of line 65, `FindClusters(seurat_obj, resolution = 0.7)`.
  Result: 2,934 cells in the same order as `sven_liver_for_paper.RDS`; 15 clusters (0 to 14); **all 2,934
  cells get the stored cluster number** (`seurat_clusters` and `RNA_snn_res.0.7` identical, adjusted Rand
  index 1); cluster sizes 0:754 1:491 2:364 3:312 4:215 5:170 6:167 7:163 8:72 9:67 10:65 11:51 12:24 13:13
  14:6. (Line 65 passes `subset()` the cells of all hashtags in the selected clusters; `subset()` keeps
  those present in `liver`.)
- **Test B1, bin20 conversion** (step 35425723.1, 20 s, peak memory 5.1 GiB; comparison step
  35425723.2). Copy of `02_Y01453J8_bin20_h5ad_to_rds.R` with the input and output folders changed
  (`B1_script_copy.diff`), run with `Rscript` on the archived `Y01453J8_tissue_bin20_raw_counts.h5ad`. The
  new Seurat object is `identical()` to the archived `Y01453J8_tissue_bin20_raw_counts.rds` that RCTD read
  (24,886 genes x 293,543 bins; counts, metadata and spatial embedding identical; SeuratObject stamp 5.4.0
  in both).
- **Test B2, RCTD** (step 35425723.3, 20 min 25 s, peak memory 15.5 GiB, 16 CPUs; the archived job:
  19 min 18 s, 15,757 MB of memory used according to `my-job-stats`, 16 CPUs). Copy of `rctd_script_NI.R` with the working directory and
  the output file changed (`B2_script_copy.diff`); the three inputs are the archived files of lines 13, 25
  and 57 (Tabula Muris droplet liver object, `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS`, the bin20
  rds). Result compared with `rctd_MNP1_2_Tcell_combined_no_sinusoids.rds` (step 35425723.4): same class
  (SpatialExperiment), dimensions (6 cell types x 289,703 bins), cell-type and bin names in the same
  order, `spatialCoords` identical; **the `weights` assay is `identical()` to the archived one (maximum
  absolute difference 0)**. The run log is the same as the archived log `3_RCTD/slurm-32372823.out` from
  the tidyverse start-up message (with its ten package versions) to the last RCTD message (reference of
  4,301 cells and 35,856 genes, DE gene counts per cell type, the three likelihood and sigma values); only
  the two memory tables printed by `gc()` differ.
- **sessionInfo:** `tests/results/rctd_seurat/sessionInfo_seurat_scripts.txt` (test D, written by the
  driver), `B1_sessionInfo_bin20_run.txt` and `B2_sessionInfo_rctd_run.txt` (written at the end of the
  `Rscript` runs by the R profile `tests/sessioninfo_at_exit.R`), and `*_compare_sessionInfo.txt` of the
  comparison steps.
- **Earlier runs** (their files are not included in this repository). Job 35423373 ran the same
  build and tests with the same results (D: 2,934 of 2,934 cells; B1 identical; B2 weights identical, RCTD
  19 min 33 s, 15.9 GiB), but did not record `sessionInfo()` of the B1 and B2 runs. Job 35424778
  ran the copies through a wrapper that
  `source()`d them and then wrote `sessionInfo()`. That changes the behaviour of the top-level
  `on.exit(h5$close_all())` in `02_Y01453J8_bin20_h5ad_to_rds.R`: inside `source()` it runs immediately,
  closes the HDF5 file, and the next read fails (`Error ...: id is invalid`). The job was cancelled and the
  wrapper replaced by the R profile, which leaves the `Rscript` run unchanged.

## dbscan_regions: build and test

- **Build.** Job 35423372 (account uom00060, 8 CPUs): 76 of 76 packages at the locked versions in 5.7 min
  (downloads 0.6 min, compilation 5.0 min). No failures.
- **Test C** (`tests/test_dbscan_regions.R`, job step 35423372.0, 2 min 1 s, peak memory 3.1 GiB).
  Copies of `liverstomics_hotspot_dbscan_NI.R` and `colocalisation_script_NI.R` in which only paths are
  changed: the two inputs, the output folder (to `/tmp`), the `setwd()` folder and, because of that, the
  RCTD file name on line 77, which was relative to it (`tests/results/dbscan_regions/C_script_copy.diff`).
  As the archive documentation says, `colocalisation_script_NI.R` is sourced first; then lines 1-40, 64-253 and 307-615
  of the copy are run in order (reading the inputs, the long weight table, bin selection, DBSCAN, region
  outlines and membership, the check figures, saving `hs_out` and `sub_df`, colocalisation per region,
  pooled, outside the regions and in 12 random windows, enrichment, figures). Not run: lines 42-61 (module
  maps and a save/reload of the RCTD object that later lines do not use) and 258-304 (cell-type maps; `ct`
  is set by hand).
  - Bins with `module_5 > 4.0`: 1,882; DBSCAN clusters with at least 20 bins: 12 (the numbers given in the archive documentation).
  - `enrichment/enrichment_of_celltypes.csv`: **identical byte for byte** to
    `5_DBSCAN_regions_enrichment/result_tables/enrichment_of_celltypes.csv`.
  - `sub_df_containsdbscan_info.rds`: `identical()` to the archived object (29,442 rows). Regions: 12, bins
    per region 1:1029 2:252 3:629 4:229 5:1063 6:230 7:418 8:268 9:164 10:153 11:215 12:257, 4,907 bins in
    all, the same as in the archived object.
  - `hs_out.rds`: equal to the archived object by `all.equal()`, not `identical()`: 23 of the 1,761,258
    module scores differ in the last bit (at most 4.4e-16); `data.table::fread()` 1.18.4 and the reader
    of the original session turned these 23 decimal strings of `hotspot_spot_scores.csv` into adjacent
    doubles. The original session therefore probably used another data.table version. The bin selection
    (`module_5 > 4`) is identical.
  - Colocalisation. Regions and random windows hold fewer than 2,000 bins, so the "random sample" of
    `nbins = 2000` bins is a permutation of all bins and the Spearman matrices do not depend on the seed.
    Of these 24 tables, 19 are identical byte for byte to the archived `rho_region_40um<k>.csv` and
    `background_regions/rho_rand_<k>.csv`; 5 (regions 1, 3, 5, 9 and `rand_3`) differ by at most 1.6e-4
    (region 9; the others at most 5.2e-6). The cause was not established; ties between neighbourhood means
    that are equal only up to rounding, which change Spearman ranks, are one possible explanation. The 12
    random windows are evidently the same as in the original run (`set.seed(42)`; 11 of their 12 tables
    are identical byte for byte). The tables that use a true random sample (`rho_background40um.csv`,
    pooled regions, pooled windows) were only checked to run.
  - The per-region enrichment (lines 588-599) is not written by the script; the driver saves it as
    `tables/enrichment_per_region_TEST_DRIVER.csv` for the record. No archived table to compare.
- **sessionInfo:** `tests/results/dbscan_regions/sessionInfo_dbscan_regions.txt`.

## What was not tested, and limits

- **Other machines.** All builds ran on Spartan (Red Hat Enterprise Linux 9.8, x86_64) with network
  access. `rctd_seurat` and `dbscan_regions` used the Spartan module R 4.5.0 (EasyBuild, GCC 11.3.0), which
  also provides the system libraries (HDF5 1.12.2, ImageMagick 7.1.0, GLPK, ...); an R 4.5.0 from another
  source (built from source, rig, or conda `r-base=4.5.0`) was not tried. `cellchat` used conda and should
  work on any x86_64 Linux machine with conda, but was tested only on Spartan.
- **Binaries are not the originals.** The lock tables fix package versions and sources, not compiled
  binaries. `cellchat` packages are compiled with the conda compiler gcc 15.1, while the original library
  was compiled with GCC 11.3.0 (and 22 packages came as conda-forge builds for R 4.3.3). `rctd_seurat`
  and `dbscan_regions` are compiled under the same module, R and compiler as the original library. All
  tested outputs are nevertheless identical (except the small colocalisation differences of part 3,
  step 5).
- **dbscan_regions is a proposal.** The original step 5 ran R 4.3.1 with a package library that was not
  found. The test shows
  that this library reproduces the enrichment table and region membership exactly; the one-bit
  differences in `hs_out.rds` show that the original session read the Hotspot table with another
  data.table version, which cannot be identified from the archive.
- **Random parts.** The colocalisation tables computed from a random sample of 2,000 bins (background,
  pooled regions, pooled random windows) were only run, not compared, because the scripts set no seed
  before them (no seed is set before line 443 of `liverstomics_hotspot_dbscan_NI.R`). The outputs with `100` in their name come from a state of the functions
  that is not archived and were not reproduced (the script copy writes them from the radius-40 result, as
  the archived script does).
- **Figures** (PDF and PNG) were written by the test runs but not compared with the archived figures.
  CellChat `.rds` objects are not archived and were not compared.
- **Part 1 (D).** Only lines 1-3, 6, 20, 65, 98 and 101 of `liver_scRNAseq_processed_NI.R` were run (read
  `GEX_HTO_Processed.rds`, select the cells, `FindClusters(resolution = 0.7)`). The plotting, annotation,
  signature-score and dot-plot sections and the scripts `sven_liverscRNAseq.R` and `Tcell_combined_ref.R`
  were only covered by loading their packages (Seurat 5.5.1, ggplot2 4.0.3, tidyverse 2.0.0, dplyr 1.2.1,
  ggrepel 0.9.8). The original part 1 objects were saved by R 4.4.0; the package versions of those
  sessions are not recorded, and the step that built `GEX_HTO_Processed.rds` is not archived.
- **Not covered by these three environments:** the Python steps (`2_bin20_matrix/01_*.py`, the Hotspot
  notebook), SAW and Cell Ranger; `2_bin20_matrix/03_Y01453J8_bin20_run_RCTD.R` (an earlier RCTD run whose
  result is not used); `2_CellChat/commands_for_two_extra_tables/make_example_LR_expression_support.sh`
  (it used the `Rscript` of the conda environment `rseurat`; it needs only Seurat and would probably run
  with the `cellchat` library, not tested); the RStudio editor buffer in `2_CellChat/earlier_analysis_2025-11/`.
- **Harmless messages.** Each `srun` step log begins with `error: couldn't chdir to '/tmp/...': No such
  file or directory: going to /tmp instead`: SLURM tries the working directory before the job's private
  `/tmp` is mounted. The test commands use absolute paths (the `rctd_seurat` and `dbscan_regions` steps
  also `cd` into their folder), so this has no effect on the outputs. The CellChat logs show `TRUE` after many
  lines; this comes from the archived scripts (`log_line()` passes `flush = TRUE` to `cat()`, which prints it).
