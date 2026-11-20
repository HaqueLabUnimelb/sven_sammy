# 3. Liver Stereo-seq v1.3: RCTD, Hotspot, DBSCAN regions and cell-type enrichment

One Stereo-seq (STOmics) T FF v1.3 chip of mouse liver, serial number `Y01453J8`. The analysis goes
through five steps:

1. SAW maps the reads and builds the expression matrix.
2. The matrix is binned to bin20 (20 x 20 chip positions per bin; 293,543 bins, 24,886 genes).
3. RCTD estimates, for every bin, the weight of each cell type of a single-cell reference.
4. Hotspot finds modules of genes with a shared spatial pattern and gives every bin a score per module.
5. Bins with a high score for module 5 are clustered with DBSCAN. Each cluster is enlarged to a region
   with a kernel-density outline (the scripts call these regions granulomas). The RCTD weights inside the
   regions are compared with the rest of the chip (cell-type enrichment), and cell-type colocalisation is
   computed inside the regions, outside them and in random windows.

All files are byte-for-byte copies. [`../MANIFEST.tsv`](../MANIFEST.tsv) lists the original path, owner,
last-modified time and md5 of each file. The scripts use absolute paths on the Spartan cluster
(University of Melbourne), mainly under `/data/projects/punim1236/liver_stomicsv1.3/`.

Lengths in the scripts (`radius`, `eps`, `h`, `R`) are in the coordinate units of the bin20 matrix, in
which neighbouring bins are 20 units apart. The cell-type names (`MNP1`, `MNP2`, `infected CD4`,
`naive CD4`, `uninfected MNP`) are the labels of the scRNA-seq reference (see
[`../1_scRNAseq_preprocessing_clustering_annotation/README.md`](../1_scRNAseq_preprocessing_clustering_annotation/README.md#objects-and-labels)).
The four files with `_NI` in their name come from the folder
`/data/projects/punim1236/liver_stomicsv1.3/scripts_for_NI/`, in which their author saved them on
19 November 2026; the folder also holds `liver_scRNAseq_processed_NI.R` of part 1. Comments in the part 1
scripts speak of the "NI revision"; the archived files do not expand `NI`.

## Files

| Step | File | What it does |
|---|---|---|
| 1. SAW | `1_SAW/run_2026-07-06_original_image/config.yaml` | Settings and the full `saw count` command that SAW 8.2.2 recorded for the run of 6 July 2026 (SLURM job 26830116): kit version `Stereo-seq T FF V1.3`, sequencing type `PE75_50+100`, organism mouse, tissue liver, the original DAPI image, 32 threads, memory 200. Its output `Y01453J8.tissue.gef` is the input of step 2. |
| | `1_SAW/run_2026-07-06_original_image/saw_liver_26830116.log` | SLURM log of that run. |
| | `1_SAW/run_2026-07-06_original_image/SAWv1.3_pipeline_modified.slurm` | SLURM script that submitted that run: `saw count` with the original DAPI image, 32 CPUs, 200 GB and `--memory 200`, output in `/data/projects/punim1236/AAA_Chen_is_here/Ollie_SAW/output_liver_stomics` (the output now sits in `/data/projects/punim1236/output_liver_stomics/`). The file was later removed from Spartan; this copy is the version kept in the edit history of code-server (the browser version of the VS Code editor) in the home folder of `czhu3869`, saved at 11:25:29 on 6 July 2026, 20 seconds before the job was submitted. The settings it prints match the log. |
| | `1_SAW/SAWv1.3_pipeline.sh` | SLURM script that runs `saw count`. This is the script of a second SAW run of the chip (10 July 2026, job 27045467). It differs from the script of the first run in the image folder (`IMAGEDIR`, a contrast-enhanced copy of the DAPI image), the working and output folder (`WORKDIR`), the SLURM request (partition `long`, 16 CPUs, 100 GB, 31 days instead of `sapphire`, 32 CPUs, 200 GB, 48 hours) and `--memory 300`. The output of the second run is not used in steps 2 to 5. |
| 2. bin20 matrix | `2_bin20_matrix/run_Y01453J8_bin20_RCTD.slurm` | SLURM wrapper that runs the three scripts below, one after the other. Each step is skipped when its output file exists. |
| | `2_bin20_matrix/01_Y01453J8_bin20_gef_to_h5ad.py` | Reads `Y01453J8.tissue.gef` of the first SAW run with stereopy (`st.io.read_gef(bin_size=20)`) and writes the raw counts as `Y01453J8_tissue_bin20_raw_counts.h5ad`, with the bin coordinates in `obs["x"]` and `obs["y"]`. This file is the input of step 4. The stereopy version of the run is not recorded; the Python environment that the wrapper names held stereopy 1.6.2 and Python 3.8.20 when it was removed on 18 November 2026. |
| | `2_bin20_matrix/02_Y01453J8_bin20_h5ad_to_rds.R` | Converts the h5ad file to a Seurat object, `Y01453J8_tissue_bin20_raw_counts.rds` (no gene or bin filter). This file is the input of step 3. |
| | `2_bin20_matrix/03_Y01453J8_bin20_run_RCTD.R` | A first RCTD run in doublet mode with an earlier reference. It was replaced by step 3 and its result is not used in steps 4 and 5. It is archived because the wrapper calls it. |
| | `2_bin20_matrix/slurm-Y01453J8_bin20_RCTD-27327937.out` | Log of job 27327937 (16 July 2026), which ran `02` and `03` (R 4.5.0, spacexr 1.2.0). `01` had written the h5ad file a few minutes earlier in job 27327611, which then failed; no log of that job was found. |
| 3. RCTD | `3_RCTD/rctd_script_NI.R` | Builds the reference: hepatocytes of the Tabula Muris droplet liver object (1,764 cells) merged with the annotated liver scRNA-seq object `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS` (2,537 cells; labels naive CD4, infected CD4, MNP1, MNP2, uninfected MNP). Runs spacexr on the bin20 object: `createRctd()` with default arguments, then `runRctd(rctd_mode = "full", max_cores = 4)`. |
| | `3_RCTD/RCTD_CLAHE_bin20_liverSTOmics_newenv_project2394.sh` | SLURM wrapper of the RCTD run: R 4.5.0 (`module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0`), the R library of step 2, 16 CPUs, 100 GB. It calls the script under its working name `RCTD_CLAHE_bin20_liverSTOmics_newenv.R`, which is the same file as `rctd_script_NI.R`. |
| | `3_RCTD/slurm-32372823.out` | Log of the RCTD run of 6 October 2026 (job 32372823, 19 minutes). It shows the reference (4,301 cells, 35,856 genes, the six cell types). The log does not name the output file; `rctd_MNP1_2_Tcell_combined_no_sinusoids.rds`, the result read in step 5, was written when this job ended. |
| 4. Hotspot | `4_Hotspot/hotspot_liverstomics_NI.ipynb` | Jupyter notebook (Python 3.9.13, hotspotsc 1.1.3, scanpy 1.10.3). Reads the bin20 h5ad file, keeps genes detected in at least 1 % of the bins (3,467 genes), runs Hotspot and writes the gene tables, module tables, per-bin module scores and module maps. |
| | `4_Hotspot/result_tables/` | Five small tables from `/data/projects/punim1236/liver_stomicsv1.3/hotspot/`, written on 20 October 2026: `hotspot_autocorrelation.csv` (C, Z, p-value and FDR per gene), `hotspot_modules.csv` and `hotspot_gene_modules.csv` (module of each gene), `hotspot_module_genes_ranked.csv` (genes of each module ranked by correlation with the module score) and `module_5_genelist.txt` (the 152 genes of module 5). |
| 5. Regions, enrichment, colocalisation | `5_DBSCAN_regions_enrichment/liverstomics_hotspot_dbscan_NI.R` | Reads the RCTD result and the per-bin Hotspot scores. Selects the bins, runs DBSCAN, defines the regions, runs the colocalisation and the enrichment and draws the figures. |
| | `5_DBSCAN_regions_enrichment/colocalisation_script_NI.R` | Three functions for the colocalisation. `coloc_analysis_replicates_fast()` is the one called by the script above. |
| | `5_DBSCAN_regions_enrichment/result_tables/enrichment_of_celltypes.csv` | The enrichment table of 27 October 2026 (the file that line 562 of the script writes). |

## Parameters

| Step | Setting, as written in the code |
|---|---|
| Hotspot genes | `sc.pp.filter_genes(adata, min_cells=int(0.01 * adata.n_obs))` |
| Hotspot model | `hotspot.Hotspot(adata, layer_key="counts", model="danb", latent_obsm_key="spatial", umi_counts_obs_key="total_counts")` |
| Hotspot neighbours | `hs.create_knn_graph(weighted_graph=False, n_neighbors=15)` in the notebook as archived. The result tables were written by a run with `n_neighbors=30`; see note 1 below. |
| Hotspot modules | 500 genes with FDR < 0.05 and the highest Z; `hs.create_modules(min_gene_threshold=20, core_only=False, fdr_threshold=0.05)`: 6 modules of 97, 62, 26, 105, 152 and 33 genes; `hs.calculate_module_scores()` |
| Bins of interest | `module_5 > 4.0`: 1,882 of 293,543 bins |
| DBSCAN | `dbscan(xy, eps = 200, minPts = 15)`; clusters with fewer than 20 bins are dropped; 12 regions |
| Region outline | For each cluster: `MASS::kde2d(d$x, d$y, h = c(200, 200), n = 400)`, density divided by its maximum; every bin of the chip with a value of at least 0.05 belongs to the region (lines 213 to 236). Lines 161 to 163 and 198 to 200 draw this outline with `stat_density_2d(contour_var = "ndensity", breaks = 0.05, h = c(200, 200))`. The cell-type maps and the zoom figures (lines 266 to 304) draw a different contour, `breaks = 0.1, h = c(150, 150)`. |
| Enrichment, all regions together | Per cell type, RCTD weights of the bins inside any region against all other bins: `wilcox.test(a, b)`, `log2FC = log2((mean(a) + 1e-3) / (mean(b) + 1e-3))`, `auc = W / (length(a) * length(b))`, p-values adjusted with `p.adjust(method = "BH")` (lines 533 to 562) |
| Enrichment, per region | Each region against all bins outside the regions: log2FC and AUC as above (lines 588 to 599) |
| Colocalisation | `coloc_analysis_replicates_fast(radius = 40, nbins = 2000)`: mean weight of every cell type over the neighbours within `radius` of each bin (`dbscan::frNN`), up to 2,000 randomly sampled bins per input table, Spearman correlation between cell types. Run per region, for the 12 regions pooled, for all bins outside the regions, and for 12 random windows. |
| Random windows | `R <- 300`, `n_regions <- 12`, `min_spots <- 50`, `min_sep <- 6 * R`, `set.seed(42)`; centres are bins more than `R` away from every region (lines 432 to 468) |

## How to reproduce

1. SAW: edit the paths at the top of `run_2026-07-06_original_image/SAWv1.3_pipeline_modified.slurm`
   (`WORKDIR`, `SAWDIR`, `REFDIR`, `FASTQDIR`, `MASKDIR`, `IMAGEDIR`) and the `--output` and `--error` lines,
   then submit it with `sbatch` (SAW 8.2.2, FASTQ files, chip mask, original DAPI image, mouse reference
   index). `config.yaml` records the full `saw count` command of the run.
2. bin20 matrix: edit the paths at the top of `run_Y01453J8_bin20_RCTD.slurm` and of the scripts `01` and
   `02`, then `sbatch run_Y01453J8_bin20_RCTD.slurm`. `01` needs Python with stereopy. Before `02`
   starts, the wrapper checks that the R packages Matrix, hdf5r, Seurat, spacexr, SummarizedExperiment,
   SpatialExperiment and S4Vectors are installed (R 4.5). The third step of the wrapper (script `03`) is
   not needed: either edit the paths in `03` as well, or expect the job to stop with an error after `02`
   has written the rds file.
3. RCTD: in `rctd_script_NI.R` edit the working directory (line 1), the three input paths (lines 13, 25,
   57) and the output path (line 59, see note 2). In the wrapper edit the `cd` line, the script path and
   `R_LIBS_USER` (an R 4.5 library with the packages loaded at the top of the script), then submit it
   with `sbatch`.
4. Hotspot: open the notebook in Jupyter and edit the input path in the second code cell and the output
   folder in the `os.chdir` cell. Run the `os.chdir` cell before the first cell that writes files. The
   tables used in step 5 come from `n_neighbors=30` (note 1), while the archived cell has 15: set the
   value to 30 to follow the archived results. Step 5 selects module 5 by its number, so compare your
   `hotspot_modules.csv` with `4_Hotspot/result_tables/hotspot_modules.csv` before going on.
5. Regions, enrichment, colocalisation: in an R session run
   `source("colocalisation_script_NI.R", encoding = "UTF-8")`, then run
   `liverstomics_hotspot_dbscan_NI.R` section by section after editing the paths. The output folder needs
   the subfolders `zoom/`, `background_regions/` and `enrichment/`.

## Notes on the archived state

1. **Hotspot neighbours.** The Hotspot cell of the archived notebook shows `n_neighbors=15` and the output
   of that run (top gene `Gm42418`, Z = 813.08). The tables written to disk on 20 October 2026, including
   `hotspot_autocorrelation.csv` (top gene `Gm42418`, Z = 1080.62) and the per-bin scores read in step 5,
   match the Hotspot cell of an earlier saved notebook,
   `/data/projects/punim1236/b_cells/Untitled3.ipynb` (20 October 2026, 17:12), which has
   `n_neighbors=30`. Its 12 other code cells and their outputs are identical to the archived notebook; it
   lacks the last two code cells, which tidy the column names and write the per-bin scores. A checkpoint
   of 20 October 2026, 18:24
   (`/data/projects/punim1236/b_cells/.ipynb_checkpoints/hotspot_liverstomics-checkpoint.ipynb`) holds
   all cells with `n_neighbors=30`.
2. **RCTD output name.** Line 59 of `rctd_script_NI.R` names the output
   `rctd_MNP1_2_Tcell_splitup_no_sinusoids.rds`. The run of 6 October 2026 used the reference on line 25
   and wrote `rctd_MNP1_2_Tcell_combined_no_sinusoids.rds`, which is the file read in step 5. The name on
   line 59 is that of a later run (job 33977712, 29 October 2026) with the "split up" reference. Run
   unchanged on Spartan, the archived script would overwrite the result of that later run. "No sinusoids"
   means that only hepatocytes are taken from the Tabula Muris reference (line 19).
3. **Region outline.** The outline is a kernel-density contour. A comment in the script (line 148) calls it
   a convex hull; the function `hull_area()` (lines 432 to 435) is defined and not called. The contour
   drawn on the cell-type maps and zoom figures (`breaks = 0.1, h = c(150, 150)`) is not the outline that
   defines the regions.
4. **Interactive scripts.** `liverstomics_hotspot_dbscan_NI.R` was run section by section. It does not load
   `colocalisation_script_NI.R` itself. Lines 260 to 265 set the cell type `ct` six times; the plot that
   follows was saved once per cell type by changing `ct` by hand. Lines 411 to 415 save the same object
   under a `40um` and a `100um` name; on Spartan the two files differ, and the outputs with `100` in their
   name come from a state of the functions with `radius = 100` that is not archived. The colocalisation
   samples bins at random and no seed is set before line 443.
5. **Length labels.** File names such as `coloc_region40um_*.pdf` use the label `um` for `radius = 40`,
   which is in coordinate units. No script states the physical length of one unit.
6. **"CLAHE" in file names.** CLAHE is the contrast enhancement of the DAPI image used for the second SAW
   run. The RCTD wrapper and its working script carry the word in their names, but they read the bin20
   matrix of the first SAW run.

## Not in this archive

- Input data and result objects: FASTQ files, chip mask, DAPI image, SAW output, the bin20 h5ad and rds
  files, the Tabula Muris object, RCTD results and the per-bin Hotspot scores
  (see the main [README](../README.md#data-that-is-not-in-this-repository)).
- Figures. The Hotspot maps are in `/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_plots/`,
  the figures and tables of step 5 in
  `/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/`.
- R version and package versions of the R session of step 5: not recorded.
- The Python environment of step 2 (`stereopy_latest_py38`): removed from Spartan on 18 November 2026.
