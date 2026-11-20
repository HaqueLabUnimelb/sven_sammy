# Part 3. Stereo-seq v1.3: RCTD, Hotspot, DBSCAN regions and cell-type enrichment

One Stereo-seq T FF v1.3 chip of mouse liver (`Y01453J8`). The expression matrix is binned to bin20
(20 x 20 chip positions per bin; 293,543 bins). RCTD estimates the weight of each reference cell type per
bin. Hotspot finds gene modules with a shared spatial pattern. Bins with a high module 5 score are
clustered with DBSCAN and enlarged to regions (called granulomas in the scripts). The RCTD weights inside
the regions are compared with the rest of the chip (enrichment) and cell-type colocalisation is computed.

| Step | Folder | Main files | Output used by the next step |
|---|---|---|---|
| 1. SAW (Stereo-seq Analysis Workflow) | `1_SAW/` | `run_2026-07-06_original_image/SAWv1.3_pipeline_modified.slurm`, `config.yaml`, log | `Y01453J8.tissue.gef` |
| 2. bin20 matrix | `2_bin20_matrix/` | `01_Y01453J8_bin20_gef_to_h5ad.py`, `02_Y01453J8_bin20_h5ad_to_rds.R`, wrapper `run_Y01453J8_bin20_RCTD.slurm` | `Y01453J8_tissue_bin20_raw_counts.h5ad` and `.rds` |
| 3. RCTD | `3_RCTD/` | `rctd_script_NI.R`, wrapper, log | `rctd_MNP1_2_Tcell_combined_no_sinusoids.rds` |
| 4. Hotspot | `4_Hotspot/` | `hotspot_liverstomics_NI.ipynb`; `result_tables/` | `hotspot_spot_scores.csv` |
| 5. Regions and enrichment | `5_DBSCAN_regions_enrichment/` | `liverstomics_hotspot_dbscan_NI.R`, `colocalisation_script_NI.R`; `result_tables/` | enrichment table, figures |

`1_SAW/SAWv1.3_pipeline.sh` (script of a second SAW run, 10 July 2026) and
`2_bin20_matrix/03_Y01453J8_bin20_run_RCTD.R` (a first RCTD run, 16 July 2026) belong to runs whose output is
not used in the later steps; they are kept because they document the processing of the chip.

## Run

Environments: [`../environment/README.md`](../environment/README.md). The scripts contain absolute paths of
the Spartan cluster; edit them first.

| Step | Environment | What to run |
|---|---|---|
| 1 | SAW 8.2.2 | edit the paths at the top of `SAWv1.3_pipeline_modified.slurm`, then `sbatch` it |
| 2 | conda `sven_sammy_stereopy`, then R `rctd_seurat` | `python 01_Y01453J8_bin20_gef_to_h5ad.py`, then `Rscript 02_Y01453J8_bin20_h5ad_to_rds.R` (edit input and output paths at the top of each) |
| 3 | R `rctd_seurat` | edit lines 1, 13, 25, 57 and 59 of `rctd_script_NI.R`, then `Rscript rctd_script_NI.R` (16 CPUs, 100 GB, about 20 minutes) |
| 4 | conda `sven_sammy_hotspot` | in the notebook, edit the input path (second code cell) and the `os.chdir` folder (eighth code cell), set `n_neighbors=30` (third code cell), run all cells. Run top to bottom, the fifth code cell writes `hotspot_plots/local_correlations.pdf` and `.png` before the `os.chdir` cell, into the starting folder |
| 5 | R `dbscan_regions` (reconstructed; the original session ran R 4.3.1) | in R, `source("colocalisation_script_NI.R", encoding = "UTF-8")`, then run `liverstomics_hotspot_dbscan_NI.R` section by section (edit paths; create the output folder with its subfolders `zoom/`, `background_regions/`, `enrichment/`, and the folder of line 59) |

## Key settings

Lengths are in the units of the bin coordinates `x` and `y`, in which neighbouring bins are 20 apart
(`eps = 200` is 10 bins).

| Step | Setting | Where |
|---|---|---|
| RCTD reference | Tabula Muris liver hepatocytes (1,764 cells) + the 2,537 annotated liver cells of part 1 | `rctd_script_NI.R` lines 13 to 41 |
| RCTD | `createRctd()` defaults (among them `UMI_min = 100`: weights for 289,703 of the 293,543 bins); `runRctd(rctd_mode = "full")` | `rctd_script_NI.R` lines 129 to 139 |
| Hotspot genes | genes detected in at least 1 % of bins (3,467) | notebook, second code cell |
| Hotspot | `model = "danb"`, `n_neighbors = 30` (the archived cell shows 15, see below); top 500 genes with FDR < 0.05; `create_modules(min_gene_threshold = 20, core_only = False, fdr_threshold = 0.05)`: 6 modules | notebook, third and fourth code cells |
| Bins of interest | `module_5 > 4.0` (1,882 bins) | `liverstomics_hotspot_dbscan_NI.R` lines 114 to 115 |
| DBSCAN | `eps = 200`, `minPts = 15`, clusters with fewer than 20 bins dropped: 12 regions | `liverstomics_hotspot_dbscan_NI.R` lines 124 to 130 |
| Region | 2-D kernel density per cluster (`MASS::kde2d`, `h = 200`), bins with normalised density at least 0.05 | `liverstomics_hotspot_dbscan_NI.R` lines 213 to 236 |
| Enrichment | RCTD weights inside versus outside the regions: Wilcoxon test, `log2((mean_in + 0.001) / (mean_out + 0.001))`, AUC, BH-adjusted p | `liverstomics_hotspot_dbscan_NI.R` lines 533 to 562 |
| Colocalisation | neighbourhood means within `radius = 40`, 2,000 random bins, Spearman correlation; regions, rest of the chip, 12 random windows (`R = 300`, `set.seed(42)`) | functions in `colocalisation_script_NI.R`; calls in `liverstomics_hotspot_dbscan_NI.R` lines 335 to 523 |

## Notes

- **Hotspot neighbours.** The tables used downstream were written by a run with `n_neighbors = 30`; the
  archived notebook cell shows 15, from a later re-run.
- **Hotspot re-runs differ slightly.** Hotspot picks tied neighbours at random and the notebook sets no
  seed, so every run writes slightly different tables. Module 5 is the same module in every run; four new
  runs gave 1,854 to 1,891 bins with `module_5 > 4` (archived run: 1,882). See
  [`../environment/README.md`](../environment/README.md#notes).
- **RCTD output name.** Line 59 of `rctd_script_NI.R` names `rctd_MNP1_2_Tcell_splitup_no_sinusoids.rds`;
  the result used in step 5 is `rctd_MNP1_2_Tcell_combined_no_sinusoids.rds`, made with the reference on
  line 25.
- **Region outline.** The region is a kernel-density contour; a comment in the script calls it a convex
  hull.

Evidence for these notes, job ids, versions and other details: [`../docs/PROVENANCE.md`](../docs/PROVENANCE.md#part-3-stereo-seq-v13).
