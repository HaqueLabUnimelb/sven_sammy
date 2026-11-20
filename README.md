# sven_sammy

Analysis scripts of the Haque Lab (University of Melbourne) for one mouse liver project, archived for
publication. The repository covers three parts:

1. single-cell RNA-seq (scRNA-seq) of liver CD4 T cells and mononuclear phagocytes (MNPs): Cell Ranger,
   clustering and cell annotation;
2. CellChat ligand-receptor analysis between MNP2 cells and CD4 T cells, with all parameters and every
   plotting script;
3. Stereo-seq (STOmics) v1.3 of a liver section: RCTD deconvolution, Hotspot gene modules, DBSCAN regions
   and cell-type enrichment.

Terms used in all parts:

- `MNP1` and `MNP2` are two clusters of MNPs in the liver scRNA-seq data. `naive CD4`, `infected CD4`,
  `cytotoxic CD4`, `effector CD4`, `proliferating CD4` and `uninfected MNP` are further cluster labels.
  Part 1 lists which cluster carries which label.
- The sample names of the scRNA-seq library distinguish infected and uninfected mice. The archived files
  do not name the infection.
- The five files with `_NI` in their name are kept by their author in one folder on Spartan,
  `/data/projects/punim1236/liver_stomicsv1.3/scripts_for_NI/`. Comments in the part 1 scripts speak of
  the "NI revision"; the archived files do not expand `NI` and do not name the manuscript.

The archived files are byte-for-byte copies of the files on the lab's storage on the Spartan cluster
(University of Melbourne), copied on 20 November 2026. They were not edited for this archive. Four files
come from records in the home folder of the account `czhu3869` instead: the contents of an RStudio editor
buffer, a SAW script from the edit history of the code-server editor (both copied byte for byte), and the
text of two commands taken from the log of the session that ran them. The folder READMEs say which files
these are.
`MANIFEST.tsv` lists, for every archived file, its original path, its owner (cluster account), its
last-modified time and its md5. Six scripts were last saved by their author on 19 November 2026, after the
results had been produced (`liver_scRNAseq_processed_NI.R`, `sven_liverscRNAseq.R`, `rctd_script_NI.R`,
`hotspot_liverstomics_NI.ipynb`, `liverstomics_hotspot_dbscan_NI.R`, `colocalisation_script_NI.R`). The
folder READMEs list the differences from the results that could be identified.

## Where the scripts are

| Part | Folder | Main files |
|---|---|---|
| 1. scRNA-seq: Cell Ranger, clustering, cell annotation | [`1_scRNAseq_preprocessing_clustering_annotation/`](1_scRNAseq_preprocessing_clustering_annotation) | `1_cellranger/` (Cell Ranger 8.0.1 input tables and command records), `2_seurat_clustering_annotation/liver_scRNAseq_processed_NI.R` (clustering, annotation, UMAPs, dot plot), `sven_liverscRNAseq.R` and `Tcell_combined_ref.R` (annotated objects for RCTD and CellChat). No script was found for the hashtag demultiplexing and first clustering. |
| 2. CellChat | [`2_CellChat/`](2_CellChat) | `run_cellchat_MNP2_infected_vs_naive_CD4_latest.R`, `run_cellchat_MNP2_cytotoxic_vs_naive_CD4_latest.R`, their SLURM wrappers and logs, and the two output folders with the run settings, result tables and figures; `commands_for_two_extra_tables/`; `earlier_analysis_2025-11/` (code and one result table of the November 2025 CellChat analysis of the same cells) |
| 3. Stereo-seq v1.3 | [`3_STOmics_liver_v1.3/`](3_STOmics_liver_v1.3) | `1_SAW/` (SAW 8.2.2), `2_bin20_matrix/` (bin20 count matrix), `3_RCTD/rctd_script_NI.R`, `4_Hotspot/hotspot_liverstomics_NI.ipynb`, `5_DBSCAN_regions_enrichment/liverstomics_hotspot_dbscan_NI.R` and `colocalisation_script_NI.R` |

Each folder has its own README with one row per file, the parameters, and notes on the archived state.

## How to reproduce

The parts depend on each other in this order:

```
scRNA-seq object GEX_HTO_Processed.rds --> part 1 R scripts --> annotated Seurat objects
annotated Seurat objects --> CellChat (part 2)
annotated Seurat objects --> RCTD reference (part 3, step 3)
Stereo-seq: SAW --> bin20 matrix --> RCTD (part 3, steps 1 to 3)
            bin20 matrix --> Hotspot (part 3, step 4)
            RCTD + Hotspot --> DBSCAN regions --> enrichment, colocalisation (part 3, step 5)
```

1. Get the input data (next section). The scripts use absolute paths on Spartan, for example
   `/data/projects/punim1236/liver_scRNAseq/`; replace them with your own paths.
2. Part 1: the R scripts in `2_seurat_clustering_annotation/` start from `GEX_HTO_Processed.rds` and are
   run section by section in R with Seurat 5. The script that built this object was not found, and the
   archived Cell Ranger records (`1_cellranger/`) belong to a later run with different counts, so the step
   from FASTQ to this object cannot be reproduced exactly.
3. Part 2: set `input_rds` and `outdir` in each CellChat script and run it, directly or through its SLURM
   wrapper (R 4.5.1, CellChat 2.2.0.9001, 8 cores, 32 GB, about 12 minutes per run).
4. Part 3: follow the five steps in [`3_STOmics_liver_v1.3/README.md`](3_STOmics_liver_v1.3/README.md)
   (SAW 8.2.2, stereopy, R 4.5.0 with spacexr 1.2.0, Python 3.9 with hotspotsc 1.1.3, R with dbscan and
   MASS).

Most R scripts and the notebook were written for interactive use and are run section by section, not from
top to bottom. The README of each folder says which lines to run and which lines need attention.

## Data that is not in this repository

No sequencing data, microscopy images or R/Python data objects are stored here. Part 2 includes its result
tables and figures; part 3 includes six small result tables. `DATA_NOT_IN_REPOSITORY.tsv` lists each input
and result object with its path on Spartan, size and md5. The main ones:

| Part | Data | Size |
|---|---|---|
| 1 | Sequencing run folder (BCL) and FASTQ files of the scRNA-seq library | 47 GB and 20 GB |
| 1 | `GEX_HTO_Processed.rds`: Seurat object with all cells, the starting point of the R scripts | 186 MB |
| 1, 2, 3 | `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS`: annotated liver cells; input of CellChat run 1 and of RCTD | 48 MB |
| 1, 2 | `sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS`: annotated liver cells with three CD4 subset labels; input of CellChat run 2 | 48 MB |
| 3 | Stereo-seq FASTQ files, chip mask and DAPI image of chip `Y01453J8` | 162 GB, 3.1 GB, 2.5 GB |
| 3 | `Y01453J8_tissue_bin20_raw_counts.h5ad` and `.rds`: bin20 count matrix | 149 MB and 730 MB |
| 3 | `droplet_Liver_seurat_tiss.Robj`: Tabula Muris droplet liver object (hepatocyte reference) | 214 MB |
| 3 | `rctd_MNP1_2_Tcell_combined_no_sinusoids.rds`: RCTD result | 83 MB |
| 3 | `hotspot_spot_scores.csv`: Hotspot module scores per bin | 43 MB |

## Open points

These items are not in the archive or do not agree with each other. Details are in the folder READMEs.

- Part 1: the R script that created `GEX_HTO_Processed.rds` and the Cell Ranger output it read were not
  found. The code that saved `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS` was not found. Four SLURM
  scripts (the bcl2fastq script and three Cell Ranger scripts) are not archived because they contain a
  personal e-mail address; their commands are given in the part 1 README or in the archived Cell Ranger
  command records.
- Part 2: the code of the November 2025 CellChat analysis was found only in its last edited state (an
  RStudio editor buffer, archived); the code of the June 2025 CellChat analyses was not found.
- Part 3: the Hotspot notebook shows `n_neighbors=15`, while the result tables used downstream were
  written by a run with `n_neighbors=30`. The RCTD script names an output file that differs from the file
  read in the next step.
