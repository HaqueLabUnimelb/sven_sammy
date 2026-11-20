# Provenance notes

This file holds the details behind the three part READMEs: where each archived file came from, how the
archived scripts relate to the results, which files were recovered from editor records, session logs or job records, and what
could not be found. The part READMEs give the short version.

Contents:

- [About the archived files](#about-the-archived-files)
- [Part 1: scRNA-seq](#part-1-scrna-seq)
- [Part 2: CellChat](#part-2-cellchat)
- [Part 3: Stereo-seq v1.3](#part-3-stereo-seq-v13)
- [Files recovered from records](#files-recovered-from-records)
- [Not found](#not-found)

## About the archived files

- The analysis scripts and results are byte-for-byte copies of the files on the lab's storage on the
  Spartan cluster (University of Melbourne), copied on 20 November 2026. They were not edited for this
  repository. Four files come from records in the home folder of the account `czhu3869` instead (see
  [Files recovered from records](#files-recovered-from-records)).
- Six scripts were last saved by their author on 19 November 2026, after the results had been produced
  (6 to 27 October 2026): `liver_scRNAseq_processed_NI.R`, `sven_liverscRNAseq.R`, `rctd_script_NI.R`,
  `hotspot_liverstomics_NI.ipynb`, `liverstomics_hotspot_dbscan_NI.R`, `colocalisation_script_NI.R`. The
  differences from the results that could be identified are listed below.
- The five files with `_NI` in their name are kept by their author in the folder
  `/data/projects/punim1236/liver_stomicsv1.3/scripts_for_NI/` on Spartan. Comments in the part 1 scripts
  speak of the "NI revision"; the archived files do not expand `NI` and do not name the manuscript.
- The cell labels (`naive CD4`, `infected CD4`, `MNP1`, `MNP2`, ...) are cluster names of the scRNA-seq
  analysis; the part 1 README has the table. The sample names of the scRNA-seq library distinguish
  infected and uninfected mice; the archived files do not name the infection.

## Part 1: scRNA-seq

### From FASTQ to the analysed object

The R scripts of this part start from one Seurat object, `GEX_HTO_Processed.rds` (9,461 cells, processed on
2 June 2025). The script that built this object and the Cell Ranger output it read were not found.
`1_cellranger/` holds the records of a later Cell Ranger run of the same sequencing run (November 2025).
The count matrix of that later run is not the one in the object: 9,368 cell barcodes are shared, and the
per-cell UMI totals in the object are about 1.5 times those of the November matrix. The step from FASTQ to
`GEX_HTO_Processed.rds` can therefore not be reproduced exactly from this archive.

### Cell Ranger runs and archived files in detail

| Step | File | What it does |
|---|---|---|
| 1. FASTQ and counting, first run (May to June 2025) | not archived here | See "First run" below the table. |
| 1. FASTQ and counting, second run (November 2025) | `1_cellranger/2025-11_mkfastq_and_count/samplesheet.csv`, `Split_Output_new/_invocation`, `Split_Output_new/_versions` | `cellranger mkfastq` 8.0.1: sample indices `SI-TT-C11` (gene expression, `GEX_C11`) and `SI-NT-B1` (hashtags, `CSP_B1`). `_invocation` is the command record written by Cell Ranger. The submission script (SLURM job 18789607) is kept in the SLURM accounting database of Spartan: it loads `CellRanger/8.0.1`, `GCC/11.3.0` and `bcl2fastq2/2.20.0` and runs `cellranger mkfastq --id=Split_Output_new --run=/data/projects/punim1236/liver_scRNAseq/250527_VH00915_604_AAH2VJHM5-10042-micheleC --csv=samplesheet.csv --localcores=24 --localmem=128`. It is not archived as a file because its header contains a personal e-mail address. |
| | `1_cellranger/2025-11_mkfastq_and_count/library.csv`, `new_cellranger/_cmdline`, `_invocation`, `_versions` | `cellranger count` 8.0.1 with both libraries (`Gene Expression` and `Antibody Capture`), reference `refdata-gex-mm10-2020-A`, `--create-bam=false`, 24 cores, 128 GB (SLURM job 18810731). `_cmdline` is the full command. The SLURM wrapper `run.slurm` that submitted it is not archived here because its header contains a personal e-mail address. |
| | `1_cellranger/2025-11_mkfastq_and_count/new_cellranger/outs/feature_reference.csv`, `outs/metrics_summary.csv` | The hashtag reference that the count run used, and its summary metrics (9,560 cells). |
| | `../2_CellChat/earlier_analysis_2025-11/rstudio_editor_buffer_2025-11-21.R`, lines 23 to 98 | Seurat processing of the November 2025 count matrix into `new_livescRNA.RDS`: cells kept with `nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10` (lines 57 to 62), `HTODemux(seu, assay = "HTO", positive.quantile = 0.99)` (lines 77 to 81), then PCA and UMAP. This object is not used by the archived analyses; the file is described in part 2. |
| 2. Seurat object, hashtag demultiplexing, first clustering | no script in this archive | Result: `GEX_HTO_Processed.rds`. The processing calls recorded inside the object are listed below. |
| 3. Clean-up, liver subset, clustering, annotation, figures | `2_seurat_clustering_annotation/liver_scRNAseq_processed_NI.R` | The script of the "NI revision" (line 254: `# Re-annotate scRNAseq dataset for NI revision and perform dotplot`). Reads `GEX_HTO_Processed.rds`. Removes contaminating clusters (saves `sven_refined.RDS`), keeps the four liver hashtags, clusters again with `FindClusters(resolution = 0.7)` (saves `sven_liver_for_paper.RDS`), draws the gene UMAPs, adds the MNP1, MNP2 and pDC signature scores, assigns the cell-type labels of the revision (line 260) and draws the annotated UMAPs and the dot plot (`Dotplot_NI_revision_v3.pdf`). |
| | `2_seurat_clustering_annotation/gene_signatures/MNP1.txt`, `MNP2.txt`, `pDCs.txt`, `DCs.txt` | Gene lists for `AddModuleScore`. The revision script reads the first three (lines 180 to 182), the working script reads all four. On Spartan they are in `/data/projects/punim1236/liver_scRNAseq/Ollie/`. |
| 4. Annotated objects for RCTD and CellChat | `2_seurat_clustering_annotation/sven_liverscRNAseq.R` | The complete working script (941 lines) from which the revision script was taken. Besides `Tcell_combined_ref.R`, it is the only archived file with code for the annotated reference objects. Line 414 labels the resolution-0.7 clusters; line 417 assigns the labels again with the MNP clusters merged into one label `MNP`; lines 427 to 430 remove the DC and pDC clusters and save `sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS`. Lines 860 to 872 read `sven_liver_Tcells_combined_forRCTD.RDS`, assign the labels `MNP1`, `MNP2`, `uninfected MNP` and the CD4 subsets, and save `sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS`. |
| | `2_seurat_clustering_annotation/Tcell_combined_ref.R` | Reads `sven_liver_new_cytotoxic_cluster_forRCTD.RDS`, merges the four CD4 labels into `Tcell` and saves `sven_liver_Tcells_combined_forRCTD.RDS`. |

First run (May to June 2025):

- Two scripts of this run are in the lab repository `HaqueLabUnimelb/Cellranger` (private; commits of
  2 June 2025). `bcl_to_fastq.slurm` runs bcl2fastq2 2.20.0:
  `bcl2fastq --runfolder-dir /data/projects/punim1236/liver_scRNAseq/250527_VH00915_604_AAH2VJHM5-10042-micheleC --output-dir /data/projects/punim1236/liver_scRNAseq/fastq_output --no-lane-splitting`.
  `run.slurm` runs Cell Ranger 8.0.1 for gene expression only (no hashtag library, no feature reference):
  `cellranger count --id=Sample1 --transcriptome=/data/projects/punim1236/liver_scRNAseq/mm10/refdata-gex-mm10-2020-A --fastqs=/data/projects/punim1236/liver_scRNAseq/fastq_output --sample=Sample1 --localcores=32 --localmem=128 --create-bam true`.
- Two attempts of 2 June 2025 are on Spartan: `Sample1_HTO/` (gene expression and hashtags together,
  failed) and `Sample1_HTOonly/` (hashtag library only, interrupted). The run that produced the matrix
  behind `GEX_HTO_Processed.rds` was not found.
- FASTQ files converted on 31 May 2025 (with the bcl2fastq folders `Reports/` and `Stats/`) are in
  `/data/projects/punim1236/liver_scRNAseq/liver_scRNAseq_raw/`, next to `feature_ref.csv` (last changed
  2 June 2025, 14:03), which is byte-identical to the hashtag reference of the November 2025 run
  (`new_cellranger/outs/feature_reference.csv`).

### Processing recorded in `GEX_HTO_Processed.rds`

Read from the object on 20 November 2026: `object@commands` (all calls time-stamped 2 June 2025, 16:30;
the file was saved on 3 June 2025) and the metadata.

| Step | Recorded call or content |
|---|---|
| Genes and cells | 32,285 genes, 9,461 cells; assays `RNA` and `HTO` (8 hashtags) |
| Hashtag normalisation | `NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")` |
| Hashtag demultiplexing | Metadata columns `HTO_maxID`, `HTO_secondID`, `HTO_margin`, `HTO_classification`, `HTO_classification.global`, `hash.ID` (the columns written by Seurat `HTODemux`): 7,425 singlets, 1,989 doublets, 47 negatives. Doublets and negatives are still in this object; the later scripts drop them by selecting hashtag names. |
| RNA normalisation | `NormalizeData(seurat_obj)` |
| Variable genes | `FindVariableFeatures(seurat_obj, nfeatures = 2000)` |
| Scaling and PCA | `ScaleData(seurat_obj)`, `RunPCA(seurat_obj, npcs = 30)` |
| Neighbours and clusters | `FindNeighbors(seurat_obj, dims = 1:20)`, `FindClusters(seurat_obj, resolution = 0.3)`: 18 clusters (0 to 17) |
| UMAP | `RunUMAP(seurat_obj, dims = 1:20)` |

The quality-control thresholds and the `HTODemux` quantile are not recorded in the object. The object was
saved by R 4.4.0 (version stamp of the `.rds` file), as were the annotated objects used in parts 2 and 3;
the package versions are not recorded.

### Objects and labels

`liver_scRNAseq_processed_NI.R` keeps clusters 0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13 and 14 of the first
clustering (8,319 cells, `sven_refined.RDS`). `sven_liver_for_paper.RDS` holds the 2,934 cells of these
clusters, except cluster 3, that carry one of the four liver hashtags. They are clustered at resolution 0.7
into 15 clusters (0 to 14). MNP1 and MNP2 are clusters 1 and 6 of this clustering. Removing the DC and pDC
clusters (4, 8, 9, 12, 13, 14) leaves the 2,537 cells of the objects used for RCTD and CellChat. Labels by
resolution-0.7 cluster, read from the saved objects:

| Cluster (cells) | Revision figures, `liver_scRNAseq_processed_NI.R` line 260 | `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS` | `sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS` |
|---|---|---|---|
| 0 (754) | naive CD4 | naive CD4 | naive CD4 |
| 1 (491) | MNP1 | MNP1 | MNP1 |
| 2 (364) | uninfected MNP | uninfected MNP | uninfected MNP |
| 3 (312), 5 (170) | effector CD4 | infected CD4 | effector CD4 |
| 6 (167) | MNP2 | MNP2 | MNP2 |
| 7 (163) | effector CD4 | infected CD4 | cytotoxic CD4 |
| 10 (65), 11 (51) | effector CD4 | infected CD4 | proliferating CD4 |
| 4, 8, 12, 13, 14 | DC | not in the object | not in the object |
| 9 | pDC | not in the object | not in the object |

In both saved objects the label of cluster 0 is stored with a leading space (`" naive CD4"`). The labels
`naive CD4` and `infected CD4` name clusters, not mice: all 754 `naive CD4` cells carry hashtag B1, and
735 of the 761 `infected CD4` cells carry B1 or B2, both from infected mice (the other 26 carry B3 or
B4).

### Running the interactive scripts

The scripts were written for interactive use and are run section by section. They select clusters by
number, so they need the saved objects of the original analysis; an object rebuilt from FASTQ cannot be
expected to carry the same cluster numbers.

1. Cell Ranger 8.0.1 (the November 2025 run): run `cellranger mkfastq` with the arguments in
   `Split_Output_new/_invocation` (`all_mkfastq_args`) and `samplesheet.csv`. Then run the command in
   `new_cellranger/_cmdline` after editing the paths in it and in `library.csv`; point `--feature-ref` to
   `new_cellranger/outs/feature_reference.csv`.
2. Clustering and annotation: open `liver_scRNAseq_processed_NI.R` in R with Seurat 5 and tidyverse and
   replace `/data/projects/punim1236/liver_scRNAseq/` (also written as `/data/gpfs/projects/...`) with
   your folder; the gene lists are read from the subfolder `Ollie/` of that folder (lines 180 to 182).
   - `sven_liver_for_paper.RDS` on Spartan (2,934 cells, no `annotation` column) corresponds to
     `sven.liver` as defined on line 65. Line 80 defines `sven.liver` again with first-clustering cluster 3
     included (2,985 cells). To obtain the 2,934-cell object, run line 65 and then lines 98 to 104 without
     running lines 80 to 95. Line 114 then cannot be run, because it plots the column `annotation`.
   - Lines 229 to 251 use the objects `immune`, `MNP1_` and `MNP2_`, which this script does not define. In
     the working script `MNP1_` and `MNP2_` are the gene lists (lines 117 and 118) and `immune` is read
     from one of the reference objects (lines 874 to 878). The violin plots compare MNP1 and MNP2, so they
     need one of the two `..._MNP1_MNP2_forRCTD.RDS` objects; the scripts do not record which one was
     used.
   - Line 256 reads `sven_liver_for_paper.RDS` again, without the signature scores. Run lines 184 to 209
     again before the dot plot (lines 275 to 284), which shows the scores `MNP1`, `MNP2` and `pDCs`.
3. Reference objects: `Tcell_combined_ref.R` starts from `sven_liver_new_cytotoxic_cluster_forRCTD.RDS`
   (2,537 cells) and saves `sven_liver_Tcells_combined_forRCTD.RDS`; lines 860 to 872 of
   `sven_liverscRNAseq.R` then assign the labels and save the "split up" object. The working script cannot
   be run from top to bottom (lines 201 and 224 are incomplete expressions), and it contains further
   analyses that are outside the three topics of this archive.

### Not in the repository

- The R script that created `GEX_HTO_Processed.rds`, and the Cell Ranger output it read (see above).
- The code that saved `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS`. No script that saves it was
  found; the labels of this object are in the table above.
- The code that saved `sven_liver_new_cytotoxic_cluster_forRCTD.RDS`, the object that
  `Tcell_combined_ref.R` reads. Its labels are those of line 414 of the working script (`MNP infected` and
  `naive MNP` kept apart), but lines 427 to 430 save a different file (`..._forRCTD_v2.RDS`, with the
  merged label `MNP`). Both files hold the same 2,537 cells.
- `sven_liver_Tcells_combined_forRCTD.RDS` on Spartan is dated 5 August 2026, before the last save of
  `Tcell_combined_ref.R` (18 September 2026), and has a metadata column `annotation_combined` that no
  archived script creates. Its `annotation` labels are the ones the script assigns.
- The SLURM wrappers of the Cell Ranger runs, which contain a personal e-mail address (see the table).
- Seurat objects, count matrices and FASTQ files (listed in [`../DATA_NOT_IN_REPOSITORY.tsv`](../DATA_NOT_IN_REPOSITORY.tsv)).
- R and Seurat versions of steps 2 to 4: not recorded in the scripts. The saved objects carry the
  SeuratObject version stamp 5.1.0.

## Part 2: CellChat

### The two runs in detail

| | Run 1: infected CD4 vs naive CD4 | Run 2: cytotoxic CD4 vs naive CD4 |
|---|---|---|
| R script | `run_cellchat_MNP2_infected_vs_naive_CD4_latest.R` | `run_cellchat_MNP2_cytotoxic_vs_naive_CD4_latest.R` |
| SLURM wrapper | `run_cellchat_MNP2_infected_vs_naive_CD4_existing_latest.sbatch` | `run_cellchat_MNP2_cytotoxic_vs_naive_CD4_latest.sbatch` |
| Input Seurat object | `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS` | `sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS` |
| Cell groups used (`annotation` column) | `MNP2` (167 cells), `infected CD4` (761), `naive CD4` (754) | `MNP2` (167 cells), `cytotoxic CD4` (163), `naive CD4` (754) |
| Output folder | `MNP2_CD4_20261023/` | `MNP2_cytotoxicCD4_vs_naiveCD4_20261023/` |
| SLURM job and log | 33696088, `slurm_logs/cellchat_mnp2_cd4_existing_latest_33696088.out/.err` | 33701062, `slurm_logs/cellchat_mnp2_cyto_cd4_latest_33701062.out/.err` |
| Run time, peak memory | 12 min 32 s, 13.5 GiB | 10 min 58 s, 13.1 GiB |
| Ligand-receptor pairs in the result, CD4 to MNP2 | 51 (37 higher with infected CD4, 14 higher with naive CD4) | 41 (23 higher with cytotoxic CD4, 18 higher with naive CD4) |
| Ligand-receptor pairs in the result, MNP2 to CD4 | 50 (39 and 11) | 42 (25 and 17) |

The two R scripts are the same code. They differ only in the input file, the output folder name and the
group label (`infected CD4` in run 1, `cytotoxic CD4` in run 2). The input objects hold 32,285 genes and
2,537 cells. They differ in how the CD4 T-cell clusters other than the naive cluster are labelled: one
label `infected CD4` in the "combined" object, three labels `cytotoxic CD4`, `effector CD4` and
`proliferating CD4` in the "split up" object.

### What each script does

1. Reads the Seurat object and uses the `RNA` assay and the `annotation` column.
2. Writes the cell counts per label (`01_input_annotation_counts.csv`), the counts per label and hashtag
   (`02_input_annotation_by_HTO.csv`) and the run settings (`00_run_config.txt`).
3. Builds two CellChat objects: `MNP2` with the infected (or cytotoxic) CD4 group, and `MNP2` with `naive CD4`.
   Each is analysed on its own with the parameters in the table below.
4. Extracts the significant ligand-receptor pairs (p-value threshold 0.05) in each direction (CD4 to MNP2,
   MNP2 to CD4) and joins the two objects pair by pair. A pair that is not significant in one object gets
   probability 0 in that object; the column `detected_in` marks these pairs. For every pair the script
   computes `log(prob + 1)` in each object. The value plotted is `+log(prob + 1)` of the infected (or
   cytotoxic) object when its probability is the higher one, and `-log(prob + 1)` of the naive object
   otherwise. Tables: `03_differential_LR_CD4_to_MNP2_*.csv` and `04_differential_LR_MNP2_to_CD4_*.csv`.
5. Draws the bar charts of step 4: all pairs (`05_`, `07_`) and the 25 largest per side (`06_`, `08_`).
6. Lifts both objects to the three common groups, merges them and draws the standard CellChat comparison
   plots: total number and strength of interactions (`09_`), differential interaction circle plots by
   strength and by number (`10_`, `11_`), differential heat maps by strength and by number (`12_`, `13_`).
7. Writes `99_sessionInfo.txt`.

Every figure is saved as PDF and as PNG (300 dpi).

### All parameters

Values are the ones written in the scripts. The CellChat objects saved by the runs store the same
`computeCommunProb` values in `object@options$parameter` (read back on 20 November 2026).

| Step | Function | Setting |
|---|---|---|
| Random seed | `set.seed(1)` | set at the start of the script |
| Parallel workers | `future::plan("multisession")` | 8 (`SLURM_CPUS_PER_TASK`) |
| Object | `createCellChat(group.by = "cellchat_group", assay = "RNA")` | two cell groups per object |
| Database | `subsetDB(CellChatDB.mouse, search = c("Secreted Signaling", "Cell-Cell Contact"))` | 1,942 interactions: 1,280 secreted signalling, 662 cell-cell contact |
| Over-expressed genes | `identifyOverExpressedGenes(do.fast = FALSE)` | other arguments at CellChat defaults |
| Over-expressed interactions | `identifyOverExpressedInteractions()` | defaults. Pairs kept: 317 (MNP2 + infected CD4), 176 (MNP2 + cytotoxic CD4), 347 (MNP2 + naive CD4) |
| Communication probability | `computeCommunProb(type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = TRUE, nboot = 100, seed.use = 1L)` | stored defaults: `Kh = 0.5`, `n = 1`, no spatial or contact-dependent options |
| Filter | `filterCommunication(min.cells = 10)` | |
| Pathway level and network | `computeCommunProbPathway()`, `aggregateNet()`, `rankNetPairwise()`, `netAnalysis_computeCentrality(slot.name = "netP")` | defaults |
| Pairs reported | `subsetCommunication(sources.use, targets.use, thresh = 0.05)` | p-value threshold 0.05 |
| Merge for comparison plots | `liftCellChat()` to `MNP2`, the infected (or cytotoxic) CD4 group and `naive CD4`; `mergeCellChat(merge.data = FALSE, cell.prefix = TRUE)` | |

### Software

- R 4.5.1, run with the `Rscript` of the conda environment `/data/gpfs/projects/punim2394/conda/rctd_env`
  (`99_sessionInfo.txt`).
- CellChat 2.2.0.9001 (`99_sessionInfo.txt`, `00_run_config.txt`), library
  `/data/gpfs/projects/punim2394/AAA_Sorry_Its_Chen_Here/Rlib_new`. The installed package records its
  source as GitHub `jinworks/CellChat`, commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`.
- Seurat 5.4.0, SeuratObject 5.3.0, ggplot2 4.0.2, patchwork 1.3.2, future 1.70.0, ComplexHeatmap 2.26.1,
  dplyr 1.2.0 (`99_sessionInfo.txt`).

### Notes on the archived runs

- **Earlier analysis in the run settings.** One line of `00_run_config.txt` refers to the earlier analysis
  (see "Earlier CellChat analysis (November 2025)" below): `database=CellChatDB.mouse subset to Secreted
  Signaling and Cell-Cell Contact, matching the existing 2025-11-19 CellChat object`. The line
  `identifyOverExpressedGenes=do.fast FALSE to match the historical scripts/templates` does not name those
  scripts; the November 2025 code calls `identifyOverExpressedGenes()` with its default arguments, and
  `do.fast = FALSE` is written in a CellChat script of another lab project
  (`/data/projects/punim1236/b_cells/CellChat/IdeasGrant/CellChat_Ideasgrant.R`, line 114). The output
  folder name ends in `matched_original_params` for these reasons. `FigS2E_style` and `screenshot` refer to
  the bar-chart layout of panel E of a Figure S2 ("Cytotoxic vs Naive CD4+ -> MNP") that served as the
  template; that figure is not in this archive.
- **`latest` and `existing` in file names.** The default CellChat library written in the R scripts,
  `Rlib_cellchat_latest`, does not exist. The wrappers point `CELLCHAT_LIB` to the library that was already
  installed (`Rlib_new`, CellChat 2.2.0.9001).
- **Output folder names.** The scripts write to `cellchat_latest_MNP2_CD4_20261023_matched_original_params/`
  and `cellchat_latest_MNP2_cytotoxicCD4_vs_naiveCD4_20261023_matched_original_params/` (see `output_dir` in
  `00_run_config.txt`). The folders were renamed afterwards to `MNP2_CD4_20261023/` and
  `MNP2_cytotoxicCD4_vs_naiveCD4_20261023/`.
- **SLURM job state.** Both jobs are recorded by SLURM as `FAILED` with exit code 139. R crashed with a
  segmentation fault while it was shutting down, after the script had written all outputs and logged
  `DONE` (end of the `.out` log). The same crash happens in this R environment when a script that only
  loads CellChat and reads the saved objects ends (checked on 20 November 2026).
- **Log format.** Each `log_line()` call of the script also prints the word `TRUE` after its line (the
  argument `flush = TRUE` is printed by `cat()`), so many lines of the `.out` logs start with `TRUE` and
  the logs end with `TRUE`. It has no meaning.
- **Two tables that the scripts did not write.** `MNP2_CD4_20261023/example_LR_expression_support.csv`
  and `MNP2_CD4_20261023/MNP2_vs_original_combinedMNP_comparison.csv` were written about one hour after
  run 1 by two separate commands; see "Commands that wrote the two extra tables" below.
- **Not archived here.** The CellChat result objects (`cellchat_*.rds`, 1.4 to 4.9 MB each) and
  `99_sessionInfo.rds` stay on Spartan; rerunning the scripts recreates them.

### Commands that wrote the two extra tables

The folder `commands_for_two_extra_tables/` holds the exact text of the two commands, recovered from the
log of the session that ran them on 23 October 2026 in `/data/scratch/projects/punim1236/AAA_Chen_is_here`.
They write into the output folder under its name at that time
(`cellchat_latest_MNP2_CD4_20261023_matched_original_params/`, now `MNP2_CD4_20261023/`).

| File | Run at | What it does |
|---|---|---|
| `make_MNP2_vs_original_combinedMNP_comparison.sh` | 11:13 | Python (pandas). Takes from `earlier_analysis_2025-11/AAA_complete_altogether.csv` the pairs sent by `cytotoxic CD4` or `naive CD4` to `MNP infected` (clusters MNP1 and MNP2 together), joins them with `03_differential_LR_CD4_to_MNP2_infected_minus_naive.csv` of run 1, and writes `MNP2_vs_original_combinedMNP_comparison.csv`. |
| `make_example_LR_expression_support.sh` | 11:33 | R, with the `Rscript` of `/data/scratch/projects/punim1236/python_envs/rseurat`. For Ifng, Ifngr1, Ifngr2, Thy1, Adgre5, Cd6 and Alcam: mean log-normalised expression (Seurat `data` layer) and percentage of cells with expression above 0 per label of `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS`; writes `example_LR_expression_support.csv`. |

On 20 November 2026 both commands were run again. For the re-run the input path of the Python command was
changed to the renamed folder and both output paths were changed to a scratch folder; nothing else was
changed. Both tables were reproduced byte for byte.

### Earlier CellChat analysis (November 2025)

`earlier_analysis_2025-11/rstudio_editor_buffer_2025-11-21.R` is a copy of the contents of an RStudio editor
buffer that RStudio kept in `~/.local/share/rstudio/sources/` of the account `czhu3869` (last changed on
21 November 2025). RStudio's record of the file name and saved state of this buffer is gone, and no saved
file with this content was found. No other code of the November 2025 CellChat analysis was found in the
home folder of `czhu3869` or in the project folders. The buffer also holds other work:

| Lines | What they do |
|---|---|
| 1 to 21 | library calls and package installation |
| 23 to 98 | Seurat processing of the November 2025 Cell Ranger run into `new_livescRNA.RDS` (see part 1; this object is not used by the archived analyses) |
| 106 to 144 | CellChat on `sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS`, cell groups = column `annotation`; database `CellChatDB <- CellChatDB.mouse` and `subsetDB(CellChatDB, search = c("Secreted Signaling", "Cell-Cell Contact"), key = "annotation")`; `identifyOverExpressedGenes(cellchat)` with default arguments; `computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, population.size = TRUE)`; `filterCommunication(cellchat, min.cells = 10)`; saves `2025_11_19_newCellchat_combined_MNP.rds` |
| 146 to 160 | reads `cellchat_humanSkin_LS.rds` (a CellChat object of June 2025), draws network plots, and writes `subsetCommunication(cellchat, slot.name = "net")` to `AAA_complete_altogether.csv` |
| 162 to 334 | bar charts of `log(prob + 1)` comparing `cytotoxic CD4` and `naive CD4` as senders to, and receivers from, `MNP` (black and grey bars, then red and blue, then black and grey); the last chart is saved as `interaction_bar.pdf` |
| 344 to 390 | further CellChat plots of an object whose groups are hashtag names |
| 398 to 525 | RCTD on a different spatial data set (Curio Seeker liver pucks); outside the three topics of this archive |

The buffer holds the last state of the code. Two CellChat objects of this analysis are on Spartan:

- `2025_11_19_newCellchat_combined_MNP.rds` (21 November 2025): groups naive CD4 (stored with a leading
  space, `" naive CD4"`), cytotoxic CD4, effector CD4, proliferating CD4 and `MNP` (1,022 cells), the labels of `..._forRCTD_v2.RDS`; this is the object
  that lines 106 to 144 write.
- `2025_11_19_newCellchat.rds` (19 November 2025): the same 2,537 cells with `MNP infected` (658) and
  `naive MNP` (364) kept apart, the labels of `sven_liver_new_cytotoxic_cluster_forRCTD.RDS`. The version of
  the code that read this file is not preserved; the buffer reads `_v2` at line 106.

Both objects store the `computeCommunProb` settings `truncatedMean`, `trim = 0.1`, `raw.use = TRUE`,
`population.size = TRUE`, `nboot = 100`, `seed.use = 1` and the 1,942-interaction database subset.
`earlier_analysis_2025-11/AAA_complete_altogether.csv` (19 November 2025) is identical, byte for byte, to
`write.csv(subsetCommunication(<2025_11_19_newCellchat.rds>, slot.name = "net"), row.names = FALSE)`
(checked on 20 November 2026). The table was therefore written from the 19 November object, in the way
lines 159 and 160 do it, and not from the object that line 146 reads. Which code drew the Figure S2E panel
shown in the template screenshot is not recorded.

Compared with the October 2026 runs, the November 2025 code reads `..._forRCTD_v2.RDS`, in which the MNP
clusters form one group `MNP` (October 2026: `MNP2` alone, with `infected CD4` or `cytotoxic CD4`), uses
all cell groups of the object in one CellChat object (October 2026: two groups per object) and calls
`identifyOverExpressedGenes()` with default arguments (October 2026: `do.fast = FALSE`). The database
subset and the stored `computeCommunProb` settings are the same.

CellChat objects of June 2025 (2,555 or 1,284 cells, earlier labels or hashtag names) are also on Spartan.
Their code was not found; an R history file of the account `czhu3869` holds only fragments, such as the
line that saved `cellchat_humanSkin_LS.rds`.

## Part 3: Stereo-seq v1.3

### Archived files in detail

| Step | File | What it does |
|---|---|---|
| 1. SAW | `1_SAW/run_2026-07-06_original_image/config.yaml` | Settings and the full `saw count` command that SAW 8.2.2 recorded for the run of 6 July 2026 (SLURM job 26830116): kit version `Stereo-seq T FF V1.3`, sequencing type `PE75_50+100`, organism mouse, tissue liver, the original DAPI image, 32 threads, memory 200. Its output `Y01453J8.tissue.gef` is the input of step 2. |
| | `1_SAW/run_2026-07-06_original_image/saw_liver_26830116.log` | SLURM log of that run. |
| | `1_SAW/run_2026-07-06_original_image/SAWv1.3_pipeline_modified.slurm` | SLURM script that submitted that run: `saw count` with the original DAPI image, 32 CPUs, 200 GB and `--memory 200`, output in `/data/projects/punim1236/AAA_Chen_is_here/Ollie_SAW/output_liver_stomics` (the output now sits in `/data/projects/punim1236/output_liver_stomics/`). The file was later removed from Spartan; this copy is the version kept in the edit history of code-server (the browser version of the VS Code editor) in the home folder of `czhu3869`, saved at 11:25:29 on 6 July 2026, 20 seconds before the job was submitted. The settings it prints match the log. |
| | `1_SAW/SAWv1.3_pipeline.sh` | SLURM script that runs `saw count`. This is the script of a second SAW run of the chip (10 July 2026, job 27045467). It differs from the script of the first run in the image folder (`IMAGEDIR`, a contrast-enhanced copy of the DAPI image), the working and output folder (`WORKDIR`), the SLURM request (partition `long`, 16 CPUs, 100 GB, 31 days instead of `sapphire`, 32 CPUs, 200 GB, 48 hours) and `--memory 300`. The output of the second run is not used in steps 2 to 5. |
| 2. bin20 matrix | `2_bin20_matrix/run_Y01453J8_bin20_RCTD.slurm` | SLURM wrapper that runs the three scripts below, one after the other. Each step is skipped when its output file exists. |
| | `2_bin20_matrix/01_Y01453J8_bin20_gef_to_h5ad.py` | Reads `Y01453J8.tissue.gef` of the first SAW run with stereopy (`st.io.read_gef(bin_size=20)`) and writes the raw counts as `Y01453J8_tissue_bin20_raw_counts.h5ad`, with the bin coordinates in `obs["x"]` and `obs["y"]`. This file is the input of step 4. The stereopy version of the run is not recorded; the Python environment that the wrapper names held stereopy 1.6.2 and Python 3.8.20 when it was removed on 18 November 2026. A rebuilt environment with these versions writes a byte-identical file ([`../environment/python/RESULTS.md`](../environment/python/RESULTS.md)). |
| | `2_bin20_matrix/02_Y01453J8_bin20_h5ad_to_rds.R` | Converts the h5ad file to a Seurat object, `Y01453J8_tissue_bin20_raw_counts.rds` (no gene or bin filter). This file is the input of step 3. |
| | `2_bin20_matrix/03_Y01453J8_bin20_run_RCTD.R` | A first RCTD run in doublet mode with an earlier reference. It was replaced by step 3 and its result is not used in steps 4 and 5. It is archived because the wrapper calls it. |
| | `2_bin20_matrix/slurm-Y01453J8_bin20_RCTD-27327937.out` | Log of job 27327937 (16 July 2026), which ran `02` and `03` (R 4.5.0, spacexr 1.2.0). `01` had written the h5ad file a few minutes earlier in job 27327611, which then failed; no log of that job was found. |
| 3. RCTD | `3_RCTD/rctd_script_NI.R` | Builds the reference: hepatocytes of the Tabula Muris droplet liver object (1,764 cells) merged with the annotated liver scRNA-seq object `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS` (2,537 cells; labels naive CD4, infected CD4, MNP1, MNP2, uninfected MNP). Runs spacexr on the bin20 object: `createRctd()` with default arguments, then `runRctd(rctd_mode = "full", max_cores = 4)`. |
| | `3_RCTD/RCTD_CLAHE_bin20_liverSTOmics_newenv_project2394.sh` | SLURM wrapper of the RCTD run: R 4.5.0 (`module load GCC/11.3.0 OpenMPI/4.1.4 R/4.5.0`), the R library of step 2, 16 CPUs, 100 GB. It calls the script under its working name `RCTD_CLAHE_bin20_liverSTOmics_newenv.R`, which is the same file as `rctd_script_NI.R`. |
| | `3_RCTD/slurm-32372823.out` | Log of the RCTD run of 6 October 2026 (job 32372823, 19 minutes). It shows the reference (4,301 cells, 35,856 genes, the six cell types). The log does not name the output file; `rctd_MNP1_2_Tcell_combined_no_sinusoids.rds`, the result read in step 5, was written when this job ended. |
| 4. Hotspot | `4_Hotspot/hotspot_liverstomics_NI.ipynb` | Jupyter notebook (Python 3.9.13, hotspotsc 1.1.3, scanpy 1.10.3). Reads the bin20 h5ad file, keeps genes detected in at least 1 % of the bins (3,467 genes), runs Hotspot and writes the gene tables, module tables, per-bin module scores and module maps. Hotspot breaks ties between equally distant neighbours at random and no seed is set, so re-runs give slightly different tables ([`../environment/python/RESULTS.md`](../environment/python/RESULTS.md)). |
| | `4_Hotspot/result_tables/` | Five small tables from `/data/projects/punim1236/liver_stomicsv1.3/hotspot/`, written on 20 October 2026: `hotspot_autocorrelation.csv` (C, Z, p-value and FDR per gene), `hotspot_modules.csv` and `hotspot_gene_modules.csv` (module of each gene), `hotspot_module_genes_ranked.csv` (genes of each module ranked by correlation with the module score) and `module_5_genelist.txt` (the 152 genes of module 5). |
| 5. Regions, enrichment, colocalisation | `5_DBSCAN_regions_enrichment/liverstomics_hotspot_dbscan_NI.R` | Reads the RCTD result and the per-bin Hotspot scores. Selects the bins, runs DBSCAN, defines the regions, runs the colocalisation and the enrichment and draws the figures. |
| | `5_DBSCAN_regions_enrichment/colocalisation_script_NI.R` | Three functions for the colocalisation. `coloc_analysis_replicates_fast()` is the one called by the script above. |
| | `5_DBSCAN_regions_enrichment/result_tables/enrichment_of_celltypes.csv` | The enrichment table of 27 October 2026 (the file that line 562 of the script writes). |

### Notes on the archived state

1. **Hotspot neighbours.** The Hotspot cell of the archived notebook shows `n_neighbors=15` and the output
   of that run (top gene `Gm42418`, Z = 813.08). The tables written to disk on 20 October 2026, including
   `hotspot_autocorrelation.csv` (top gene `Gm42418`, Z = 1080.62) and the per-bin scores read in step 5,
   match the Hotspot cell of an earlier saved notebook,
   `/data/projects/punim1236/b_cells/Untitled3.ipynb` (20 October 2026, 17:12), which has
   `n_neighbors=30`. Its 12 other code cells and their outputs are identical to the archived notebook; it
   lacks the last two code cells, which tidy the column names and write the per-bin scores. A checkpoint
   of 20 October 2026, 18:24
   (`/data/projects/punim1236/b_cells/.ipynb_checkpoints/hotspot_liverstomics-checkpoint.ipynb`) holds
   all cells; its code differs from the archived notebook only in `n_neighbors=30` (the stored output of
   its Hotspot cell ends with a `KeyboardInterrupt`).
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
   name come from a state of the functions with `radius = 100` that is not archived. A file with such a
   state is on Spartan, `/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/colocalisation_for_new_stomics.R`
   (19 October 2026): its first 236 lines equal `colocalisation_script_NI.R` except for `radius = 100` on
   lines 13, 78 and 141. Which file was sourced for those outputs is not recorded. The colocalisation
   samples bins at random and no seed is set before line 443.
5. **Length labels.** File names such as `coloc_region40um_*.pdf` use the label `um` for `radius = 40`,
   which is in coordinate units. No script states the physical length of one unit.
6. **"CLAHE" in file names.** CLAHE is the contrast enhancement of the DAPI image used for the second SAW
   run. The RCTD wrapper and its working script carry the word in their names, but they read the bin20
   matrix of the first SAW run.

### Not in the repository

- Input data and result objects: FASTQ files, chip mask, DAPI image, SAW output, the bin20 h5ad and rds
  files, the Tabula Muris object, RCTD results and the per-bin Hotspot scores
  (listed in [`../DATA_NOT_IN_REPOSITORY.tsv`](../DATA_NOT_IN_REPOSITORY.tsv)).
- Figures. The Hotspot maps are in `/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_plots/`,
  the figures and tables of step 5 in
  `/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/`.
- The R package library of the R session of step 5. The session ran R 4.3.1 (the `Producer` field of its
  PDF files and the version stamp of the `.rds` files it saved); its package versions are not recorded.
  The reconstructed library `dbscan_regions` (R 4.5.0) reproduces the enrichment table byte for byte
  ([`../environment/R/RESULTS.md`](../environment/R/RESULTS.md)).
- The Python environment of step 2 (`stereopy_latest_py38`): removed from Spartan on 18 November 2026;
  rebuilt as `sven_sammy_stereopy` ([`../environment/README.md`](../environment/README.md)).

## Files recovered from records

Four archived files were not taken from the analysis folders but from records kept in the home folder of
the account `czhu3869` on Spartan. `../MANIFEST.tsv` gives the exact origin of each.

| Repository file | Origin | Check |
|---|---|---|
| `3_STOmics_liver_v1.3/1_SAW/run_2026-07-06_original_image/SAWv1.3_pipeline_modified.slurm` | Edit history of the code-server editor: the version of `/data/projects/punim1236/AAA_Chen_is_here/Ollie_SAW/SAWv1.3_pipeline_modified.slurm` saved at 11:25:29 on 6 July 2026 (the file was removed later) | SLURM job 26830116 was submitted 20 seconds later; the settings the script prints equal the job log and `config.yaml` |
| `2_CellChat/earlier_analysis_2025-11/rstudio_editor_buffer_2025-11-21.R` | Contents of an RStudio editor buffer in `~/.local/share/rstudio/sources/` (last changed 21 November 2025) | `AAA_complete_altogether.csv` equals the communication table of `2025_11_19_newCellchat.rds` byte for byte; the stored CellChat settings of both November 2025 objects equal the code |
| `2_CellChat/commands_for_two_extra_tables/make_MNP2_vs_original_combinedMNP_comparison.sh` | Exact command text from the log of the session that ran it (23 October 2026, 11:13) | A re-run on 20 November 2026 reproduced the table byte for byte |
| `2_CellChat/commands_for_two_extra_tables/make_example_LR_expression_support.sh` | Exact command text from the log of the session that ran it (23 October 2026, 11:33) | A re-run on 20 November 2026 reproduced the table byte for byte |

The submission script of the `cellranger mkfastq` run of November 2025 (SLURM job 18789607) was read from
the SLURM accounting database; its command is quoted in the part 1 notes above. It is not archived as a
file because its header contains a personal e-mail address.

## Not found

These items were searched for on all storage that the account `czhu3869` can read (project and scratch
folders, home folder, archives, editor and job records) and were not found:

- the R script that created `GEX_HTO_Processed.rds` (June 2025) and the Cell Ranger output it read;
- the code that saved `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS`, `sven_liver_new_cytotoxic_cluster_forRCTD.RDS` and the column `annotation_combined`;
- the code of the June 2025 CellChat analyses;
- the package library of the R session of part 3, step 5 (R 4.3.1).
