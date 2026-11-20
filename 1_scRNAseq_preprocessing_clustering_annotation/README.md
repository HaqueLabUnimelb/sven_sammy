# 1. Liver scRNA-seq: preprocessing, clustering and cell annotation

One 10x Genomics single-cell capture with two libraries, gene expression and hashtag oligos (HTO), of CD4
T cells and mononuclear phagocytes (MNPs) from mouse liver and spleen. The run was sequenced on
27 May 2025 (NextSeq 2000, flowcell `AAH2VJHM5`, read lengths 28 + 10 + 10 + 90). The liver cells are the
ones used in the CellChat analysis ([`../2_CellChat`](../2_CellChat)) and as the reference for RCTD on the
Stereo-seq section ([`../3_STOmics_liver_v1.3`](../3_STOmics_liver_v1.3)).

All files are byte-for-byte copies. [`../MANIFEST.tsv`](../MANIFEST.tsv) lists the original path, owner,
last-modified time and md5 of each file. Paths inside the scripts point to
`/data/projects/punim1236/liver_scRNAseq/` on the Spartan cluster (University of Melbourne).

Hashtags (from `feature_reference.csv`). The names distinguish infected and uninfected mice; the archived
files do not name the infection.

| Hashtag | Sample |
|---|---|
| HTO1 | B1: Naive CD4s Liver (infected mice) |
| HTO2 | B2: Activated CD4s Liver (infected mice) |
| HTO3 | B3: MNPs Liver (infected mice) |
| HTO4 | B4: Naive MNPs Liver (uninfected mice) |
| HTO7 | B7: Naive CD4s Spleen (infected mice) |
| HTO8 | B8: Activated CD4s Spleen (infected mice) |
| HTO9 | B9: MNPs Spleen (infected mice) |
| HTO10 | B10: Naive MNPs Spleen (uninfected mice) |

## From FASTQ to the analysed object: what is and is not archived

The R scripts of this part start from one Seurat object, `GEX_HTO_Processed.rds` (9,461 cells, processed on
2 June 2025). The script that built this object and the Cell Ranger output it read were not found.
`1_cellranger/` holds the records of a later Cell Ranger run of the same sequencing run (November 2025).
The count matrix of that later run is not the one in the object: 9,368 cell barcodes are shared, and the
per-cell UMI totals in the object are about 1.5 times those of the November matrix. The step from FASTQ to
`GEX_HTO_Processed.rds` can therefore not be reproduced exactly from this archive.

## Files

| Step | File | What it does |
|---|---|---|
| 1. FASTQ and counting, first run (May to June 2025) | not archived here | Two scripts of this run are in the lab repository `HaqueLabUnimelb/Cellranger` (private; commits of 2 June 2025). `bcl_to_fastq.slurm` runs bcl2fastq2 2.20.0: `bcl2fastq --runfolder-dir /data/projects/punim1236/liver_scRNAseq/250527_VH00915_604_AAH2VJHM5-10042-micheleC --output-dir /data/projects/punim1236/liver_scRNAseq/fastq_output --no-lane-splitting`. `run.slurm` runs Cell Ranger 8.0.1 for gene expression only (no hashtag library, no feature reference): `cellranger count --id=Sample1 --transcriptome=/data/projects/punim1236/liver_scRNAseq/mm10/refdata-gex-mm10-2020-A --fastqs=/data/projects/punim1236/liver_scRNAseq/fastq_output --sample=Sample1 --localcores=32 --localmem=128 --create-bam true`. Two attempts of 2 June 2025 are on Spartan: `Sample1_HTO/` (gene expression and hashtags together, failed) and `Sample1_HTOonly/` (hashtag library only, interrupted). The run that produced the matrix behind `GEX_HTO_Processed.rds` was not found. |
| 1. FASTQ and counting, second run (November 2025) | `1_cellranger/2025-11_mkfastq_and_count/samplesheet.csv`, `Split_Output_new/_invocation`, `Split_Output_new/_versions` | `cellranger mkfastq` 8.0.1: sample indices `SI-TT-C11` (gene expression, `GEX_C11`) and `SI-NT-B1` (hashtags, `CSP_B1`). `_invocation` is the command record written by Cell Ranger. The submission script (SLURM job 18789607) is kept in the SLURM accounting database of Spartan: it loads `CellRanger/8.0.1`, `GCC/11.3.0` and `bcl2fastq2/2.20.0` and runs `cellranger mkfastq --id=Split_Output_new --run=/data/projects/punim1236/liver_scRNAseq/250527_VH00915_604_AAH2VJHM5-10042-micheleC --csv=samplesheet.csv --localcores=24 --localmem=128`. It is not archived as a file because its header contains a personal e-mail address. |
| | `1_cellranger/2025-11_mkfastq_and_count/library.csv`, `new_cellranger/_cmdline`, `_invocation`, `_versions` | `cellranger count` 8.0.1 with both libraries (`Gene Expression` and `Antibody Capture`), reference `refdata-gex-mm10-2020-A`, `--create-bam=false`, 24 cores, 128 GB (SLURM job 18810731). `_cmdline` is the full command. The SLURM wrapper `run.slurm` that submitted it is not archived here because its header contains a personal e-mail address. |
| | `1_cellranger/2025-11_mkfastq_and_count/new_cellranger/outs/feature_reference.csv`, `outs/metrics_summary.csv` | The hashtag reference that the count run used, and its summary metrics (9,560 cells). |
| | `../2_CellChat/earlier_analysis_2025-11/rstudio_editor_buffer_2025-11-21.R`, lines 23 to 98 | Seurat processing of the November 2025 count matrix into `new_livescRNA.RDS`: cells kept with `nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10` (lines 57 to 62), `HTODemux(seu, assay = "HTO", positive.quantile = 0.99)` (lines 77 to 81), then PCA and UMAP. This object is not used by the archived analyses; the file is described in part 2. |
| 2. Seurat object, hashtag demultiplexing, first clustering | no script in this archive | Result: `GEX_HTO_Processed.rds`. The processing calls recorded inside the object are listed below. |
| 3. Clean-up, liver subset, clustering, annotation, figures | `2_seurat_clustering_annotation/liver_scRNAseq_processed_NI.R` | The script of the "NI revision" (line 254: `# Re-annotate scRNAseq dataset for NI revision and perform dotplot`). Reads `GEX_HTO_Processed.rds`. Removes contaminating clusters (saves `sven_refined.RDS`), keeps the four liver hashtags, clusters again with `FindClusters(resolution = 0.7)` (saves `sven_liver_for_paper.RDS`), draws the gene UMAPs, adds the MNP1, MNP2 and pDC signature scores, assigns the cell-type labels of the revision (line 260) and draws the annotated UMAPs and the dot plot (`Dotplot_NI_revision_v3.pdf`). |
| | `2_seurat_clustering_annotation/gene_signatures/MNP1.txt`, `MNP2.txt`, `pDCs.txt`, `DCs.txt` | Gene lists for `AddModuleScore`. The revision script reads the first three (lines 180 to 182), the working script reads all four. On Spartan they are in `/data/projects/punim1236/liver_scRNAseq/Ollie/`. |
| 4. Annotated objects for RCTD and CellChat | `2_seurat_clustering_annotation/sven_liverscRNAseq.R` | The complete working script (941 lines) from which the revision script was taken. Besides `Tcell_combined_ref.R`, it is the only archived file with code for the annotated reference objects. Line 414 labels the resolution-0.7 clusters; line 417 assigns the labels again with the MNP clusters merged into one label `MNP`; lines 427 to 430 remove the DC and pDC clusters and save `sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS`. Lines 860 to 872 read `sven_liver_Tcells_combined_forRCTD.RDS`, assign the labels `MNP1`, `MNP2`, `uninfected MNP` and the CD4 subsets, and save `sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS`. |
| | `2_seurat_clustering_annotation/Tcell_combined_ref.R` | Reads `sven_liver_new_cytotoxic_cluster_forRCTD.RDS`, merges the four CD4 labels into `Tcell` and saves `sven_liver_Tcells_combined_forRCTD.RDS`. |

## Processing recorded in `GEX_HTO_Processed.rds`

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

The quality-control thresholds and the `HTODemux` quantile are not recorded in the object.

## Objects and labels

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

## How to reproduce

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

## Not in this archive

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
- Seurat objects, count matrices and FASTQ files (see the main [README](../README.md#data-that-is-not-in-this-repository)).
- R and Seurat versions of steps 2 to 4: not recorded in the scripts. The saved objects carry the
  SeuratObject version stamp 5.1.0.
