# Part 1. scRNA-seq: preprocessing, clustering and cell annotation

One 10x Genomics capture with two libraries, gene expression and 8 hashtag oligos (HTO), of CD4 T cells and
mononuclear phagocytes (MNPs) from the liver and spleen of infected and uninfected mice (sequenced on
27 May 2025). These scripts produce the annotated liver cells that part 2 (CellChat) uses and that part 3
uses as the RCTD reference.

## Files

| File | Purpose |
|---|---|
| `1_cellranger/2025-11_mkfastq_and_count/` | Records of a later Cell Ranger 8.0.1 run (November 2025): `samplesheet.csv`, `library.csv`, hashtag reference, command records (`Split_Output_new/_invocation`, `new_cellranger/_cmdline`) and metrics |
| `2_seurat_clustering_annotation/liver_scRNAseq_processed_NI.R` | Main script: removes contaminating clusters, keeps the liver cells, clusters them (resolution 0.7), assigns cell labels, draws UMAPs and the dot plot |
| `2_seurat_clustering_annotation/sven_liverscRNAseq.R` | Full working script of the same analysis. Lines 860 to 872 save `sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS` (CellChat run 2); lines 413 to 430 save `sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS` (read by the November 2025 CellChat code). No archived script saves the object used by CellChat run 1 and RCTD |
| `2_seurat_clustering_annotation/Tcell_combined_ref.R` | Merges the four CD4 labels into `Tcell` (intermediate reference object) |
| `2_seurat_clustering_annotation/gene_signatures/` | Gene lists for the MNP1, MNP2, pDC and DC module scores |

## Run

R library: `rctd_seurat` (R 4.5.0, Seurat 5.5.1; [`../environment/README.md`](../environment/README.md)). The
package versions of the original sessions are not recorded (the objects were saved by R 4.4.0).

1. Cell Ranger 8.0.1 (records of the November 2025 run): run `cellranger mkfastq` with the arguments in
   `Split_Output_new/_invocation`, then the command in `new_cellranger/_cmdline` (edit the paths in it and in
   `library.csv`; for `--feature-ref` use `new_cellranger/outs/feature_reference.csv`). This run is not the
   one whose matrix is in `GEX_HTO_Processed.rds` (9,368 shared barcodes; per-cell UMI totals about 1.5 times
   higher in the object); that run and the code from its output to `GEX_HTO_Processed.rds` were not found.
2. In R, set the data folder (the scripts use `/data/projects/punim1236/liver_scRNAseq/`) and run
   `liver_scRNAseq_processed_NI.R` section by section. It starts from `GEX_HTO_Processed.rds`, the Seurat
   object with all cells (not in this repository). The gene lists of lines 180 to 182 are the files in
   `gene_signatures/`.
3. Annotated objects: `Tcell_combined_ref.R` (it reads `sven_liver_new_cytotoxic_cluster_forRCTD.RDS`, whose
   code was not found), then lines 860 to 872 of `sven_liverscRNAseq.R`, which save the "split up" object of
   CellChat run 2. No archived code saves the "combined" object of CellChat run 1 and RCTD.

The scripts were run interactively and select clusters by number, so they need the saved objects of the
original analysis. Which lines to run in which order: [`../docs/PROVENANCE.md`](../docs/PROVENANCE.md#running-the-interactive-scripts).

## Key settings

| Step | Setting | Where |
|---|---|---|
| Hashtag demultiplexing, normalisation, first clustering | Recorded in `GEX_HTO_Processed.rds` (`@commands`): CLR for hashtags, 2,000 variable genes, 30 PCs, 20 dimensions, resolution 0.3. The hashtag calls are in metadata columns of the kind Seurat `HTODemux` writes; the call itself is not recorded | object `@commands` and metadata |
| Remove contaminating clusters | keep first-clustering clusters 0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14 | line 67 |
| Liver cells | hashtags B1 to B4 (liver) and the clusters of line 67 except cluster 3: 2,934 cells | lines 20 and 65 (line 80 includes cluster 3; skip lines 80 to 95) |
| Clustering of liver cells | `FindClusters(resolution = 0.7)`: 15 clusters | line 101 |
| Module scores | `AddModuleScore` with the gene lists | lines 184 to 209 |
| Cell labels | see the table below | line 260; `sven_liverscRNAseq.R` lines 414, 866 |

## Cell labels

Labels of the resolution-0.7 clusters (cells per cluster). DC and pDC clusters are removed from the objects
used in parts 2 and 3 (2,537 cells). Column 3 was read from the saved object (no archived code assigns these
labels); column 4 is line 866 of `sven_liverscRNAseq.R`.

| Cluster (cells) | Figures (line 260) | `..._combined_MNP1_MNP2_forRCTD.RDS` (CellChat run 1, RCTD) | `..._split_up_MNP1_MNP2_forRCTD.RDS` (CellChat run 2) |
|---|---|---|---|
| 0 (754) | naive CD4 | naive CD4 | naive CD4 |
| 1 (491) | MNP1 | MNP1 | MNP1 |
| 2 (364) | uninfected MNP | uninfected MNP | uninfected MNP |
| 3 (312), 5 (170) | effector CD4 | infected CD4 | effector CD4 |
| 6 (167) | MNP2 | MNP2 | MNP2 |
| 7 (163) | effector CD4 | infected CD4 | cytotoxic CD4 |
| 10 (65), 11 (51) | effector CD4 | infected CD4 | proliferating CD4 |
| 4, 8, 12, 13, 14 | DC | removed | removed |
| 9 | pDC | removed | removed |

The labels name clusters, not mice: the `naive CD4` cells and most `infected CD4` cells carry hashtags of
infected mice (B1, B2).

## Not in this repository

- The script that created `GEX_HTO_Processed.rds` (June 2025) and the Cell Ranger output it read.
- The code that saved `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS`.
- The code that saved `sven_liver_new_cytotoxic_cluster_forRCTD.RDS`, the input of `Tcell_combined_ref.R`.

All three were searched for and not found; details in [`../docs/PROVENANCE.md`](../docs/PROVENANCE.md#part-1-scrna-seq).
