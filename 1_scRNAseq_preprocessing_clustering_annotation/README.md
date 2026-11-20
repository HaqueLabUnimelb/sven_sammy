# Part 1. scRNA-seq: preprocessing, clustering and cell annotation

One 10x Genomics capture with two libraries, gene expression and 8 hashtag oligos (HTO), of CD4 T cells and
mononuclear phagocytes (MNPs) from the liver and spleen of infected and uninfected mice (sequenced on
27 May 2025). These scripts produce the annotated liver cells that part 2 (CellChat) uses and that part 3
uses as the RCTD reference.

## Files

| File | Purpose |
|---|---|
| `1_cellranger/run_cellranger.sh` | Standard Cell Ranger 8.0.1 workflow (`mkfastq`, then `count` of the gene expression and hashtag libraries), written for this repository with the recorded commands |
| `1_cellranger/2025-11_mkfastq_and_count/` | Records of the Cell Ranger run of November 2025: `samplesheet.csv`, `library.csv`, hashtag reference, command records (`Split_Output_new/_invocation`, `new_cellranger/_cmdline`) and metrics |
| `2_seurat_clustering_annotation/0_make_GEX_HTO_Processed.R` | Standard Seurat processing into `GEX_HTO_Processed.rds` (hashtag demultiplexing, cell filter, normalisation, first clustering), written for this repository with the settings given in the Methods of the manuscript |
| `2_seurat_clustering_annotation/liver_scRNAseq_processed_NI.R` | Main script: removes contaminating clusters, keeps the liver cells, clusters them (resolution 0.7), assigns cell labels, draws UMAPs and the dot plot |
| `2_seurat_clustering_annotation/sven_liverscRNAseq.R` | Full working script of the same analysis. Lines 860 to 872 save `sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS` (CellChat run 2); lines 413 to 430 save `sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS` (read by the November 2025 CellChat code) |
| `2_seurat_clustering_annotation/Tcell_combined_ref.R` | Merges the four CD4 labels into `Tcell` (intermediate reference object) |
| `2_seurat_clustering_annotation/gene_signatures/` | Gene lists for the MNP1, MNP2, pDC and DC module scores |

## Run

R library: `rctd_seurat` (R 4.5.0, Seurat 5.5.1; [`../environment/README.md`](../environment/README.md)).

1. Cell Ranger 8.0.1: edit the two paths at the top of `1_cellranger/run_cellranger.sh` and run it in an
   empty folder.
2. `Rscript 0_make_GEX_HTO_Processed.R new_cellranger/outs/filtered_feature_bc_matrix GEX_HTO_Processed.rds`
3. In R, set the data folder (the scripts use `/data/projects/punim1236/liver_scRNAseq/`) and run
   `liver_scRNAseq_processed_NI.R` section by section, starting from `GEX_HTO_Processed.rds`. The gene lists of
   lines 180 to 182 are the files in `gene_signatures/`.
4. Annotated objects: `Tcell_combined_ref.R`, then lines 860 to 872 of `sven_liverscRNAseq.R`.

The scripts were run interactively and select clusters by number, so they need the saved objects of the
original analysis (listed in [`../DATA_NOT_IN_REPOSITORY.tsv`](../DATA_NOT_IN_REPOSITORY.tsv)). Which lines to
run in which order: [`../docs/PROVENANCE.md`](../docs/PROVENANCE.md#running-the-interactive-scripts).

## Key settings

| Step | Setting | Where |
|---|---|---|
| Cell Ranger | `cellranger count` with both libraries, reference `refdata-gex-mm10-2020-A` (mouse GRCm38, GENCODE vM23), `--create-bam=false` | `1_cellranger/run_cellranger.sh` |
| Hashtag demultiplexing, cell filter, normalisation, first clustering | CLR for hashtags, `HTODemux(positive.quantile = 0.95)`; cells with more than 20 % mitochondrial counts, more than 60,000 UMIs or more than 7,500 detected genes removed; 2,000 variable genes, 30 PCs, 20 dimensions, resolution 0.3 | `0_make_GEX_HTO_Processed.R` (settings given in the Methods of the manuscript) |
| Remove contaminating clusters | keep first-clustering clusters 0, 1, 2, 3, 4, 6, 7, 8, 9, 12, 13, 14 | line 67 |
| Liver cells | hashtags B1 to B4 (liver) and the clusters of line 67 except cluster 3: 2,934 cells | lines 20 and 65 (line 80 includes cluster 3; skip lines 80 to 95) |
| Clustering of liver cells | `FindClusters(resolution = 0.7)`: 15 clusters | line 101 |
| Module scores | `AddModuleScore` with the gene lists | lines 184 to 209 |
| Cell labels | see the table below | line 260; `sven_liverscRNAseq.R` lines 414, 866 |

## Cell labels

Labels of the resolution-0.7 clusters (cells per cluster). DC and pDC clusters are removed from the objects
used in parts 2 and 3 (2,537 cells). Column 3 is read from the saved object; column 4 is line 866 of
`sven_liverscRNAseq.R`.

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
