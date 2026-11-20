# Part 2. CellChat: MNP2 cells and CD4 T cells

Two CellChat runs (23 October 2026). Each infers ligand-receptor communication between the `MNP2` cells
and one group of CD4 T cells of the annotated liver scRNA-seq object (part 1), and compares it with the
`naive CD4` cells. Cell labels are defined in [`../1_scRNAseq_preprocessing_clustering_annotation/README.md`](../1_scRNAseq_preprocessing_clustering_annotation/README.md#cell-labels).

| Run | Script | Input object | Groups (cells) | Output folder |
|---|---|---|---|---|
| 1 | `run_cellchat_MNP2_infected_vs_naive_CD4_latest.R` | `sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS` | MNP2 (167), infected CD4 (761), naive CD4 (754) | `MNP2_CD4_20261023/` |
| 2 | `run_cellchat_MNP2_cytotoxic_vs_naive_CD4_latest.R` | `sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS` | MNP2 (167), cytotoxic CD4 (163), naive CD4 (754) | `MNP2_cytotoxicCD4_vs_naiveCD4_20261023/` |

Each output folder holds the run settings (`00_run_config.txt`), cell counts (`01`, `02`), the
differential ligand-receptor tables CD4 to MNP2 (`03`) and MNP2 to CD4 (`04`), bar charts (`05` to `08`),
the standard CellChat comparison plots (`09` to `13`) and `99_sessionInfo.txt`. SLURM wrappers sit next to
the scripts; their logs are in `slurm_logs/`.

## Run

R library: `cellchat` ([`../environment/README.md`](../environment/README.md)).

1. In each script set `input_rds` (line 19) and `outdir` (line 20).
2. Activate the R 4.5.1 conda environment of the library and point the script to the library:
   `export CELLCHAT_LIB=<library> SEURAT_FALLBACK_LIB=<library>` (read on lines 6 and 7).
3. `Rscript run_cellchat_MNP2_infected_vs_naive_CD4_latest.R` and
   `Rscript run_cellchat_MNP2_cytotoxic_vs_naive_CD4_latest.R` (8 cores: outside SLURM set
   `SLURM_CPUS_PER_TASK=8`, read on line 23; 32 GB; about 12 minutes each).

## Key settings

| Step | Setting |
|---|---|
| Design | per run, two separate CellChat objects with two groups each (MNP2 + infected or cytotoxic CD4; MNP2 + naive CD4); their pairs are joined by name, and a pair absent from one object gets probability 0 (lines 79 to 165, 225 to 226) |
| Database | `CellChatDB.mouse`, categories Secreted Signaling and Cell-Cell Contact (1,942 interactions) |
| Over-expressed genes | `identifyOverExpressedGenes(do.fast = FALSE)` |
| Communication probability | `computeCommunProb(type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = TRUE, nboot = 100, seed.use = 1)` |
| Filter | `filterCommunication(min.cells = 10)`; pairs with p < 0.05 are reported |
| Bar charts | `log(prob + 1)` of the group with the higher probability, positive for infected or cytotoxic CD4, negative for naive CD4 |
| Software | R 4.5.1, CellChat 2.2.0.9001 (GitHub jinworks/CellChat, commit 75253cd), Seurat 5.4.0 |

## Other folders

- `commands_for_two_extra_tables/`: the commands that wrote the two extra tables in `MNP2_CD4_20261023/`
  (`MNP2_vs_original_combinedMNP_comparison.csv`, `example_LR_expression_support.csv`).
- `earlier_analysis_2025-11/`: code and one result table of the November 2025 CellChat analysis of the same
  cells, in which the MNP clusters form one group.

Details, including why SLURM marks both finished jobs as FAILED: [`../docs/PROVENANCE.md`](../docs/PROVENANCE.md#part-2-cellchat).
