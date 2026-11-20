# 2. CellChat: ligand-receptor communication between MNP2 cells and CD4 T cells

This folder holds the two CellChat runs of 23 October 2026, with their scripts, SLURM wrappers, logs, recorded
parameters, result tables and every figure the scripts drew. Each run infers ligand-receptor communication
between the cells labelled `MNP2` and one group of CD4 T cells in the annotated liver scRNA-seq object, and
compares that group with the cells labelled `naive CD4`. The labels are cluster names of the scRNA-seq
analysis (see the table "Objects and labels" in
[`../1_scRNAseq_preprocessing_clustering_annotation/README.md`](../1_scRNAseq_preprocessing_clustering_annotation/README.md#objects-and-labels)).
MNPs are mononuclear phagocytes; MNP2 is one of their clusters.

The files of the two runs are byte-for-byte copies of the files in
`/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/` on the Spartan cluster
(University of Melbourne). Two further folders hold the code and one result table of an earlier CellChat
analysis of the same cells (`earlier_analysis_2025-11/`) and the commands that wrote two extra tables
(`commands_for_two_extra_tables/`); see the sections at the end. [`../MANIFEST.tsv`](../MANIFEST.tsv) lists
the original path and md5 of each file.

## The two runs

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

## What each script does

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

## Parameters

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

## Software

- R 4.5.1, run with the `Rscript` of the conda environment `/data/gpfs/projects/punim2394/conda/rctd_env`
  (`99_sessionInfo.txt`).
- CellChat 2.2.0.9001 (`99_sessionInfo.txt`, `00_run_config.txt`), library
  `/data/gpfs/projects/punim2394/AAA_Sorry_Its_Chen_Here/Rlib_new`. The installed package records its
  source as GitHub `jinworks/CellChat`, commit `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`.
- Seurat 5.4.0, SeuratObject 5.3.0, ggplot2 4.0.2, patchwork 1.3.2, future 1.70.0, ComplexHeatmap 2.26.1,
  dplyr 1.2.0 (`99_sessionInfo.txt`).

## How to reproduce

1. Get the two input Seurat objects (not in this repository, see the main [README](../README.md#data-that-is-not-in-this-repository)).
2. In each R script edit two lines: `input_rds` (line 19) and `outdir` (line 20).
3. Install R 4.5 with Seurat 5 and CellChat at the commit above. Point the script to the libraries with the
   environment variables `CELLCHAT_LIB` and `SEURAT_FALLBACK_LIB`. The default of `CELLCHAT_LIB` written in
   the R script is a folder that does not exist, so the variable has to be set.
4. Run `Rscript <script>.R` on a machine with 8 cores and 32 GB of memory, or adapt the SLURM wrapper and
   submit it with `sbatch`. In the wrapper all paths are Spartan paths: the `Rscript` program, the script,
   `CELLCHAT_LIB`, `SEURAT_FALLBACK_LIB`, `LD_LIBRARY_PATH`, `TMPDIR` and the `--output` and `--error`
   files. One run takes about 12 minutes.

## Notes on the archived runs

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

## Commands that wrote the two extra tables

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

## Earlier CellChat analysis (November 2025)

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
