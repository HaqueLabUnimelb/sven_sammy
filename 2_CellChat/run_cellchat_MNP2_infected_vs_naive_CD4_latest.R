#!/usr/bin/env Rscript

options(width = 180, stringsAsFactors = FALSE, future.globals.maxSize = 8 * 1024^3)
set.seed(1)

project_lib <- Sys.getenv("CELLCHAT_LIB", "/data/scratch/projects/punim1236/AAA_Chen_is_here/Rlib_cellchat_latest")
base_lib <- Sys.getenv("SEURAT_FALLBACK_LIB", "/data/scratch/projects/punim1236/python_envs/rseurat/lib/R/library")
.libPaths(c(project_lib, base_lib))

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(future)
})

input_rds <- "/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS"
outdir <- "/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_latest_MNP2_CD4_20261023_matched_original_params"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
future::plan("multisession", workers = n_workers)

log_line <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "", flush = TRUE)
}

write_text <- function(lines, filename) {
  writeLines(lines, con = file.path(outdir, filename), useBytes = TRUE)
}

save_gg <- function(plot, filename_stub, width, height) {
  ggsave(file.path(outdir, paste0(filename_stub, ".pdf")), plot = plot, width = width, height = height, units = "in", limitsize = FALSE)
  ggsave(file.path(outdir, paste0(filename_stub, ".png")), plot = plot, width = width, height = height, units = "in", dpi = 300, limitsize = FALSE, bg = "white")
}

log_line("STEP 1: load input Seurat object")
seu <- readRDS(input_rds)
stopifnot(inherits(seu, "Seurat"))
DefaultAssay(seu) <- "RNA"
seu$celltype_clean <- trimws(as.character(seu$annotation))

required_groups <- c("MNP2", "infected CD4", "naive CD4")
missing_groups <- setdiff(required_groups, unique(seu$celltype_clean))
if (length(missing_groups) > 0) {
  stop("Missing required annotation groups: ", paste(missing_groups, collapse = ", "))
}

log_line("Input dimensions: ", paste(dim(seu), collapse = " x "))
log_line("CellChat version: ", as.character(packageVersion("CellChat")))
log_line("CellChat library: ", find.package("CellChat"))
log_line("Seurat version: ", as.character(packageVersion("Seurat")))

annotation_counts <- as.data.frame(table(annotation_clean = seu$celltype_clean, useNA = "ifany"))
write.csv(annotation_counts, file.path(outdir, "01_input_annotation_counts.csv"), row.names = FALSE)
annotation_hto <- as.data.frame(table(annotation_clean = seu$celltype_clean, HTO_classification = seu$HTO_classification, useNA = "ifany"))
write.csv(annotation_hto, file.path(outdir, "02_input_annotation_by_HTO.csv"), row.names = FALSE)

write_text(c(
  paste0("input_rds=", input_rds),
  paste0("output_dir=", outdir),
  paste0("R_version=", R.version.string),
  paste0("CellChat_version=", packageVersion("CellChat")),
  paste0("CellChat_path=", find.package("CellChat")),
  paste0("Seurat_version=", packageVersion("Seurat")),
  paste0("n_workers=", n_workers),
  "cell_groups=MNP2; infected CD4; naive CD4",
  "database=CellChatDB.mouse subset to Secreted Signaling and Cell-Cell Contact, matching the existing 2025-11-19 CellChat object",
  "identifyOverExpressedGenes=do.fast FALSE to match the historical scripts/templates",
  "computeCommunProb=type truncatedMean, trim 0.1, raw.use TRUE, population.size TRUE, nboot 100, seed.use 1",
  "comparison=signed log1p probability; positive bars use infected-CD4 probability when higher, negative bars use -naive-CD4 probability when higher",
  "main_direction_matching_screenshot=infected/naive CD4 -> MNP2",
  "literal_reverse_direction=MNP2 -> infected/naive CD4",
  "MNP2_selection=all cells annotated as MNP2 in the input RDS"
), "00_run_config.txt")

run_cellchat <- function(seurat_obj, cd4_group, dataset_name) {
  groups <- c("MNP2", cd4_group)
  cells_use <- colnames(seurat_obj)[seurat_obj$celltype_clean %in% groups]
  log_line("Create ", dataset_name, " object: ", length(cells_use), " cells (", paste(groups, collapse = ", "), ")")
  obj <- subset(seurat_obj, cells = cells_use)
  obj$cellchat_group <- factor(obj$celltype_clean, levels = groups)
  obj <- subset(obj, cells = colnames(obj)[!is.na(obj$cellchat_group)])
  obj$cellchat_group <- droplevels(obj$cellchat_group)
  Idents(obj) <- "cellchat_group"

  cc <- createCellChat(object = obj, group.by = "cellchat_group", assay = "RNA")
  cc@DB <- subsetDB(CellChatDB.mouse, search = c("Secreted Signaling", "Cell-Cell Contact"))
  cc <- subsetData(cc)

  log_line(dataset_name, ": identify over-expressed genes")
  cc <- identifyOverExpressedGenes(cc, do.fast = FALSE)
  log_line(dataset_name, ": identify over-expressed interactions")
  cc <- identifyOverExpressedInteractions(cc)
  log_line(dataset_name, ": compute communication probabilities")
  cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1, raw.use = TRUE, population.size = TRUE, nboot = 100, seed.use = 1L)
  log_line(dataset_name, ": filter and aggregate network")
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- rankNetPairwise(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")

  saveRDS(cc, file.path(outdir, paste0("cellchat_", gsub("[^A-Za-z0-9]+", "_", dataset_name), ".rds")), compress = "xz")
  log_line(dataset_name, ": saved CellChat object")
  cc
}

extract_pair <- function(cc, source_group, target_group, dataset_name) {
  df <- subsetCommunication(cc, sources.use = source_group, targets.use = target_group, thresh = 0.05)
  if (nrow(df) == 0) {
    warning("No significant interactions for ", dataset_name, ": ", source_group, " -> ", target_group)
    df <- data.frame(
      source = character(), target = character(), interaction_name = character(), pathway_name = character(),
      ligand = character(), receptor = character(), annotation = character(), interaction_name_2 = character(),
      prob = numeric(), pval = numeric(), stringsAsFactors = FALSE
    )
  }
  df$dataset <- dataset_name
  df
}

safe_first <- function(x, y) {
  ifelse(!is.na(x) & x != "", x, y)
}

make_diff_table <- function(df_infected, df_naive, direction_label) {
  keep_cols <- c("interaction_name", "interaction_name_2", "pathway_name", "ligand", "receptor", "annotation", "prob", "pval")
  a <- df_infected[, intersect(keep_cols, colnames(df_infected)), drop = FALSE]
  b <- df_naive[, intersect(keep_cols, colnames(df_naive)), drop = FALSE]
  colnames(a)[colnames(a) != "interaction_name"] <- paste0(colnames(a)[colnames(a) != "interaction_name"], "_infected")
  colnames(b)[colnames(b) != "interaction_name"] <- paste0(colnames(b)[colnames(b) != "interaction_name"], "_naive")
  out <- merge(a, b, by = "interaction_name", all = TRUE, sort = FALSE)

  if (!"interaction_name_2_infected" %in% colnames(out)) out$interaction_name_2_infected <- NA_character_
  if (!"interaction_name_2_naive" %in% colnames(out)) out$interaction_name_2_naive <- NA_character_
  if (!"pathway_name_infected" %in% colnames(out)) out$pathway_name_infected <- NA_character_
  if (!"pathway_name_naive" %in% colnames(out)) out$pathway_name_naive <- NA_character_
  if (!"annotation_infected" %in% colnames(out)) out$annotation_infected <- NA_character_
  if (!"annotation_naive" %in% colnames(out)) out$annotation_naive <- NA_character_

  out$prob_infected[is.na(out$prob_infected)] <- 0
  out$prob_naive[is.na(out$prob_naive)] <- 0
  out$pval_infected[is.na(out$pval_infected)] <- 1
  out$pval_naive[is.na(out$pval_naive)] <- 1
  out$interaction_label <- safe_first(out$interaction_name_2_infected, out$interaction_name_2_naive)
  out$interaction_label[is.na(out$interaction_label) | out$interaction_label == ""] <- out$interaction_name[is.na(out$interaction_label) | out$interaction_label == ""]
  out$interaction_label <- gsub(" - ", " -> ", out$interaction_label, fixed = TRUE)
  out$interaction_label <- toupper(out$interaction_label)
  out$pathway_name <- safe_first(out$pathway_name_infected, out$pathway_name_naive)
  out$annotation <- safe_first(out$annotation_infected, out$annotation_naive)
  out$detected_in <- ifelse(out$prob_infected > 0 & out$prob_naive > 0, "both", ifelse(out$prob_infected > 0, "infected_CD4_only", "naive_CD4_only"))
  out$log_prob_infected <- log1p(out$prob_infected)
  out$log_prob_naive <- log1p(out$prob_naive)
  out$diff_log_prob <- out$log_prob_infected - out$log_prob_naive
  out$higher_condition <- ifelse(out$prob_infected >= out$prob_naive, "Infected CD4", "Naive CD4")
  out$signed_log_prob <- ifelse(out$higher_condition == "Infected CD4", out$log_prob_infected, -out$log_prob_naive)
  out$direction <- direction_label
  out <- out[out$signed_log_prob != 0, , drop = FALSE]
  out <- out[order(-out$signed_log_prob, out$interaction_label), , drop = FALSE]
  rownames(out) <- NULL
  out
}

plot_diff <- function(df, title_text, prefix, max_each_side = Inf) {
  plot_df <- df[!is.na(df$signed_log_prob) & df$signed_log_prob != 0, , drop = FALSE]
  if (is.finite(max_each_side)) {
    pos <- plot_df[plot_df$signed_log_prob > 0, , drop = FALSE]
    neg <- plot_df[plot_df$signed_log_prob < 0, , drop = FALSE]
    pos <- head(pos, max_each_side)
    neg <- neg[order(neg$signed_log_prob, neg$interaction_label), , drop = FALSE]
    neg <- head(neg, max_each_side)
    plot_df <- rbind(pos, neg)
    plot_df <- plot_df[order(-plot_df$signed_log_prob, plot_df$interaction_label), , drop = FALSE]
  }
  if (nrow(plot_df) == 0) {
    warning("No non-zero differential interactions for ", prefix)
    return(invisible(NULL))
  }

  plot_df$interaction_label <- factor(plot_df$interaction_label, levels = plot_df$interaction_label)
  n_pos <- sum(plot_df$signed_log_prob > 0)
  n_neg <- sum(plot_df$signed_log_prob < 0)
  y_max <- max(plot_df$signed_log_prob, 0)
  y_min <- min(plot_df$signed_log_prob, 0)
  y_pad <- max(0.0002, 0.08 * (y_max - y_min))

  p <- ggplot(plot_df, aes(x = interaction_label, y = signed_log_prob, fill = higher_condition)) +
    geom_col(width = 0.82, colour = NA) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.45) +
    scale_fill_manual(values = c("Infected CD4" = "black", "Naive CD4" = "grey72"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
    labs(title = title_text, x = NULL, y = "log(prob + 1)") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "plain"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7, colour = "black"),
      axis.text.y = element_text(size = 11, colour = "black"),
      axis.title.y = element_text(size = 14),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor = element_line(colour = "grey94", linewidth = 0.25),
      plot.margin = margin(12, 18, 12, 12)
    )

  if (n_pos > 0 && n_neg > 0) {
    p <- p + geom_vline(xintercept = n_pos + 0.5, colour = "black", linewidth = 0.7)
  }
  if (n_pos > 0) {
    p <- p + annotate("text", x = max(1, n_pos * 0.58), y = y_max + y_pad * 0.55, label = "Infected CD4+", size = 5.2, colour = "black", hjust = 0.5)
  }
  if (n_neg > 0) {
    p <- p + annotate("text", x = n_pos + max(1, n_neg * 0.50), y = y_min - y_pad * 0.40, label = "Naive CD4+", size = 5.2, colour = "grey50", hjust = 0.5)
  }

  width <- max(9, min(26, 4.5 + 0.14 * nrow(plot_df)))
  height <- 7.5
  save_gg(p, prefix, width, height)
  invisible(p)
}

log_line("STEP 2: run CellChat separately for infected-CD4 and naive-CD4 comparisons")
cc_infected <- run_cellchat(seu, "infected CD4", "MNP2_with_infected_CD4")
cc_naive <- run_cellchat(seu, "naive CD4", "MNP2_with_naive_CD4")

log_line("STEP 3: extract source-target ligand-receptor communications")
cd4_to_mnp2 <- make_diff_table(
  extract_pair(cc_infected, "infected CD4", "MNP2", "MNP2_with_infected_CD4"),
  extract_pair(cc_naive, "naive CD4", "MNP2", "MNP2_with_naive_CD4"),
  "CD4 -> MNP2"
)
mnp2_to_cd4 <- make_diff_table(
  extract_pair(cc_infected, "MNP2", "infected CD4", "MNP2_with_infected_CD4"),
  extract_pair(cc_naive, "MNP2", "naive CD4", "MNP2_with_naive_CD4"),
  "MNP2 -> CD4"
)

write.csv(cd4_to_mnp2, file.path(outdir, "03_differential_LR_CD4_to_MNP2_infected_minus_naive.csv"), row.names = FALSE)
write.csv(mnp2_to_cd4, file.path(outdir, "04_differential_LR_MNP2_to_CD4_infected_minus_naive.csv"), row.names = FALSE)

log_line("CD4 -> MNP2 differential LR pairs: ", nrow(cd4_to_mnp2))
log_line("MNP2 -> CD4 differential LR pairs: ", nrow(mnp2_to_cd4))

log_line("STEP 4: plot screenshot-style differential LR bars")
plot_diff(cd4_to_mnp2, "Infected vs Naive CD4  →  MNP2", "05_FigS2E_style_CD4_to_MNP2_all")
plot_diff(cd4_to_mnp2, "Infected vs Naive CD4  →  MNP2 (top 25 per side)", "06_FigS2E_style_CD4_to_MNP2_top25", max_each_side = 25)
plot_diff(mnp2_to_cd4, "MNP2  →  Infected vs Naive CD4", "07_reverse_direction_MNP2_to_CD4_all")
plot_diff(mnp2_to_cd4, "MNP2  →  Infected vs Naive CD4 (top 25 per side)", "08_reverse_direction_MNP2_to_CD4_top25", max_each_side = 25)

log_line("STEP 5: create lifted merged object and standard CellChat comparison plots")
common_groups <- c("MNP2", "infected CD4", "naive CD4")
cc_naive_lifted <- liftCellChat(cc_naive, common_groups)
cc_infected_lifted <- liftCellChat(cc_infected, common_groups)
cc_merged <- mergeCellChat(
  list(`Naive CD4` = cc_naive_lifted, `Infected CD4` = cc_infected_lifted),
  add.names = c("Naive CD4", "Infected CD4"),
  merge.data = FALSE,
  cell.prefix = TRUE
)
saveRDS(cc_merged, file.path(outdir, "cellchat_merged_lifted_MNP2_CD4.rds"), compress = "xz")

gg_count <- compareInteractions(cc_merged, show.legend = FALSE, measure = "count", title.name = "Total inferred interactions")
gg_weight <- compareInteractions(cc_merged, show.legend = FALSE, measure = "weight", title.name = "Total interaction strength")
save_gg(gg_count + gg_weight, "09_total_interactions_count_and_weight", 8, 4)

pdf(file.path(outdir, "10_netVisual_diffInteraction_weight.pdf"), width = 7, height = 7)
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "weight", title.name = "Differential interaction strength (red: higher in Infected CD4)")
dev.off()
png(file.path(outdir, "10_netVisual_diffInteraction_weight.png"), width = 2100, height = 2100, res = 300, bg = "white")
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "weight", title.name = "Differential interaction strength (red: higher in Infected CD4)")
dev.off()

pdf(file.path(outdir, "11_netVisual_diffInteraction_count.pdf"), width = 7, height = 7)
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "count", title.name = "Differential interaction number (red: higher in Infected CD4)")
dev.off()
png(file.path(outdir, "11_netVisual_diffInteraction_count.png"), width = 2100, height = 2100, res = 300, bg = "white")
netVisual_diffInteraction(cc_merged, weight.scale = TRUE, measure = "count", title.name = "Differential interaction number (red: higher in Infected CD4)")
dev.off()

ht_weight <- netVisual_heatmap(cc_merged, measure = "weight", title.name = "Differential interaction strength (red: higher in Infected CD4)")
pdf(file.path(outdir, "12_netVisual_heatmap_weight.pdf"), width = 7, height = 5)
ComplexHeatmap::draw(ht_weight)
dev.off()
png(file.path(outdir, "12_netVisual_heatmap_weight.png"), width = 2100, height = 1500, res = 300, bg = "white")
ComplexHeatmap::draw(ht_weight)
dev.off()

ht_count <- netVisual_heatmap(cc_merged, measure = "count", title.name = "Differential interaction number (red: higher in Infected CD4)")
pdf(file.path(outdir, "13_netVisual_heatmap_count.pdf"), width = 7, height = 5)
ComplexHeatmap::draw(ht_count)
dev.off()
png(file.path(outdir, "13_netVisual_heatmap_count.png"), width = 2100, height = 1500, res = 300, bg = "white")
ComplexHeatmap::draw(ht_count)
dev.off()

write_text(capture.output(sessionInfo()), "99_sessionInfo.txt")
saveRDS(sessionInfo(), file.path(outdir, "99_sessionInfo.rds"))
log_line("DONE")
