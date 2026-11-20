library(dplyr)
library(ggplot2)
library(CellChat)
library(patchwork)
library(Seurat)
devtools::install_github("jinworks/CellChat")



install.packages('dplyr')
#一定要跑
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jokergoo/circlize")
devtools::install_github('immunogenomics/presto')


library(Seurat)
library(dplyr)
library(Matrix)



data_dir <- "/data/projects/punim1236/liver_scRNAseq/new_cellranger/outs/filtered_feature_bc_matrix"

counts <- Read10X(data.dir = data_dir)


rna_counts <- counts$`Gene Expression`
hto_counts <- counts$`Antibody Capture`
common_cells <- intersect(colnames(rna_counts), colnames(hto_counts))

length(common_cells)     # 看看交集多少
rna_counts <- rna_counts[, common_cells]
hto_counts <- hto_counts[, common_cells]

seu <- CreateSeuratObject(
  counts = rna_counts,
  project = "Liver_scRNA_hashing")

# 正常加入 HTO assay
seu@assays[["HTO"]] <- CreateAssayObject(counts = hto_counts)

DefaultAssay(seu) <- "RNA"

seu


# 计算线粒体比例
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")


seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt    < 10
)



# 切换到 HTO assay，做 CLR 标准化
DefaultAssay(seu) <- "HTO"

seu <- NormalizeData(seu)


# 基于 HTO 的最近邻图（用于 demux）
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)

# 运行 HTODemux（Seurat 自带，取阳性阈值可调）
seu <- HTODemux(
  seu,
  assay            = "HTO",
  positive.quantile = 0.99
)


table(seu$hash.ID)

# 看一下每个 hash.ID 在 UMAP 上的分布（先用 RNA 跑个粗略 UMAP）
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, reduction = "pca")

DimPlot(seu, reduction = "umap", group.by = "hash.ID")

FeaturePlot(seu, features = "Zeb2" )

saveRDS(seu,file="/data/projects/punim1236/liver_scRNAseq/new_cellranger/new_livescRNA.RDS")


### CELL CHAT ### CELL CHAT ### CELL CHAT ### CELL CHAT ### CELL CHAT ### CELL CHAT




cc_seurat <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS")

DimPlot(cc_seurat, group = "annotation")

data.input <- cc_seurat[["RNA"]]$data

Idents(cc_seurat) <- cc_seurat$hash.ID

labels <- cc_seurat@meta.data[["annotation"]]

meta <- data.frame(labels = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

CellChatDB <- CellChatDB.mouse  
showDatabaseCategory(CellChatDB)

CellChatDB.use <- subsetDB(CellChatDB, 
                           search = c("Secreted Signaling", "Cell-Cell Contact"), 
                           key = "annotation")

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 16) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

#cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, population.size = TRUE)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


saveRDS(cellchat, file = "/data/projects/punim1236/liver_scRNAseq/2025_11_19_newCellchat_combined_MNP.rds")

cellchat <- readRDS("/data/projects/punim1236/liver_scRNAseq/cellchat_humanSkin_LS.rds")

pathways.show <- c("TNF") 

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

cellchat@netP$pathways
df.net <- subsetCommunication(cellchat, slot.name = "net")
write.csv(df.net, file = "/data/projects/punim1236/liver_scRNAseq/AAA_complete_altogether.csv", row.names = FALSE)



df.net <- df.net %>%
  mutate(`log(prob+1)` = log(prob + 1),
         highlight = ifelse(pathway_name %in% highlight_pathways, pathway_name, "Other"))

exp <- df.net %>%
  filter(source == "MNP", target == "cytotoxic CD4") %>%
  mutate(
    log_prob = log(prob + 1),
    group = "Cytotoxic CD4",
    highlight = ifelse(pathway_name %in% highlight_pathways, pathway_name, "Other")
  ) %>%
  select(interaction_name, log_prob, group, highlight)

ctrl <- df.net %>%
  filter(source == "MNP", target == " naive CD4") %>%
  mutate(
    log_prob = log(prob + 1),
    group = "Naive CD4",
    highlight = ifelse(pathway_name %in% highlight_pathways, pathway_name, "Other")
  ) %>%
  select(interaction_name, log_prob, group, highlight)


all_interactions <- union(exp$interaction_name, ctrl$interaction_name)

template <- expand.grid(
  interaction_name = all_interactions,
  group = c("Cytotoxic CD4", "Naive CD4"),
  stringsAsFactors = FALSE
)

library(dplyr)
library(ggplot2)

df_long <- template %>%
  left_join(bind_rows(exp, ctrl), by = c("interaction_name", "group")) %>%
  mutate(
    log_prob = ifelse(is.na(log_prob), 0, log_prob),
    highlight = ifelse(is.na(highlight), "Other", highlight),
    color_group = ifelse(group == "Cytotoxic CD4", highlight, "Other"),
    interaction_name = factor(interaction_name, levels = unique(interaction_name[order(-log_prob)]))
  )


order_levels <- df_long %>%
  filter(group == "Cytotoxic CD4") %>%
  arrange(desc(log_prob)) %>%
  pull(interaction_name) %>%
  unique()

df_long <- df_long %>%
  mutate(
    interaction_name = factor(interaction_name, levels = rev(order_levels))
  )

library(ggplot2)

ggplot(df_long, aes(y = interaction_name, x = log_prob, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c("Cytotoxic CD4" = "black", "Naive CD4" = "gray80"),
    name = "Group"
  ) +
  labs(
    x = "log(prob + 1)",
    y = "Interaction",
    title = "MNP ➝ Cytotoxic CD4 vs Naive CD4"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )





cellchat@netP$pathways
# 获取通信数据
df.net <- subsetCommunication(cellchat, slot.name = "net")

# 设置实验组和对照组的 log(prob+1)
exp <- df.net %>%
  filter(source == "cytotoxic CD4", target == "MNP") %>%
  mutate(log_prob_exp = log(prob + 1)) %>%
  select(interaction_name, log_prob_exp)

ctrl <- df.net %>%
  filter(source == " naive CD4", target == "MNP") %>%
  mutate(log_prob_ctrl = log(prob + 1)) %>%
  select(interaction_name, log_prob_ctrl)

# 合并并计算差值
df_wide <- full_join(exp, ctrl, by = "interaction_name") %>%
  mutate(
    log_prob_exp = ifelse(is.na(log_prob_exp), 0, log_prob_exp),
    log_prob_ctrl = ifelse(is.na(log_prob_ctrl), 0, log_prob_ctrl),
    log_prob_diff = log_prob_exp - log_prob_ctrl
  )

# 按照差值排序绘图顺序
df_wide <- df_wide %>%
  arrange(log_prob_diff) %>%
  mutate(interaction_name = factor(interaction_name, levels = interaction_name))


ggplot(df_wide, aes(x = log_prob_diff, y = interaction_name, fill = log_prob_diff > 0)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(
    values = c("TRUE" = "red", "FALSE" = "blue"),
    labels = c("Naive CD4 ➝ MNP", "Cytotoxic CD4 ➝ MNP"),
    name = "Group"
  ) +
  labs(
    x = "log(prob + 1) (Cytotoxic - Naive)",
    y = "Interaction",
    title = "Cytotoxic CD4 vs Naive CD4 ➝ MNP"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top"
  )



df_wide$interaction_name <- factor(df_wide$interaction_name, levels = rev(levels(df_wide$interaction_name)))



p <- ggplot(df_wide, aes(x = log_prob_diff, y = interaction_name, fill = log_prob_diff > 0)) +
  geom_bar(stat = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "gray80"),
    labels = c("Naive CD4 ➝ MNP", "Cytotoxic CD4 ➝ MNP"),
    name = "Group"
  ) +
  labs(
    x = "log(prob + 1) (Cytotoxic - Naive)",
    y = "Interaction",
    title = "Cytotoxic CD4 vs Naive CD4 ➝ MNP"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "top",
    # 关键：适当加大边距（单位可以用 "pt" 或 "cm"）
    plot.margin = margin(t = 9, r = 30, b = 20, l = 18, unit = "pt")
  ) +
  coord_flip()

ggsave(
  "interaction_bar.pdf",
  plot   = p,
  width  = 16,   # 单位 inch，适当加大
  height = 10,
  units  = "in",
  dpi    = 1200
)









netVisual_bubble(
  object = cellchat,
  sources.use = c("B2: Activated CD4s Liver (infected mice)", "B3: MNPs Liver (infected mice)"),
  targets.use = c("B2: Activated CD4s Liver (infected mice)", "B3: MNPs Liver (infected mice)"),
  remove.isolate = FALSE,
  color.heatmap = 'viridis'
)



pathways.show <- c("FASLG", "IFN-II")

netVisual_bubble(cellchat, signaling = pathways.show, sources.use = "B2: Activated CD4s Liver (infected mice)", targets.use = "B3: MNPs Liver (infected mice)")


pathways.show <- c("MHC-II") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

sources.use <- "B2: Activated CD4s Liver (infected mice)"
targets.use <- "B3: MNPs Liver (infected mice)"
netVisual_bubble(cellchat, 
                 sources.use = sources.use, 
                 targets.use = targets.use, 
                 signaling = c("FASLG", "IFN-II", "LIGHT", "APRIL"), 
                 remove.isolate = FALSE)


ptm = Sys.time()
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

plotGeneExpression(cellchat, signaling = "APRIL", enriched.only = TRUE, type = "violin")







### RCTD

devtools::install_github("dmcable/spacexr")

setwd('/data/gpfs/projects/punim1236/liver_curio')
install.packages("tidyverse")
library(cowplot)
library(tidyverse)
library(spacexr)
library(Seurat)
library(spacexr)
library(Matrix)
library(Rcpp)

# Set working directory
setwd('/data/gpfs/projects/punim1236/liver_curio')

# Load and update reference datasets
load(file = 'tabula_muris/droplet_Liver_seurat_tiss.Robj')
tiss <- UpdateSeuratObject(tiss)

# Subset tiss to keep desired cell types
tiss_v2 <- subset(tiss, cells = tiss@meta.data %>% 
                    filter(cell_ontology_class %in% c('duct epithelial cell', 
                                                      'endothelial cell of hepatic sinusoid', 
                                                      'hepatocyte')) %>% 
                    rownames())
tiss_v2 <- SetIdent(tiss_v2, value = tiss_v2$cell_ontology_class)

# Verify cell types in tiss_v2
table(tiss_v2@meta.data$cell_ontology_class)

# Load and subset immune dataset
immune <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS")
immune <- SetIdent(immune, value = immune$annotation)

table(immune@meta.data$annotation)

tiss <- merge(tiss_v2,immune)

# Join layers to combine counts into a single layer
tiss <- JoinLayers(tiss)

# Verify the layers in the merged object
Layers(tiss[["RNA"]])

# Assign cell types to a new metadata column
tiss$cell_type <- tiss@active.ident

# Verify the cell types in the merged object
table(tiss$cell_type)

# Clean up
rm(tiss_v2, immune)
gc()




### do RCTD ###

## import data
puck1_loc <- '/data/projects/punim1236/liver_curio/data/liver1_seurat_binned_50um.rds'
puck2_loc <- '/data/projects/punim1236/liver_curio/data/liver2_seurat_binned_50um.rds'
rctd_save_loc <- "/data/gpfs/projects/punim1236/liver_curio/rctd_revision_ollie_new_ref_v2.rds"

## import puck
puck1 <- readRDS(file=puck1_loc)
puck2 <- readRDS(file=puck2_loc)


## get list of pucks
pucks <- list(puck1,puck2)

for(i in 1:2){
  
  pucks[[i]] <- UpdateSeuratObject(pucks[[i]])
  
  ## get materials for converting puck to SpatialRNA
  counts <- pucks[[i]]@assays$RNA@counts
  coords <- as.data.frame(Embeddings(pucks[[i]],reduction="map"))
  colnames(coords) <- c("x","y")
  
  nUMI <- as.numeric(pucks[[i]]@meta.data$nCount_RNA)
  names(nUMI) <- colnames(counts)
  
  
  ## convert puck to SpatialRNA format for RCTD
  pucks[[i]] <- new(
    "SpatialRNA",
    counts = counts,
    coords = coords,
    nUMI = nUMI)
}

## Convert reference to class 'Reference'
cell_types <- as.factor(tiss@meta.data$cell_type)
names(cell_types) <- rownames(tiss@meta.data)

# Extract counts matrix from Seurat v5 assay
counts <- LayerData(tiss, assay = "RNA", layer = "counts")
# Ensure counts is a dgCMatrix (it usually is, but this is a safety check)
if (!is(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

# Create Reference object
reference <- Reference(
  counts = counts,
  cell_types = cell_types
)

## Create RCTD object
rctd.obj <- create.RCTD.replicates(
  pucks, reference,
  replicate_names = c('puck1','puck2'),
  max_cores = 1,
  CELL_MIN_INSTANCE = 10)

## remove
rm(puck1,puck2,pucks,reference)
gc()



## run.RCTD - this is the step that runs RCTD, takes around 1 day for a typical puck of 30-40k spots.
rctd.obj <- run.RCTD.replicates(rctd.obj, doublet_mode = 'doublet')
saveRDS(rctd.obj,file=rctd_save_loc)











