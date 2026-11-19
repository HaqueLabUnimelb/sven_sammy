library(Seurat)
library(ggplot2)
library(tidyverse)

#install.packages('tidyverse')
sven <- readRDS("/data/projects/punim1236/liver_scRNAseq/GEX_HTO_Processed.rds")

DimPlot(sven, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()

table(sven$hash.ID)

# Change order of Hash.ID
sven$hash.ID <- factor(sven$hash.ID, levels = c('B4: Naive MNPs Liver (uninfected mice)', 'B10: Naive MNPs Spleen (uninfected mice)', 'B3: MNPs Liver (infected mice)', 'B9: MNPs Spleen (infected mice)', 'B1: Naive CD4s Liver (infected mice)', 'B7: Naive CD4s Spleen (infected mice)', 'B2: Activated CD4s Liver (infected mice)', 'B8: Activated CD4s Spleen (infected mice)', 'Doublet', 'Negative'))
DimPlot(sven, split.by = "hash.ID", order = F, pt.size = 0.1, label = F)+NoAxes()+ coord_fixed()+ NoLegend()

MNPs.hash <- subset(sven, cells = sven@meta.data %>% filter(hash.ID %in% c('B4: Naive MNPs Liver (uninfected mice)', 'B10: Naive MNPs Spleen (uninfected mice)', 'B3: MNPs Liver (infected mice)', 'B9: MNPs Spleen (infected mice)')) %>% rownames())
DimPlot(MNPs.hash, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
DimPlot(MNPs.hash, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

liver <- subset(sven, cells = sven@meta.data %>% filter(hash.ID %in% c('B4: Naive MNPs Liver (uninfected mice)', 'B3: MNPs Liver (infected mice)', 'B1: Naive CD4s Liver (infected mice)', 'B2: Activated CD4s Liver (infected mice)')) %>% rownames())
DimPlot(liver, split.by = "hash.ID", order = F, pt.size = 0.1, label = F)+NoAxes()+ coord_fixed()+ NoLegend()


Tcells.hash <- subset(sven, cells = sven@meta.data %>% filter(hash.ID %in% c('B1: Naive CD4s Liver (infected mice)', 'B7: Naive CD4s Spleen (infected mice)', 'B2: Activated CD4s Liver (infected mice)', 'B8: Activated CD4s Spleen (infected mice)')) %>% rownames())
DimPlot(Tcells.hash, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

DimPlot(sven, group.by = "hash.ID", order = F, pt.size = 0.1, label = F, cols = c('#f4a2b6', '#f4a2b6', '#d44094', '#d44094', '#a1dab4', '#a1dab4', '#225ea8', '#225ea8', 'white', 'white'))+NoAxes()+ coord_fixed()+ NoLegend()



DimPlot(sven, split.by = "HTO_maxID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

DimPlot(sven)+NoAxes()+ coord_fixed()+ NoLegend()

# B cell marker
FeaturePlot(sven, "Ms4a1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cd79a", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Ighm", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

# Erythrocyte marker
FeaturePlot(sven, "Hbb-bt", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Hepatocyte marker or just dying cells
FeaturePlot(sven, "Alb", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "percent.mt", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Aldob", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Cluster 15 Neutrophil markers
FeaturePlot(sven, "Mpo", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cxcr2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Ngp", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Cluster 17 Keratinocytes/epithelial cells
FeaturePlot(sven, "Stfa2l1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cstdc4", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


sven.markers <- FindAllMarkers(sven, only.pos = TRUE)
sven.markers <- read.csv(file = "/data/projects/punim1236/liver_scRNAseq/sven_cluster_markers.csv")

sven.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
DoHeatmap(sven, features = top10$gene) + NoLegend()

sven.liver <- subset(liver, cells = sven@meta.data %>% filter(seurat_clusters %in% c('0', '1', '2', '4', '6', '7', '8', '9', '12', '13', '14')) %>% rownames())

sven.refined <- subset(sven, cells = sven@meta.data %>% filter(seurat_clusters %in% c('0', '1', '2', '3', '4', '6', '7', '8', '9', '12', '13', '14')) %>% rownames())


Tcells.hash <- subset(sven.refined, cells = sven.refined@meta.data %>% filter(hash.ID %in% c('B1: Naive CD4s Liver (infected mice)', 'B7: Naive CD4s Spleen (infected mice)', 'B2: Activated CD4s Liver (infected mice)', 'B8: Activated CD4s Spleen (infected mice)')) %>% rownames())
DimPlot(Tcells.hash, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

MNPs.hash <- subset(sven.refined, cells = sven.refined@meta.data %>% filter(hash.ID %in% c('B4: Naive MNPs Liver (uninfected mice)', 'B10: Naive MNPs Spleen (uninfected mice)', 'B3: MNPs Liver (infected mice)', 'B9: MNPs Spleen (infected mice)')) %>% rownames())
DimPlot(MNPs.hash, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

saveRDS(sven.refined, file = "/data/projects/punim1236/liver_scRNAseq/sven_refined.RDS" )

sven.refined <- readRDS(file = "/data/projects/punim1236/liver_scRNAseq/sven_refined.RDS")

sven.liver <- subset(sven.refined, cells = sven.refined@meta.data %>% filter(hash.ID %in% c('B4: Naive MNPs Liver (uninfected mice)', 'B3: MNPs Liver (infected mice)', 'B1: Naive CD4s Liver (infected mice)', 'B2: Activated CD4s Liver (infected mice)')) %>% rownames())


DimPlot(sven.liver, group.by = "hash.ID", order = F, pt.size = 0.1, label = F)+NoAxes()+ coord_fixed()
DimPlot(sven.liver, split.by = "hash.ID", order = F, pt.size = 0.1, label = F)+NoAxes()+ coord_fixed()
DimPlot(sven.liver, group.by = "seurat_clusters", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()

# Assign cluster IDs in prep for cell-type deconvolution with RCTD
sven.liver@meta.data <- sven.liver@meta.data %>%
  mutate(annotation = str_replace_all(seurat_clusters, c("\\b0\\b"="liver cytotoxic CD4", "\\b1\\b"="liver MNP infected", "\\b2\\b"="liver naive CD4", "\\b3\\b"="other", "\\b4\\b"="liver MNP naive", "\\b6\\b"="liver cytotoxic CD4", "\\b7\\b"="liver cytotoxic CD4", "\\b8\\b"="other", "\\b9\\b"="other", "\\b12\\b"="other", "\\b13\\b"="other", "\\b14\\b"="other"), ))

DimPlot(sven.liver, group.by = "annotation", order = F, pt.size = 0.1, label = F)+NoAxes()+ coord_fixed()
DimPlot(sven.liver, group.by = "seurat_clusters", order = F, pt.size = 0.1, label = F)+NoAxes()+ coord_fixed()

DimPlot(sven.liver, group.by = "annotation", order = F, pt.size = 0.1, label = F, cols = c('#225ea8', '#d44094', 'darkgrey', 'darkgrey', 'grey'))+NoAxes()+ coord_fixed()
table(sven.liver@meta.data$annotation)


seurat_obj <- sven.liver

# Split up clusters (into cytotoxic and effector)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.7)

# Removed cluster 3, small population in T cell data set, base don too few cells
saveRDS(seurat_obj, file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_for_paper.RDS")

seurat_obj <- readRDS(file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_for_paper.RDS")


DimPlot(seurat_obj, reduction = "umap", label = T, group.by = "seurat_clusters")+NoAxes()+ coord_fixed()+ NoLegend()
ggsave(filename = "liver_scRNAseq_UMAP_SEURAT_CLUSTERS.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")

DimPlot(seurat_obj, reduction = "umap", label = T, split.by = "hash.ID")+NoAxes()+ coord_fixed()+ NoLegend()

DimPlot(seurat_obj, reduction = "umap", group.by = "annotation", label = F, order = F)+NoAxes()+ coord_fixed()+ NoLegend()
ggsave(filename = "liver_scRNAseq_UMAP_annotated.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")

DimPlot(seurat_obj, reduction = "umap", group.by = "hash.ID", label = F, order = F)+NoAxes()+ coord_fixed()+ NoLegend()
ggsave(filename = "liver_scRNAseq_UMAP_hashtag.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")



FeaturePlot(seurat_obj, "Cd4", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Cd4.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Cd44", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Cd44.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Tbx21", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Tbx21.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Cxcr6", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Cxcr6.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Ifng", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Ifng.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")

FeaturePlot(seurat_obj, "S1pr5", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "S1pr5.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Klrg1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Klrg1.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Cx3cr1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Cx3cr1.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Zeb2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Zeb2.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Mki67", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Mki67.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")

FeaturePlot(seurat_obj, "Lyz2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Lyz2.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Itgax", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Cd11c.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Itgam", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Cd11b.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")


FeaturePlot(seurat_obj, "Fas", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Fas.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Fasl", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Fasl.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Ifngr1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Ifngr1.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Ifngr2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Ifngr2.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Tnf", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Tnf.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Tnfrsf1a", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Tnfrsf1a.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Tnfrsf1b", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Tnfrsf1b.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")


FeaturePlot(seurat_obj, "Gzma", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Gzma.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Prf1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Prf1.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Gzmb", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Gzmb.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")


# Plotting genes onto liver data, dotplot, for revision

# Add signatures

MNP2 <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/MNP2.txt")
MNP1 <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/MNP1.txt")
pDCs <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/pDCs.txt")

seurat_obj <- AddModuleScore(seurat_obj,
                             features = MNP1,
                             name="MNP1")

# Rename directly
colnames(seurat_obj@meta.data)[
  colnames(seurat_obj@meta.data) == "MNP11"
] <- "MNP1"

seurat_obj <- AddModuleScore(seurat_obj,
                             features = MNP2,
                             name="MNP2")

# Rename directly
colnames(seurat_obj@meta.data)[
  colnames(seurat_obj@meta.data) == "MNP21"
] <- "MNP2"

seurat_obj <- AddModuleScore(seurat_obj,
                             features = pDCs,
                             name="pDCs")

# Rename directly
colnames(seurat_obj@meta.data)[
  colnames(seurat_obj@meta.data) == "pDCs1"
] <- "pDCs"

seurat_obj$hash.ID <- factor(seurat_obj$hash.ID, levels = c('B4: Naive MNPs Liver (uninfected mice)', 'B3: MNPs Liver (infected mice)', 'B2: Activated CD4s Liver (infected mice)', 'B1: Naive CD4s Liver (infected mice)'))

DimPlot(seurat_obj, reduction = "umap", label = T, group.by = "hash.ID", order = )+NoAxes()+ coord_fixed()+ NoLegend()
DimPlot(seurat_obj, reduction = "umap", label = T, split.by = "hash.ID")+NoAxes()+ coord_fixed()

FeaturePlot(seurat_obj, "pDCs", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "MNP1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "MNP2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

FeaturePlot(seurat_obj, "Cxcr3", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Il23r", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Pecam1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Btla", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Tnfrsf4", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# VlnPlot comparing expression of genes between MNP1 and MNP2: CD11c, MHC-II, FAS

VlnPlot(immune, "Itgax", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()
VlnPlot(immune, "H2-Aa", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()
VlnPlot(immune, "H2-Ab1", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()
VlnPlot(immune, "Fas", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()
VlnPlot(immune, "Ccr2", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()
VlnPlot(immune, "Adgre1", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()
VlnPlot(immune, "Fcgr1", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()


immune <- AddModuleScore(immune, features = MNP1_, name="mnp1_")
immune <- AddModuleScore(immune, features = MNP2_, name="mnp2_")

VlnPlot(immune, "mnp1_1", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()
VlnPlot(immune, "mnp2_1", group.by = "annotation") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()


# Re-annotate scRNAseq dataset for NI revision and perform dotplot

seurat_obj <- readRDS(file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_for_paper.RDS")

DimPlot(seurat_obj, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(annotation = str_replace_all(seurat_clusters, c("\\b0\\b"="naive CD4", "\\b1\\b"="MNP1", "\\b2\\b"="uninfected MNP", "\\b3\\b"="effector CD4", "\\b4\\b"="DC", "\\b5\\b"="effector CD4", "\\b6\\b"="MNP2", "\\b7\\b"="effector CD4", "\\b8\\b"="DC", "\\b9\\b"="pDC", "\\b10\\b"="effector CD4", "\\b11\\b"="effector CD4", "\\b12\\b"="DC", "\\b13\\b"="DC", "\\b14\\b"="DC")))
DimPlot(seurat_obj, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()
DimPlot(seurat_obj, order = F, pt.size = 0.1, label = F, group.by = "annotation")+NoAxes()+ coord_fixed()

# Re-order annotation

seurat_obj$annotation <- factor(seurat_obj$annotation, levels = c('uninfected MNP', 'MNP1', 'MNP2', 'DC', 'pDC', 'naive CD4', 'effector CD4'))
DimPlot(seurat_obj, order = F, pt.size = 0.1, label = F, group.by = "annotation", cols = c('#feb24c', '#fc4e2a', '#b10026', '#c994c7', '#8c6bb1', '#a1dab4', '#225ea8'))+NoAxes()+ coord_fixed()

ggsave(filename = "liver_scRNAseq_UMAP_annotated_NI_revision.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")

DimPlot(seurat_obj, order = T, pt.size = 0.1, label = F, group.by = "hash.ID", cols = c('#feb24c', '#fc4e2a', '#a1dab4', '#225ea8'))+NoAxes()+ coord_fixed()
ggsave(filename = "liver_scRNAseq_UMAP_annotated_NI_revision_hashtag.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 16, height = 16, units = "cm")


features = c("Mki67", "Prf1", "Gzmb", "Gzma", "Fasl", "Fas", "S1pr5", "Klrg1", "Cx3cr1", "Zeb2", "Cd44", "Ifng", "Il23r", "Cxcr3", "Cxcr6", "Tbx21", "Btla", "Pecam1", "Cd4", "pDCs", "MNP2", "MNP1", "Cd8a", "H2-Ab1", "Flt3", "Itgax", "Itgam", "Lyz2", "Adgre1", "Ly6c2", "Ccr2")

Idents(seurat_obj) <- "annotation"

DotPlot(seurat_obj, features = features, dot.min = 0, scale = T, scale.by = 'size')& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1)& RotatedAxis()& coord_flip()
a <- DotPlot(seurat_obj, features = features, dot.min = 0, scale = T, scale.by = 'size')& coord_fixed()& RotatedAxis()& coord_flip()
b <- a + scale_colour_gradient2(low = "#253494", mid = "white", high = "#d7191c") +
  scale_x_discrete(limits=rev)
b
ggsave(filename = "Dotplot_NI_revision_v3.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 15, height = 19, units = "cm")
