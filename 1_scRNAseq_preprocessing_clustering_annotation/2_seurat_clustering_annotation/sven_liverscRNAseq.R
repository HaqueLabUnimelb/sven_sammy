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
# T cell marker
FeaturePlot(sven, "Cd3e", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cd4", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cd8a", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Ifng", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

FeaturePlot(sven, "Cd44", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Zeb2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cx3cr1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cxcr6", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Mki67", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# Cytotoxic T cell - Cluster 0
FeaturePlot(sven, "Gzma", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

# Cluster 3 -
FeaturePlot(sven, "Cxcr3", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# Cluster 7 - Nk1.1
FeaturePlot(sven, "Klrb1c", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Zbtb16", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# CD122
FeaturePlot(sven, "Il2rb", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Tox", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# Cluster 6 - Treg
FeaturePlot(sven, "percent.mt", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Foxp3", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# Cluster 7 - Tcrg - GD T cells
FeaturePlot(sven, "Tcrg-C", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# Prolifertaing T cells - Cluster 6
FeaturePlot(sven, "Mki67", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

# Cluster 10
FeaturePlot(sven, "Ms4a2", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Fcer1a", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cd63", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Cd9", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# MNP marker
FeaturePlot(sven, "Lyz2", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# MNP marker (Cd11b)
FeaturePlot(sven, "Itgam", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# MNP marker (Cd11c)
FeaturePlot(sven, "Itgax", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

#MNP markers from Pham et al
#MNP2
FeaturePlot(sven, "Clec4n", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Ly6i", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
#MNP1
FeaturePlot(sven, "F13a1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Ly6i", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

FeaturePlot(sven, "Nrg1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Ly6c2", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven, "Ms4a3", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


#pDCs
FeaturePlot(sven.refined, "Siglech", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Klk1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Cox6a2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

#DCs
FeaturePlot(sven.refined, "Cst3", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Tbc1d4", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Flt3", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Wdfy4", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Naaa", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# CSR signature? Maybe a bit in cluster 7?
MNP2 = c('Adgre4', 'Clec4n', 'Ly6i', 'Pla2g7', 'Treml4')

MNP2_ <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/MNP2.txt")
MNP1_ <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/MNP1.txt")
pDCs_ <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/pDCs.txt")
DCs_ <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/DCs.txt")

sven.refined <- AddModuleScore(sven.refined,
                          features = DCs_,
                          name="DCs_")


FeaturePlot(sven.refined, "DCs_1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)



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

DimPlot(sven, order = F, pt.size = 0.1, label = F, cols = c('#225ea8', '#d44094', '#a1dab4', '#225ea8', '#f4a2b6', 'black', '#225ea8', '#225ea8', '#d44094', '#d44094', 'black', 'black', '#d44094', '#d44094', '#d44094', 'black', 'black', 'black'))+NoAxes()+ coord_fixed()+ NoLegend()
DimPlot(sven.refined, group.by = "hash.ID", order = F, pt.size = 0.1, label = F, cols = c('#f4a2b6', '#f4a2b6', '#d44094', '#d44094', '#a1dab4', '#a1dab4', '#225ea8', '#225ea8', 'white', 'white'))+NoAxes()+ coord_fixed()+ NoLegend()

DimPlot(sven.refined, order = F, pt.size = 0.1, label = T, split.by = 'hash.ID')+NoAxes()+ coord_fixed()

Tcells.hash <- subset(sven.refined, cells = sven.refined@meta.data %>% filter(hash.ID %in% c('B1: Naive CD4s Liver (infected mice)', 'B7: Naive CD4s Spleen (infected mice)', 'B2: Activated CD4s Liver (infected mice)', 'B8: Activated CD4s Spleen (infected mice)')) %>% rownames())
DimPlot(Tcells.hash, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

MNPs.hash <- subset(sven.refined, cells = sven.refined@meta.data %>% filter(hash.ID %in% c('B4: Naive MNPs Liver (uninfected mice)', 'B10: Naive MNPs Spleen (uninfected mice)', 'B3: MNPs Liver (infected mice)', 'B9: MNPs Spleen (infected mice)')) %>% rownames())
DimPlot(MNPs.hash, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

saveRDS(sven.refined, file = "/data/projects/punim1236/liver_scRNAseq/sven_refined.RDS" )

sven.liver <- subset(sven.refined, cells = sven.refined@meta.data %>% filter(hash.ID %in% c('B4: Naive MNPs Liver (uninfected mice)', 'B3: MNPs Liver (infected mice)', 'B1: Naive CD4s Liver (infected mice)', 'B2: Activated CD4s Liver (infected mice)')) %>% rownames())

Infected.liver.hash <- subset(sven.refined, cells = sven.refined@meta.data %>% filter(hash.ID %in% c('B3: MNPs Liver (infected mice)', 'B2: Activated CD4s Liver (infected mice)')) %>% rownames())

DimPlot(Infected.liver.hash, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
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

saveRDS(sven.liver, file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_ref_for_RCTD.RDS" )
sven.liver <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_ref_for_RCTD.RDS")

# Remove the 'other cells' to leave behind an enriched population of liver CD4s and MNPs for spatial deconvolution onto Curio slice
Infected.liver.hash <- subset(sven.liver, cells = sven.liver@meta.data %>% filter( %in% c('B3: MNPs Liver (infected mice)', 'B2: Activated CD4s Liver (infected mice)')) %>% rownames())



Infected.liver.hash.refined <- subset(Infected.liver.hash, cells = Infected.liver.hash@meta.data %>% filter(seurat_clusters %in% c('1', '0', '6', '7')) %>% rownames())

DimPlot(Infected.liver.hash.refined, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

table(Infected.liver.hash.refined@meta.data$hash.ID)

# MNPs - 618 cells, activated CD4s - 666 cells

DimPlot(Infected.liver.hash.refined, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
DimPlot(Infected.liver.hash.refined, group.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
FeaturePlot(Infected.liver.hash.refined, "Ifng", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

saveRDS(Infected.liver.hash.refined, file = "/data/projects/punim1236/liver_scRNAseq/sven_refined_infected_liver.RDS" )
Infected.liver.hash.refined <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_refined_infected_liver.RDS" )

# Sammy requested UMAPs of the following genes

sven.refined <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_refined.RDS" )

Infected.liver.hash.refined <- subset(sven.refined, cells = sven.refined@meta.data %>% filter( %in% c('1', '0', '6', '7')) %>% rownames())


# Fas
FeaturePlot(sven.refined, "Fas", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Fas ligand
FeaturePlot(sven.refined, "Fasl", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# TNF
FeaturePlot(sven.refined, "Tnf", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# TNFR1
FeaturePlot(sven.refined, "Tnfrsf1a", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# TNFR2
FeaturePlot(sven.refined, "Tnfrsf1b", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# TRAIL
FeaturePlot(sven.refined, "Tnfsf10", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# DR4/TRAILR - doesn't exist?
FeaturePlot(sven.refined, "Tnfrsf10a", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# DR5/TRAILR
FeaturePlot(sven.refined, "Tnfrsf10b", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# MHC-II
FeaturePlot(sven.refined, "H2-Ab1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# T cell receptor alpha constant
FeaturePlot(sven.refined, "Trac", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# IFNG
FeaturePlot(sven.refined, "Ifng", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# IFNGR
FeaturePlot(sven.refined, "Ifngr1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Ifngr2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Perforin
FeaturePlot(sven.refined, "Prf1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Granzymes
FeaturePlot(sven.refined, "Gzma", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Gzmb", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Zeb2 
FeaturePlot(sven.refined, "Zeb2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# CX3CR1
FeaturePlot(sven.refined, "Cx3cr1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# CXCR6
FeaturePlot(sven.refined, "Cxcr6", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(sven.refined, "Klrg1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


FeaturePlot(sven.refined, "Klrg1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

# Jan Schroeder list of DE genes from bulkRNAseq in IFNg+ CD4+
de.list <- readRDS(file = "/data/projects/punim1236/liver_scRNAseq/BulkRNAseq/sven_dge2.Rds")

Infected.liver.hash.refined <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_refined_infected_liver.RDS")
DimPlot(Infected.liver.hash.refined, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
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


seurat_obj <- subset(sven.liver, cells = sven.liver@meta.data %>% filter(seurat_clusters %in% c('0', '1', '2', '4', '6', '7', '8', '9', '12', '13', '14')) %>% rownames())

#seurat_obj <- subset(seurat_obj, cells = seurat_obj@meta.data %>% filter(seurat_clusters %in% c('0', '3', '5', '6', '7', '8')) %>% rownames())

cytotoxic.de.cluster6.markers <- FindMarkers(seurat_obj, ident.1 = "6", ident.2 = NULL)
write.csv(cytotoxic.de.cluster6.markers, file = "/data/projects/punim1236/liver_scRNAseq/cytotoxic.de.cluster6.markers.csv")

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

seurat_obj <- AddModuleScore(seurat_obj,
                               features = MNP1_,
                               name="mnp1_")
seurat_obj <- AddModuleScore(seurat_obj,
                             features = MNP2_,
                             name="mnp2_")
seurat_obj <- AddModuleScore(seurat_obj,
                             features = pDCs_,
                             name="pdc_")


FeaturePlot(seurat_obj, "mnp1_1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "MNP1.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "mnp2_1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "MNP2.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "pdc_1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "pDC.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "Flt3", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Flt3.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")
FeaturePlot(seurat_obj, "H2-Ab1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "H2-Ab1.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")


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

FeaturePlot(seurat_obj, "Mki67", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Fas
FeaturePlot(seurat_obj, "Fas", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Fas ligand
FeaturePlot(seurat_obj, "Fasl", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# TNF
FeaturePlot(seurat_obj, "Tnf", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# TNFR1
FeaturePlot(seurat_obj, "Tnfrsf1a", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# TNFR2
FeaturePlot(seurat_obj, "Tnfrsf1b", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# TRAIL
FeaturePlot(seurat_obj, "Tnfsf10", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# DR4/TRAILR - doesn't exist?
FeaturePlot(seurat_obj, "Tnfrsf10a", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# DR5/TRAILR
FeaturePlot(seurat_obj, "Tnfrsf10b", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# MHC-II
FeaturePlot(seurat_obj, "H2-Ab1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# T cell receptor alpha constant
FeaturePlot(seurat_obj, "Cd4", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# IFNG
FeaturePlot(seurat_obj, "Ifng", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# IFNGR
FeaturePlot(seurat_obj, "Ifngr1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Ifngr2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

FeaturePlot(seurat_obj, "Maf", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Jun", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Zeb2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# Assign cluster IDs in prep for cell-type deconvolution with RCTD
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(annotation = str_replace_all(seurat_clusters, c("\\b0\\b"=" naive CD4", "\\b1\\b"="MNP infected", "\\b2\\b"="naive MNP", "\\b3\\b"="effector CD4", "\\b4\\b"="DCs", "\\b5\\b"="effector CD4", "\\b6\\b"="MNP infected", "\\b7\\b"="cytotoxic CD4", "\\b8\\b"="DCs", "\\b9\\b"="pDCs", "\\b10\\b"="proliferating CD4", "\\b11\\b"="proliferating CD4", "\\b12\\b"="DCs", "\\b13\\b"="DCs", "\\b14\\b"="DCs"), ))

seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(annotation = str_replace_all(seurat_clusters, c("\\b0\\b"=" naive CD4", "\\b1\\b"="MNP", "\\b2\\b"="MNP", "\\b3\\b"="effector CD4", "\\b4\\b"="DCs", "\\b5\\b"="effector CD4", "\\b6\\b"="MNP", "\\b7\\b"="cytotoxic CD4", "\\b8\\b"="DCs", "\\b9\\b"="pDCs", "\\b10\\b"="proliferating CD4", "\\b11\\b"="proliferating CD4", "\\b12\\b"="DCs", "\\b13\\b"="DCs", "\\b14\\b"="DCs"), ))

DimPlot(seurat_obj, reduction = "umap", label = F, group.by = "annotation")+NoLegend() +NoAxes()+ coord_fixed()
ggsave(filename = "liver_scRNAseq_UMAP_SEURAT_CLUSTERS_updated.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")

saveRDS(seurat_obj, file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_new_cytotoxic_cluster.RDS")
seurat_obj <- readRDS(file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_new_cytotoxic_cluster.RDS")

table(seurat_obj$annotation)

seurat_obj <- subset(seurat_obj, cells = seurat_obj@meta.data %>% filter(annotation %in% c(' naive CD4', 'effector CD4', 'cytotoxic CD4', 'MNP', 'proliferating CD4')) %>% rownames())
DimPlot(seurat_obj, reduction = "umap", label = T, group.by = "annotation")+NoLegend() +NoAxes()+ coord_fixed()

saveRDS(seurat_obj, file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS")
seurat_obj <- readRDS(file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_new_cytotoxic_cluster_forRCTD_v2.RDS")


# Take just CD4s and run DEG on cytotoxic (can't use MNPs for this, confounds analysis)
seurat_obj2 <- subset(seurat_obj, cells = seurat_obj@meta.data %>% filter(annotation %in% c(' naive CD4', 'effector CD4', 'cytotoxic CD4')) %>% rownames())
DimPlot(seurat_obj2, reduction = "umap", label = T, group.by = "annotation")+NoAxes()+ coord_fixed()
table(seurat_obj$annotation)

Idents(seurat_obj2) <- "annotation"
cytotoxic.de.markers <- FindMarkers(seurat_obj2, ident.1 = "cytotoxic CD4", ident.2 = NULL)

# Klrg1

write.csv(cytotoxic.de.markers, file = "/data/projects/punim1236/liver_scRNAseq/cytotoxic.de.markers.csv")


immune <- readRDS("/data/projects/punim1236/spleen_st/Seurat_objects/Manuscript_reference_D7.rds")
table(immune$cell_type)
Idents(immune) <- "cell_type"
immune2 <- subset(immune, cells = immune@meta.data %>% filter(cell_type %in% c('CD4+ iNKT', 'CD4+ T bystander', 'CD4+ T proliferating', 'CD4+ Tcmp', 'CD4+ Tfh', 'CD4+ Th1', 'CD4+ Treg')) %>% rownames())
table(immune2$cell_type)
Th1.de.markers <- FindMarkers(immune2, ident.1 = "CD4+ Th1", ident.2 = NULL)




write.csv(Th1.de.markers, file = "/data/projects/punim1236/liver_scRNAseq/Th1.de.markers.csv")


# Take just the T cells used in the initial RCTD
immune3 <- subset(immune, cells = immune@meta.data %>% filter(cell_type %in% c('CD4+ iNKT', 'CD4+ T bystander', 'CD4+ Th1', 'CD8+ T bystander')) %>% rownames())
table(immune3$cell_type)
Th1.de.markers.v2 <- FindMarkers(immune3, ident.1 = "CD4+ Th1", ident.2 = NULL)

write.csv(Th1.de.markers.v2, file = "/data/projects/punim1236/liver_scRNAseq/Th1.de.markers.v2.csv")


# Add gene names as a column
Th1.de.markers.v2$gene <- rownames(Th1.de.markers.v2)

library(ggrepel)
genes_to_label <- c("Cx3cr1", "Zeb2", "Klrg1", "S1pr5")

# Add a column to indicate significance (customize thresholds as needed)
Th1.de.markers.v2$gene_color <- "Not Significant"  # default

Th1.de.markers.v2$gene_color[Th1.de.markers.v2$p_val_adj < 0.05 & abs(Th1.de.markers.v2$avg_log2FC) > 0.5] <- "Significant"
Th1.de.markers.v2$gene_color[Th1.de.markers.v2$gene %in% genes_to_label] <- "Labelled"




library(ggplot2)
library(ggrepel)

ggplot(Th1.de.markers.v2, aes(x = avg_log2FC, y = -log10(p_val_adj), color = gene_color)) +
  geom_point(alpha = 0.8, size = 1.0) +
  scale_color_manual(values = c("Not Significant" = "gray",
                                "Significant" = "red",
                                "Labelled" = "black")) +
  geom_text_repel(
    data = subset(Th1.de.markers.v2, gene %in% genes_to_label),
    aes(label = gene),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf
  ) +
  theme_minimal() +
  labs(title = "Volcano Plot with Highlighted Genes",
       x = "Log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme(plot.title = element_text(hjust = 0.5))


# Jan Schroeder DGE list from bulkRNAseq of act T cells comapred to WT T cells

  
markers <- read.table(file = "/data/projects/punim1236/liver_scRNAseq/Act_Tcell_sig_bulkRNAseq.txt")

sven.liver <- AddModuleScore(sven.liver,
                          features = markers,
                          name="Act_Tcell_markers")
FeaturePlot(sven.liver, "Act_Tcell_markers1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
  








# MNP subsetting for Justin

seurat_obj <- subset(MNPs.hash, cells = MNPs.hash@meta.data %>% filter(seurat_clusters %in% c('1', '4', '8', '9', '12', '13', '14')) %>% rownames())

DimPlot(seurat_obj, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
DimPlot(seurat_obj, split.by = "hash.ID", order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

# Recluster these populations

seurat_obj <- SCTransform(object = seurat_obj, verbose = FALSE)

# Dimensional reduction - elbow plot and dimheatmap to help visualise number of dimensions to use for UMAP
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
ElbowPlot(seurat_obj, ndims = 50)
DimHeatmap(seurat_obj, dims = 1:20, cells = 200, balanced = TRUE)

# Determine percent of variation associated with each PC
pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# RunUMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:17)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:17)
seurat_obj2 <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj3 <- FindClusters(seurat_obj, resolution = 0.7)

DimPlot(seurat_obj3, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

DimPlot(seurat_obj2, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
DimPlot(seurat_obj2, order = F, pt.size = 0.1, label = T, split.by = "hash.ID")+NoAxes()+ coord_fixed()+ NoLegend()


FeaturePlot(seurat_obj3, "Lyz2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

MNP.markers <- FindAllMarkers(seurat_obj3, only.pos = TRUE)
write.csv(MNP.markers, file = "/data/projects/punim1236/liver_scRNAseq/MNP.markers.csv")

MNP.markers <- read.csv(file = "/data/projects/punim1236/liver_scRNAseq/MNP.markers.csv")

# Cluster 0
FeaturePlot(seurat_obj3, "Gbp2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Cluster 1
FeaturePlot(seurat_obj3, "Fn1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj3, "Apoe", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

saveRDS(seurat_obj3, file = "/data/projects/punim1236/liver_scRNAseq/MNPs/reclustered_MNPs.RDS")

DimPlot(seurat_obj, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

seurat_obj3 <- FindClusters(seurat_obj, resolution = 0.7)

DimPlot(seurat_obj, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
MNP.markers_same_UMAP <- FindAllMarkers(seurat_obj, only.pos = TRUE)

write.csv(MNP.markers_same_UMAP, file = "/data/projects/punim1236/liver_scRNAseq/MNP.markers_same_UMAP.csv")

MNP.markers_same_UMAP <- read.csv(file = "/data/projects/punim1236/liver_scRNAseq/MNP.markers_same_UMAP.csv")

# Do the split by hash.ID for each cluster, may inform liver and spleen specific or naive vs activated MNPs

clustercounts <- table(seurat_obj$hash.ID, seurat_obj$seurat_clusters)
clustercounts
write.csv(clustercounts, file = "/data/projects/punim1236/liver_scRNAseq/MNPs/clustercounts_MNPs.csv")


# Markers for Justin
FeaturePlot(MNPs, "Ly6c2", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Itgam", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Itgax", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Adgre1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Mertk", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Timd4", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Spn", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Fcgr1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Ly6g", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Cd9", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "H2-Ab1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Fas", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Tnfrsf1a", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Vsig4", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Clec1b", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Cx3cr1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

#pDCs
FeaturePlot(MNPs, "Siglech", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Klk1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Cox6a2", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

#DCs
FeaturePlot(MNPs, "Cst3", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Tbc1d4", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Flt3", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Wdfy4", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(MNPs, "Naaa", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)


# CSR signature? Maybe a bit in cluster 7?
MNP2 = c('Adgre4', 'Clec4n', 'Ly6i', 'Pla2g7', 'Treml4')

MNP2_ <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/MNP2.txt")
MNP1_ <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/MNP1.txt")
pDCs_ <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/pDCs.txt")
DCs_ <- read.table(file = "/data/gpfs/projects/punim1236/liver_scRNAseq/Ollie/DCs.txt")

seurat_obj3 <- AddModuleScore(seurat_obj3,
                               features = MNP2_,
                               name="MNP2_")

FeaturePlot(seurat_obj3, "DCs_1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj3, "pDCs_1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj3, "MNP1_1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj3, "MNP2_1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

DimPlot(seurat_obj3, order = F, pt.size = 0.1, label = F, cells.highlight=WhichCells(seurat_obj3, idents = c("4")), cols.highlight = c("#00BD5C"), cols= "grey")+NoAxes()+ coord_fixed()+ NoLegend()

MNPs <- subset(seurat_obj3, cells = seurat_obj3@meta.data %>% filter(seurat_clusters %in% c('0', '1', '2', '3', '12')) %>% rownames())
DimPlot(MNPs, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
DimPlot(MNPs, order = F, pt.size = 0.1, label = T, split.by = 'hash.ID')+NoAxes()+ coord_fixed()+ NoLegend()

MNP.markers <- FindAllMarkers(MNPs, only.pos = TRUE)
write.csv(MNP.markers, file = "/data/projects/punim1236/liver_scRNAseq/MNPs/subsetted_MNP.markers.csv")

# Signatures for Justin

# Cd11c
VlnPlot(MNPs, "Itgax") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()

# Cd11c and Cd11b
Cd11b_c = c('Itgax', 'Itgam')
MNPs <- AddModuleScore(MNPs,
                              features = list(Cd11b_c),
                              name="Cd11b_c_real_one")

FeaturePlot(MNPs, "Cd11b_c_real_one1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
VlnPlot(MNPs, "Cd11b_c_real_one1") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()

# Recruited macrophages
recruited_macs = c('Ly6c2', 'Itgax', 'Itgam', 'Fcgr1', 'Adgre1')
MNPs <- AddModuleScore(MNPs,
                       features = list(recruited_macs),
                       name="recruited_macs")

FeaturePlot(MNPs, "recruited_macs1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
VlnPlot(MNPs, "recruited_macs1") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()

# Kupffer cells
Kuppfer_cells = c('Vsig4', 'Clec1b', 'Timd4', 'Mertk', 'Adgre1', 'Fcgr1')
MNPs <- AddModuleScore(MNPs,
                       features = list(Kuppfer_cells),
                       name="Kuppfer_cells")

FeaturePlot(MNPs, "Kuppfer_cells1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
VlnPlot(MNPs, "Kuppfer_cells1") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()


# Splenic resident macs
splenic_res_macs = c('Timd4', 'Mertk', 'Adgre1', 'Fcgr1')
MNPs <- AddModuleScore(MNPs,
                       features = list(splenic_res_macs),
                       name="splenic_res_macs")

FeaturePlot(MNPs, "splenic_res_macs1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
VlnPlot(MNPs, "splenic_res_macs1") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()

FeaturePlot(MNPs, "Mertk", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

saveRDS(MNPs, file = "/data/projects/punim1236/liver_scRNAseq/MNPs/MNPs_for_Justin.RDS")

MNPs <- readRDS(file = "/data/projects/punim1236/liver_scRNAseq/MNPs/MNPs_for_Justin.RDS")


# Re-cluster MNPs - we thinkw e cna see a specific cluster in cluster 3 at the bottom right, that resembles resident macrophages
# After doing this at res of 2.0 we eventually get a cluster that seems to encompass the resident macs genes, but there are so many other clusters, looks a bit noisy
# Probably best just to leave the clustering as it is, don't think it really matters, can just say this region of UMAP resembles resident macs
DimPlot(MNPs, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()
MNPs_v2 <- FindClusters(MNPs, resolution = 2.0)
DimPlot(MNPs_v2, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()+ NoLegend()

# Next, add cell death markers to recrutied and resident macs signature

# Cell death markers signature
cell_death = c('H2-Ab1', 'Fas', 'Tnfrsf1a')
MNPs <- AddModuleScore(MNPs,
                       features = list(cell_death),
                       name="cell_death")

FeaturePlot(MNPs, "cell_death1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
VlnPlot(MNPs, "cell_death1") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()

# Recruited macrophages
recruited_macs_death = c('Ly6c2', 'Itgax', 'Itgam', 'Fcgr1', 'Adgre1', 'H2-Ab1', 'Fas', 'Tnfrsf1a')
MNPs <- AddModuleScore(MNPs,
                       features = list(recruited_macs_death),
                       name="recruited_macs_death")

FeaturePlot(MNPs, "recruited_macs_death1", order = F)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
VlnPlot(MNPs, "recruited_macs_death1") +
  stat_summary(fun.y = mean, geom='point', size = 25, colour = "black", shape = 95)+ NoLegend()

MNPs.markers.Justin <- FindAllMarkers(MNPs, only.pos = TRUE)
write.csv(MNPs.markers.Justin, file = "/data/projects/punim1236/liver_scRNAseq/MNPs.markers.Justin.csv")


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

features = c("Lyz2", "Itgam", "MNP1", "MNP2", "Itgax", "Flt3", "H2-Ab1", "pDCs", "Cd8a", "Cd4", "Tbx21", "Cxcr6", "Cxcr3", "Il23r", "Pecam1", "Btla", "Tnfrsf14", "Ifng", "Zeb2", "Cx3cr1", "Klrg1", "S1pr5", "Fasl", "Gzma", "Gzmb", "Prf1", "Mki67")

Idents(seurat_obj) <- "seurat_clusters"

DotPlot(seurat_obj, features = features, dot.min = 0, scale = T, scale.by = 'size')& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1)& RotatedAxis()& coord_flip()
a <- DotPlot(seurat_obj, features = features, dot.min = 0, scale = T, scale.by = 'size')& coord_fixed()& RotatedAxis()& coord_flip()
b <- a + scale_colour_gradient2(low = "#253494", mid = "white", high = "#d7191c") +
    scale_x_discrete(limits=rev)
b
ggsave(filename = "Dotplot.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 15, height = 15, units = "cm")

FeaturePlot(seurat_obj, "Gzma", order = T, split.by = 'hash.ID')&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

# Sven wants to see if T-bet, Eomes, Blimp-1 upreg in effetcro T cells
FeaturePlot(seurat_obj, "Cd3d", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Tbx21", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Eomes", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Prdm1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)

seurat_obj$annotation <- factor(seurat_obj$annotation,
                          levels = c(" naive CD4","effector CD4","cytotoxic CD4","proliferating CD4", "MNP"))


table(seurat_obj$annotation)
VlnPlot(seurat_obj, "Tbx21", group.by = "annotation")& coord_fixed()
VlnPlot(seurat_obj, "Eomes", group.by = "annotation")& coord_fixed()
VlnPlot(seurat_obj, "Prdm1", group.by = "annotation")& coord_fixed()

DimPlot(seurat_obj, order = F, pt.size = 0.1, label = F, group.by = "annotation")+NoAxes()+ coord_fixed()+ NoLegend()


# Checking for classicla monocytes in scRNAseq data

seurat_obj <- readRDS(file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_for_paper.RDS")

# Plot by hash.ID, naive vs infected MNPs
DimPlot(seurat_obj, order = F, pt.size = 0.1, label = F, group.by = "hash.ID")+NoAxes()+ coord_fixed()
DimPlot(seurat_obj, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()



# Classical monocyte markers
FeaturePlot(seurat_obj, "Ly6c2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Ly6c2.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")

FeaturePlot(seurat_obj, "Ccr2", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Ccr2.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")

FeaturePlot(seurat_obj, "Plac8", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Chil3", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "F13a1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Vcan", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Upon infection classical monocytes should upregulate the following markers
FeaturePlot(seurat_obj, "Fas", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "H2-Ab1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "H2-Aa", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Cd74", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Fcgr1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
# Other markers (resident macrophages)
FeaturePlot(seurat_obj, "Itgax", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
FeaturePlot(seurat_obj, "Adgre1", order = T)&NoAxes()& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95)
ggsave(filename = "Adgre1.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 12, height = 12, units = "cm")



## load immune data and subset into MNP1 and MNP2 (MNP2 maybe more closely associated with T cells, stronger MHCII)
immune <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_forRCTD.RDS")
DimPlot(immune, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()
table(immune@meta.data$annotation)
DimPlot(immune, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()

immune@meta.data <- immune@meta.data %>%
  mutate(annotation = str_replace_all(seurat_clusters, c("\\b0\\b"=" naive CD4", "\\b1\\b"="MNP1", "\\b2\\b"="uninfected MNP", "\\b3\\b"="effector CD4", "\\b5\\b"="effector CD4", "\\b6\\b"="MNP2", "\\b7\\b"="cytotoxic CD4", "\\b10\\b"="proliferating CD4", "\\b11\\b"="proliferating CD4"), ))

table(immune@meta.data$annotation)
DimPlot(immune, order = F, pt.size = 0.1, label = T, group.by = "annotation")+NoAxes()+ coord_fixed()


saveRDS(immune, file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS")

immune <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS")

immune <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_split_up_MNP1_MNP2_forRCTD.RDS")
immune <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS")
immune <- readRDS("/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_forRCTD.RDS")

DimPlot(immune, order = F, pt.size = 0.1, label = F, group.by = "annotation")+NoAxes()+ coord_fixed()

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


features = c("Lyz2", "Itgam", "Itgax", "Flt3", "H2-Ab1", "Cd8a", "MNP1.1", "MNP2.1", "pDCs", "Cd4", "Pecam1", "Btla", "Tbx21", "Cxcr6", "Cxcr3", "Il23r", "Ifng", "Klrg1", "S1pr5", "Fasl", "Gzma", "Gzmb", "Prf1", "Mki67", "Adreg1", "Ly6C", "Ccr2")

features = c("Mki67", "Prf1", "Gzmb", "Gzma", "Fasl", "Fas", "S1pr5", "Klrg1", "Cx3cr1", "Zeb2", "Cd44", "Ifng", "Il23r", "Cxcr3", "Cxcr6", "Tbx21", "Btla", "Pecam1", "Cd4", "pDCs", "MNP2", "MNP1", "Cd8a", "H2-Ab1", "Flt3", "Itgax", "Itgam", "Lyz2", "Adgre1", "Ly6c2", "Ccr2")

Idents(seurat_obj) <- "annotation"

DotPlot(seurat_obj, features = features, dot.min = 0, scale = T, scale.by = 'size')& coord_fixed()& scale_colour_viridis_c(option='magma', direction=-1, begin=0.1)& RotatedAxis()& coord_flip()
a <- DotPlot(seurat_obj, features = features, dot.min = 0, scale = T, scale.by = 'size')& coord_fixed()& RotatedAxis()& coord_flip()
b <- a + scale_colour_gradient2(low = "#253494", mid = "white", high = "#d7191c") +
  scale_x_discrete(limits=rev)
b
ggsave(filename = "Dotplot_NI_revision_v3.pdf", device = "pdf", path = "/data/projects/punim1236/liver_scRNAseq/", width = 15, height = 19, units = "cm")
