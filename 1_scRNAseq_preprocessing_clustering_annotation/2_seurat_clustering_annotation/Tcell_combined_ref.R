library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)

#install.packages('tidyverse')
sven <- readRDS(  "/data/projects/punim1236/liver_scRNAseq/sven_liver_new_cytotoxic_cluster_forRCTD.RDS")

DimPlot(sven, order = F, pt.size = 0.1, label = T)+NoAxes()+ coord_fixed()

table(sven$annotation)

sven$annotation <- recode(
  sven$annotation,
  " naive CD4" = "Tcell",
  "cytotoxic CD4" = "Tcell",
  "effector CD4" = "Tcell",
  "proliferating CD4" = "Tcell",
  "MNP infected" = "MNP infected",
  "naive MNP" = "naive MNP"
)

DimPlot(sven, order = F, pt.size = 0.1, label = T, group.by = "annotation")+NoAxes()+ coord_fixed()
saveRDS(sven, file = "/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_forRCTD.RDS")

