RSCR=/data/scratch/projects/punim1236/python_envs/rseurat/bin/Rscript; "$RSCR" - <<'RS'
suppressPackageStartupMessages(library(Seurat))
x <- readRDS('/data/projects/punim1236/liver_scRNAseq/sven_liver_Tcells_combined_MNP1_MNP2_forRCTD.RDS')
x$celltype_clean <- trimws(as.character(x$annotation))
DefaultAssay(x) <- 'RNA'
d <- Seurat::GetAssayData(x, assay='RNA', layer='data')
genes <- c('Ifng','Ifngr1','Ifngr2','Thy1','Adgre5','Cd6','Alcam')
groups <- c('naive CD4','infected CD4','MNP1','MNP2','uninfected MNP')
res <- do.call(rbind, lapply(genes, function(g) {
  do.call(rbind, lapply(groups, function(gp) {
    cells <- colnames(x)[x$celltype_clean == gp]
    v <- d[g, cells, drop=TRUE]
    data.frame(gene=g, group=gp, mean=mean(v), pct_expressed=mean(v>0)*100, n=length(cells))
  }))
}))
write.csv(res, '/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/cellchat_latest_MNP2_CD4_20261023_matched_original_params/example_LR_expression_support.csv', row.names=FALSE)
print(res, row.names=FALSE, digits=4)
RS