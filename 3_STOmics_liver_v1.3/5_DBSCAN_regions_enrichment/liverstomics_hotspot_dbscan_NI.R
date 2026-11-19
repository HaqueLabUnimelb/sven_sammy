library(data.table); library(dplyr)

#load seurat object with cell weights

rctd.res <- readRDS(file = "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/rctd_MNP1_2_Tcell_combined_no_sinusoids.rds")

hs_out <- fread("/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_results/hotspot_spot_scores.csv", colClasses = c(barcode = "character")) |>
  as.data.frame()
rownames(hs_out) <- hs_out$barcode
mod_cols <- grep("^module_", colnames(hs_out), value = TRUE)

# check the IDs match before adding anything
sum(hs_out$barcode %in% colnames(rctd.res))     # spots shared with your SpatialExperiment
length(hs_out$barcode); ncol(rctd.res)

# Add to your R object
library(SummarizedExperiment)
colData(rctd.res)[, mod_cols] <- hs_out[colnames(rctd.res), mod_cols]

library(ggplot2); library(dplyr)

mod_cols <- grep("^module_\\d+$", colnames(hs_out), value = TRUE)


plot_module <- function(m, xlim = NULL, ylim = NULL) {
  ggplot(hs_out %>% arrange(.data[[m]]),
         aes(x, y, colour = .data[[m]])) +
    geom_point(size = 0.2) +
    scale_colour_viridis_c(option = "viridis", name = "Score", oob = scales::squish) +
    coord_fixed(xlim = xlim, ylim = ylim) +          # NULL = show everything
    theme_void() +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title       = element_text(colour = "white", hjust = 0.5),
      legend.text      = element_text(colour = "black"),
      legend.title     = element_text(colour = "black")
    ) +
    ggtitle(m)
}

plot_module("module_1")      
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module1_hotspot.png"), width = 6, height = 6, dpi = 300, bg = "white")# whole array
plot_module("module_2")      
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module2_hotspot.png"), width = 6, height = 6, dpi = 300, bg = "white")# whole array
plot_module("module_3")      
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module3_hotspot.png"), width = 6, height = 6, dpi = 300, bg = "white")# whole array
plot_module("module_4")      
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module4_hotspot.png"), width = 6, height = 6, dpi = 300, bg = "white")# whole array
plot_module("module_5")      
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module5_hotspot.png"), width = 6, height = 6, dpi = 300, bg = "white")# whole array
plot_module("module_6")      
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module6_hotspot.png"), width = 6, height = 6, dpi = 300, bg = "white")# whole array


plot_module("module_5", xlim = c(5500, 7000), ylim = c(13000, 14500))   # zoomed
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom.pdf"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom.png"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
saveRDS(rctd.res, file = "/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_results/rctd_res_hotspot_modules.RDS")

rctd.res <- readRDS(file = "/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_results/rctd_res_hotspot_modules.RDS")


#Now essentially need to perform the DBscan but this time on the module 5 score

suppressMessages({
  library(SpatialExperiment)
  library(ggplot2)
  library(viridisLite)
  library(patchwork)
  library(Matrix)
})

### Load in the RCTD output of liver stomics
setwd('/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/')

rctd_loc <- "rctd_MNP1_2_Tcell_combined_no_sinusoids.rds"
out_dir  <- "."

res <- readRDS(rctd_loc)

coords <- as.data.frame(spatialCoords(res))
colnames(coords) <- c("x", "y")

weights <- assay(res, "weights")          # cell_type x spot matrix, colSums ~1
cell_types <- rownames(weights)

## build one long data frame: x, y, cell_type, weight
df <- do.call(rbind, lapply(cell_types, function(ct) {
  data.frame(
    x = coords$x,
    y = coords$y,
    cell_type = ct,
    weight = as.numeric(weights[ct, ])
  )
}))

## keep a consistent, readable facet order (by total contribution, descending)
ct_order <- names(sort(Matrix::rowSums(weights), decreasing = TRUE))
df$cell_type <- factor(df$cell_type, levels = ct_order)

## y is flipped to match image/array orientation (row 0 at top)
y_range <- range(df$y)


### Claude DBscan script - inputting the spatial deconvolution results

library(dbscan)
library(ggplot2)

# res = data.frame with x, y and one column per cell type (as built earlier)
# 1) Pick the spots of interest, e.g. spots with a high B cell weight
library(dplyr)
sel <- hs_out %>%
  filter(module_5 > 4.0)         # use your own cell type and threshold
xy <- sel[, c("x", "y")]


# 2) Choose eps from the k-nearest-neighbour distance plot
kNNdistplot(xy, k = 15)             # look for the "elbow"
abline(h = 200, lty = 2)           # set h at the elbow, then use it as eps

# 3) Run DBSCAN
db <- dbscan(xy, eps = 200, minPts = 15)
xy$cluster <- factor(db$cluster)   # 0 = noise

# 4) Optionally drop small clusters
sizes <- table(xy$cluster)
keep  <- names(sizes[sizes >= 20 & names(sizes) != "0"])
xy$cluster[!xy$cluster %in% keep] <- "0"

# 5) Plot

# Capture all background spots in grey
bg <- unique(hs_out[, c("x", "y")])

a <- ggplot() +
  geom_point(data = bg, aes(x, y), colour = "grey90", size = 0.3) +
  geom_point(data = xy, aes(x, y, colour = cluster), size = 0.3) +
  coord_fixed() +
  theme_void()
a

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module5_4.0_k15_eps200.png"),
       a, width = 15, height = 12, dpi = 300, bg = "white")


# now wish to to do a convex hulll to loosely capture grnauloma region to then perform co-loclaisation analysis

library(dplyr)
library(ggplot2)



# xy has x, y, cluster (from dbscan); drop noise
clust <- xy %>% filter(cluster != "0")

b <- ggplot() +
  geom_point(data = bg, aes(x, y), colour = "grey85", size = 0.3) +
  geom_point(data = clust, aes(x, y, colour = cluster), size = 0.5) +
  stat_density_2d(data = clust, aes(x, y, group = cluster, colour = cluster),
                  geom = "path", contour_var = "ndensity", breaks = 0.05,
                  h = c(200, 200)) +
  coord_fixed() + theme_void()

b

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module5_4.0_k15_eps200_border.png"),   b, width = 15, height = 12, dpi = 300, bg = "white")



# fill in dbscan regions

library(dplyr)
library(ggplot2)

clust <- xy %>% filter(cluster != "0")

# one colour per cluster, reused for fill, outline and points
cl_names <- sort(unique(as.character(clust$cluster)))
cols <- setNames(scales::hue_pal()(length(cl_names)), cl_names)

# filled layer for each cluster
fills <- lapply(cl_names, function(cl) {
  stat_density_2d_filled(
    data = filter(clust, cluster == cl),
    aes(x, y), inherit.aes = FALSE,
    contour_var = "ndensity", breaks = c(0.05, 1.05),  # same 10% level as your outline
    h = c(200, 200),                                  # same smoothing as before
    fill = cols[[cl]], alpha = 1, colour = NA
  )
})

c <- ggplot() +
  geom_point(data = bg, aes(x, y), colour = "grey85", size = 0.3) +
  fills +                                              # filled clusters
  geom_point(data = clust, aes(x, y, colour = cluster), size = 0.5) +
  stat_density_2d(data = clust, aes(x, y, group = cluster, colour = cluster),
                  geom = "path", contour_var = "ndensity", breaks = 0.05,
                  h = c(200, 200)) +                   # outline on top
  scale_colour_manual(values = cols) +
  coord_fixed() + theme_void()

c

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/module5_4.0_k15_eps200_border_fill.png"),   c, width = 15, height = 12, dpi = 300, bg = "white")


### Now to extratc the DBscanned regions (the expanded region) as dataframe

library(MASS); library(dplyr)

h_val <- 200      # same as your plot
level <- 0.05      # same as your breaks
pad   <- 4 * h_val

# every unique spot position in the full dataset
pos <- hs_out %>% distinct(x, y)

hit <- lapply(cl_names, function(cl) {
  d <- filter(clust, cluster == cl)
  
  k <- MASS::kde2d(d$x, d$y, h = c(h_val, h_val), n = 400,
                   lims = c(range(d$x) + c(-pad, pad),
                            range(d$y) + c(-pad, pad)))
  k$z <- k$z / max(k$z)                        # ndensity
  
  # look up the density grid cell for every spot
  ix <- findInterval(pos$x, k$x, all.inside = TRUE)
  iy <- findInterval(pos$y, k$y, all.inside = TRUE)
  in_grid <- pos$x >= min(k$x) & pos$x <= max(k$x) &
    pos$y >= min(k$y) & pos$y <= max(k$y)
  
  keep <- in_grid & k$z[cbind(ix, iy)] >= level
  pos[keep, ] %>% mutate(cluster = cl)
}) %>% bind_rows()

# attach back to the long data: all cell type weights for those spots
sub_df <- df %>% inner_join(hit, by = c("x", "y"))


# chekc this subsetting has worked

library(ggplot2); library(dplyr)

d <- ggplot() +
  geom_point(data = bg, aes(x, y), colour = "grey85", size = 0.3) +
  geom_point(data = hit, aes(x, y, colour = cluster), size = 0.6) +        # all spots kept
  geom_point(data = clust, aes(x, y), colour = "black", size = 0.3) +      # original DBSCAN spots
  scale_colour_manual(values = cols) +
  coord_fixed() + theme_void()
d
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/afterdbscan_subsetted_check.png"),   d, width = 15, height = 12, dpi = 300, bg = "white")




library(ggplot2); library(dplyr)

ct <- "infected CD4"   # change to any cell type
ct <- "MNP2"
ct <- "hepatocyte"
ct <- "uninfected MNP"
ct <- "MNP1"
ct <- " naive CD4"
d <- ggplot() +
  geom_point(data = filter(df, cell_type == ct) %>% arrange(weight),   # high weights drawn last
             aes(x, y, colour = weight), size = 0.1) +
  stat_density_2d(data = clust, aes(x, y, group = cluster),
                  geom = "path", contour_var = "ndensity", breaks = 0.1,
                  h = c(150, 150), colour = "black", linewidth = 0.5) +
  scale_colour_viridis_c(option='magma', direction=-1, begin=0.1, end=0.95) +
  coord_fixed() + theme_void()

d

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/afterdbscan_subsetted_check_infectedCD4.png"),   d, width = 15, height = 12, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/afterdbscan_subsetted_check_MNP2.png"),   d, width = 15, height = 12, dpi = 300, bg = "white")

xl <- c(5500, 7000)
yl <- c(13000, 14500)

d <- ggplot() +
  geom_point(data = filter(df, cell_type == ct) %>% arrange(weight),   # high weights drawn last
             aes(x, y, colour = weight), size = 0.2) +
  stat_density_2d(data = clust, aes(x, y, group = cluster),
                  geom = "path", contour_var = "ndensity", breaks = 0.1,
                  h = c(150, 150), colour = "black", linewidth = 0.5) +
  scale_colour_viridis_c(option = "magma", direction = -1, begin = 0.1, end = 0.95) +
  coord_fixed(xlim = xl, ylim = yl) + theme_void()

d

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_infected_CD4.pdf"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_MNP2.pdf"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_hepatocyte.pdf"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_uninfectedMNP.pdf"), width = 2.1, height = 2.1, dpi = 300, bg = "white")

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_infected_CD4.png"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_MNP2.png"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_hepatocyte.png"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_uninfectedMNP.png"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_MNP1.png"), width = 2.1, height = 2.1, dpi = 300, bg = "white")
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/zoom/module5_hotspot_zoom_naiveCD4.png"), width = 2.1, height = 2.1, dpi = 300, bg = "white")


saveRDS(hs_out, file = "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/hs_out.rds")
saveRDS(sub_df, file = "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/sub_df_containsdbscan_info.rds")

#~~~~~~~~~~~~
#~~~~~~~~~~~~
#~~~~~~~~~~~~

### co-localisation analysis on grnauloma regions, indiviudally and in comaprison to non grnauloma regions

#~~~~~~~~~~~~
#~~~~~~~~~~~~
#~~~~~~~~~~~~

library(dplyr); library(tidyr)


sub_wide <- sub_df %>%
  dplyr::select(x, y, cluster, cell_type, weight) %>%
  pivot_wider(names_from = cell_type, values_from = weight, values_fill = 0)

cell_cols <- setdiff(colnames(sub_wide), c("x", "y", "cluster"))

region_list <- sub_wide %>%
  split(.$cluster) %>%
  lapply(function(d) as.data.frame(d[, c("x", "y", cell_cols)]))

sapply(region_list, nrow)

outs_regions <- lapply(region_list, function(d)
  coloc_analysis_replicates_fast(list(d)))

outs_regions[["7"]]$plot

for (nm in names(outs_regions)) {
  ggsave(
    filename = file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/", paste0("coloc_region40um_", nm, ".pdf")),
    plot     = outs_regions[[nm]]$plot + ggtitle(paste("Region", nm)),
    width = 5, height = 5
  )
}


# r = 40
outs_pooled <- coloc_analysis_replicates_fast(unname(region_list))
outs_pooled$plot
ggsave(filename = "outs_40_all_dbscanned_regions.pdf", device = "pdf", path = "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/", width = 15, height = 15, units = "cm")

library(dplyr)

for (nm in names(outs_regions)) {
  write.csv(
    outs_regions[[nm]]$rho_mat,
    file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/", paste0("rho_region_40um", nm, ".csv"))
  )
}


# Plot below allows you to see the siz eof the neighbourhoods

library(ggplot2); library(dplyr)

r  <- 40
bg <- distinct(df, x, y)

circle <- function(cx, cy, r, id, n = 100) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t), id = id)
}

ex      <- hit %>% slice_sample(n = 5)
circles <- bind_rows(lapply(seq_len(nrow(ex)), function(i) circle(ex$x[i], ex$y[i], r, i)))

e <- ggplot() +
  geom_point(data = bg, aes(x, y), colour = "grey85", size = 0.2) +
  geom_point(data = hit, aes(x, y, colour = cluster), size = 0.4) +
  geom_path(data = circles, aes(x, y, group = id), colour = "red", linewidth = 0.6) +
  coord_fixed() + theme_void()

e

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/afterdbscan_subsetted_check_radius40.png"),   e, width = 15, height = 12, dpi = 300, bg = "white")


### Do co-localisation of entire array minus the dbscanned regions, i.e. not granulomas



map <- df %>%
  dplyr::select(x, y, cell_type, weight) %>%
  pivot_wider(names_from = cell_type, values_from = weight, values_fill = 0) %>%
  as.data.frame()

dim(map)         # spots x (2 + number of cell types)
colnames(map) 

region_spots <- hit %>% distinct(x, y)
bg_wide <- map %>% anti_join(region_spots, by = c("x", "y"))

nrow(map); nrow(region_spots); nrow(bg_wide)   # the last should equal the first minus the second

# co-localisation on the rest of the array
outs_background <- coloc_analysis_replicates_fast(list(as.data.frame(bg_wide)))
f <- outs_background$plot
f
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/bakcground_minus_dbscan_rgeions_coloc40um.png"),   f, width = 5, height = 5, dpi = 300, bg = "white")
write.csv(outs_background$rho_mat, "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/rho_background40um.csv")

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/bakcground_minus_dbscan_rgeions_coloc100um.png"),   f, width = 5, height = 5, dpi = 300, bg = "white")
write.csv(outs_background$rho_mat, "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/rho_background100um.csv")



# Or check the enriuchment of cells in grnauloma vs non? Co-localisation might not be greta?

### Assign 11 random regions in the rest of the array, as comparator


library(dplyr); library(dbscan); library(ggplot2)

## ---- 0) setup ---------------------------------------------------------
cell_cols    <- setdiff(colnames(map), c("x", "y"))
region_spots <- hit %>% distinct(x, y)
bg_wide      <- map %>% anti_join(region_spots, by = c("x", "y"))

## ---- 1) window size matched to your granulomas -----------------------
hull_area <- function(d) {
  h <- d[chull(d$x, d$y), c("x", "y")]
  abs(sum(h$x * c(h$y[-1], h$y[1]) - c(h$x[-1], h$x[1]) * h$y)) / 2
}
R <- 300

## ---- 2) settings you can tweak ---------------------------------------
n_regions <- 12
min_spots <- 50        # reject windows over sparse tissue or the edge
min_sep   <- 6 * R     # minimum distance between window centres (raise for more spread)
max_tries <- 10000
set.seed(42)

## ---- 3) candidate centres: background spots at least R from any granuloma
d_gran <- dbscan::kNN(as.matrix(region_spots[, c("x", "y")]), k = 1,
                      query = as.matrix(bg_wide[, c("x", "y")]))$dist[, 1]
cand <- bg_wide[d_gran > R, ]

## ---- 4) pick dispersed random windows --------------------------------
centres   <- data.frame(x = numeric(), y = numeric())
rand_list <- list()
tries     <- 0

while (length(rand_list) < n_regions && tries < max_tries) {
  tries <- tries + 1
  cx <- cand[sample(nrow(cand), 1), ]
  
  # enforce the minimum gap between windows
  if (nrow(centres) > 0 &&
      any(sqrt((centres$x - cx$x)^2 + (centres$y - cx$y)^2) < min_sep)) next
  
  in_win <- sqrt((bg_wide$x - cx$x)^2 + (bg_wide$y - cx$y)^2) <= R
  if (sum(in_win) < min_spots) next
  
  centres <- rbind(centres, data.frame(x = cx$x, y = cx$y))
  rand_list[[length(rand_list) + 1]] <- as.data.frame(bg_wide[in_win, c("x", "y", cell_cols)])
}

if (length(rand_list) < n_regions)
  warning("Only found ", length(rand_list), " windows. Lower min_sep or R, or lower min_spots.")

names(rand_list) <- paste0("rand_", seq_along(rand_list))
sapply(rand_list, nrow)      # compare with sapply(region_list, nrow)

## ---- 5) check where the windows landed -------------------------------
circle <- function(cx, cy, r, id, n = 100) {
  t <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(t), y = cy + r * sin(t), id = id)
}
circ <- bind_rows(lapply(seq_len(nrow(centres)),
                         function(i) circle(centres$x[i], centres$y[i], R, i)))

f <- ggplot() +
  geom_point(data = distinct(map, x, y), aes(x, y), colour = "grey88", size = 0.4) +
  geom_point(data = hit, aes(x, y), colour = "firebrick", size = 0.4) +   # granulomas
  geom_path(data = circ, aes(x, y, group = id), colour = "black") +       # random windows
  coord_fixed() + theme_void()

f

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/background_regions/randomwindows.png"),   f, width = 12, height = 12, dpi = 300, bg = "white")

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/background_regions/randomwindows.pdf"),   f, width = 6, height = 6, dpi = 300, bg = "white")


outs_rand <- lapply(rand_list, function(d)
  coloc_analysis_replicates_fast(list(d), radius = 40, nbins = 2000))


f <- outs_rand$plot
f

names(outs_rand)                 # e.g. "rand_1" ... "rand_11"

f <- outs_rand[["rand_1"]]$plot  # or outs_rand[[1]]$plot
f

for (nm in names(outs_rand)) {
  ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/background_regions/", paste0("coloc_", nm, ".pdf")),
         outs_rand[[nm]]$plot + ggtitle(nm), width = 6, height = 5)
  write.csv(outs_rand[[nm]]$rho_mat, file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/background_regions/", paste0("rho_", nm, ".csv")))
}





# one pooled heatmap across all 11 windows
outs_rand_pooled <- coloc_analysis_replicates_fast(unname(rand_list), radius = 40, nbins = 2000)
f <- outs_rand_pooled$plot
f
ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/background_regions/background_minus_dbscan_rgeions_coloc40um.png"),   f, width = 5, height = 5, dpi = 300, bg = "white")





# Enrichment of clel type sin db scan vs non dbscan - analysis below is done for eahc spot in or out of dbscanned rgeions, not done on a by region basis

library(dplyr); library(tidyr); library(ggplot2)

cell_cols    <- setdiff(colnames(map), c("x", "y"))
region_spots <- hit %>% distinct(x, y) %>% mutate(inside = TRUE)

dat <- map %>%
  left_join(region_spots, by = c("x", "y")) %>%
  mutate(inside = ifelse(is.na(inside), FALSE, TRUE))

table(dat$inside)    # spots outside vs inside

# Enrichment per celltype

enrich <- lapply(cell_cols, function(ct) {
  a <- dat[[ct]][dat$inside]
  b <- dat[[ct]][!dat$inside]
  w <- wilcox.test(a, b)
  data.frame(
    celltype  = ct,
    mean_in   = mean(a),
    mean_out  = mean(b),
    log2FC    = log2((mean(a) + 1e-3) / (mean(b) + 1e-3)),   # small pseudocount avoids dividing by 0
    auc       = unname(w$statistic) / (length(a) * length(b)), # 0.5 = no difference, >0.5 = higher inside
    p         = w$p.value
  )
}) %>% bind_rows() %>%
  mutate(padj = p.adjust(p, method = "BH")) %>%
  arrange(desc(log2FC))

enrich

write.csv(enrich, "/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/enrichment/enrichment_of_celltypes.csv")


# plot

ggplot(enrich, aes(reorder(celltype, log2FC), log2FC, fill = log2FC > 0)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("royalblue3", "red3"), guide = "none") +
  labs(x = NULL, y = "log2 fold change (inside vs outside regions)") +
  theme_minimal()

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/enrichment/foldchange.pdf"), width = 5, height = 5, dpi = 300, bg = "white")

# distributions
dat %>%
  pivot_longer(all_of(cell_cols), names_to = "celltype", values_to = "weight") %>%
  ggplot(aes(celltype, weight, fill = inside)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("grey80", "firebrick"), labels = c("outside", "inside"), name = NULL) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/enrichment/distribution.png"), width = 5, height = 5, dpi = 300, bg = "white")


# one regions

reg_levels <- sort(unique(hit$cluster))
outside    <- dat[!dat$inside, ]

enrich_reg <- lapply(reg_levels, function(cl) {
  inreg <- map %>% semi_join(filter(hit, cluster == cl), by = c("x", "y"))
  lapply(cell_cols, function(ct) {
    a <- inreg[[ct]]; b <- outside[[ct]]
    data.frame(region = cl, celltype = ct,
               log2FC = log2((mean(a) + 1e-3) / (mean(b) + 1e-3)),
               auc    = unname(wilcox.test(a, b)$statistic) / (length(a) * length(b)))
  }) %>% bind_rows()
}) %>% bind_rows()


ct_order  <- c("MNP2", "infected CD4", " naive CD4", "MNP1", "uninfected MNP", "hepatocyte")      # your order, spelled exactly as in the data
reg_order <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")                              # region labels, as character

enrich_reg <- enrich_reg %>%
  mutate(celltype = factor(celltype, levels = ct_order),
         region   = factor(region,   levels = reg_order))

ggplot(enrich_reg, aes(celltype, region, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = "royalblue3", mid = "white", high = "red3", midpoint = 0) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(file.path("/data/projects/punim1236/liver_stomicsv1.3/OllieRCTDs/MNP1_2_Tcellcombined_no_sinusoids_fullmode_RCTDs/DBscan_on_module5_hotspot/enrichment/enrichment_foldchnage_heatmap.pdf"), width = 5, height = 5, dpi = 300, bg = "white")

