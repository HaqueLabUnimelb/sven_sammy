"""One replicate of the notebook's autocorrelation step (cells 2 and 4).

Usage: python hotspot_autocorr_replicate.py IN_H5AD OUT_CSV N_NEIGHBORS JOBS
"""
import sys

import numpy as np
import scanpy as sc
import hotspot

in_h5ad, out_csv, k, jobs = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

adata = sc.read_h5ad(in_h5ad)
adata.layers["counts"] = adata.X.copy()
adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()
adata.obs["total_counts"] = np.asarray(adata.layers["counts"].sum(1)).ravel()
sc.pp.filter_genes(adata, min_cells=int(0.01 * adata.n_obs))

hs = hotspot.Hotspot(adata, layer_key="counts", model="danb",
                     latent_obsm_key="spatial", umi_counts_obs_key="total_counts")
hs.create_knn_graph(weighted_graph=False, n_neighbors=k)
hs_results = hs.compute_autocorrelations(jobs=jobs)
hs_results.to_csv(out_csv)
print(out_csv, hs_results.loc["Gm42418", ["C", "Z"]].tolist())
