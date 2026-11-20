"""Test a: run the Hotspot steps of hotspot_liverstomics_NI.ipynb as a script.

The code below is the code of the notebook's cells, in top-to-bottom order, without the
cells and lines that only draw figures, and with two changes:
  * n_neighbors=30 (the value of the run that wrote the archived tables) instead of 15;
  * jobs=JOBS instead of jobs=8 (number of worker processes only; every gene and every
    gene pair is computed independently, so the value does not change the results).
It writes the same tables, with the same code, under OUTDIR/hotspot_results and
OUTDIR/hotspot_plots. adata_with_hotspot.h5ad (1.3 GB) is not written.

It then repeats only the autocorrelation step with n_neighbors=15, the value shown in the
archived notebook, and writes OUTDIR/n15/hotspot_autocorrelation.csv.

If the environment variable HOTSPOT_TEST_SEED is set, numpy.random.seed(HOTSPOT_TEST_SEED)
is called first (seeded test: the neighbour graph of Hotspot depends on numpy's global
random state, see hotspot_knn_check.py). The archived notebook sets no seed.

Usage: python hotspot_test_a.py IN_H5AD OUTDIR JOBS
"""
import os
import re
import sys
import time

import numpy as np

if os.environ.get("HOTSPOT_TEST_SEED"):
    np.random.seed(int(os.environ["HOTSPOT_TEST_SEED"]))
    print("numpy.random.seed", os.environ["HOTSPOT_TEST_SEED"])
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import hotspot

IN_H5AD, OUTDIR, JOBS = sys.argv[1], sys.argv[2], int(sys.argv[3])
os.makedirs(OUTDIR, exist_ok=True)
t0 = time.time()


def stamp(msg):
    print(f"[{time.time() - t0:8.1f} s] {msg}", flush=True)


for mod in ("numpy", "pandas", "scipy", "sklearn", "numba", "anndata", "scanpy", "hotspot", "statsmodels", "h5py"):
    m = __import__(mod)
    print(mod, getattr(m, "__version__", "?"), getattr(m, "__file__", "?"))

# ---- notebook cell 2 -------------------------------------------------------------
adata = sc.read_h5ad(IN_H5AD)
adata.layers["counts"] = adata.X.copy()
adata.obsm["spatial"] = adata.obs[["x", "y"]].to_numpy()
adata.obs["total_counts"] = np.asarray(adata.layers["counts"].sum(1)).ravel()

# keep expressed genes only; this matters a lot for speed
sc.pp.filter_genes(adata, min_cells=int(0.01 * adata.n_obs))
stamp(f"read and filtered: {adata.n_obs} bins x {adata.n_vars} genes")

# ---- notebook cell 4 (n_neighbors=30) --------------------------------------------
hs = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model="danb",
    latent_obsm_key="spatial",
    umi_counts_obs_key="total_counts",
)

hs.create_knn_graph(weighted_graph=False, n_neighbors=30)

hs_results = hs.compute_autocorrelations(jobs=JOBS)
print(hs_results.sort_values("Z", ascending=False).head(20))
stamp("autocorrelations done (n_neighbors=30)")

# ---- notebook cell 5 -------------------------------------------------------------
top_genes = hs_results.loc[hs_results.FDR < 0.05].sort_values("Z", ascending=False).head(500).index

hs.compute_local_correlations(top_genes, jobs=JOBS)
stamp("local correlations done")

modules = hs.create_modules(min_gene_threshold=20, core_only=False, fdr_threshold=0.05)
print(modules.value_counts())

module_scores = hs.calculate_module_scores()
adata.obs[module_scores.columns] = module_scores
stamp("module scores done")

# ---- notebook cell 7 (without the plotting loop) ---------------------------------
adata.obs = adata.obs.drop(columns=[c for c in adata.obs.columns if not isinstance(c, str)])
module_scores.columns = [f"module_{m}" for m in module_scores.columns]
adata.obs[module_scores.columns] = module_scores

# ---- notebook cell 9: os.chdir to the output folder ------------------------------
os.chdir(OUTDIR)

# ---- notebook cell 10 (without the plotting loop) --------------------------------
outdir = "hotspot_plots"
os.makedirs(outdir, exist_ok=True)
adata.obs = adata.obs.drop(columns=[c for c in adata.obs.columns if not isinstance(c, str)])
module_scores.columns = [f"module_{m}" for m in module_scores.columns]
adata.obs[module_scores.columns] = module_scores

# ---- notebook cell 12 ------------------------------------------------------------
outdir = "hotspot_plots"
os.makedirs(outdir, exist_ok=True)

mod_df = pd.DataFrame({"gene": modules.index, "module": modules.values})
mod_df = mod_df[mod_df["module"] != -1]
mod_df = mod_df.merge(hs_results[["C", "Z", "FDR"]],
                      left_on="gene", right_index=True, how="left")
mod_df = mod_df.sort_values(["module", "Z"], ascending=[True, False])
mod_df.to_csv(os.path.join(outdir, "hotspot_module_genes.csv"), index=False)

depth = adata.obs["total_counts"].values
rows = []

for m in sorted(mod_df["module"].unique()):
    genes = mod_df.loc[mod_df["module"] == m, "gene"].tolist()

    X = adata[:, genes].layers["counts"]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    X = np.log1p(X / depth[:, None] * 1e4).astype("float32")

    s = adata.obs[f"module_{m}"].values.astype("float32")

    Xz = (X - X.mean(0)) / X.std(0)
    sz = (s - s.mean()) / s.std()
    r = Xz.T @ sz / len(s)

    rows.append(pd.DataFrame({"gene": genes, "module": m, "corr_with_module": r}))

mod_df2 = mod_df.merge(pd.concat(rows), on=["gene", "module"]).sort_values(
    ["module", "corr_with_module"], ascending=[True, False])

mod_df2.to_csv(os.path.join(outdir, "hotspot_module_genes_ranked.csv"), index=False)
mod_df2.groupby("module").head(20).to_csv(
    os.path.join(outdir, "hotspot_module_top20_ranked.csv"), index=False)

for m, g in mod_df2.groupby("module"):
    g.to_csv(os.path.join(outdir, f"module_{m}_genes.csv"), index=False)
    g["gene"].to_csv(os.path.join(outdir, f"module_{m}_genelist.txt"),
                     index=False, header=False)

print(mod_df2.groupby("module").head(10).groupby("module")["gene"].apply(list))

# ---- notebook cell 13 ------------------------------------------------------------
outdir = "hotspot_results"
os.makedirs(outdir, exist_ok=True)

hs_results.to_csv(os.path.join(outdir, "hotspot_autocorrelation.csv"))
modules.to_csv(os.path.join(outdir, "hotspot_modules.csv"))
hs.local_correlation_z.to_csv(os.path.join(outdir, "hotspot_local_correlation_z.csv"))
module_scores.to_csv(os.path.join(outdir, "hotspot_module_scores.csv"))

# ---- notebook cell 14: writes adata_with_hotspot.h5ad; not repeated here ---------

# ---- notebook cell 16 ------------------------------------------------------------
module_scores.columns = ["module_" + re.sub(r"^(module_)+", "", str(c)) for c in module_scores.columns]
old = [c for c in adata.obs.columns if re.fullmatch(r"(module_)+\d+", str(c))]
adata.obs = adata.obs.drop(columns=old)
adata.obs[module_scores.columns] = module_scores
print([c for c in adata.obs.columns if "module" in c])

# ---- notebook cell 17 ------------------------------------------------------------
outdir = "hotspot_results"
os.makedirs(outdir, exist_ok=True)

mod_cols = [c for c in adata.obs.columns if c.startswith("module_")]

out = pd.DataFrame(adata.obsm["spatial"], index=adata.obs_names.astype(str), columns=["x", "y"])
out = pd.concat([out, adata.obs.loc[:, mod_cols]], axis=1)
out.index.name = "barcode"
out.to_csv(os.path.join(outdir, "hotspot_spot_scores.csv"))

gene_tab = pd.DataFrame({"gene": modules.index, "module": modules.values}).merge(
    hs_results[["C", "Z", "FDR"]], left_on="gene", right_index=True, how="left")
gene_tab.to_csv(os.path.join(outdir, "hotspot_gene_modules.csv"), index=False)
stamp("n_neighbors=30 tables written")

# ---- autocorrelation only, n_neighbors=15 (the value in the archived notebook) ---
hs15 = hotspot.Hotspot(
    adata,
    layer_key="counts",
    model="danb",
    latent_obsm_key="spatial",
    umi_counts_obs_key="total_counts",
)
hs15.create_knn_graph(weighted_graph=False, n_neighbors=15)
hs15_results = hs15.compute_autocorrelations(jobs=JOBS)
print(hs15_results.sort_values("Z", ascending=False).head(20))
os.makedirs("n15", exist_ok=True)
hs15_results.to_csv(os.path.join("n15", "hotspot_autocorrelation.csv"))
stamp("n_neighbors=15 autocorrelations written")
