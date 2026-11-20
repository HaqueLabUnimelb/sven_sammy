"""Show where Hotspot's run-to-run differences come from.

hotspot 1.1.3, Hotspot.create_knn_graph(), default approx_neighbors=True, calls
    pynndescent.NNDescent(coords, n_neighbors=n_neighbors + 1)
without a random_state, so the neighbour graph depends on numpy's global random state.
This script builds the graph of the notebook (n_neighbors=30 and 15) twice in one process
and counts the bins whose neighbour sets differ, in one of three modes:
  unseeded   as in the notebook
  seeded     numpy.random.seed(0) before each of the two builds
It also builds the exact graph (scikit-learn ball tree, the approx_neighbors=False path)
and counts the bins whose NNDescent neighbour set differs from the exact one, and among
those, the bins where the k-th neighbour distance differs (that is, where NNDescent did
not return k nearest neighbours; the other bins differ only in how ties at the k-th
distance were broken).

Usage: python hotspot_knn_check.py IN_H5AD MODE    (run with NUMBA_NUM_THREADS set or not)
"""
import os
import sys
import time

import numpy as np
import pandas as pd
import anndata as ad
from pynndescent import NNDescent
from sklearn.neighbors import NearestNeighbors

in_h5ad, mode = sys.argv[1], sys.argv[2]
adata = ad.read_h5ad(in_h5ad, backed="r")
coords = pd.DataFrame(adata.obs[["x", "y"]].to_numpy(), index=adata.obs_names).values
threads = os.environ.get("NUMBA_NUM_THREADS", "default")


def nnd(k):
    if mode == "seeded":
        np.random.seed(0)
    t = time.time()
    ind, dist = NNDescent(coords, n_neighbors=k + 1).neighbor_graph
    return ind[:, 1:], dist[:, 1:], time.time() - t


def n_rows_differ(a, b):
    return int((np.sort(a, axis=1) != np.sort(b, axis=1)).any(axis=1).sum())


for k in (30, 15):
    i1, d1, t1 = nnd(k)
    i2, d2, t2 = nnd(k)
    ex_d, ex_i = NearestNeighbors(n_neighbors=k, algorithm="ball_tree").fit(coords).kneighbors()
    kth_ok = np.isclose(d1.max(axis=1), ex_d.max(axis=1))
    print(f"KNNCHECK mode={mode} numba_threads={threads} n_neighbors={k} bins={coords.shape[0]} "
          f"rows_differ_between_two_NNDescent_builds={n_rows_differ(i1, i2)} "
          f"rows_differ_NNDescent_vs_exact={n_rows_differ(i1, ex_i)} "
          f"rows_with_wrong_kth_distance={int((~kth_ok).sum())} "
          f"build_seconds={t1:.1f},{t2:.1f}", flush=True)
