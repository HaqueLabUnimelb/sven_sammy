"""Compare the h5ad written by the test run of 01_Y01453J8_bin20_gef_to_h5ad.py with the
archived one.

Part 1 checks the fields named in the task: X (shape, data, indices, indptr), obs_names,
var_names, obs x, y and nCount_RNA. Part 2 walks both HDF5 files and compares every
dataset (dtype, shape, values) and every attribute.

Usage: python compare_h5ad.py NEW.h5ad ARCHIVED.h5ad OUT.json
"""
import hashlib
import json
import sys

import h5py
import numpy as np
import anndata as ad

new_p, old_p, out_json = sys.argv[1:4]
R = {}


def say(k, v):
    R[k] = v
    print(f"{k}: {v}", flush=True)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


say("anndata_version", ad.__version__)
say("md5_new", md5(new_p))
say("md5_archived", md5(old_p))
say("file_bytes_identical", R["md5_new"] == R["md5_archived"])

# ---- part 1: named fields --------------------------------------------------------
n = ad.read_h5ad(new_p)
o = ad.read_h5ad(old_p)
say("shape_new", list(n.shape))
say("shape_archived", list(o.shape))
say("X_type_new", type(n.X).__name__)
say("X_type_archived", type(o.X).__name__)
say("X_dtype_new", str(n.X.dtype))
say("X_dtype_archived", str(o.X.dtype))
say("X_nnz_new", int(n.X.nnz))
say("X_nnz_archived", int(o.X.nnz))
say("X_data_identical", bool(np.array_equal(n.X.data, o.X.data)))
say("X_indices_identical", bool(np.array_equal(n.X.indices, o.X.indices)))
say("X_indptr_identical", bool(np.array_equal(n.X.indptr, o.X.indptr)))
say("obs_names_identical", bool(list(n.obs_names) == list(o.obs_names)))
say("var_names_identical", bool(list(n.var_names) == list(o.var_names)))
for c in ["x", "y", "nCount_RNA", "nFeature_RNA"]:
    same = bool(np.array_equal(n.obs[c].values, o.obs[c].values)) and str(n.obs[c].dtype) == str(o.obs[c].dtype)
    say(f"obs_{c}_identical", same)
    say(f"obs_{c}_dtype", [str(n.obs[c].dtype), str(o.obs[c].dtype)])
say("obs_columns", [list(n.obs.columns), list(o.obs.columns)])
say("var_columns", [list(n.var.columns), list(o.var.columns)])
say("obs_index_name", [n.obs.index.name, o.obs.index.name])
say("var_index_name", [n.var.index.name, o.var.index.name])
del n, o

# ---- part 2: every dataset and attribute -----------------------------------------
diffs = []
storage = []
n_ds = 0


def attrs_equal(x, y):
    if set(x.attrs.keys()) != set(y.attrs.keys()):
        return False
    for k in x.attrs.keys():
        vx, vy = x.attrs[k], y.attrs[k]
        if isinstance(vx, np.ndarray) or isinstance(vy, np.ndarray):
            if not np.array_equal(np.asarray(vx), np.asarray(vy)):
                return False
        elif vx != vy:
            return False
    return True


with h5py.File(new_p, "r") as fn, h5py.File(old_p, "r") as fo:
    names_n, names_o = [], []
    fn.visit(names_n.append)
    fo.visit(names_o.append)
    if sorted(names_n) != sorted(names_o):
        diffs.append(f"object names differ: only new {sorted(set(names_n) - set(names_o))}, only archived {sorted(set(names_o) - set(names_n))}")
    if not attrs_equal(fn, fo):
        diffs.append("root attributes differ")
    for name in sorted(set(names_n) & set(names_o)):
        x, y = fn[name], fo[name]
        if not attrs_equal(x, y):
            diffs.append(f"{name}: attributes differ: {dict(x.attrs)} vs {dict(y.attrs)}")
        if isinstance(x, h5py.Dataset):
            n_ds += 1
            if x.dtype != y.dtype or x.shape != y.shape:
                diffs.append(f"{name}: dtype/shape {x.dtype}{x.shape} vs {y.dtype}{y.shape}")
                continue
            if not np.array_equal(x[()], y[()]):
                diffs.append(f"{name}: values differ")
            if x.compression != y.compression or x.chunks != y.chunks:
                storage.append(f"{name}: compression {x.compression} vs {y.compression}, chunks {x.chunks} vs {y.chunks}")
say("hdf5_datasets_compared", n_ds)
say("hdf5_differences", diffs)
say("hdf5_storage_layout_differences", storage)
say("hdf5_all_datasets_and_attributes_identical", len(diffs) == 0)

with open(out_json, "w") as fh:
    json.dump(R, fh, indent=1)
