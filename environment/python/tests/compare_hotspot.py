"""Compare Hotspot tables of a test run with the archived tables.

Usage:
  python compare_hotspot.py --new NEWDIR --archived ARCHDIR --repo-tables REPODIR \
      --notebook ARCHIVED_NOTEBOOK [--untitled3 NOTEBOOK_N30] --json OUT.json

NEWDIR holds hotspot_results/ and hotspot_plots/ of the test run (and n15/ for test a).
ARCHDIR is the folder of the archived run (with hotspot_results/ and hotspot_plots/).
REPODIR is 3_STOmics_liver_v1.3/4_Hotspot/result_tables of the repository.
"""
import argparse
import hashlib
import json
import os
import re

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--new", required=True)
ap.add_argument("--archived", required=True)
ap.add_argument("--repo-tables", required=True)
ap.add_argument("--notebook", required=True)
ap.add_argument("--untitled3", default=None)
ap.add_argument("--json", required=True)
a = ap.parse_args()

R = {}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def same_bytes(p, q):
    return os.path.exists(p) and os.path.exists(q) and md5(p) == md5(q)


def maxabs(x, y):
    return float(np.nanmax(np.abs(np.asarray(x, dtype=float) - np.asarray(y, dtype=float))))


def say(key, value):
    R[key] = value
    print(f"{key}: {value}")


new_res = os.path.join(a.new, "hotspot_results")
new_plt = os.path.join(a.new, "hotspot_plots")
arc_res = os.path.join(a.archived, "hotspot_results")
arc_plt = os.path.join(a.archived, "hotspot_plots")

# ---- 1. per-gene autocorrelation -------------------------------------------------
n_ac = pd.read_csv(os.path.join(new_res, "hotspot_autocorrelation.csv"), index_col=0)
o_ac = pd.read_csv(os.path.join(arc_res, "hotspot_autocorrelation.csv"), index_col=0)
say("autocorr.n_genes_new", int(n_ac.shape[0]))
say("autocorr.n_genes_archived", int(o_ac.shape[0]))
say("autocorr.same_gene_order", bool(list(n_ac.index) == list(o_ac.index)))
say("autocorr.same_gene_set", bool(set(n_ac.index) == set(o_ac.index)))
oa = o_ac.reindex(n_ac.index)
for c in ["C", "Z", "Pval", "FDR"]:
    say(f"autocorr.max_abs_diff_{c}", maxabs(n_ac[c], oa[c]))
say("autocorr.max_rel_diff_Z", float(np.nanmax(np.abs(n_ac["Z"] - oa["Z"]) / np.abs(oa["Z"]).clip(lower=1e-300))))
say("autocorr.file_identical_to_archived",
    same_bytes(os.path.join(new_res, "hotspot_autocorrelation.csv"), os.path.join(arc_res, "hotspot_autocorrelation.csv")))
say("autocorr.file_identical_to_repo",
    same_bytes(os.path.join(new_res, "hotspot_autocorrelation.csv"), os.path.join(a.repo_tables, "hotspot_autocorrelation.csv")))

# ---- 2. the 500 genes ------------------------------------------------------------
def top500(r):
    return list(r.loc[r.FDR < 0.05].sort_values("Z", ascending=False).head(500).index)


t_new = top500(n_ac)
t_old = top500(o_ac)
o_lcz = pd.read_csv(os.path.join(arc_res, "hotspot_local_correlation_z.csv"), index_col=0)
n_lcz = pd.read_csv(os.path.join(new_res, "hotspot_local_correlation_z.csv"), index_col=0)
say("top500.n_new", len(t_new))
say("top500.same_set_as_archived_autocorr", bool(set(t_new) == set(t_old)))
say("top500.same_order_as_archived_autocorr", bool(t_new == t_old))
say("top500.same_order_as_archived_local_corr_matrix", bool(t_new == list(o_lcz.index)))
say("top500.n_in_common_with_archived", len(set(t_new) & set(o_lcz.index)))

# ---- 3. local correlation Z matrix -----------------------------------------------
say("localcorrz.same_labels", bool(list(n_lcz.index) == list(o_lcz.index) and list(n_lcz.columns) == list(o_lcz.columns)))
say("localcorrz.max_abs_diff", maxabs(n_lcz.values, o_lcz.reindex(index=n_lcz.index, columns=n_lcz.columns).values))
say("localcorrz.file_identical_to_archived",
    same_bytes(os.path.join(new_res, "hotspot_local_correlation_z.csv"), os.path.join(arc_res, "hotspot_local_correlation_z.csv")))

# ---- 4. module of every gene -----------------------------------------------------
n_mod = pd.read_csv(os.path.join(new_res, "hotspot_modules.csv"), index_col=0).iloc[:, 0]
o_mod = pd.read_csv(os.path.join(arc_res, "hotspot_modules.csv"), index_col=0).iloc[:, 0]
both = n_mod.index.intersection(o_mod.index)
say("modules.n_genes_new", int(n_mod.size))
say("modules.n_genes_archived", int(o_mod.size))
say("modules.n_genes_same_module", int((n_mod.reindex(both) == o_mod.reindex(both)).sum()))
say("modules.n_genes_different_module", int((n_mod.reindex(both) != o_mod.reindex(both)).sum()))
say("modules.sizes_new", {int(k): int(v) for k, v in n_mod.value_counts().sort_index().items()})
say("modules.sizes_archived", {int(k): int(v) for k, v in o_mod.value_counts().sort_index().items()})
say("modules.file_identical_to_archived",
    same_bytes(os.path.join(new_res, "hotspot_modules.csv"), os.path.join(arc_res, "hotspot_modules.csv")))
say("modules.file_identical_to_repo",
    same_bytes(os.path.join(new_res, "hotspot_modules.csv"), os.path.join(a.repo_tables, "hotspot_modules.csv")))

n_gm = pd.read_csv(os.path.join(new_res, "hotspot_gene_modules.csv"))
o_gm = pd.read_csv(os.path.join(arc_res, "hotspot_gene_modules.csv"))
say("gene_modules.same_gene_and_module_columns", bool(n_gm[["gene", "module"]].equals(o_gm[["gene", "module"]])))
say("gene_modules.file_identical_to_archived",
    same_bytes(os.path.join(new_res, "hotspot_gene_modules.csv"), os.path.join(arc_res, "hotspot_gene_modules.csv")))
say("gene_modules.file_identical_to_repo",
    same_bytes(os.path.join(new_res, "hotspot_gene_modules.csv"), os.path.join(a.repo_tables, "hotspot_gene_modules.csv")))

# ---- 5. per-bin module scores ----------------------------------------------------
n_sp = pd.read_csv(os.path.join(new_res, "hotspot_spot_scores.csv"), index_col=0)
o_sp = pd.read_csv(os.path.join(arc_res, "hotspot_spot_scores.csv"), index_col=0)
say("spot_scores.n_bins_new", int(n_sp.shape[0]))
say("spot_scores.n_bins_archived", int(o_sp.shape[0]))
say("spot_scores.same_barcodes_in_order", bool(list(n_sp.index) == list(o_sp.index)))
say("spot_scores.same_columns", bool(list(n_sp.columns) == list(o_sp.columns)))
say("spot_scores.xy_identical", bool(n_sp[["x", "y"]].equals(o_sp[["x", "y"]])))
mods = [c for c in n_sp.columns if c.startswith("module_")]
per = {c: maxabs(n_sp[c], o_sp[c]) for c in mods}
say("spot_scores.max_abs_diff_per_module", per)
say("spot_scores.max_abs_diff_all_modules", max(per.values()))
say("spot_scores.n_values_not_identical", int(sum((n_sp[c].values != o_sp[c].values).sum() for c in mods)))
say("spot_scores.module5_bins_above_4_new", int((n_sp["module_5"] > 4.0).sum()))
say("spot_scores.module5_bins_above_4_archived", int((o_sp["module_5"] > 4.0).sum()))
say("spot_scores.file_identical_to_archived",
    same_bytes(os.path.join(new_res, "hotspot_spot_scores.csv"), os.path.join(arc_res, "hotspot_spot_scores.csv")))

n_ms = pd.read_csv(os.path.join(new_res, "hotspot_module_scores.csv"), index_col=0)
o_ms = pd.read_csv(os.path.join(arc_res, "hotspot_module_scores.csv"), index_col=0)
say("module_scores_csv.header_new", list(n_ms.columns))
say("module_scores_csv.same_index", bool(list(n_ms.index) == list(o_ms.index)))
say("module_scores_csv.max_abs_diff_values", maxabs(n_ms.values, o_ms.values))

# ---- 6. ranked module gene tables (hotspot_plots) --------------------------------
n_rk = pd.read_csv(os.path.join(new_plt, "hotspot_module_genes_ranked.csv"))
o_rk = pd.read_csv(os.path.join(a.repo_tables, "hotspot_module_genes_ranked.csv"))
say("ranked.n_rows_new_archived", [int(n_rk.shape[0]), int(o_rk.shape[0])])
say("ranked.same_gene_module_order", bool(n_rk[["gene", "module"]].equals(o_rk[["gene", "module"]])))
mk = n_rk.merge(o_rk, on=["gene", "module"], suffixes=("_new", "_old"))
say("ranked.n_gene_module_pairs_in_both", int(mk.shape[0]))
say("ranked.max_abs_diff_corr_with_module_common", maxabs(mk["corr_with_module_new"], mk["corr_with_module_old"]) if len(mk) else None)
say("ranked.file_identical_to_repo",
    same_bytes(os.path.join(new_plt, "hotspot_module_genes_ranked.csv"), os.path.join(a.repo_tables, "hotspot_module_genes_ranked.csv")))
say("module5_genelist.identical_to_repo",
    same_bytes(os.path.join(new_plt, "module_5_genelist.txt"), os.path.join(a.repo_tables, "module_5_genelist.txt")))
plot_csvs = sorted(f for f in os.listdir(arc_plt) if f.endswith((".csv", ".txt")))
say("hotspot_plots.tables_identical_to_archived",
    {f: same_bytes(os.path.join(new_plt, f), os.path.join(arc_plt, f)) for f in plot_csvs})


# ---- 7. printed tables in the notebooks ------------------------------------------
def printed_top20(nb_path, marker):
    nb = json.load(open(nb_path, encoding="utf-8"))
    for cell in nb["cells"]:
        if cell["cell_type"] == "code" and marker in "".join(cell["source"]):
            for o in cell.get("outputs", []):
                txt = "".join(o.get("data", {}).get("text/plain", ""))
                rows = re.findall(r"^(\S+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.e]+)\s+([-0-9.e]+)$", txt, flags=re.M)
                if rows:
                    return rows
    return []


def check_printed(rows, table, label):
    ok = 0
    bad = []
    for g, c, z, _, _ in rows:
        if g in table.index and f"{table.loc[g, 'C']:.6f}" == c and f"{table.loc[g, 'Z']:.6f}" == z:
            ok += 1
        else:
            bad.append(g)
    say(f"{label}.rows_printed", len(rows))
    say(f"{label}.rows_matching_to_6_decimals", ok)
    if bad:
        say(f"{label}.rows_not_matching", bad)


n15p = os.path.join(a.new, "n15", "hotspot_autocorrelation.csv")
if os.path.exists(n15p):
    n15 = pd.read_csv(n15p, index_col=0)
    say("n15.Gm42418_C", float(n15.loc["Gm42418", "C"]))
    say("n15.Gm42418_Z", float(n15.loc["Gm42418", "Z"]))
    say("n15.Gm42418_matches_notebook_C_0.477215_Z_813.077882",
        f"{n15.loc['Gm42418', 'C']:.6f}" == "0.477215" and f"{n15.loc['Gm42418', 'Z']:.6f}" == "813.077882")
    check_printed(printed_top20(a.notebook, "n_neighbors=15"), n15, "n15.archived_notebook_top20")
if a.untitled3:
    check_printed(printed_top20(a.untitled3, "n_neighbors=30"), n_ac, "n30.untitled3_top20")

with open(a.json, "w") as fh:
    json.dump(R, fh, indent=1)
