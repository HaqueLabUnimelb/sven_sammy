"""Compare the archived Hotspot run with independent replicate runs in the new environment.

Hotspot's neighbour graph is random (pynndescent without a seed), so no run reproduces
another one exactly. The question is whether the archived run looks like one more
replicate: are its differences from the new runs as large as the differences between new
runs, and does it fall inside their spread?

Usage:
  python compare_replicates.py --archived ARCHDIR --notebook ARCHIVED_NOTEBOOK \
      --untitled3 NOTEBOOK_N30 --auto30 CSV... --auto15 CSV... --full DIR... --json OUT.json
  --auto30/--auto15: hotspot_autocorrelation.csv files of replicate runs
  --full: output folders of full replicate runs (each with hotspot_results/)
"""
import argparse
import itertools
import json
import os
import re

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

ap = argparse.ArgumentParser()
ap.add_argument("--archived", required=True)
ap.add_argument("--notebook", required=True)
ap.add_argument("--untitled3", required=True)
ap.add_argument("--auto30", nargs="+", required=True)
ap.add_argument("--auto15", nargs="+", required=True)
ap.add_argument("--full", nargs="+", required=True)
ap.add_argument("--json", required=True)
a = ap.parse_args()
R = {}


def say(k, v):
    R[k] = v
    print(f"{k}: {v}", flush=True)


def summary(values):
    v = np.asarray(values, dtype=float)
    return {"min": float(v.min()), "median": float(np.median(v)), "max": float(v.max()), "n": int(v.size)}


def printed_top20(nb_path, marker):
    nb = json.load(open(nb_path, encoding="utf-8"))
    for cell in nb["cells"]:
        if cell["cell_type"] == "code" and marker in "".join(cell["source"]):
            for o in cell.get("outputs", []):
                txt = "".join(o.get("data", {}).get("text/plain", ""))
                rows = re.findall(r"^(\S+)\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.e]+)\s+([-0-9.e]+)$", txt, flags=re.M)
                if rows:
                    return pd.DataFrame([(g, float(c), float(z)) for g, c, z, _, _ in rows], columns=["gene", "C", "Z"]).set_index("gene")
    return None


def autocorr_block(label, files, archived):
    reps = [pd.read_csv(f, index_col=0) for f in files]
    genes = reps[0].index
    for r in reps:
        assert set(r.index) == set(genes)
    K = len(reps)
    say(f"{label}.n_replicates", K)
    for col in ("C", "Z"):
        M = np.column_stack([r.loc[genes, col].values for r in reps])
        pair = [float(np.max(np.abs(M[:, i] - M[:, j]))) for i, j in itertools.combinations(range(K), 2)]
        say(f"{label}.{col}.max_abs_diff_between_replicates", summary(pair))
        if archived is not None:
            av = archived.loc[genes, col].values
            arc = [float(np.max(np.abs(av - M[:, i]))) for i in range(K)]
            say(f"{label}.{col}.max_abs_diff_archived_vs_replicate", summary(arc))
            below = (av < M.min(axis=1)).sum()
            above = (av > M.max(axis=1)).sum()
            say(f"{label}.{col}.genes_archived_outside_replicate_range",
                {"below_min": int(below), "above_max": int(above), "of": int(len(genes)),
                 "fraction": float((below + above) / len(genes)),
                 "expected_fraction_if_archived_is_one_more_replicate": float(2 / (K + 1))})
            # one replicate held out in the same way, for comparison
            held = []
            for i in range(K):
                rest = np.delete(M, i, axis=1)
                held.append(float(((M[:, i] < rest.min(axis=1)) | (M[:, i] > rest.max(axis=1))).mean()))
            say(f"{label}.{col}.fraction_outside_range_for_each_held_out_replicate_vs_other_{K - 1}", summary(held))
            mean = M.mean(axis=1)
            sd = M.std(axis=1, ddof=1)
            z = (av - mean) / np.where(sd > 0, sd, np.nan)
            say(f"{label}.{col}.archived_minus_replicate_mean_in_replicate_sd", {
                "median_abs": float(np.nanmedian(np.abs(z))), "max_abs": float(np.nanmax(np.abs(z))),
                "mean_signed": float(np.nanmean(z))})
    return reps


def top500(r):
    return list(r.loc[r.FDR < 0.05].sort_values("Z", ascending=False).head(500).index)


arc_res = os.path.join(a.archived, "hotspot_results")
arch30 = pd.read_csv(os.path.join(arc_res, "hotspot_autocorrelation.csv"), index_col=0)

# ---- autocorrelation, n_neighbors=30, against the archived table ----------------
reps30 = autocorr_block("n30", a.auto30, arch30)
t_arc = set(top500(arch30))
t_rep = [set(top500(r)) for r in reps30]
say("n30.top500_overlap_archived_vs_replicate", summary([len(t_arc & t) for t in t_rep]))
say("n30.top500_overlap_between_replicates", summary([len(x & y) for x, y in itertools.combinations(t_rep, 2)]))

# printed top-20 tables in the notebooks (6 decimals): inside the replicate range?
def printed_check(label, printed, reps):
    if printed is None:
        say(f"{label}.printed_rows", 0)
        return
    out = {}
    for col in ("C", "Z"):
        M = np.column_stack([r.loc[printed.index, col].round(6).values for r in reps])
        p = printed[col].values
        inside = ((p >= M.min(axis=1)) & (p <= M.max(axis=1)))
        out[col] = {"rows": int(len(p)), "inside_replicate_range": int(inside.sum()),
                    "max_abs_diff_to_replicate_mean": float(np.max(np.abs(p - M.mean(axis=1))))}
    say(f"{label}.printed_top20_vs_replicates", out)


printed_check("n30.untitled3", printed_top20(a.untitled3, "n_neighbors=30"), reps30)

# ---- autocorrelation, n_neighbors=15, against the printed notebook values -------
reps15 = autocorr_block("n15", a.auto15, None)
p15 = printed_top20(a.notebook, "n_neighbors=15")
printed_check("n15.archived_notebook", p15, reps15)
g = np.array([[r.loc["Gm42418", "C"], r.loc["Gm42418", "Z"]] for r in reps15])
say("n15.Gm42418_replicates_C", summary(g[:, 0]))
say("n15.Gm42418_replicates_Z", summary(g[:, 1]))
say("n15.Gm42418_printed", {"C": 0.477215, "Z": 813.077882})

# ---- full runs: modules and per-bin scores ----------------------------------------
runs = {"archived": a.archived}
for i, d in enumerate(a.full, 1):
    runs[f"rep{i}"] = d
say("full.runs", runs)
mods, spots = {}, {}
for name, d in runs.items():
    res = os.path.join(d, "hotspot_results")
    mods[name] = pd.read_csv(os.path.join(res, "hotspot_modules.csv"), index_col=0).iloc[:, 0]
    spots[name] = pd.read_csv(os.path.join(res, "hotspot_spot_scores.csv"), index_col=0)
say("full.module_sizes", {n: {int(k): int(v) for k, v in m.value_counts().sort_index().items()} for n, m in mods.items()})
marker = ["Tmsb4x", "Lyz2", "H2-Aa", "Cd74", "Alb", "Glul", "Gm42418", "Apoa2", "Gm26917"]
say("full.module_of_marker_genes", {n: {g: (int(m[g]) if g in m.index else None) for g in marker} for n, m in mods.items()})


def pair_stats(x, y):
    mx, my = mods[x], mods[y]
    common = mx.index.intersection(my.index)
    same = int((mx[common] == my[common]).sum())
    ari = float(adjusted_rand_score(mx[common], my[common]))
    m5x = set(mx[mx == 5].index)
    m5y = set(my[my == 5].index)
    sx, sy = spots[x], spots[y]
    assert list(sx.index) == list(sy.index)
    r = {c: float(np.corrcoef(sx[c], sy[c])[0, 1]) for c in sx.columns if c.startswith("module_")}
    d = {c: float(np.max(np.abs(sx[c] - sy[c]))) for c in sx.columns if c.startswith("module_")}
    bx = set(sx.index[sx["module_5"] > 4.0])
    by = set(sy.index[sy["module_5"] > 4.0])
    return {"genes_in_both_top500": int(len(common)), "same_module_label": same, "adjusted_rand_index": ari,
            "module5_genes": [len(m5x), len(m5y)], "module5_jaccard": len(m5x & m5y) / len(m5x | m5y),
            "spot_score_pearson_min_over_modules": min(r.values()), "spot_score_pearson": r,
            "spot_score_max_abs_diff": d,
            "module5_gt4_bins": [len(bx), len(by)], "module5_gt4_bins_common": len(bx & by),
            "module5_gt4_jaccard": len(bx & by) / len(bx | by)}


names = list(runs)
pairs = {}
for x, y in itertools.combinations(names, 2):
    pairs[f"{x}~{y}"] = pair_stats(x, y)
    print(f"{x}~{y}: " + json.dumps({k: v for k, v in pairs[f'{x}~{y}'].items() if k not in ('spot_score_pearson', 'spot_score_max_abs_diff')}))
R["full.pairs"] = pairs
for key in ("same_module_label", "adjusted_rand_index", "module5_jaccard", "spot_score_pearson_min_over_modules", "module5_gt4_jaccard"):
    arc = [v[key] for k, v in pairs.items() if k.startswith("archived~")]
    rep = [v[key] for k, v in pairs.items() if not k.startswith("archived~")]
    say(f"full.{key}.archived_vs_replicates", summary(arc))
    say(f"full.{key}.between_replicates", summary(rep))
for key in ("spot_score_max_abs_diff",):
    arc = [max(v[key].values()) for k, v in pairs.items() if k.startswith("archived~")]
    rep = [max(v[key].values()) for k, v in pairs.items() if not k.startswith("archived~")]
    say(f"full.{key}_all_modules.archived_vs_replicates", summary(arc))
    say(f"full.{key}_all_modules.between_replicates", summary(rep))

with open(a.json, "w") as fh:
    json.dump(R, fh, indent=1)
