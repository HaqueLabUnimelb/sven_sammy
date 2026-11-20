#!/usr/bin/env python3
"""Write cellchat_r.yml from the conda-meta folder of the conda environment that ran CellChat.

usage: python3 make_conda_yml.py /data/gpfs/projects/punim2394/conda/rctd_env > ../cellchat_r.yml
"""
import glob, json, os, sys

E = sys.argv[1]
rows = []
for f in sorted(glob.glob(os.path.join(E, "conda-meta", "*.json"))):
    d = json.load(open(f))
    rows.append((d["name"], d["version"], d["build"], d.get("channel", "")))
chans = sorted(set(r[3] for r in rows))
print("# Conda environment that provides R 4.5.1 for the CellChat step (part 2).")
print("#")
print("# Lines 1 to %d of the dependency list are the exact packages (name=version=build) of the conda" % len(rows))
print("# environment that ran the CellChat scripts on 23 October 2026 (/data/gpfs/projects/punim2394/conda/rctd_env,")
print("# created with `conda create -n rctd_env r-base=4.5.1 -c conda-forge`; read from its conda-meta folder).")
print("# Two of them are changed: sysroot_linux-64 2.12 and kernel-headers_linux-64 2.6.32 (C library")
print("# headers of CentOS 6 used by the conda compilers) are replaced by 2.17 and 3.10, because the C++")
print("# library of gcc 15 does not compile against the 2.12 headers (error: 'timespec_get' has not been")
print("# declared, e.g. in Rcpp). They are used only when packages are compiled; R itself is unchanged.")
print("# The packages after them are not in that environment; they are added so that R packages can be")
print("# compiled from source inside this environment (headers and build tools).")
print("#")
print("# channels recorded in conda-meta: %s" % ", ".join(chans))
print("#")
print("# usage:  conda env create -p <prefix> -f cellchat_r.yml")
print("#         <prefix>/bin/Rscript install_cellchat_env.R <library_dir>")
print("name: cellchat_r451")
print("channels:")
print("  - conda-forge")
print("  - nodefaults")
print("dependencies:")
REPLACE = {"sysroot_linux-64": "2.17", "kernel-headers_linux-64": "3.10"}
for n, v, b, c in rows:
    if n in REPLACE:
        print("  - %s=%s   # original environment: %s=%s=%s" % (n, REPLACE[n], n, v, b))
    else:
        print("  - %s=%s=%s" % (n, v, b))
print("  # added for building R packages from source (not in the original environment)")
for extra in ["pkg-config", "cmake", "zlib", "libxml2", "glpk", "gmp", "nlopt"]:
    print("  - %s" % extra)
