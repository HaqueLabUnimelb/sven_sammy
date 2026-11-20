#!/usr/bin/env python3
"""Build a lock table (package, version, source, origin) from original R libraries.

usage: make_lock.py --roots <file> --libs <lib1,lib2,...> --labels <lab1,lab2,...> [--si <sessionInfo pkgs file>]
                    [--extra name=version=source=origin ...] > lock.tsv

roots: one package name per line (or name_version tokens).  The closure over Depends, Imports and
LinkingTo is resolved through the libraries in the given order (first hit wins, like .libPaths()).
Base packages (Priority: base) are dropped.  When --si is given, every package listed there must
resolve to the same version, else the script stops.
"""
import argparse, os, re, sys

BASE = {"base", "compiler", "datasets", "graphics", "grDevices", "grid", "methods", "parallel",
        "splines", "stats", "stats4", "tcltk", "tools", "utils", "R", "translations"}


def parse_dcf(path):
    d, key = {}, None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line[:1] in (" ", "\t") and key:
                d[key] += " " + line.strip()
            elif ":" in line:
                key, val = line.split(":", 1)
                key = key.strip()
                d[key] = val.strip()
    return d


def deps(d):
    out = []
    for fld in ("Depends", "Imports", "LinkingTo"):
        for tok in d.get(fld, "").split(","):
            name = re.sub(r"\(.*?\)", "", tok).strip()
            if name and name not in BASE:
                out.append(name)
    return out


def source_of(d):
    if d.get("RemoteType") == "github" or d.get("GithubRepo"):
        user = d.get("RemoteUsername", d.get("GithubUsername"))
        repo = d.get("RemoteRepo", d.get("GithubRepo"))
        sha = d.get("RemoteSha", d.get("GithubSHA1"))
        return "GitHub %s/%s@%s" % (user, repo, sha)
    repo = d.get("Repository", "")
    if repo.startswith("Bioconductor") or "biocViews" in d:
        return repo if repo.startswith("Bioconductor") else "Bioconductor"
    if repo == "CRAN" or d.get("Priority") == "recommended":
        return "CRAN"
    if repo.startswith("RSPM") or repo.startswith("PPM"):
        return "CRAN"
    return repo or "UNKNOWN"


ap = argparse.ArgumentParser()
ap.add_argument("--roots", required=True)
ap.add_argument("--libs", required=True)
ap.add_argument("--labels", required=True)
ap.add_argument("--si")
ap.add_argument("--extra", nargs="*", default=[])
a = ap.parse_args()
libs = a.libs.split(",")
labels = a.labels.split(",")


def tok_name(t):
    m = re.match(r"^(.+?)_([0-9][0-9.\-]*)$", t)
    return (m.group(1), m.group(2)) if m else (t, None)


roots = [tok_name(l.strip())[0] for l in open(a.roots) if l.strip()]
resolved = {}
missing = []
stack = [r for r in roots if r not in BASE]
while stack:
    p = stack.pop()
    if p in resolved or p in BASE:
        continue
    hit = None
    for lib, lab in zip(libs, labels):
        dp = os.path.join(lib, p, "DESCRIPTION")
        if os.path.isfile(dp):
            hit = (lab, parse_dcf(dp))
            break
    if hit is None:
        missing.append(p)
        resolved[p] = None
        continue
    if hit[1].get("Priority") == "base":
        continue
    resolved[p] = hit
    stack.extend(deps(hit[1]))

if missing:
    sys.exit("unresolved packages: " + ", ".join(sorted(missing)))

if a.si:
    for line in open(a.si):
        n, v = tok_name(line.strip())
        if not n or n in BASE:
            continue
        if n not in resolved:
            sys.exit("sessionInfo package %s not in closure" % n)
        if v and resolved[n][1]["Version"] != v:
            sys.exit("version mismatch %s: sessionInfo %s, library %s" % (n, v, resolved[n][1]["Version"]))

rows = []
for p, (lab, d) in resolved.items():
    rows.append((p, d["Version"], source_of(d), lab))
for e in a.extra:
    n, v, s, o = e.split("=")
    rows.append((n, v, s, o))
rows.sort(key=lambda r: r[0].lower())
print("package\tversion\tsource\torigin")
for r in rows:
    print("\t".join(r))
