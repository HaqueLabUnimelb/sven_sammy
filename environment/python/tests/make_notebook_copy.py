"""Test b, step 1: write a copy of hotspot_liverstomics_NI.ipynb with three changes.

  1. the input path in sc.read_h5ad(...)            -> IN_H5AD
  2. the target of os.chdir(...)                    -> OUTDIR
  3. hs.create_knn_graph(..., n_neighbors=15)       -> n_neighbors=30

Each original string must occur exactly once in the notebook. The script prints a diff
of all code cells (original against copy) so that the log shows that nothing else
changed. Outputs stored in the notebook are left as they are; nbconvert --execute
replaces them.

Usage: python make_notebook_copy.py ARCHIVED_NOTEBOOK COPY_PATH IN_H5AD OUTDIR
"""
import difflib
import json
import sys

src, dst, in_h5ad, outdir = sys.argv[1:5]

EDITS = [
    ('sc.read_h5ad("/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.h5ad")',
     f'sc.read_h5ad("{in_h5ad}")'),
    ('os.chdir("/data/projects/punim1236/liver_stomicsv1.3/hotspot/")',
     f'os.chdir("{outdir}")'),
    ('hs.create_knn_graph(weighted_graph=False, n_neighbors=15)',
     'hs.create_knn_graph(weighted_graph=False, n_neighbors=30)'),
]

with open(src, encoding="utf-8") as fh:
    nb = json.load(fh)

orig_code = []
new_code = []
for cell in nb["cells"]:
    if cell["cell_type"] != "code":
        continue
    text = "".join(cell["source"])
    orig_code.append(text)
    for old, new in EDITS:
        text = text.replace(old, new)
    new_code.append(text)
    cell["source"] = text.splitlines(keepends=True)

joined = "\n".join(orig_code)
for old, _ in EDITS:
    n = joined.count(old)
    if n != 1:
        sys.exit(f"ERROR: expected exactly one occurrence of {old!r}, found {n}")

with open(dst, "w", encoding="utf-8") as fh:
    json.dump(nb, fh, indent=1, ensure_ascii=False)
    fh.write("\n")

diff = list(difflib.unified_diff(
    "\n\n".join(orig_code).splitlines(), "\n\n".join(new_code).splitlines(),
    fromfile="archived notebook (code cells)", tofile="test copy (code cells)", lineterm="", n=0))
print("\n".join(diff))
changed = sum(1 for line in diff if line.startswith("-") and not line.startswith("---"))
expected = sum(1 for old, new in EDITS if old != new)
print(f"code lines changed: {changed} (expected {expected})")
if changed != expected:
    sys.exit("ERROR: unexpected number of changed lines")
