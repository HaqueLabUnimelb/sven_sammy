# Python environments for the sven_sammy archive: build and test results

Tests run on 20 November 2026 on Spartan (University of Melbourne), partition `sapphire`
(Intel Xeon Gold 6448H). Every environment was built inside a SLURM job on the node-local
disk `/tmp` with Miniforge3 24.7.1 (conda 24.7.1, libmamba solver), from conda-forge and
PyPI only, with a fresh `HOME`, `PYTHONNOUSERSITE=1` and no `PYTHONPATH`, and removed at the
end of the job.

## Summary

| Test | Result |
|---|---|
| 2. `01_Y01453J8_bin20_gef_to_h5ad.py` in `sven_sammy_stereopy` | The h5ad file it writes is byte-identical to the archived file (same md5). |
| 1c. pandas command of `2_CellChat/commands_for_two_extra_tables/` in `sven_sammy_hotspot` | Output CSV byte-identical to the archived CSV. |
| 1b. notebook executed with `jupyter nbconvert --execute` | Runs without error. With the same random seed, the notebook and the script (test a) write byte-identical tables (21 of 21 files). |
| 1a/1b. Hotspot tables against the archived tables | Not identical, and no run can reproduce them exactly: Hotspot builds its neighbour graph with a random tie-break that the notebook does not seed, so two runs in the same environment also differ. Compared with 12 new runs of the autocorrelation step and 4 new full runs, the archived run behaves like one more run: its differences from the new runs are as large as the differences between new runs (max abs difference of C 3.8e-4 to 9.6e-4 against 3.8e-4 to 1.0e-3 between new runs; 479 to 483 of about 490 shared genes in the same module against 474 to 483). |
| Lock files | Both environments rebuilt from `*.lock.txt` + `*.pip-freeze.txt` give the same `pip list`. |

## Environment files

| File | What it is |
|---|---|
| `sven_sammy_hotspot.yml` | Python 3.9.13. The 25 packages that the notebook's first cell installed with pip, at the versions pip reported, plus `fonttools` 4.25.0, all from PyPI; the other packages that the notebook used from the Anaconda3/2022.10 module (scipy 1.9.1, h5py 3.7.0, statsmodels 0.13.2, numexpr 2.8.3, bottleneck 1.3.5, networkx 2.8.4, tqdm 4.64.1, ...) and Jupyter (ipykernel 6.15.2, nbconvert 6.4.4, ...) at the versions of that module, from conda-forge. numpy and pandas are listed twice: conda installs the Anaconda 2022.10 versions (1.21.5, 1.4.4) and pip replaces them by the versions the notebook used (1.24.4, 2.3.3), as the original user-site installation did. |
| `sven_sammy_stereopy.yml` | Python 3.8.20 (the same conda-forge build as the removed environment `stereopy_latest_py38`), stereopy 1.6.2, gefpy 1.3.2, the five other packages the script imports and the 168 further packages that `pip install stereopy==1.6.2 ipython` needs, each pinned to the version the original run could import (see "Notes"). |
| `*.lock.txt` | `conda list --explicit --md5` of the built environment (conda packages, with URLs and md5). |
| `*.pip-freeze.txt` | `pip list --format=freeze` of the built environment: `name==version` of every Python distribution, conda- or pip-installed. |
| `tests/logs/*.conda-env-export.yml` | `conda env export --no-builds` (the `prefix:` line removed). |
| `tests/logs/*.conda-list.txt`, `*.pip-freeze-raw.txt`, `*.python-version.txt` | `conda list`, raw `pip freeze`, `sys.version`. |

`sys.version`:
- `sven_sammy_hotspot`: `3.9.13 | packaged by conda-forge | (main, May 27 2022, 16:58:50) [GCC 10.3.0]`
- `sven_sammy_stereopy`: `3.8.20 | packaged by conda-forge | (default, Sep 30 2024, 17:52:49) [GCC 13.3.0]`

## How a reader creates and uses the environments

Packages in `~/.local/lib/python3.X` take precedence over the environment, so set
`export PYTHONNOUSERSITE=1` (and unset `PYTHONPATH`) before creating or using either environment:
`conda env create` runs pip, which otherwise skips packages that are already in `~/.local` and can
uninstall other versions from there.

Create from the YAML files (Linux x86_64, conda or mamba with access to conda-forge and PyPI):

```
export PYTHONNOUSERSITE=1
unset PYTHONPATH
conda env create -f sven_sammy_hotspot.yml
conda env create -f sven_sammy_stereopy.yml
```

Or rebuild the exact tested environments from the lock files (with the same two settings):

```
conda create -n sven_sammy_hotspot --file sven_sammy_hotspot.lock.txt
conda run -n sven_sammy_hotspot python -m pip install --no-deps -r sven_sammy_hotspot.pip-freeze.txt
conda create -n sven_sammy_stereopy --file sven_sammy_stereopy.lock.txt
conda run -n sven_sammy_stereopy python -m pip install --no-deps -r sven_sammy_stereopy.pip-freeze.txt
```

bin20 matrix (step 2 of `3_STOmics_liver_v1.3`; edit `gef_file` and `output_dir` first):

```
conda activate sven_sammy_stereopy
python 01_Y01453J8_bin20_gef_to_h5ad.py
```

Hotspot notebook (step 4; edit the input path and the `os.chdir` target, and set
`n_neighbors=30` to follow the archived tables). Run it top to bottom without opening it:

```
conda activate sven_sammy_hotspot
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 \
    --output hotspot_liverstomics_NI_executed.ipynb hotspot_liverstomics_NI.ipynb
```

or register the environment as a kernel of an existing Jupyter installation. A kernel started by
another Jupyter installation does not inherit `PYTHONNOUSERSITE` from the shell, so the command writes it
into the kernel definition:

```
conda activate sven_sammy_hotspot
python -m ipykernel install --user --name sven_sammy_hotspot --display-name "Python 3.9 (sven_sammy_hotspot)" \
    --env PYTHONNOUSERSITE 1
```

The pandas command of `2_CellChat/commands_for_two_extra_tables/` (edit the paths first):

```
conda activate sven_sammy_hotspot
bash make_MNP2_vs_original_combinedMNP_comparison.sh
```

## Test 2: bin20 matrix with stereopy (`sven_sammy_stereopy`)

Job 35423815 (`tests/stereopy_env_test.sbatch`, log `tests/logs/stereopy_env_test_35423815.log`).
A copy of `01_Y01453J8_bin20_gef_to_h5ad.py` in which only `output_dir` points to `/tmp`
(`tests/logs/stereopy_35423815_script_diff.txt`) read the archived `Y01453J8.tissue.gef`.
`tests/compare_h5ad.py` compared the h5ad it wrote with the archived
`/data/projects/punim1236/liver_curio/2026714RCTDnew/bin20/Y01453J8_tissue_bin20_raw_counts.h5ad`
(`tests/logs/stereopy_35423815_compare.txt` and `.json`):

| Check | Result |
|---|---|
| whole file | byte-identical: md5 `000f5f0fef13bd5aa5971c3d232b3570` for both (148,652,398 bytes) |
| `X` | `csr_matrix`, int64, 293,543 x 24,886, 57,212,888 non-zeros; `data`, `indices`, `indptr` identical |
| `obs_names`, `var_names` | identical (index names `bin_id`, `gene`) |
| `obs['x']`, `obs['y']` (uint32), `obs['nCount_RNA']`, `obs['nFeature_RNA']` (int64) | identical, same dtypes |
| all HDF5 objects | 9 datasets compared; values, dtypes, attributes, compression and chunking identical |

The script printed the same summary as the archived run's downstream log (mean 356.92 and
median 351 UMI per bin, mean 194.90 and median 192 genes per bin).

Which stereopy version wrote the archived file: the test cannot tell. Side check, job
35424101 (`sbatch stereopy_env_test.sbatch sven_sammy_stereopy check_stereopy_1.5.1_gefpy_1.1.20 "stereopy==1.5.1 gefpy==1.1.20"`):
the same environment with stereopy 1.5.1 and gefpy 1.1.20 (the versions of the other two
stereopy environments on Spartan) installed over it also wrote a byte-identical file
(same md5; `tests/logs/stereopy_35424101_compare.txt`). Stereopy 1.6.2 is supported by
the package list of the removed environment, and none of the 148 Python package versions
in that list was released on PyPI after 6 July 2026; the environment was created on
7 July 2026 (its Jupyter kernel was registered at 22:18 that day), nine days before the
run. Version "1.2.2" (named in `/data/projects/punim1236/liver_stomicsv1.3/2026714RCTDnew/README.md`)
does not exist: PyPI and the `stereopy` conda channel have no stereopy 1.2.2, and
`stereo.__version__` of 1.6.2 is "1.6.2".

Run time and memory (GNU time): build 1 min 40 s; script 41 s, peak RSS 8.0 GiB; comparison
10 s, 2.4 GiB. Whole job 3 min 33 s, 8 CPUs, 96 GB requested, sacct MaxRSS 11.3 GiB.

## Test 1c: pandas command in `sven_sammy_hotspot`

Job 35423407 (`tests/hotspot_env_test.sbatch`). A copy of
`2_CellChat/commands_for_two_extra_tables/make_MNP2_vs_original_combinedMNP_comparison.sh`
with two path changes (`tests/logs/hotspot_35423407_test_c_command_diff.txt`): `newp` points
to the renamed folder `/data/projects/punim1236/liver_scRNAseq/Oct22newcellchat/MNP2_CD4_20261023/`,
and `out` to `/tmp`. Run with the environment's `python`: 0.4 s, peak RSS 110 MB.
The output is byte-identical to `2_CellChat/MNP2_CD4_20261023/MNP2_vs_original_combinedMNP_comparison.csv`
(`cmp`; md5 `78fcc574d4a5b1ff26e053dbff977380` for both). Printed output:
`tests/logs/hotspot_35423407_test_c_stdout.txt`.

## Tests 1a and 1b: Hotspot in `sven_sammy_hotspot`

Archived tables: `/data/projects/punim1236/liver_stomicsv1.3/hotspot/hotspot_results/` and
`hotspot_plots/` (20 October 2026), and the repository copies in
`3_STOmics_liver_v1.3/4_Hotspot/result_tables/`. All runs read the archived
`Y01453J8_tissue_bin20_raw_counts.h5ad` and kept 3,467 genes, as the notebook did.

### Why no run reproduces the archived tables exactly

`Hotspot.create_knn_graph()` of hotspotsc 1.1.3 has `approx_neighbors=True` by default and
then calls `pynndescent.NNDescent(coords, n_neighbors=n_neighbors + 1)` without a
`random_state`; the notebook sets no seed. On the bin20 grid many neighbours are equally
far away: for a bin surrounded by tissue, the 30 nearest bins are the 28 bins within 3 bin
widths plus 2 of the 8 bins at sqrt(10) bin widths (for 15: the 12 bins within 2 bin widths
plus 3 of the 8 at sqrt(5)). Which of the tied bins enter the graph is random.
`tests/hotspot_knn_check.py` (job 35427214, `tests/logs/hotspot_35427214_knn_n30_knn_check.txt`)
shows this directly:

| Two graph builds in one process | n_neighbors=30 | n_neighbors=15 |
|---|---|---|
| unseeded, as in the notebook: bins whose neighbour set differs | 276,609 of 293,543 | 280,445 of 293,543 |
| `numpy.random.seed(0)` before each build (default numba threads, and 1 thread) | 0 | 0 |
| bins where the k-th neighbour distance of NNDescent differs from the exact graph (scikit-learn ball tree) | 0 | 0 |

So NNDescent always returns neighbours at the correct distances, and the runs differ only
in how ties are broken. Every later step is deterministic (see the seeded pair below).
The random state of the archived run was never recorded, so its graph, and therefore its
exact tables, cannot be re-created in any environment. The tests below therefore check two
things: that the notebook and the script compute the same thing in this environment
(seeded pair), and that the archived tables are what an unseeded run in this environment
produces (replicates).

### Test 1a: the notebook's Hotspot code as a script, n_neighbors=30

`tests/hotspot_test_a.py`: the notebook's code cells top to bottom without the plotting
lines, with `n_neighbors=30` and `jobs=JOBS` (worker processes only). Job 35423407, 23
workers: 15 min 11 s (autocorrelation 4 min, local correlations 6 min 42 s, module scores
17 s, the extra n_neighbors=15 autocorrelation 4 min), peak RSS 3.8 GiB. Comparison with the
archived tables (`tests/logs/hotspot_35423407_test_a_vs_archived_rerun.txt`; the comparison
in job 35423407 itself stopped at the ranked-gene table, which has a different number of
rows, and was re-run in job 35425958 with `tests/compare_rerun.sbatch`):

| Table | Result |
|---|---|
| `hotspot_autocorrelation.csv` | same 3,467 genes; max abs difference C 4.90e-4, Z 2.13, Pval 0.26, FDR 0.26 (the large Pval and FDR differences are for genes with Z near 0); gene order (by Z) not identical |
| 500 genes for the local correlations | 490 in common with the archived 500 |
| `hotspot_modules.csv` | 483 of the 490 shared genes in the same module; module sizes 95, 62, 25, 102, 160, 34 (22 unassigned) against archived 97, 62, 26, 105, 152, 33 (25 unassigned) |
| `hotspot_local_correlation_z.csv` | max abs difference 0.95 over the shared genes |
| `hotspot_spot_scores.csv` | same 293,543 bins in the same order, x and y identical; max abs difference of the module scores 0.36, 0.46, 0.70, 0.29, 0.67, 0.86 (modules 1 to 6); bins with `module_5 > 4.0`: 1,891 against 1,882, 1,862 in common |
| `hotspot_module_genes_ranked.csv` | 461 gene-module pairs in common; max abs difference of `corr_with_module` 0.0022 |
| files byte-identical to the archived ones | none (one gene list, `module_2_genelist.txt`, is identical) |

n_neighbors=15, autocorrelation only: Gm42418 C 0.477286, Z 813.153 against the notebook's
printed C 0.477215, Z 813.077882. Twelve new runs with n_neighbors=15 give Gm42418 C from
0.476998 to 0.477523 and Z from 812.814 to 813.642, so the printed values lie inside the
spread of new runs; for the 20 genes the notebook prints, C lies inside the range of the 12
new runs for 19 genes and Z for 20 genes (`tests/logs/compare_replicates_35427218.txt`).

### Test 1b: the notebook executed with nbconvert

`tests/make_notebook_copy.py` writes a copy with three changes: the input path of
`sc.read_h5ad` (here the same file, so this line is unchanged), the `os.chdir` target
(a folder in `/tmp`) and `n_neighbors=15` to 30. The diff of all code cells is in
`tests/logs/hotspot_35423407_test_b_notebook_diff.txt`. Command:
`jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=python3 --output ..._executed.ipynb ...`.

- Job 35423407 (run at the same time as test a): finished without error in 18 min 45 s,
  peak RSS 3.9 GiB; 14 of the 15 code cells executed (the 15th is empty), no error output.
  The first cell printed "Requirement already satisfied" for hotspotsc and scanpy. As
  expected from the random graph, all 21 tables it wrote differ from those of test a
  (`tests/logs/hotspot_35423407_test_b_vs_test_a.txt`); its printed top gene was
  Gm42418 C 0.455423, Z 1080.583453 (archived 0.455460, 1080.620850).
- Seeded pair, job 35427216: the script (`HOTSPOT_TEST_SEED=0`, which calls
  `numpy.random.seed(0)` before anything else) and the notebook copy (seeded by an IPython
  startup file `00-seed.py` in the job's private `IPYTHONDIR`, so the notebook itself has
  only the three changes), both with `NUMBA_NUM_THREADS=1`, ran side by side with 8
  workers each (26 min 5 s and 26 min 30 s, peak RSS 3.8 and 3.9 GiB). All 21 tables were
  byte-identical (`tests/logs/hotspot_35427216_seeded_b_seeded_vs_a_seeded.txt`): the 6
  files of `hotspot_results/` (autocorrelation, modules, gene modules, local correlation Z,
  module scores, spot scores) and the 15 tables and gene lists of `hotspot_plots/`. The
  notebook run had no error output.

### The archived run against new runs

`tests/compare_replicates.py` (job 35427218, `tests/logs/compare_replicates_35427218.txt`
and `.json`) compared the archived tables with 12 new autocorrelation runs per neighbour
count (8 runs of the autocorrelation step alone, jobs 35427214 and 35427215, and the
autocorrelation of the 4 full runs) and with 4 new full runs (test a of job 35423407, two
unseeded runs `a2` and `a3` of job 35427217, and the seeded script run of job 35427216).
"Between new runs" uses all pairs of new runs.

| n_neighbors=30 | archived against new runs | between new runs |
|---|---|---|
| max abs difference of C (per pair) | 3.8e-4 to 9.6e-4 (median 5.7e-4) | 3.8e-4 to 1.03e-3 (median 5.3e-4) |
| max abs difference of Z (per pair) | 0.80 to 3.04 (median 2.08) | 0.90 to 3.07 (median 1.53) |
| genes whose archived C (Z) lies outside the range of the 12 new runs | 13.9% (13.8%) | 15.4% to 17.8% for one new run against the other 11 (15.4% expected for an exchangeable run) |
| mean signed difference of the archived C and Z from the mean of the new runs | -0.02 SD of the new runs | |
| overlap of the 500 genes | 488 to 492 (median 490) | 485 to 492 (median 489) |
| top 20 genes printed in `Untitled3.ipynb` (n_neighbors=30, 6 decimals) inside the range of the new runs | C 17 of 20, Z 17 of 20 | |

| Full runs | archived against new runs (4 pairs) | between new runs (6 pairs) |
|---|---|---|
| genes in both 500-gene lists with the same module | 479 to 483 (median 481) | 474 to 483 (median 481.5) |
| adjusted Rand index of the module labels | 0.945 to 0.970 (median 0.952) | 0.917 to 0.967 (median 0.953) |
| Jaccard index of the module 5 genes | 0.894 to 0.938 (median 0.918) | 0.877 to 0.923 (median 0.908) |
| Jaccard index of the bins with `module_5 > 4.0` | 0.967 to 0.974 (median 0.969) | 0.967 to 0.972 (median 0.969) |
| Pearson r of the per-bin scores, modules 1 to 5 | 0.9964 to 0.9993 | 0.9946 to 0.9993 |
| Pearson r of the per-bin scores, module 6 | 0.9940 to 0.9947 | 0.9954 to 0.9958 |
| max abs difference of the per-bin scores, all modules | 0.86 to 1.17 | 0.79 to 1.20 |

Module numbers mean the same in all runs: Alb is in module 1, Glul in 2, Gm42418 in 3,
Apoa2 in 4, Tmsb4x, Lyz2, H2-Aa and Cd74 in 5 and Gm26917 in 6, in the archived run and in
all new runs, so step 5 (which selects module 5 by its number) selects the same module.
Module 5 has 149 to 160 genes in the new runs (archived 152), and 1,854 to 1,891 bins with
`module_5 > 4.0` (archived 1,882).

The only measure in which the archived run agrees less well with the new runs than the new
runs agree with each other is the per-bin score of module 6. Its gene list differs from all four new runs in four genes near the
cut-offs: Prkg1 and Slco3a1 are in module 6 in the archived run and in module 5 in every new
run, and Lncpint and Pard3 are in module 6 of every new run but just outside the archived
run's 500 genes (Z ranks 503 and 517; ranks 470 to 524 and 434 to 510 in the 12 new runs).

Conclusion: in this environment the notebook and the script compute the same thing, and
the archived tables are within the run-to-run spread of new runs. The environment can
therefore be said to reproduce the archived Hotspot results up to the randomness of the
neighbour graph, but not byte for byte. A reader who re-runs the notebook will get tables
that differ from the archived ones by about as much as the "between new runs" column.

## Records and lock check

- `sven_sammy_hotspot`: 133 conda packages in the lock file, 90 Python distributions in the
  pip list. `sven_sammy_stereopy`: 26 conda packages, 178 Python distributions; its pip list
  equals the YAML pins.
- Lock check (end of jobs 35423407 and 35423815, `tests/env_common.sh: lock_check`): a
  second environment created with `conda create --file <name>.lock.txt` followed by
  `pip install --no-deps -r <name>.pip-freeze.txt` has the same `pip list` as the
  environment built from the YAML (`tests/logs/<name>.lockcheck.log`; 19 s and 46 s).
- Rebuild check: the four jobs of `hotspot_replicates.sbatch` built `sven_sammy_hotspot`
  again from the YAML; each time the lock file and the pip list were identical to the
  recorded ones.

## Jobs

All jobs on partition `sapphire` (Intel Xeon Gold 6448H). Build times and per-step times
are from GNU time (`tests/logs/*_times.txt`); MaxRSS is sacct's peak for the whole batch
step (several processes at once in most jobs).

| Job | sbatch file | Account | CPUs / memory requested | Elapsed | MaxRSS | Content |
|---|---|---|---|---|---|---|
| 35423407 | `hotspot_env_test.sbatch` | uom00060 | 32 / 128 GB | 20 min 27 s | 12.3 GiB | build 56 s; test c 0.4 s; tests a (15 min 11 s) and b (18 min 45 s) at the same time; comparisons; lock check |
| 35423815 | `stereopy_env_test.sbatch` | punim1445 | 8 / 96 GB | 3 min 33 s | 11.3 GiB | build 1 min 40 s; script 41 s (peak RSS 8.0 GiB); comparison 10 s; lock check 46 s |
| 35424101 | `stereopy_env_test.sbatch sven_sammy_stereopy check_stereopy_1.5.1_gefpy_1.1.20 "stereopy==1.5.1 gefpy==1.1.20"` | uom00103 | 8 / 96 GB | 3 min 17 s | 6.0 GiB | side check with stereopy 1.5.1 (script 25 s, peak RSS 2.2 GiB) |
| 35425958 | `compare_rerun.sbatch` | uom00060 | 1 / 16 GB | 12 s | | comparison of test a of job 35423407, re-run |
| 35427214 | `hotspot_replicates.sbatch knn_n30` | punim2394 | 12 / 48 GB | 9 min 41 s | 23.1 GiB | neighbour-graph check; 8 autocorrelation runs, n_neighbors=30 (6 min 5 s) |
| 35427215 | `hotspot_replicates.sbatch n15` | uom00103 | 10 / 40 GB | 6 min 50 s | 23.7 GiB | 8 autocorrelation runs, n_neighbors=15 (5 min 41 s) |
| 35427216 | `hotspot_replicates.sbatch seeded` | uom00060 | 16 / 64 GB | 27 min 36 s | 10.8 GiB | seeded script and seeded notebook, 8 workers each (26 min) |
| 35427217 | `hotspot_replicates.sbatch unseeded` | punim1445 | 16 / 64 GB | 30 min 54 s | 11.1 GiB | two unseeded script runs, 7 workers each (29 min) |
| 35427218 | `compare_replicates.sbatch` | uom00060 | 2 / 24 GB | 9 s | 0.3 GiB | `compare_replicates.py` |

Earlier submissions that are not part of the results:
35423408 (`sven_sammy_stereopy.yml` with `pip=25.0.1`, which conda-forge does not have for
Python 3.8; log not kept), 35423537 (stereopy environment without IPython; log in
`tests/logs/superseded/`) and 35423937 (a separate YAML for stereopy 1.5.1 that pip could
not build; replaced by job 35424101). Jobs 35424563 and 35426837 (first versions of the
replicate job, 48 and 32 CPUs) were cancelled before they started.

Run time of Hotspot: in the notebook of 20 October 2026 the local correlations took
1 h 2 min with 8 workers; here they took 17 min with 8 workers (job 35427216), 21 min with
7 workers (job 35427217) and 6 min 42 s with 23 workers (job 35423407). The
autocorrelation step takes about 4 min whatever the number of workers.

## Notes on the environment files

`sven_sammy_hotspot`
- The notebook's kernel was the Python of the module Anaconda3/2022.10 with the user site of
  the account that ran it. The 25 packages of the `pip install hotspotsc scanpy` cell came
  from PyPI wheels; the YAML installs the same wheels. All other packages it used came from
  the Anaconda base; the YAML installs conda-forge builds of the same versions (read from
  `/apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/conda-meta`).
  conda-forge has no build of fonttools 4.25.0, so it comes from PyPI.
- The builds differ from Anaconda's: the Anaconda base linked scipy and numexpr to MKL,
  conda-forge links them to OpenBLAS 0.3.25, and h5py uses HDF5 1.12.2 instead of 1.10.6.
  numpy itself is the PyPI wheel 1.24.4 with its own OpenBLAS, as in the original.
- On import, pandas warns that numexpr 2.8.3 and bottleneck 1.3.5 are older than it
  supports; the archived notebook shows the same two warnings.
- In the notebook copy, the first cell (`pip install hotspotsc scanpy`) only reports
  "Requirement already satisfied".

`sven_sammy_stereopy`
- The archived run (job 27327611, 16 July 2026) used
  `.../conda_envs/stereopy_latest_py38/bin/python3.8`, and its wrapper
  `run_Y01453J8_bin20_RCTD.slurm` does not set `PYTHONNOUSERSITE`, so Python also imported
  from the user site `~/.local/lib/python3.8/site-packages` of `czhu3869`.
  That user site still exists: 52 of its packages were installed on 7 July 2026, while the
  environment was set up (panel, holoviews, datashader, hvplot, param, SQLAlchemy, urllib3,
  pillow 10.4.0, tqdm 4.65.0, scanpy 1.9.6, pytz 2026.2, ...); they are missing from the
  environment's own package list. The YAML pins every package that
  `pip install stereopy==1.6.2 ipython` installs to the version of the environment list or,
  where the user site had one, of the user site (the user site comes first on `sys.path`).
  Left out because the script does not use them: the Jupyter kernel packages of the
  environment (ipykernel, jupyter_client, ...), and the user-site packages that are not
  stereopy dependencies (sphinx and its extensions, joypy, pkgconfig, and six packages
  installed in March to May 2026).
- stereopy imports IPython but does not declare it. The first build without IPython failed
  at `import stereo` (`tests/logs/superseded/stereopy_env_test_35423537.log`); that build
  also showed that pip, unpinned, now picks newer versions than the original for 14
  packages (for example plotly 7.1.0, python-dateutil 3.0.0, scanpy 1.9.8).
- conda-forge has no pip 25.0.1 for Python 3.8: conda installs pip 24.3.1 and pip then
  upgrades itself to 25.0.1, the version of the original environment.
- The built environment's `pip list` equals the YAML pins exactly (no package added by the
  resolver).

## Files in this folder

| File | Content |
|---|---|
| `RESULTS.md` | this file |
| `sven_sammy_hotspot.yml`, `sven_sammy_stereopy.yml` | environment files |
| `sven_sammy_hotspot.lock.txt`, `sven_sammy_stereopy.lock.txt` | explicit conda lock files |
| `sven_sammy_hotspot.pip-freeze.txt`, `sven_sammy_stereopy.pip-freeze.txt` | pip lists (`name==version`) |
| `tests/env_common.sh` | shell functions: isolation of HOME, TMPDIR and caches, timing, records, lock check |
| `tests/hotspot_env_test.sbatch` | build of `sven_sammy_hotspot`, tests a, b and c, lock check |
| `tests/hotspot_test_a.py` | test a: the notebook's Hotspot code as a script (optional seed) |
| `tests/make_notebook_copy.py` | copy of the notebook with the three permitted changes, with a diff |
| `tests/make_command_copy.py` | copy of a script with literal path replacements, with a diff (tests c and 2) |
| `tests/compare_hotspot.py` | one Hotspot run against the archived tables |
| `tests/compare_rerun.sbatch` | re-run of `compare_hotspot.py` for test a of job 35423407 |
| `tests/hotspot_replicates.sbatch` | neighbour-graph check, seeded pair, replicate runs (four parts) |
| `tests/hotspot_knn_check.py` | neighbour-graph check |
| `tests/hotspot_autocorr_replicate.py` | one run of the autocorrelation step |
| `tests/compare_replicates.sbatch`, `tests/compare_replicates.py` | the archived run against all new runs |
| `tests/stereopy_env_test.sbatch` | build of `sven_sammy_stereopy`, test 2, lock check, side checks of other versions |
| `tests/compare_h5ad.py` | h5ad comparison (named fields, then every HDF5 dataset and attribute) |
| `tests/logs/` | job logs, per-step times, diffs of the copied scripts, comparison outputs (`.txt` and `.json`), environment records; `check_stereopy_1.5.1_gefpy_1.1.20/` holds the records of the side check, `superseded/` the log of the failed first stereopy build |

The sbatch files are submitted from the `tests` folder (`cd tests && sbatch ...`); they find
the YAML files in the parent folder and the root of this repository two folders above it
(`../..`). The data paths and the folder for kept outputs (`/path/to/kept_outputs`, about
450 MB) are variables at the top of each file; the `#SBATCH` lines are for Spartan.

Changes made to the test files after the runs, none of which changes what a run does:
the literal paths of the staging folder and of the repository clone in the sbatch files were
replaced by these relative forms (the tests ran in a folder outside the repository, next to a
local clone of it); the folder for kept outputs, a temporary folder on Spartan, was replaced
by `/path/to/kept_outputs`; the `--account` line of `stereopy_env_test.sbatch` was
changed between submissions (accounts were rotated, and some jobs set `--account` on the
command line, as the Jobs table shows); `stereopy_env_test.sbatch` got its optional side-check arguments after job 35423815 (used
by job 35424101); `hotspot_test_a.py` got the optional `HOTSPOT_TEST_SEED` after job
35423407 (without it the code is the same); the ranked-table step of `compare_hotspot.py`
was made robust to tables of different length after job 35423407 (hence job 35425958).

The logs are program output with three changes: the local folder in which the tests ran is
shown as `<build_dir>`, the folder for kept outputs as `<kept_outputs_dir>`, and the tqdm
progress bars were removed from the four large Hotspot run logs
(`hotspot_35423407_test_a.log`, `hotspot_35427216_seeded_a_seeded.log`,
`hotspot_35427217_unseeded_a2.log`, `hotspot_35427217_unseeded_a3.log`; about 99 % of their
lines). The tables written by the new Hotspot runs (441 MB) are not in this repository.
