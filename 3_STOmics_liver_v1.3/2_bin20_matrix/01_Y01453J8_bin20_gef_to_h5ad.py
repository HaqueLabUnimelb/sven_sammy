#!/usr/bin/env python3
"""
Convert a Stereo-seq tissue.gef file into a bin20 raw-count H5AD file
for downstream Seurat and RCTD analysis.

Sample: Y01453J8
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import stereo as st
from scipy import sparse


# =============================================================================
# Input and output
# =============================================================================

gef_file = (
    "/data/projects/punim1236/output_liver_stomics/"
    "Y01453J8/outs/feature_expression/Y01453J8.tissue.gef"
)

output_dir = Path(
    "/data/gpfs/projects/punim1236/"
    "liver_curio/2026714RCTDnew/bin20"
)

output_h5ad = output_dir / "Y01453J8_tissue_bin20_raw_counts.h5ad"

output_dir.mkdir(parents=True, exist_ok=True)


# =============================================================================
# Read bin20 GEF
# =============================================================================

print("=" * 70)
print("Reading Stereo-seq tissue.gef as bin20")
print("=" * 70)

print(f"Input GEF: {gef_file}")

data = st.io.read_gef(
    file_path=gef_file,
    bin_size=20
)

print(f"Data shape: {data.shape}")
print(f"Number of bins: {data.n_cells}")
print(f"Number of genes: {data.n_genes}")


# =============================================================================
# Extract raw-count matrix
# =============================================================================

expression_matrix = data.exp_matrix

if not sparse.issparse(expression_matrix):
    expression_matrix = sparse.csr_matrix(expression_matrix)
else:
    expression_matrix = expression_matrix.tocsr()

print(f"\nExpression matrix type: {type(expression_matrix)}")
print(f"Expression matrix shape: {expression_matrix.shape}")
print(f"Number of non-zero entries: {expression_matrix.nnz:,}")


# =============================================================================
# Check whether X contains raw integer counts
# =============================================================================

print("\nChecking raw counts...")

if expression_matrix.nnz > 0:
    sample_values = expression_matrix.data[: min(1000, expression_matrix.nnz)]

    print(f"Sample non-zero values: {sample_values[:10]}")

    is_integer = np.allclose(
        sample_values,
        np.round(sample_values)
    )

    print(f"Are sampled values integers? {is_integer}")

    if not is_integer:
        raise ValueError(
            "Expression matrix contains non-integer values. "
            "The matrix may not contain raw counts."
        )

    expression_matrix.data = np.round(
        expression_matrix.data
    ).astype(np.int64)


# =============================================================================
# Extract bin IDs and coordinates
# =============================================================================

cell_ids = np.asarray(data.cell_names).astype(str)
coords = np.asarray(data.position)

if coords.ndim != 2 or coords.shape[1] < 2:
    raise ValueError(
        f"Unexpected coordinate shape: {coords.shape}"
    )

if expression_matrix.shape[0] != len(cell_ids):
    raise ValueError(
        "Number of rows in the expression matrix does not match "
        "the number of bin IDs."
    )

if coords.shape[0] != len(cell_ids):
    raise ValueError(
        "Number of coordinate rows does not match the number of bin IDs."
    )

print(f"\nCoordinate shape: {coords.shape}")
print(f"Example bin IDs: {cell_ids[:5]}")
print(f"Example coordinates:\n{coords[:5, :2]}")


# =============================================================================
# Extract gene names
# =============================================================================

genes_df = data.genes.to_df()
fallback_gene_names = np.asarray(data.gene_names).astype(str)

print(f"\nGene metadata columns: {genes_df.columns.tolist()}")

if "real_gene_name" in genes_df.columns:
    real_gene_names = genes_df["real_gene_name"]

    valid_real_names = (
        real_gene_names.notna()
        & (real_gene_names.astype(str) != "")
        & (real_gene_names.astype(str).str.lower() != "nan")
    )

    gene_names = fallback_gene_names.copy()
    gene_names[valid_real_names.to_numpy()] = (
        real_gene_names.loc[valid_real_names]
        .astype(str)
        .to_numpy()
    )

    print("Using real_gene_name where available.")
else:
    gene_names = fallback_gene_names
    print("real_gene_name not found; using data.gene_names.")

if expression_matrix.shape[1] != len(gene_names):
    raise ValueError(
        "Number of matrix columns does not match the number of genes."
    )

print(f"Example gene names: {gene_names[:10]}")


# =============================================================================
# Create AnnData
# =============================================================================

print("\nCreating AnnData object...")

obs = pd.DataFrame(index=pd.Index(cell_ids, name="bin_id"))
var = pd.DataFrame(index=pd.Index(gene_names, name="gene"))

adata = ad.AnnData(
    X=expression_matrix,
    obs=obs,
    var=var
)

adata.obs_names_make_unique()
adata.var_names_make_unique()

adata.obs["x"] = coords[:, 0]
adata.obs["y"] = coords[:, 1]

adata.obs["nCount_RNA"] = np.asarray(
    adata.X.sum(axis=1)
).ravel()

adata.obs["nFeature_RNA"] = np.asarray(
    (adata.X > 0).sum(axis=1)
).ravel()


# =============================================================================
# Final checks
# =============================================================================

print("\nAnnData summary:")
print(adata)

print(f"Shape: {adata.shape}")
print(
    f"Mean UMI per bin: "
    f"{adata.obs['nCount_RNA'].mean():.2f}"
)
print(
    f"Median UMI per bin: "
    f"{adata.obs['nCount_RNA'].median():.2f}"
)
print(
    f"Mean genes per bin: "
    f"{adata.obs['nFeature_RNA'].mean():.2f}"
)
print(
    f"Median genes per bin: "
    f"{adata.obs['nFeature_RNA'].median():.2f}"
)


# =============================================================================
# Save
# =============================================================================

print(f"\nSaving H5AD:\n{output_h5ad}")

adata.write_h5ad(
    output_h5ad,
    compression="gzip"
)

print("\nDone.")