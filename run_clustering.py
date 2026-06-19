"""CLUSTERING: a genome-tile matrix + spectral embedding, not genes + PCA.

scATAC has no gene-count matrix. The feature matrix is the genome binned into fixed TILES (here 500bp),
giving a cell x ~6-million-tile matrix that is extremely sparse and near-binary (a tile is open or not
in a given cell). Two consequences shape the pipeline:
  - You SELECT the most informative tiles (top ~250k by accessibility) before embedding.
  - You don't run PCA on log-counts. The standard ATAC dimensionality reduction is a SPECTRAL
    embedding (a graph-Laplacian / LSI-style method on a cell-cell similarity matrix) that is
    appropriate for sparse binary data. From there it's the usual kNN graph -> Leiden -> UMAP.

This script runs (or loads the cached) full pipeline and reports the structure it recovers.

  python run_clustering.py
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from src import data


def run():
    adata = data.get_processed()       # builds + caches on first call, then loads instantly
    obs = data.obs_df(adata)           # backed obs -> pandas
    var = adata.var[:]
    var = var.to_pandas() if hasattr(var, "to_pandas") else var

    n_selected = int(var["selected"].sum()) if "selected" in var.columns else data.N_FEATURES
    print(f"cells: {adata.n_obs}")
    print(f"genome tiles ({data.BIN_SIZE}bp bins): {adata.n_vars:,}   (sparse, near-binary)")
    print(f"selected features for embedding       : {n_selected:,}")
    print(f"spectral embedding                    : {adata.obsm['X_spectral'].shape[1]} components "
          f"(the ATAC analogue of PCA, built for sparse binary data)")
    if "X_umap" in list(adata.obsm.keys()):
        print(f"UMAP                                  : {adata.obsm['X_umap'].shape}")

    vc = pd.Series(np.asarray(obs["leiden"]).astype(str)).value_counts()
    print(f"\nLeiden clusters: {len(vc)}")
    for c in sorted(vc.index, key=int):
        print(f"  cluster {c:>2s}: {vc[c]:4d} cells")

    print("\nTakeaway: the ATAC feature space is millions of genome tiles (sparse, near-binary), not\n"
          "genes — so you select informative tiles and use a SPECTRAL embedding rather than PCA on\n"
          "log-counts. Clusters are open-chromatin states; run_gene_activity.py turns them into cell\n"
          "types via gene-activity scores.")
    return len(vc)


if __name__ == "__main__":
    run()
