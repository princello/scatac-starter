"""GENE ACTIVITY: turning open chromatin into cell types via marker genes.

ATAC clusters are open-chromatin states, not labeled cell types. To annotate them you bridge to known
biology: sum a cell's accessibility over each gene's body + promoter to get a GENE-ACTIVITY score (a
proxy for expression), then score canonical marker-gene sets per cluster. This is how you connect a
purely epigenomic readout to the cell-type vocabulary everyone (and every RNA atlas) speaks.

  python run_gene_activity.py

(Gene activity is a proxy, not RNA — it's good enough to assign broad lineages; for fine states you'd
co-embed with paired/transferred scRNA. The point here is the ATAC->cell-type bridge.)
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import scanpy as sc
import snapatac2 as snap
from src import data

# canonical PBMC marker-gene sets (broad lineages)
MARKERS = {
    "T cell":     ["CD3D", "CD3E", "CD3G", "IL7R", "CD8A"],
    "B cell":     ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "Monocyte":   ["CD14", "LYZ", "S100A8", "FCN1", "CST3"],
    "NK cell":    ["NKG7", "GNLY", "KLRD1", "NCAM1", "KLRF1"],
    "Dendritic":  ["FCER1A", "CLEC9A", "CLEC10A", "IRF8"],
}


def run():
    adata = data.get_processed()
    leiden = np.asarray(data.obs_df(adata)["leiden"]).astype(str)

    # build the cell x gene activity matrix from accessibility over gene body + promoter
    gm = snap.pp.make_gene_matrix(adata, snap.genome.hg38)
    gm = gm.to_memory() if hasattr(gm, "to_memory") else gm
    gm.obs["leiden"] = pd.Categorical(leiden)
    print(f"gene-activity matrix: {gm.shape[0]} cells x {gm.shape[1]} genes")

    # standard normalization, then background-corrected marker scores per cell
    sc.pp.normalize_total(gm, target_sum=1e4)
    sc.pp.log1p(gm)
    present = {}
    for ct, genes in MARKERS.items():
        g = [x for x in genes if x in set(gm.var_names)]
        present[ct] = g
        sc.tl.score_genes(gm, g, score_name=ct)

    # mean marker score per (cluster, cell type) -> z-score across clusters -> assign argmax
    scores = gm.obs.groupby("leiden", observed=True)[list(MARKERS)].mean()
    z = (scores - scores.mean()) / scores.std(ddof=0)
    sizes = pd.Series(leiden).value_counts()

    print("\ncluster -> cell type (by strongest marker-set activity):\n")
    print(f"  {'clust':5s} {'n':>5s}  {'label':10s} {'z':>5s}   runner-up")
    rows = []
    for c in sorted(z.index, key=int):
        ranked = z.loc[c].sort_values(ascending=False)
        best, second = ranked.index[0], ranked.index[1]
        rows.append((c, best))
        print(f"  {c:>5s} {int(sizes[c]):5d}  {best:10s} {ranked.iloc[0]:5.1f}   "
              f"{second} ({ranked.iloc[1]:.1f})")

    label_counts = pd.Series([r[1] for r in rows]).value_counts()
    print("\nclusters per assigned lineage:", label_counts.to_dict())
    missing = {ct: [x for x in MARKERS[ct] if x not in present[ct]] for ct in MARKERS}
    missing = {k: v for k, v in missing.items() if v}
    if missing:
        print("(markers absent from the gene-activity matrix:", missing, ")")
    print("\nTakeaway: gene-activity scores bridge open chromatin to the marker-gene language of cell\n"
          "types, so epigenomic clusters become T / B / Monocyte / NK / DC lineages. It's a proxy for\n"
          "expression (not RNA); for fine states you'd transfer labels from paired scRNA.")
    return rows


if __name__ == "__main__":
    run()
