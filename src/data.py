"""Load + process the 10x PBMC 5k scATAC-seq dataset with snapATAC2.

snapATAC2 ships the 10x Genomics 5k PBMC scATAC-seq run as a fragment file (one BED-like row per
sequenced fragment: chrom, start, end, cell-barcode, count). The whole analysis is built FROM that
file — there is no gene-count matrix to start from. The pipeline here is the standard snapATAC2 flow:

  fragments -> per-cell QC (TSS enrichment) -> filter real cells -> genome TILE matrix (500bp bins)
            -> select features -> SPECTRAL embedding -> kNN graph -> Leiden clusters

The expensive steps (import + tile matrix) are run once and the result is cached to a backed .h5ad,
so the downstream scripts (clustering report, gene activity, peak calling) reuse it instantly.
"""
from __future__ import annotations
from pathlib import Path

GENOME = None  # set lazily to snap.genome.hg38 (the dataset is aligned to hg38/GRCh38)
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PROCESSED = DATA_DIR / "pbmc5k_processed.h5ad"

# QC thresholds (the two ATAC-specific gates): minimum unique fragments and minimum TSS enrichment.
MIN_COUNTS = 1000
MIN_TSSE = 5.0
BIN_SIZE = 500          # tile width in bp
N_FEATURES = 250_000    # most-accessible tiles kept for the embedding
N_COMPS = 30            # spectral components
N_NEIGHBORS = 50
RESOLUTION = 1.0


def genome():
    import snapatac2 as snap
    return snap.genome.hg38


def fragment_file():
    """Path to the 10x PBMC-5k fragment file (snapATAC2 downloads + caches it on first call)."""
    import snapatac2 as snap
    return snap.datasets.pbmc5k("fragment")


def import_with_qc(min_num_fragments: int = 100, file=None):
    """Import all barcodes (>= min_num_fragments) and compute TSS enrichment. NOT yet cell-filtered —
    this is the raw view used by run_qc.py to show the cell-vs-debris separation. The 10x fragment
    file is coordinate-sorted, so sorted_by_barcode=False is required (else snapATAC2 panics)."""
    import snapatac2 as snap
    g = snap.genome.hg38
    data = snap.pp.import_fragments(fragment_file(), chrom_sizes=g,
                                    min_num_fragments=min_num_fragments,
                                    sorted_by_barcode=False, file=file)
    snap.metrics.tsse(data, g)      # adds data.obs['tsse']
    return data


def build_processed(path: Path = PROCESSED):
    """Run the full pipeline from fragments and cache the result to a backed .h5ad."""
    import snapatac2 as snap
    g = snap.genome.hg38
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    data = snap.pp.import_fragments(fragment_file(), chrom_sizes=g, min_num_fragments=100,
                                    sorted_by_barcode=False, file=str(path))
    snap.metrics.tsse(data, g)
    snap.pp.filter_cells(data, min_counts=MIN_COUNTS, min_tsse=MIN_TSSE)   # real cells only
    snap.pp.add_tile_matrix(data, bin_size=BIN_SIZE)                       # cell x 500bp-tile matrix
    snap.pp.select_features(data, n_features=N_FEATURES, verbose=False)
    snap.tl.spectral(data, n_comps=N_COMPS)        # ATAC's LSI-style DR (NOT PCA on log counts)
    snap.pp.knn(data, n_neighbors=N_NEIGHBORS, use_rep="X_spectral")
    snap.tl.leiden(data, resolution=RESOLUTION)
    snap.tl.umap(data, use_rep="X_spectral")       # for visualization; cached so all scripts can read it
    return data


def get_processed(path: Path = PROCESSED):
    """Load the cached processed dataset, building (and caching) it on first call."""
    import snapatac2 as snap
    if path.exists():
        return snap.read(str(path))
    return build_processed(path)


def obs_df(adata):
    """Return a pandas DataFrame of .obs (snapATAC2 backed AnnData materializes obs as polars)."""
    o = adata.obs[:]
    return o.to_pandas() if hasattr(o, "to_pandas") else o
