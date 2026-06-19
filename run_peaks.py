"""PEAKS: the ATAC-specific payoff — call peaks, find differentially-accessible regions, read off TFs.

Genome tiles are fine for clustering but blunt for interpretation. The real regulatory units are
PEAKS: contiguous stretches of open chromatin. The standard workflow:
  1. call peaks PER CLUSTER with MACS3 (pseudobulk the cells in each cluster, then peak-call), so
     cell-type-specific elements aren't washed out by calling on all cells at once;
  2. MERGE the per-cluster peak sets into one non-overlapping reference;
  3. build a cell x PEAK matrix;
  4. find MARKER PEAKS per cluster = differentially-accessible regions (DARs);
  5. (optional) MOTIF ENRICHMENT in each cluster's marker peaks -> the transcription factors that
     likely drive that cell type's regulatory program (e.g. CEBPA in monocytes, LEF1/TCF7 in T cells,
     EOMES in NK, EBF1 in B cells).

Peak calling (step 1) is the slow part (~8 min, single-process on macOS); the merged peak set is
cached to data/ so re-runs are instant. Steps 4-5 are the interpretation.

  python run_peaks.py                # peaks + differential accessibility
  python run_peaks.py --motifs       # also run TF-motif enrichment (downloads hg38 FASTA + cis-BP)
"""
from __future__ import annotations
import argparse
import numpy as np
import snapatac2 as snap
from src import data

PEAKS_FILE = data.DATA_DIR / "pbmc5k_merged_peaks.txt"


def merged_peak_list(adata):
    """Per-cluster MACS3 peak calling + merge, cached to PEAKS_FILE (the slow ~8-min step)."""
    if PEAKS_FILE.exists():
        peaks = [ln.strip() for ln in PEAKS_FILE.read_text().splitlines() if ln.strip()]
        print(f"loaded {len(peaks):,} merged peaks from cache")
        return peaks
    print("calling peaks per cluster with MACS3 (single-process; ~8 min on macOS) ...")
    per_cluster = snap.tl.macs3(adata, groupby="leiden", inplace=False, n_jobs=1)
    total = sum(df.shape[0] for df in per_cluster.values())
    print(f"  MACS3: {len(per_cluster)} clusters, {total:,} peaks before merge "
          f"(median {int(np.median([df.shape[0] for df in per_cluster.values()]))}/cluster)")
    merged = snap.tl.merge_peaks(per_cluster, snap.genome.hg38)
    peaks = merged["Peaks"].to_list()
    data.DATA_DIR.mkdir(parents=True, exist_ok=True)
    PEAKS_FILE.write_text("\n".join(peaks) + "\n")
    print(f"  merged -> {len(peaks):,} non-overlapping peaks (cached)")
    return peaks


def run(motifs: bool = False, top_n: int = 8):
    adata = data.get_processed()
    leiden = np.asarray(data.obs_df(adata)["leiden"]).astype(str)

    peaks = merged_peak_list(adata)
    pmat = snap.pp.make_peak_matrix(adata, use_rep=peaks)
    pmat = pmat.to_memory() if hasattr(pmat, "to_memory") else pmat
    pmat.obs["leiden"] = leiden
    print(f"\npeak matrix: {pmat.shape[0]} cells x {pmat.shape[1]:,} peaks")

    # differential accessibility: marker peaks per cluster
    mr = snap.tl.marker_regions(pmat, groupby="leiden", pvalue=0.01)
    print("\nmarker (differentially-accessible) peaks per cluster:")
    for c in sorted(mr, key=int):
        print(f"  cluster {c:>2s}: {len(mr[c]):5d} DARs")
    print("  (hundreds-to-thousands of interpretable regions per cluster, vs ~200k noisy raw tiles)")

    if motifs:
        print("\nrunning TF-motif enrichment on each cluster's marker peaks "
              "(downloads hg38 FASTA + cis-BP motifs on first use) ...")
        motif_db = snap.datasets.cis_bp(unique=True)
        enr = snap.tl.motif_enrichment(motifs=motif_db, regions=mr,
                                       genome_fasta=snap.genome.hg38, method="hypergeometric")
        print("\ntop enriched TF motifs per cluster (marker peaks):")
        for c in sorted(enr, key=int):
            df = enr[c]
            col = "name" if "name" in df.columns else df.columns[0]
            sort_col = "adjusted p-value" if "adjusted p-value" in df.columns else (
                "p-value" if "p-value" in df.columns else df.columns[-1])
            try:
                top = df.sort(sort_col).head(top_n)[col].to_list()
            except Exception:
                top = df.head(top_n)[col].to_list()
            print(f"  cluster {c:>2s}: {', '.join(str(x) for x in top)}")
        print("\n(enriched TFs name the regulators of each state — here CEBPA/AP-1 in monocytes,\n"
              " LEF1/TCF7/RUNX1 in T cells, EOMES/T-box in NK, EBF1/POU2F in B cells — off the chromatin.)")

    print("\nTakeaway: peaks (not tiles) are the regulatory units; call them per cluster so cell-type\n"
          "elements survive, then differential accessibility + motif enrichment turn open chromatin\n"
          "into a mechanistic story (which regions, driven by which TFs) — the analysis RNA can't give.")
    return mr


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--motifs", action="store_true", help="also run TF-motif enrichment (heavy downloads)")
    a = ap.parse_args()
    run(motifs=a.motifs)
