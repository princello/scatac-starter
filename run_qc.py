"""QC: telling real cells from empty droplets — with ATAC metrics, not RNA ones.

scATAC QC is different from scRNA QC. There are no genes/counts/MT% to threshold; instead the two
canonical gates are:
  - number of unique fragments per barcode (sequencing depth / a real cell vs an empty droplet), and
  - TSS enrichment (TSSe): how concentrated a barcode's reads are at transcription start sites
    relative to background. Real cells have open, accessible promoters -> high TSSe; ambient/debris
    barcodes have flat, genome-uniform signal -> low TSSe.
Crucially, depth ALONE isn't enough: a barcode can have many fragments but low TSSe (ambient DNA /
debris). You need both gates. We also show the fragment-size distribution, whose nucleosome banding
(a nucleosome-free peak <100bp + a mono-nucleosome bump ~200bp) is itself a library-quality readout.

  python run_qc.py
"""
from __future__ import annotations
import numpy as np
from src import data


def run():
    # import ALL barcodes (no cell filter yet) so we can see the cell-vs-debris separation
    adata = data.import_with_qc(min_num_fragments=100, file=None)
    obs = adata.obs
    n_frag = obs["n_fragment"].to_numpy()
    tsse = obs["tsse"].to_numpy()
    n0 = adata.n_obs
    print(f"barcodes imported (>=100 fragments): {n0}")
    print(f"  fragments/barcode: median {np.median(n_frag):.0f}, 90th pct {np.percentile(n_frag,90):.0f}, "
          f"max {n_frag.max():.0f}")
    print(f"  TSS enrichment   : median {np.median(tsse):.1f}, range {tsse.min():.1f}-{tsse.max():.1f}")

    deep = n_frag >= data.MIN_COUNTS
    passing = deep & (tsse >= data.MIN_TSSE)
    # the key point: depth alone admits debris that the TSSe gate removes
    deep_lowtsse = int((deep & (tsse < data.MIN_TSSE)).sum())
    print(f"\nGate 1 — depth only (>= {data.MIN_COUNTS} fragments)     : {int(deep.sum())} barcodes")
    print(f"Gate 2 — depth AND TSSe (>= {data.MIN_TSSE})             : {int(passing.sum())} barcodes")
    print(f"  -> {deep_lowtsse} deep barcodes are REJECTED by the TSSe gate (high depth, low promoter\n"
          f"     enrichment = ambient/debris). Depth alone would have kept them.")
    print(f"\nfinal: {n0} barcodes -> {int(passing.sum())} real cells "
          f"({passing.mean():.0%} kept; ~{1-passing.mean():.0%} are empty droplets / debris)")

    # fragment-size distribution: nucleosome banding as a library-quality readout
    try:
        import snapatac2 as snap
        fsd = snap.metrics.frag_size_distr(adata, max_recorded_size=1000, inplace=False)
        fsd = np.asarray(fsd, dtype=float)
        nfr = fsd[:100].sum() / fsd.sum()          # nucleosome-free (<100 bp)
        mono = fsd[147:294].sum() / fsd.sum()      # mono-nucleosome band (~1 nucleosome + linker)
        peak = int(np.argmax(fsd[10:]) + 10)
        print(f"\nfragment-size distribution: nucleosome-free (<100bp) {nfr:.0%}, "
              f"mono-nucleosome (147-294bp) {mono:.0%}, modal size {peak}bp")
        print("  (a clear sub-100bp peak + a ~200bp shoulder = well-resolved nucleosome banding = good library)")
    except Exception as e:
        print(f"\n(fragment-size distribution skipped: {e})")

    print("\nTakeaway: ATAC QC gates on fragment depth AND TSS enrichment, not genes/MT%. TSSe is what\n"
          "separates real open-chromatin cells from deep-but-flat ambient barcodes — depth alone lies.")
    return int(passing.sum())


if __name__ == "__main__":
    run()
