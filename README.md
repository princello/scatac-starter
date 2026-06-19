# scatac-starter

A runnable, honest walkthrough of the **core single-cell ATAC-seq pipeline** on a real 10x dataset —
from raw fragments to cell types, differentially-accessible regions, and the transcription factors that
drive them — using [`snapATAC2`](https://kzhang.org/SnapATAC2/). Each script is one self-contained
step that prints its result and the methodological lesson behind it.

The point isn't the dataset; it's the **decisions that are specific to ATAC** and trip up people coming
from scRNA: QC on TSS enrichment instead of genes, a genome-tile feature space instead of a count
matrix, a spectral embedding instead of PCA, gene-*activity* scores to recover cell types, and peak
calling + motif enrichment to read regulation straight off the chromatin.

## Dataset

10x Genomics **5k PBMC scATAC-seq** (`snapatac2.datasets.pbmc5k`), aligned to hg38. snapATAC2 ships it
as a **fragment file** (one BED-like row per sequenced fragment) — there is no gene-count matrix to
start from; the whole analysis is built from fragments. The first run downloads it (~1 GB) to the
snapATAC2 cache.

PBMC is the right harness for a sanity check: we know the answer should be T cells, monocytes, B cells,
and NK cells, with monocyte- and lymphocyte-specific transcription factors driving their open chromatin.

## The four steps

| Script | Question | Method |
|---|---|---|
| `run_qc.py` | Which barcodes are real cells, not empty droplets? | fragment depth **and** TSS enrichment gates; fragment-size / nucleosome banding |
| `run_clustering.py` | What are the cell states? | genome **tile** matrix (500bp) → select features → **spectral** embedding → kNN → Leiden → UMAP |
| `run_gene_activity.py` | What cell *types* are those clusters? | gene-activity scores (accessibility over gene body+promoter) → marker-gene scoring per cluster |
| `run_peaks.py` | Which regions, driven by which TFs? | MACS3 peaks per cluster → merge → peak matrix → marker peaks (DARs) → TF-motif enrichment |

## What each step actually shows (verified output)

**1. QC — depth alone lies.** scATAC has no genes/counts/MT% to threshold. The two ATAC gates are
fragment depth and **TSS enrichment** (how concentrated a barcode's reads are at promoters). Of
**34,384 barcodes → 5,173 real cells (~15% kept; ~85% empty droplets/debris)**. The key point: a depth-
only gate keeps 5,543 barcodes, but **370 of those are deep yet have low TSSe** (ambient DNA / debris) —
the TSSe gate is what removes them. The fragment-size distribution shows clean nucleosome banding (41%
nucleosome-free <100bp, 35% mono-nucleosome), itself a library-quality readout. *Lesson: gate on depth
AND TSS enrichment, not genes/MT%.*

**2. Clustering — tiles + spectral, not genes + PCA.** The feature matrix is the genome binned into
**6,062,095 500bp tiles** — sparse and near-binary (a tile is open or not in a cell). You select the
~250,000 most-accessible tiles, then use a **spectral embedding** (a graph-Laplacian / LSI-style method
for sparse binary data), *not* PCA on log-counts. kNN → Leiden recovers **14 clusters** (131–799 cells).
*Lesson: the ATAC feature space and dimensionality reduction are fundamentally different from RNA's.*

**3. Gene activity — bridging chromatin to cell types.** ATAC clusters are open-chromatin states, not
labels. Summing accessibility over each gene's body+promoter gives a **gene-activity** matrix (5,173 ×
60,606 genes), a proxy for expression. Scoring canonical PBMC markers per cluster recovers the expected
composition — **5 Monocyte, 5 T-cell, 2 B-cell, 2 NK-cell** clusters (the largest cluster, 799 cells, is
monocytes). *Lesson: gene-activity scores connect a purely epigenomic readout to the marker-gene
vocabulary every RNA atlas speaks; it's a proxy, so for fine states you'd transfer labels from paired
scRNA.*

**4. Peaks — the ATAC-specific payoff.** Tiles are blunt; the real regulatory units are **peaks**. Call
peaks **per cluster** with MACS3 (so cell-type-specific elements survive — 992,321 peaks total, median
~66k/cluster), **merge** into **186,498** non-overlapping reference peaks, build a cell × peak matrix,
and find **marker peaks** = differentially-accessible regions (DARs): now 505–11,680 *interpretable*
regions per cluster, versus the ~200k noise you'd get from raw tiles. With `--motifs`, TF-motif
enrichment on each cluster's marker peaks names the regulators — and they line up with the cell types:

| cluster lineage | top enriched motifs | known role |
|---|---|---|
| Monocyte | **CEBPA**, FOS, JUNB, BACH | CEBPA = master myeloid TF; AP-1 |
| T cell | **LEF1, TCF7, TCF7L1/2, RUNX1** | canonical (naive) T-cell program |
| NK / effector | **EOMES, T-box (TBR1/TBX)** | EOMES = cytotoxic-lineage master |
| B cell | **EBF1, POU2F, POU2AF1** | EBF1 = master B-cell TF |

*Lesson: peak calling + differential accessibility + motif enrichment turn open chromatin into a
mechanistic, TF-level story — the analysis scRNA can't give you.*

The steps are deliberately consistent: the lineages from gene activity (step 3) are exactly the groups
whose marker peaks and motifs you interpret in step 4 — and where they disagree, the motifs win. The one
cluster gene activity couldn't confidently call (near-zero marker z) is enriched for **TAL1, LYL1, MYF5**
— hematopoietic-progenitor/erythroid factors — i.e. not a mature lymphoid type at all, which is exactly
why the lineage-marker scores were flat. A coherent pipeline tells one story.

## Setup

The snapATAC2/scanpy stack is happiest on **Python 3.11**.

```bash
cd scatac-starter
python3.11 -m venv .venv        # or: conda create -p .venv python=3.11
.venv/bin/pip install -r requirements.txt
```

## Run

```bash
.venv/bin/python run_qc.py
.venv/bin/python run_clustering.py      # first run builds + caches the processed file (~2-3 min)
.venv/bin/python run_gene_activity.py
.venv/bin/python run_peaks.py           # first run calls peaks with MACS3 (~8 min), then caches them
.venv/bin/python run_peaks.py --motifs  # optional: TF-motif enrichment (downloads hg38 FASTA + cis-BP)
```

The first `run_clustering.py` triggers the ~1 GB fragment download and caches the processed dataset to
`./data/` (gitignored); `run_peaks.py` caches the merged peak set so re-runs are instant.

## Layout

```
src/data.py            # load fragments; full pipeline (QC -> tiles -> spectral -> Leiden -> UMAP); cached
run_qc.py
run_clustering.py
run_gene_activity.py
run_peaks.py
requirements.txt
```

## Notes / honest caveats

- **MACS3 runs single-process here (`n_jobs=1`)**: on macOS the multi-process peak caller crashes
  (forked workers can't share the backed HDF5 file), so peak calling takes ~8 min. On Linux you can
  raise `n_jobs`. The merged peaks are cached either way.
- The 10x fragment file is coordinate-sorted, so import uses `sorted_by_barcode=False` (snapATAC2 panics
  otherwise) — a small but real gotcha.
- Gene activity is a proxy for expression, not RNA; it's good for broad lineages, not fine states.
- This is the snapATAC2 path; ArchR (R) covers the same workflow and is the other common toolkit.

## Related projects

Part of a small single-cell / computational-biology portfolio — each a runnable, honestly-documented pipeline:

- **[perturbseq-starter](https://github.com/princello/perturbseq-starter)** — Perturb-seq / CRISPR-screen analysis: guide assignment → Mixscape → pseudobulk-vs-per-cell DE → E-distance.
- **[scatac-starter](https://github.com/princello/scatac-starter)** — single-cell ATAC: TSS-enrichment QC → tile/spectral clustering → gene activity → MACS3 peaks + TF-motif enrichment. *(this repo)*
- **[scrna-workflow](https://github.com/princello/scrna-workflow)** — the same multi-sample scRNA pipeline as a reproducible DAG in both Snakemake and Nextflow.
