# scVelo Pipeline

[![GEO](https://img.shields.io/badge/NCBI%20GEO-GSE249285-blue)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249285)

RNA velocity analysis pipeline using scVelo, implementing **both** analysis
modes from the official tutorials:

1. **RNA Velocity Basics** — stochastic / steady-state model
2. **Dynamical Modeling** — likelihood-based EM with latent time

Supports Seurat object integration.  Tested on macOS with Apple Silicon
(M1/M2/M4 Pro).

## Demo Data

The **LungOrganoid** demo dataset is included in this repository under `LungOrganoid/results/` (output figures, tables, and HTML report). The raw sequencing data and processed input files (`.h5ad`, `.rda`) are deposited at **NCBI GEO** with accession number **[GSE249285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE249285)**.

To reproduce the demo analysis, download the input files from GEO and place them in the `LungOrganoid/` directory:
- `LOBLMday14WTA_scVelo.h5ad` — h5ad with spliced/unspliced layers
- `LOBLMday14WTA_Seurat.rda` — Seurat object with metadata and embeddings

## Directory Structure

```
scvelo_pipeline/
├── README.md
├── report.qmd                      # Branded Quarto report template
├── references.bib                  # BibTeX references (Bergen, La Manno, Wolf)
├── scripts/
│   ├── run_scvelo_pipeline.py      # Main pipeline (stochastic + dynamical)
│   ├── build_anndata_from_seurat.R # R helper: export Seurat metadata/embeddings
│   └── requirements.txt            # Python dependencies
└── <case_dir>/                     # One folder per case (e.g. Yamada_CD8/)
    ├── *.h5ad                      # Input: h5ad with spliced/unspliced layers
    ├── *_Seurat.rda / .qs          # Input: Seurat object (.rda or .qs)
    ├── spliced/                    # Input (optional): 10X spliced matrix
    ├── unspliced/                  # Input (optional): 10X unspliced matrix
    └── results/                    # Output directory
        ├── <case>_scVelo.html      # Branded Quarto HTML report (self-contained)
        ├── figures/
        │   ├── stochastic/         # Part 1 plots (QC + stochastic velocity)
        │   └── dynamical/          # Part 2 plots (dynamical modeling)
        ├── tables/                 # TSV and styled HTML tables
        └── data/                   # h5ad checkpoints and final files
```

## Environment Setup

### Create conda environment

```bash
conda create -n scvelo_env python=3.9 -y
conda activate scvelo_env
pip install -r scripts/requirements.txt
```

### Installed dependencies

| Package      | Version | Purpose                        |
|--------------|---------|--------------------------------|
| scvelo       | 0.3.3   | RNA velocity                   |
| scanpy       | 1.9.3   | Single-cell analysis           |
| anndata      | 0.8.0   | Annotated data objects         |
| numpy        | 1.24.4  | Numerical computing            |
| scipy        | 1.13.1  | Scientific computing           |
| pandas       | 2.3.3   | Data frames                    |
| scikit-learn | 1.6.1   | Machine learning               |
| matplotlib   | 3.9.4   | Plotting                       |
| igraph       | 0.11.9  | Graph analysis (PAGA)          |
| louvain      | 0.8.2   | Louvain clustering             |
| leidenalg    | 0.10.2  | Leiden clustering              |
| hnswlib      | 0.8.0   | Fast neighbor search           |
| numba        | 0.60.0  | JIT compilation                |
| umap-learn   | 0.5.11  | UMAP embeddings                |
| h5py         | 3.14.0  | HDF5 file I/O                  |
| loompy       | 3.0.8   | Loom file support              |
| beautifulsoup4 | —     | HTML table parsing (PAGA)      |

### Additional requirements

- **R 4.0+** with `Seurat`, `Matrix`, and (if using `.qs` files) `qs` packages; `Rscript` on PATH
- **Quarto ≥ 1.4** — used to render the branded HTML report
  (`report.qmd`).  Install from <https://quarto.org/docs/get-started/>.
  The base Python (used by Quarto internally) must have `pyyaml`,
  `nbformat`, and `nbclient` installed.

## Usage

The pipeline supports two input modes.  All commands should be run from the
`scvelo_pipeline/` directory with the conda environment activated.

### Mode A: h5ad + Seurat metadata (recommended)

Use when you have an h5ad file with spliced/unspliced layers and a Seurat object
with metadata (clusters, embeddings).  The Seurat file can be `.rda` or `.qs`.
For Seurat-based annotations, `scripts/build_anndata_from_seurat.R` loads the
Seurat object as `mBC` and exports `orig.ident` directly from
`mBC@meta.data$orig.ident`.

```bash
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 \
       VECLIB_MAXIMUM_THREADS=6 NUMEXPR_MAX_THREADS=6 \
       NUMBA_DISABLE_JIT=1 MPLCONFIGDIR=/tmp/mpl_scvelo
conda activate scvelo_env

python scripts/run_scvelo_pipeline.py \
  --h5ad_file LungOrganoid/LOBLMday14WTA_scVelo.h5ad \
  --seurat_file LungOrganoid/LOBLMday14WTA_Seurat.rda \
  --cell_prefix LOBLMday14WTA_ \
  --meta_cols orig.ident,TagIDs,seurat_clusters \
  --case_dir LungOrganoid \
  --n_pcs 30 --n_neighbors 30 --n_jobs 6
```

### Mode B: Preprocessed h5ad (optional)

```bash
python scripts/run_scvelo_pipeline.py \
  --preprocessed your_data.h5ad \
  --case_dir <case> \
  --n_pcs 30 --n_neighbors 30 --n_jobs 6
```

## Pipeline Workflow

The pipeline runs **two complete velocity analyses** sequentially, then
automatically renames all output files to report-matching names and renders
a branded Quarto HTML report.

### Preprocessing

| Step | Function | Details |
|------|----------|---------|
| Gene filtering | Species-specific removal | Human: mitochondrial `MT-*`, ribosomal `RPL/RPS/RPLP`, pseudogene-like `RPL*P*`, `RPS*P*`, `-ps`; Mouse: mitochondrial `mt-*`, ribosomal `Rpl/Rps/Rplp`, pseudogene-like `-ps`; for both: remove features annotated as `pseudogene` |
| Filter & normalize | `scv.pp.filter_and_normalize(min_shared_counts=20, n_top_genes=2000)` | Selects top 2,000 highly variable genes (matches official scVelo tutorial) |
| Neighbors & moments | `scv.pp.moments(n_pcs=30, n_neighbors=30)` | Computes first/second order moments and the neighborhood graph used downstream |

### QC (before both parts)

| Step | Function | Output |
|------|----------|--------|
| Spliced/unspliced proportions | `scv.pl.proportions()` | `figures/stochastic/fig-01-proportions.png` |

### Part 1 — Stochastic Velocity (RNA Velocity Basics)

| Step | Function | Output |
|------|----------|--------|
| 1a. Stochastic velocity | `scv.tl.velocity()` + `velocity_graph()` | Velocity layer |
| 1b. Stream plots | `velocity_embedding_stream()` | `figures/stochastic/fig-{02-05}-stoch-stream-{basis}-{cvar}.png` |
| 1b′. Arrow plots | `velocity_embedding()` | `figures/stochastic/stoch-arrows-{basis}-{cvar}.png` |
| 1c. Phase portraits | `scv.pl.velocity()` | `figures/stochastic/fig-06-phase-portraits.png` |
| 1d. Top genes per cluster | `rank_velocity_genes()` + `scv.pl.scatter()` | `figures/stochastic/stoch-top-genes-cluster-{NN}.png`, `tables/rank-velocity-genes.tsv` |
| 1e. Speed & coherence | `velocity_confidence()` | `figures/stochastic/fig-07-velocity-length-{basis}.png`, `fig-08-velocity-confidence-{basis}.png`, `fig-09-velocity-confidence-combined.png`, `tables/tbl-01-velocity-confidence-by-cluster.html` |
| 1f. Velocity graph | `scv.pl.velocity_graph()` | `figures/stochastic/fig-10-velocity-graph-{basis}.png` |
| 1g. Cell transitions | `scv.utils.get_cell_transitions()` | `figures/stochastic/fig-11-cell-transitions-{basis}.png` |
| 1h. Velocity pseudotime | `velocity_pseudotime()` | `figures/stochastic/fig-{12,13}-pseudotime-{basis}.png` |
| 1i. PAGA velocity graph | `scv.tl.paga()` + `scv.pl.paga()` | `figures/stochastic/fig-14-paga-{basis}.png`, `tables/tbl-02-paga-transitions.html` |
| Checkpoint | `adata.write()` | `data/adata_scvelo_stochastic.h5ad` |

### Part 2 — Dynamical Modeling

| Step | Function | Output |
|------|----------|--------|
| 2a. Recover dynamics | `scv.tl.recover_dynamics()` | Fitted kinetic params in `adata.var` |
| 2b. Dynamical velocity | `scv.tl.velocity(mode='dynamical')` + `velocity_graph()` | Velocity layer (overwritten) |
| 2c. Stream plots | `velocity_embedding_stream()` | `figures/dynamical/fig-{16-19}-dyn-stream-{basis}-{cvar}.png` |
| 2c′. Arrow plots | `velocity_embedding()` | `figures/dynamical/dyn-arrows-{basis}-{cvar}.png` |
| 2d. Kinetic rate histograms | Histograms of α, β, γ | `figures/dynamical/fig-15-kinetic-rates.png` |
| 2e. Latent time | `scv.tl.latent_time()` | `figures/dynamical/fig-{20,21}-latent-time-{basis}.png`, `fig-20{a-c}-latent-time-{basis}-{cvar}.png` |
| 2f. Top genes heatmap | `scv.pl.heatmap()` (top 150 genes) | `figures/dynamical/fig-25-heatmap-top-genes.png` |
| 2g. Top-likelihood portraits | `scv.pl.scatter(basis=top15)` | `figures/dynamical/fig-22-top-likelihood-genes.png` |
| 2h. Expression vs latent time | `scv.pl.scatter(x='latent_time')` | `figures/dynamical/fig-23-gene-expression-phase.png`, `fig-24-gene-expression-latent-time.png` |
| 2i. Dynamical genes / cluster | `rank_dynamical_genes()` + scatter | `figures/dynamical/fig-26-dyn-top-genes-{cluster}.png`, `tables/dynamical-genes-by-cluster.tsv` |
| Final save | `adata.write()` | `data/<case>_dynamical.h5ad` |

### Post-processing (automated)

| Step | Description |
|------|-------------|
| Rename outputs | All 60+ output files renamed to report-matching convention (`fig-NN-`, `tbl-NN-`, clean kebab-case) |
| Update report template | `report.qmd` auto-updated with case directory paths, human-readable label (e.g. "Yamada CD8"), and cell count from `adata.n_obs` |
| Render report | Quarto renders `report.qmd` → `<case>_scVelo.html` (self-contained branded HTML) |
| Move & clean up | Report moved to `results/`; legacy `report.html` removed |

## Custom Colour Palette

All categorical plots use a 104-colour palette combining four sets:

| Source                    | Count | Description                         |
|---------------------------|-------|-------------------------------------|
| `colors_dutch`            | 20    | Flat UI Dutch palette               |
| `colors_spanish`          | 20    | Flat UI Spanish palette             |
| `colors_custom`           | 15    | Paired/qualitative research palette |
| `pal_igv("default")(51)`  | 49*   | IGV genomics viewer palette         |

\* 49 unique after de-duplication (104 total).

The palette is defined as `CUSTOM_PALETTE` in `run_scvelo_pipeline.py` and
assigned to every categorical variable (`seurat_clusters`, `orig.ident`,
`TagIDs`, etc.) via `adata.uns['{var}_colors']` before plotting.

## Multi-Variable Plotting

Stream, arrow, and latent-time plots are generated for **each combination** of
embedding × colour variable:

| Embedding | Colour variables                      | Plot types               |
|-----------|---------------------------------------|--------------------------|
| UMAP      | seurat_clusters, orig.ident, TagIDs   | stream, arrows, latent time |
| FItSNE    | seurat_clusters, orig.ident, TagIDs   | stream, arrows, latent time |

Both stochastic and dynamical results get their own set (prefixed
`stoch-` and `dyn-` respectively).

TagIDs plots use `legend_loc="right margin"` to keep labels readable;
all other categorical variables use on-data annotation.

## Quarto Report

The pipeline produces a branded self-contained HTML report via Quarto
(`report.qmd`).  Features:

- **ImmunoGene Teqs. Inc** branding with teal colour scheme and company logo
- 20 scientific description boxes and 15 "How to Read" interpretation guides
- All figures and tables embedded inline with numbered cross-references
- Mathematical formulas for RNA velocity ODEs and EM model
- Interactive table of contents (left sidebar)
- BibTeX citations (Bergen 2020, La Manno 2018, Wolf 2018, Wolf 2019)
- Cluster tables numerically sorted (0, 1, 2, …, 10 — not lexicographic)
- Session information section with Python/R runtime and package versions

The report renders automatically at the end of the pipeline.  To re-render
manually:

```bash
quarto render report.qmd --execute --output <case>_scVelo.html
mv <case>_scVelo.html <case>/results/
```

## All Command-Line Options

```
--h5ad_file       h5ad file with spliced/unspliced layers
--seurat_file     Seurat object file (.rda or .qs) for metadata and embeddings
                  (--seurat_rda is accepted as a backward-compatible alias)
--preprocessed    Preprocessed h5ad (standalone mode)
--cell_prefix     Prefix to add to h5ad barcodes for Seurat matching
--meta_cols       Comma-separated metadata columns (default: orig.ident,TagIDs,seurat_clusters)
--case_dir        Case directory; outputs go to <case_dir>/results/ organised into
                  figures/, tables/, and data/ subdirectories
--outdir          Output directory (alternative to --case_dir)
--n_pcs           Number of principal components (default: 30)
--n_neighbors     Number of neighbors for kNN/moments graph (default: 30)
--n_jobs          Parallel jobs for dynamics recovery (default: 6)
```

## Output Files

All outputs are written under `<case_dir>/results/`, organised into
subdirectories:

### `<case>_scVelo.html`

Self-contained branded Quarto report with all figures, tables, scientific
descriptions, and interpretation guides.  Open in any web browser.

### `data/`

| File | Description |
|------|-------------|
| `<case>_dynamical.h5ad` | Final AnnData (dynamical velocity, latent time, gene rankings) |
| `adata_scvelo_stochastic.h5ad` | Checkpoint after Part 1 (stochastic velocity) |

### `tables/`

| File | Description |
|------|-------------|
| `rank-velocity-genes.tsv` | Stochastic velocity gene ranking by cluster |
| `dynamical-genes-by-cluster.tsv` | Dynamical gene ranking by cluster |
| `tbl-01-velocity-confidence-by-cluster.html` | Styled HTML table (coolwarm gradient) — mean velocity length & confidence per cluster |
| `tbl-02-paga-transitions.html` | Styled HTML table (Blues gradient) — PAGA transition confidence matrix |

### `figures/stochastic/` — QC & Part 1

| File pattern | Description |
|--------------|-------------|
| `fig-01-proportions.png` | Spliced / unspliced / ambiguous proportions per cluster |
| `fig-{02-05}-stoch-stream-{basis}-{cvar}.png` | Velocity streamlines (report figures) |
| `stoch-stream-{basis}-{cvar}.png` | Velocity streamlines (supplementary) |
| `stoch-arrows-{basis}-{cvar}.png` | Velocity arrows |
| `fig-06-phase-portraits.png` | Phase portraits for top 8 velocity genes |
| `stoch-top-genes-cluster-{NN}.png` | Top 5 velocity genes per cluster (coloured) |
| `fig-07-velocity-length-{basis}.png` | Differentiation speed |
| `fig-08-velocity-confidence-{basis}.png` | Velocity coherence |
| `fig-09-velocity-confidence-combined.png` | Combined velocity length & confidence scatter |
| `fig-10-velocity-graph-{basis}.png` | Cell-cell velocity transition network |
| `fig-11-cell-transitions-{basis}.png` | Cell transition trajectories overlaid on velocity graph |
| `fig-{12,13}-pseudotime-{basis}.png` | Velocity pseudotime |
| `fig-14-paga-{basis}.png` | PAGA directed cluster graph (legend on right margin) |

### `figures/dynamical/` — Part 2

| File pattern | Description |
|--------------|-------------|
| `fig-{16-19}-dyn-stream-{basis}-{cvar}.png` | Velocity streamlines (report figures) |
| `dyn-stream-{basis}-{cvar}.png` | Velocity streamlines (supplementary) |
| `dyn-arrows-{basis}-{cvar}.png` | Velocity arrows |
| `fig-15-kinetic-rates.png` | Transcription / splicing / degradation rate histograms |
| `fig-{20,21}-latent-time-{basis}.png` | Latent time (continuous gnuplot colourmap) |
| `fig-20{a-c}-latent-time-{basis}-{cvar}.png` | Embedding coloured by cluster variable |
| `latent-time-{basis}-{cvar}.png` | Latent time supplementary (FItSNE embeddings) |
| `fig-22-top-likelihood-genes.png` | Phase portraits for 15 highest-likelihood genes |
| `fig-23-gene-expression-phase.png` | Top 4 driver genes — spliced vs unspliced |
| `fig-24-gene-expression-latent-time.png` | Top 4 driver genes — expression over latent time |
| `fig-25-heatmap-top-genes.png` | Top 150 genes sorted by latent time |
| `fig-26-dyn-top-genes-{cluster}.png` | Top 5 dynamical genes per cluster (coloured) |

Where `{basis}` = `umap` or `fitsne`, `{cvar}` = `seurat-clusters`,
`orig-ident`, or `TagIDs`, `{cluster}` = cluster name,
`{NN}` = zero-padded cluster number (00, 01, …, 10).

## Memory & Performance (24 GB RAM, Apple M4 Pro)

Set thread limits to avoid OOM on macOS:

```bash
export OMP_NUM_THREADS=6
export OPENBLAS_NUM_THREADS=6
export VECLIB_MAXIMUM_THREADS=6
export NUMEXPR_MAX_THREADS=6
```

| Dataset size | Estimated time | Recommended n_jobs |
|-------------|----------------|--------------------|
| < 5,000 cells | 5–15 min | 6 |
| 5,000–15,000 cells | 15–60 min | 6 |
| > 15,000 cells | 1–3 hours | 4–6 |


## References

- scVelo: https://scvelo.readthedocs.io/
  - [RNA Velocity Basics](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html)
  - [Dynamical Modeling](https://scvelo.readthedocs.io/en/stable/DynamicalModeling.html)
- scanpy: https://scanpy.readthedocs.io/
- Seurat: https://satijalab.org/seurat/
- Quarto: https://quarto.org/
