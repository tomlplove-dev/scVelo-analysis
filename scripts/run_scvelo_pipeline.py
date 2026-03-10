#!/usr/bin/env python3
"""
scVelo Pipeline -- Stochastic (Basics) + Dynamical Modeling
===========================================================
Implements both RNA Velocity analyses from the scVelo tutorials:
  Part 1. RNA Velocity Basics  (stochastic / steady-state)
  Part 2. Dynamical Modeling   (likelihood-based EM)

Two input modes:
  A. --h5ad_file + --seurat_rda   (h5ad layers + Seurat metadata)
  B. --preprocessed               (ready-made h5ad)
"""
import argparse, os, sys, warnings
from pathlib import Path
import re
import subprocess
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Enforce high-resolution output for all saved figures.
SAVE_DPI = 300


def p(*a, **kw):
    print(*a, flush=True, **kw)


def ensure_dir(d):
    Path(d).mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Custom colour palette (104 unique colours)
# Replicates: unique(c(colors_dutch, colors_spanish, colors_custom,
#                       pal_igv("default")(51)))
# ---------------------------------------------------------------------------
CUSTOM_PALETTE = [
    # --- colors_dutch (20) ---
    "#FFC312", "#C4E538", "#12CBC4", "#FDA7DF", "#ED4C67",
    "#F79F1F", "#A3CB38", "#1289A7", "#D980FA", "#B53471",
    "#EE5A24", "#009432", "#0652DD", "#9980FA", "#833471",
    "#EA2027", "#006266", "#1B1464", "#5758BB", "#6F1E51",
    # --- colors_spanish (20) ---
    "#40407A", "#706FD3", "#F7F1E3", "#34ACE0", "#33D9B2",
    "#2C2C54", "#474787", "#AAA69D", "#227093", "#218C74",
    "#FF5252", "#FF793F", "#D1CCC0", "#FFB142", "#FFDA79",
    "#B33939", "#CD6133", "#84817A", "#CC8E35", "#CCAE62",
    # --- colors_custom (15) ---
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
    "#FFFF99", "#B15928", "#49BEAA", "#611C35", "#2708A0",
    # --- pal_igv("default")(51) unique additions ---
    "#5050FF", "#CE3D32", "#749B58", "#F0E685", "#466983", "#BA6338",
    "#5DB1DD", "#802268", "#6BD76B", "#D595A7", "#924822", "#837B8D",
    "#C75127", "#D58F5C", "#7A65A5", "#E4AF69", "#3B1B53", "#CDDEB7",
    "#612A79", "#AE1F63", "#E7C76F", "#5A655E", "#CC9900", "#99CC00",
    "#A9A9A9", "#33CC00", "#00CC33", "#00CC99",
    "#0099CC", "#0A47FF", "#4775FF", "#FFC20A", "#FFD147", "#990033",
    "#991A00", "#996600", "#809900", "#339900", "#00991A", "#009966",
    "#008099", "#003399", "#1A0099", "#660099", "#990080", "#D60047",
    "#FF1463", "#00D68F", "#14FFB1",
]


def _assign_palette(adata, color_vars):
    """Assign CUSTOM_PALETTE colours to every categorical colour variable."""
    import numpy as np
    for cvar in color_vars:
        if cvar not in adata.obs.columns:
            continue
        if not hasattr(adata.obs[cvar].dtype, "categories"):
            continue
        n = len(adata.obs[cvar].cat.categories)
        pal = [CUSTOM_PALETTE[i % len(CUSTOM_PALETTE)] for i in range(n)]
        adata.uns[f"{cvar}_colors"] = np.array(pal)
        p(f"  Assigned {n} colours -> {cvar}")


def _ensure_clusters_key(adata, groupby):
    """Create/standardize `clusters` and keep its palette aligned to groupby."""
    import numpy as np
    import pandas as pd

    if groupby and groupby in adata.obs.columns:
        if not hasattr(adata.obs[groupby].dtype, "categories"):
            adata.obs[groupby] = adata.obs[groupby].astype(str).astype("category")

        src = adata.obs[groupby]
        src_categories = [str(x) for x in src.cat.categories]
        adata.obs["clusters"] = pd.Categorical(
            src.astype(str),
            categories=src_categories,
            ordered=bool(src.cat.ordered),
        )

        src_colors_key = f"{groupby}_colors"
        if src_colors_key in adata.uns and len(adata.uns[src_colors_key]) == len(src_categories):
            adata.uns["clusters_colors"] = np.array(list(adata.uns[src_colors_key]))
            p(f"  Synced clusters palette from {groupby}")
        else:
            pal = [CUSTOM_PALETTE[i % len(CUSTOM_PALETTE)] for i in range(len(src_categories))]
            adata.uns["clusters_colors"] = np.array(pal)
            p("  Assigned fallback palette to clusters")
        return

    if "clusters" in adata.obs.columns:
        if not hasattr(adata.obs["clusters"].dtype, "categories"):
            adata.obs["clusters"] = adata.obs["clusters"].astype(str).astype("category")
        if "clusters_colors" not in adata.uns:
            n = len(adata.obs["clusters"].cat.categories)
            pal = [CUSTOM_PALETTE[i % len(CUSTOM_PALETTE)] for i in range(n)]
            adata.uns["clusters_colors"] = np.array(pal)
            p("  Assigned palette to existing clusters key")


def _bases(adata):
    """Return list of (basis_key, display_name) for available embeddings."""
    return [(k, n) for k, n in [("umap", "UMAP"), ("FItSNE", "FItSNE")]
            if f"X_{k}" in adata.obsm]


def _save_figure(fig, path, dpi=SAVE_DPI, bbox_inches="tight"):
    """Save figure as PNG only (report output; 300 DPI minimum)."""
    out = Path(path).with_suffix(".png")
    fig.savefig(out, dpi=max(dpi, SAVE_DPI), bbox_inches=bbox_inches)
    return str(out)


def _save_current_figure(path, dpi=SAVE_DPI, bbox_inches="tight"):
    """Save current matplotlib figure as PNG."""
    import matplotlib.pyplot as plt

    return _save_figure(plt.gcf(), path, dpi=dpi, bbox_inches=bbox_inches)


def _sanitize_filename_token(text):
    """Return a filesystem-safe token for dynamic figure names."""
    token = str(text).strip()
    token = token.replace("|", "-").replace(" ", "-")
    token = re.sub(r"[^A-Za-z0-9._-]+", "-", token)
    token = re.sub(r"-{2,}", "-", token).strip("-_.")
    return token or "NA"


def _load_orig_ident_levels(levels_path):
    """Load orig.ident levels (in Seurat order) from a text file."""
    path = Path(levels_path)
    if not path.exists():
        return None
    levels = [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]
    return levels if levels else None


def _apply_orig_ident_order(adata, expected_order=None):
    """Set orig.ident categories using Seurat-exported label order."""
    import pandas as pd

    if "orig.ident" not in adata.obs.columns:
        return
    vals = adata.obs["orig.ident"].astype(str)
    present = set(vals.unique())

    if expected_order:
        expected = [str(x) for x in expected_order]
        ordered = [x for x in expected if x in present]
        extras = sorted(present - set(expected))
        if extras:
            p(f"  NOTE: orig.ident has {len(extras)} unexpected labels: {extras}")
            ordered.extend(extras)
        missing = [x for x in expected if x not in present]
        if missing:
            p(f"  NOTE: orig.ident missing expected labels: {missing}")
    else:
        # Keep order of appearance if no external order file is available.
        ordered = vals.drop_duplicates().tolist()

    adata.obs["orig.ident"] = pd.Categorical(vals, categories=ordered, ordered=True)
    p(f"  Set orig.ident category order ({len(ordered)} labels)")


def _format_gene_panel_figure(fig=None):
    """Improve readability of dense multi-panel gene figures."""
    import matplotlib.pyplot as plt

    fig = fig or plt.gcf()
    fig.patch.set_facecolor("white")
    n_axes = len(fig.axes)
    # Some scVelo panel plots can inflate canvas size unexpectedly.
    # Clamp to readable report-friendly dimensions before saving.
    if n_axes >= 12:
        fig.set_size_inches(24, 18, forward=True)
    elif n_axes >= 8:
        fig.set_size_inches(22, 14, forward=True)
    else:
        fig.set_size_inches(20, 10, forward=True)
    title_size = 20 if n_axes >= 10 else 22
    for ax in fig.axes:
        ax.set_facecolor("white")
        title = ax.get_title()
        if title:
            ax.set_title(title, fontsize=title_size, fontweight="bold")
        ax.xaxis.label.set_size(15)
        ax.yaxis.label.set_size(15)
        ax.tick_params(labelsize=12)
        leg = ax.get_legend()
        if leg is not None:
            for txt in leg.get_texts():
                txt.set_fontsize(11)

    # Compact layout without runaway tight-bbox expansion in large grids.
    try:
        fig.tight_layout()
    except Exception:
        pass


def _save_gene_panel_figure(path, dpi=SAVE_DPI):
    """Save current multi-panel gene figure with readability tuning."""
    import matplotlib.pyplot as plt

    fig = plt.gcf()
    _format_gene_panel_figure(fig)
    return _save_figure(fig, path, dpi=dpi, bbox_inches=None)


def _format_heatmap_figure(fig=None):
    """Improve readability for large latent-time heatmaps."""
    import matplotlib.pyplot as plt

    fig = fig or plt.gcf()
    fig.patch.set_facecolor("white")
    for ax in fig.axes:
        ax.set_facecolor("white")
        ax.tick_params(axis="x", labelsize=12)
        ax.tick_params(axis="y", labelsize=7)
    try:
        fig.tight_layout()
    except Exception:
        pass


def _stream_arrows(adata, bases, cvars, outdir, tag, images):
    """Generate stream + arrow velocity plots for every basis x colour var."""
    import matplotlib.pyplot as plt
    import scvelo as scv
    for bk, bn in bases:
        for cv in cvars:
            s = f"{bk}_{cv}"
            for kind, func, kw in [
                ("stream", scv.pl.velocity_embedding_stream,
                 dict(size=50, alpha=0.3)),
                ("arrows", scv.pl.velocity_embedding,
                 dict(arrow_length=3, arrow_size=3)),
            ]:
                try:
                    path = str(outdir / f"{tag}_{kind}_{s}.png")
                    fig, ax = plt.subplots(figsize=(12, 9))
                    leg = "right margin" if cv == "TagIDs" else "on data"
                    func(adata, basis=bk, color=cv, legend_loc=leg, ax=ax, show=False, **kw)
                    ax.set_title(f"{tag.title()} {kind} ({bn}) -- {cv}", fontsize=20, fontweight="bold")
                    png_path = _save_figure(fig, path, dpi=SAVE_DPI, bbox_inches="tight")
                    plt.close(fig)
                    images.append((f"{tag.title()} {kind} ({bn}) -- {cv}", png_path))
                    p(f"    saved {tag} {kind}: {s}")
                except Exception as exc:
                    p(f"    skip {tag} {kind} {s}: {exc}")



# ---------------------------------------------------------------------------
# Post-pipeline: rename outputs to match report figure numbers + render report
# ---------------------------------------------------------------------------
_RENAME_MAP = {
    # ═══ STOCHASTIC — figures in report ═══
    "figures/stochastic/proportions.png":
        "figures/stochastic/fig-01-proportions.png",
    "figures/stochastic/stochastic_stream_umap_seurat_clusters.png":
        "figures/stochastic/fig-02-stoch-stream-umap-seurat-clusters.png",
    "figures/stochastic/stochastic_stream_umap_TagIDs.png":
        "figures/stochastic/fig-03-stoch-stream-umap-TagIDs.png",
    "figures/stochastic/stochastic_stream_FItSNE_seurat_clusters.png":
        "figures/stochastic/fig-04-stoch-stream-fitsne-seurat-clusters.png",
    "figures/stochastic/stochastic_stream_FItSNE_TagIDs.png":
        "figures/stochastic/fig-05-stoch-stream-fitsne-TagIDs.png",
    "figures/stochastic/stochastic_phase_portraits.png":
        "figures/stochastic/fig-06-phase-portraits.png",
    "figures/stochastic/velocity_length_umap.png":
        "figures/stochastic/fig-07-velocity-length-umap.png",
    "figures/stochastic/velocity_confidence_umap.png":
        "figures/stochastic/fig-08-velocity-confidence-umap.png",
    "figures/stochastic/velocity_confidence_combined.png":
        "figures/stochastic/fig-09-velocity-confidence-combined.png",
    "figures/stochastic/velocity_graph_umap.png":
        "figures/stochastic/fig-10-velocity-graph-umap.png",
    "figures/stochastic/cell_transitions_umap.png":
        "figures/stochastic/fig-11-cell-transitions-umap.png",
    "figures/stochastic/velocity_pseudotime_umap.png":
        "figures/stochastic/fig-12-pseudotime-umap.png",
    "figures/stochastic/velocity_pseudotime_FItSNE.png":
        "figures/stochastic/fig-13-pseudotime-fitsne.png",
    "figures/stochastic/paga_umap.png":
        "figures/stochastic/fig-14-paga-umap.png",
    # ═══ STOCHASTIC — supplementary (not in report) ═══
    "figures/stochastic/stochastic_stream_umap_orig.ident.png":
        "figures/stochastic/stoch-stream-umap-orig-ident.png",
    "figures/stochastic/stochastic_stream_FItSNE_orig.ident.png":
        "figures/stochastic/stoch-stream-fitsne-orig-ident.png",
    "figures/stochastic/stochastic_arrows_umap_seurat_clusters.png":
        "figures/stochastic/stoch-arrows-umap-seurat-clusters.png",
    "figures/stochastic/stochastic_arrows_umap_orig.ident.png":
        "figures/stochastic/stoch-arrows-umap-orig-ident.png",
    "figures/stochastic/stochastic_arrows_umap_TagIDs.png":
        "figures/stochastic/stoch-arrows-umap-TagIDs.png",
    "figures/stochastic/stochastic_arrows_FItSNE_seurat_clusters.png":
        "figures/stochastic/stoch-arrows-fitsne-seurat-clusters.png",
    "figures/stochastic/stochastic_arrows_FItSNE_orig.ident.png":
        "figures/stochastic/stoch-arrows-fitsne-orig-ident.png",
    "figures/stochastic/stochastic_arrows_FItSNE_TagIDs.png":
        "figures/stochastic/stoch-arrows-fitsne-TagIDs.png",
    "figures/stochastic/velocity_length_FItSNE.png":
        "figures/stochastic/velocity-length-fitsne.png",
    "figures/stochastic/velocity_confidence_FItSNE.png":
        "figures/stochastic/velocity-confidence-fitsne.png",
    "figures/stochastic/velocity_graph_FItSNE.png":
        "figures/stochastic/velocity-graph-fitsne.png",
    "figures/stochastic/cell_transitions_FItSNE.png":
        "figures/stochastic/cell-transitions-fitsne.png",
    # ═══ DYNAMICAL — figures in report ═══
    "figures/dynamical/kinetic_rates.png":
        "figures/dynamical/fig-15-kinetic-rates.png",
    "figures/dynamical/dynamical_stream_umap_seurat_clusters.png":
        "figures/dynamical/fig-16-dyn-stream-umap-seurat-clusters.png",
    "figures/dynamical/dynamical_stream_umap_TagIDs.png":
        "figures/dynamical/fig-17-dyn-stream-umap-TagIDs.png",
    "figures/dynamical/dynamical_stream_FItSNE_seurat_clusters.png":
        "figures/dynamical/fig-18-dyn-stream-fitsne-seurat-clusters.png",
    "figures/dynamical/dynamical_stream_FItSNE_TagIDs.png":
        "figures/dynamical/fig-19-dyn-stream-fitsne-TagIDs.png",
    "figures/dynamical/latent_time_umap.png":
        "figures/dynamical/fig-20-latent-time-umap.png",
    "figures/dynamical/latent_time_FItSNE.png":
        "figures/dynamical/fig-21-latent-time-fitsne.png",
    "figures/dynamical/latent_time_umap_seurat_clusters.png":
        "figures/dynamical/fig-20a-latent-time-umap-seurat-clusters.png",
    "figures/dynamical/latent_time_umap_orig.ident.png":
        "figures/dynamical/fig-20b-latent-time-umap-orig-ident.png",
    "figures/dynamical/latent_time_umap_TagIDs.png":
        "figures/dynamical/fig-20c-latent-time-umap-TagIDs.png",
    "figures/dynamical/top_likelihood_genes.png":
        "figures/dynamical/fig-22-top-likelihood-genes.png",
    "figures/dynamical/gene_expression_phase.png":
        "figures/dynamical/fig-23-gene-expression-phase.png",
    "figures/dynamical/gene_expression_latent_time.png":
        "figures/dynamical/fig-24-gene-expression-latent-time.png",
    "figures/dynamical/heatmap_top_genes.png":
        "figures/dynamical/fig-25-heatmap-top-genes.png",
    # ═══ DYNAMICAL — supplementary (not in report) ═══
    "figures/dynamical/dynamical_stream_umap_orig.ident.png":
        "figures/dynamical/dyn-stream-umap-orig-ident.png",
    "figures/dynamical/dynamical_stream_FItSNE_orig.ident.png":
        "figures/dynamical/dyn-stream-fitsne-orig-ident.png",
    "figures/dynamical/dynamical_arrows_umap_seurat_clusters.png":
        "figures/dynamical/dyn-arrows-umap-seurat-clusters.png",
    "figures/dynamical/dynamical_arrows_umap_orig.ident.png":
        "figures/dynamical/dyn-arrows-umap-orig-ident.png",
    "figures/dynamical/dynamical_arrows_umap_TagIDs.png":
        "figures/dynamical/dyn-arrows-umap-TagIDs.png",
    "figures/dynamical/dynamical_arrows_FItSNE_seurat_clusters.png":
        "figures/dynamical/dyn-arrows-fitsne-seurat-clusters.png",
    "figures/dynamical/dynamical_arrows_FItSNE_orig.ident.png":
        "figures/dynamical/dyn-arrows-fitsne-orig-ident.png",
    "figures/dynamical/dynamical_arrows_FItSNE_TagIDs.png":
        "figures/dynamical/dyn-arrows-fitsne-TagIDs.png",
    "figures/dynamical/latent_time_FItSNE_seurat_clusters.png":
        "figures/dynamical/latent-time-fitsne-seurat-clusters.png",
    "figures/dynamical/latent_time_FItSNE_orig.ident.png":
        "figures/dynamical/latent-time-fitsne-orig-ident.png",
    "figures/dynamical/latent_time_FItSNE_TagIDs.png":
        "figures/dynamical/latent-time-fitsne-TagIDs.png",
    # ═══ TABLES ═══
    "tables/paga_transitions.html":
        "tables/tbl-02-paga-transitions.html",
    "tables/velocity_confidence_by_cluster.html":
        "tables/tbl-01-velocity-confidence-by-cluster.html",
    "tables/rank_velocity_genes.tsv":
        "tables/rank-velocity-genes.tsv",
    "tables/dynamical_genes_by_cluster.tsv":
        "tables/dynamical-genes-by-cluster.tsv",
}


def _rename_if_exists(src, dst):
    """Rename file if present and return 1 when renamed."""
    src = Path(src)
    dst = Path(dst)
    if src.exists():
        src.rename(dst)
        return 1
    return 0


def rename_outputs(results_dir):
    """Rename pipeline outputs to report-matching names (fig-NN-…, tbl-NN-…).

    Also renames:
    stochastic_top_genes_{N}.png → stoch-top-genes-cluster-{NN}.png
    dynamical_top_genes_{name}.png → fig-26-dyn-top-genes-{name}.png
    """
    results = Path(results_dir)
    renamed = 0

    # Static rename map
    for old_rel, new_rel in _RENAME_MAP.items():
        renamed += _rename_if_exists(results / old_rel, results / new_rel)

    # Stochastic top genes per cluster (variable names)
    for f in sorted((results / "figures" / "stochastic").glob("stochastic_top_genes_*.png")):
        cluster = f.stem.replace("stochastic_top_genes_", "")
        # Zero-pad numeric clusters
        try:
            cluster_str = f"{int(cluster):02d}"
        except ValueError:
            cluster_str = cluster
        new_name = f"stoch-top-genes-cluster-{cluster_str}.png"
        renamed += _rename_if_exists(f, f.parent / new_name)

    # Dynamical top genes per cluster (variable names)
    for f in sorted((results / "figures" / "dynamical").glob("dynamical_top_genes_*.png")):
        cluster = f.stem.replace("dynamical_top_genes_", "")
        new_name = f"fig-26-dyn-top-genes-{cluster}.png"
        renamed += _rename_if_exists(f, f.parent / new_name)

    p(f"  Renamed {renamed} output files to report-matching names")
    return renamed


def render_quarto_report(pipeline_dir, results_dir, case_name="output",
                         n_cells=None):
    """Render the branded Quarto report and move it to results_dir.

    Produces {case_name}_scVelo.html via Quarto, moves it into results_dir,
    and removes the legacy report.html generated by generate_html().
    Automatically replaces hardcoded case directory paths *and* human-readable
    case labels (e.g. "Yamada CD8" ↔ "Yamada CD4") in report.qmd so the same
    template works for any case.
    """
    report_qmd = Path(pipeline_dir) / "report.qmd"
    if not report_qmd.exists():
        p("  report.qmd not found — skipping Quarto render")
        return

    # --- Update case-specific paths & text in report.qmd ------------------
    qmd_text = report_qmd.read_text()
    changed = False

    # 1) Replace directory-style case name (e.g. Yamada_CD8 → Yamada_CD4)
    m = re.search(r'RESULTS\s*=\s*Path\("([^"]+)/results"\)', qmd_text)
    old_case = m.group(1) if m else None
    if old_case and old_case != case_name:
        qmd_text = qmd_text.replace(old_case, case_name)
        p(f"  Updated paths: {old_case} → {case_name}")
        changed = True

    # 2) Replace human-readable case name (underscores → spaces).
    #    e.g. "Yamada CD8" → "Yamada CD4"
    old_label = old_case.replace("_", " ") if old_case else None
    new_label = case_name.replace("_", " ")
    if old_label and old_label != new_label and old_label in qmd_text:
        qmd_text = qmd_text.replace(old_label, new_label)
        p(f"  Updated labels: {old_label} → {new_label}")
        changed = True

    # 3) Replace hardcoded cell count with actual count
    if n_cells is not None:
        m2 = re.search(r'\(([\d,]+) cells\)', qmd_text)
        if m2:
            old_count = m2.group(1)
            new_count = f"{n_cells:,}"
            if old_count != new_count:
                qmd_text = qmd_text.replace(
                    f"({old_count} cells)", f"({new_count} cells)")
                p(f"  Updated cell count: {old_count} → {new_count}")
                changed = True

    if changed:
        report_qmd.write_text(qmd_text)

    # Find quarto executable
    quarto = None
    for candidate in ["/Applications/quarto/bin/quarto", "quarto"]:
        try:
            subprocess.run([candidate, "--version"],
                           capture_output=True, check=True)
            quarto = candidate
            break
        except (FileNotFoundError, subprocess.CalledProcessError):
            continue

    if quarto is None:
        p("  Quarto not found — skipping report render")
        return

    out_name = f"{case_name}_scVelo.html"
    p(f"  Rendering Quarto report with {quarto} → {out_name} ...")
    result = subprocess.run(
        [quarto, "render", str(report_qmd), "--execute",
         "--output", out_name],
        cwd=str(pipeline_dir),
        capture_output=True, text=True
    )
    if result.returncode == 0:
        # Move rendered report into results directory
        src = Path(pipeline_dir) / out_name
        dst = Path(results_dir) / out_name
        if src.exists():
            src.rename(dst)
            p(f"  Quarto report → {dst}")
        # Remove legacy report.html
        legacy = Path(results_dir) / "report.html"
        if legacy.exists():
            legacy.unlink()
            p("  Removed legacy report.html")
    else:
        p("  Quarto render failed (exit code", result.returncode, ")")
        if result.stderr:
            for line in result.stderr.strip().split("\n")[-5:]:
                p("    ", line)


def generate_html(outdir, images, tables):
    """Write self-contained HTML gallery report."""
    import pandas as pd
    h = ['<html><head><meta charset="utf-8"><title>scVelo Report</title></head><body>',
         "<h1>scVelo Pipeline Report</h1>"]
    for desc, img in images:
        h.append(f"<h2>{desc}</h2>")
        img_path = Path(img)
        rel = os.path.relpath(str(img_path), str(outdir))
        h.append(
            f'<img src="{rel}" style="max-width:100%;height:auto;">'
        )
    for desc, tab in tables:
        h.append(f"<h2>{desc}</h2>")
        try:
            if str(tab).endswith(".html"):
                h.append("<div style='overflow-x:auto;'>")
                h.append(Path(tab).read_text())
                h.append("</div>")
            else:
                df = pd.read_csv(tab, sep="\t", index_col=0)
                h.append("<div style='overflow-x:auto;'><pre>")
                h.append(df.head(5).to_csv(sep="\t"))
                h.append("</pre></div>")
        except Exception:
            h.append(f"<p>{os.path.basename(str(tab))}</p>")
    h.append("</body></html>")
    out = Path(outdir) / "report.html"
    out.write_text("\n".join(h))
    p("Wrote report:", out)


# ======================================================================== #
#  MAIN                                                                      #
# ======================================================================== #
def main():
    import scvelo as scv
    import scanpy as sc
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    scv.settings.verbosity = 3
    scv.settings.presenter_view = True
    scv.settings.set_figure_params("scvelo", dpi_save=SAVE_DPI)
    plt.rcParams.update({
        "figure.figsize": (12, 9),
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.edgecolor": "white",
        "text.color": "black",
        "axes.labelcolor": "black",
        "axes.edgecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "savefig.dpi": SAVE_DPI,
        "font.size": 14,
        "axes.titlesize": 18,
        "axes.labelsize": 15,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
    })

    # -- CLI ---------------------------------------------------------------
    ap = argparse.ArgumentParser(
        description="scVelo: stochastic + dynamical RNA velocity pipeline")
    ap.add_argument("--preprocessed", help="Pre-processed AnnData (.h5ad)")
    ap.add_argument("--h5ad_file", help="h5ad with spliced/unspliced layers")
    ap.add_argument("--seurat_file", "--seurat_rda",
                        help="Seurat object file (.rda or .qs) for metadata")
    ap.add_argument("--cell_prefix", default="", help="Barcode prefix")
    ap.add_argument("--meta_cols",
                    default="orig.ident,TagIDs,seurat_clusters",
                    help="Comma-separated metadata columns")
    ap.add_argument("--outdir", help="Output directory (legacy)")
    ap.add_argument("--case_dir", help="Case directory")
    ap.add_argument("--n_pcs", type=int, default=30)
    ap.add_argument("--n_neighbors", type=int, default=30)
    ap.add_argument("--n_jobs", type=int, default=6)
    args = ap.parse_args()

    outdir = (Path(args.case_dir) / "results" if args.case_dir
              else Path(args.outdir) if args.outdir else None)
    if outdir is None:
        p("ERROR: supply --outdir or --case_dir"); sys.exit(1)
    ensure_dir(outdir)

    # -- Create output subdirectories
    fig_stoch = outdir / "figures" / "stochastic"
    fig_dyn   = outdir / "figures" / "dynamical"
    tab_dir   = outdir / "tables"
    data_dir  = outdir / "data"
    for d in [fig_stoch, fig_dyn, tab_dir, data_dir]:
        ensure_dir(d)

    adata = None
    color_by = None      # primary cluster key
    orig_ident_levels = None

    # === MODE A: h5ad + Seurat metadata ===================================
    if args.h5ad_file and args.seurat_file:
        p("=== Mode A: h5ad + Seurat metadata ===")
        h5ad_path = Path(args.h5ad_file)
        seurat_path = Path(args.seurat_file)
        meta_dir = outdir / "preprocessed_meta"
        ensure_dir(meta_dir)
        r_script = Path(__file__).parent / "build_anndata_from_seurat.R"
        cmd = ["Rscript", str(r_script), str(seurat_path),
               str(meta_dir), "none", "none"]
        p("Exporting Seurat metadata:", " ".join(cmd))
        res = subprocess.run(cmd)
        if res.returncode != 0:
            p("Rscript failed with code", res.returncode); sys.exit(1)
        p("Reading h5ad:", h5ad_path)
        adata = sc.read_h5ad(str(h5ad_path))
        p("  h5ad shape:", adata.shape)
        p("  h5ad layers:", list(adata.layers.keys()))
        prefix = args.cell_prefix
        if prefix:
            p("  Adding prefix to barcodes:", prefix)
            adata.obs_names = [prefix + bc for bc in adata.obs_names]
        obs_path = meta_dir / "obs.tsv"
        obs_df = pd.read_csv(obs_path, sep="\t", index_col=0)
        p("  Seurat metadata shape:", obs_df.shape)
        orig_ident_levels = _load_orig_ident_levels(meta_dir / "orig_ident_levels.tsv")
        if orig_ident_levels:
            p(f"  Loaded orig.ident levels from Seurat ({len(orig_ident_levels)} labels)")
        common = sorted(set(adata.obs_names) & set(obs_df.index))
        p("  Common cells:", len(common), "/", adata.n_obs,
          "(h5ad) /", len(obs_df), "(Seurat)")
        if len(common) == 0:
            p("ERROR: No matching cells. Check --cell_prefix.")
            p("  h5ad barcode sample:", adata.obs_names[:3].tolist())
            p("  Seurat barcode sample:", obs_df.index[:3].tolist())
            sys.exit(1)
        adata = adata[common].copy()
        meta_cols = [c.strip() for c in args.meta_cols.split(",") if c.strip()]
        # Always keep core Seurat annotations from mBC@meta.data.
        # `orig.ident` is the canonical annotation source requested by user.
        for required in ["orig.ident", "seurat_clusters"]:
            if required not in meta_cols:
                meta_cols.append(required)
        for col in meta_cols:
            if col in obs_df.columns:
                adata.obs[col] = obs_df.loc[adata.obs_names, col].values
                p(f"  Added metadata: {col}"
                  f" ({adata.obs[col].nunique()} unique)")
            else:
                p(f"  WARNING: {col} not in Seurat metadata")
        if "seurat_clusters" in adata.obs.columns:
            adata.obs["seurat_clusters"] = (
                adata.obs["seurat_clusters"].astype(str).astype("category"))
            color_by = "seurat_clusters"
        for emb_file, obsm_key in [("emb_FItSNE.tsv", "X_FItSNE"),
                                    ("emb_umap.tsv", "X_umap"),
                                    ("emb_pca.tsv", "X_pca")]:
            emb_path = meta_dir / emb_file
            if emb_path.exists():
                emb_df = pd.read_csv(emb_path, sep="\t", index_col=0)
                matched = [c for c in adata.obs_names if c in emb_df.index]
                if len(matched) > 0:
                    adata.obsm[obsm_key] = emb_df.loc[adata.obs_names].values
                    p(f"  Loaded embedding: {obsm_key}"
                      f" ({emb_df.shape[1]} dims)")

    # === MODE B: preprocessed h5ad ========================================
    elif args.preprocessed:
        p("=== Mode B: preprocessed h5ad ===")
        adata = sc.read_h5ad(args.preprocessed)
        if "seurat_clusters" in adata.obs.columns:
            if not hasattr(adata.obs["seurat_clusters"].dtype, "categories"):
                adata.obs["seurat_clusters"] = (
                    adata.obs["seurat_clusters"].astype(str).astype("category"))
            color_by = "seurat_clusters"
            p("Using seurat_clusters as primary grouping variable")
    else:
        p("ERROR: need --h5ad_file+--seurat_rda, or --preprocessed")
        sys.exit(1)

    p("\nAnnData:", adata.shape)
    p("Layers:", list(adata.layers.keys()))
    p("Obs columns:", list(adata.obs.columns))

    if "spliced" in adata.layers:
        adata.X = adata.layers["spliced"].copy()
        p("Set X from spliced layer")

    # -- Remove mitochondrial, ribosomal & pseudogenes ---------------------
    # Gene name conventions differ by species:
    #   Human: MT-CO1, RPL30, RPS27
    #   Mouse: mt-Co1, Rpl30, Rps27
    p("\n=== Gene filtering (MT / ribosomal / pseudogenes) ===")
    n_before = adata.n_vars
    vn = adata.var_names
    # Mitochondrial genes (human MT-, mouse mt-)
    mito = vn.str.startswith("MT-") | vn.str.startswith("mt-")
    # Ribosomal protein genes (human RPL/RPS, mouse Rpl/Rps)
    ribo = (vn.str.match(r"^RPL\d")  | vn.str.match(r"^RPS\d")  |
            vn.str.match(r"^Rpl\d")  | vn.str.match(r"^Rps\d")  |
            vn.str.match(r"^RPLP\d") | vn.str.match(r"^Rplp\d"))
    # Pseudogenes from gene symbols and optional feature annotations
    pseudo = (vn.str.match(r"(?i)^RPL\d+P\d+$") |
              vn.str.match(r"(?i)^RPS\d+P\d+$") |
              vn.str.match(r"(?i).*-ps\d*$"))
    for bc in ["gene_biotype", "biotype", "gene_type", "feature_biotype"]:
        if bc in adata.var.columns:
            pseudo = pseudo | adata.var[bc].astype(str).str.contains(
                "pseudogene", case=False, na=False).to_numpy()
    keep = ~(mito | ribo | pseudo)
    adata = adata[:, keep].copy()
    p(f"  Removed {mito.sum()} mitochondrial, {ribo.sum()} ribosomal, "
      f"{pseudo.sum()} pseudogenes")
    p(f"  Genes: {n_before} -> {adata.n_vars}")
    _apply_orig_ident_order(adata, expected_order=orig_ident_levels)

    # -- Colour variable setup ---------------------------------------------
    color_vars = []
    if color_by and color_by in adata.obs.columns:
        color_vars.append(color_by)
    for ex in ["orig.ident", "TagIDs"]:
        if ex in adata.obs.columns and ex not in color_vars:
            if not hasattr(adata.obs[ex].dtype, "categories"):
                adata.obs[ex] = adata.obs[ex].astype(str).astype("category")
            color_vars.append(ex)
    p("Colour variables:", color_vars)
    _assign_palette(adata, color_vars)
    groupby = color_by if color_by and color_by in adata.obs.columns else None
    # scVelo tutorial convention: use `clusters` as the grouping key.
    # Keep palette synchronized so cluster colours match across all plots.
    _ensure_clusters_key(adata, groupby)

    images, tables = [], []

    # ================================================================== #
    #  QC -- spliced / unspliced proportions                              #
    # ================================================================== #
    try:
        p("\n=== QC: spliced / unspliced proportions ===")
        out = str(fig_stoch / "proportions.png")
        plt.close("all")
        scv.pl.proportions(
            adata, groupby=groupby, figsize=(14, 8), show=False
        )
        out_png = _save_current_figure(
            out, dpi=SAVE_DPI, bbox_inches="tight"
        )
        plt.close("all")
        images.append(("Spliced / unspliced proportions", out_png))
        p("  Saved proportions")
    except Exception as exc:
        p("  proportions failed:", exc)

    # -- Preprocessing -----------------------------------------------------
    p("\n=== Preprocessing ===")
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

    p("\n=== Computing neighbors & moments ===")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    bases = _bases(adata)
    p("Embeddings:", [b[1] for b in bases])

    # ================================================================== #
    #  PART 1 -- Stochastic Velocity  (RNA Velocity Basics)               #
    # ================================================================== #
    p("\n" + "=" * 60)
    p("  PART 1 -- Stochastic Velocity (RNA Velocity Basics)")
    p("=" * 60)

    # 1a. Estimate stochastic velocity -------------------------------------
    p("\n--- 1a. Stochastic velocity ---")
    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)

    # 1b. Stream & arrow plots ---------------------------------------------
    p("\n--- 1b. Stream & arrow plots ---")
    _stream_arrows(adata, bases, color_vars, fig_stoch, "stochastic", images)

    # 1c. Phase portraits for top velocity genes ---------------------------
    p("\n--- 1c. Phase portraits (top velocity genes) ---")
    try:
        rank_groupby = "clusters" if "clusters" in adata.obs.columns else groupby
        if rank_groupby:
            scv.tl.rank_velocity_genes(
                adata, groupby=rank_groupby, min_corr=0.3, n_genes=200)
            df_rv = (scv.DataFrame(adata.uns["rank_velocity_genes"]["names"])
                     if hasattr(scv, "DataFrame")
                     else pd.DataFrame(adata.uns["rank_velocity_genes"]["names"]))
            tp = str(tab_dir / "rank_velocity_genes.tsv")
            df_rv.to_csv(tp, sep="\t", index=False)
            tables.append(("Ranked velocity genes (stochastic)", tp))
            p("  Saved rank_velocity_genes.tsv")
            # pick top 8 unique genes across clusters
            seen, top = set(), []
            for col in df_rv.columns:
                for g in df_rv[col]:
                    if g not in seen and g in adata.var_names:
                        top.append(g); seen.add(g)
                    if len(top) >= 8:
                        break
                if len(top) >= 8:
                    break
            if top:
                out = str(fig_stoch / "stochastic_phase_portraits.png")
                plt.close("all")
                scv.pl.velocity(
                    adata, top, ncols=2, color=rank_groupby,
                    figsize=(20, 14), show=False
                )
                out_png = _save_gene_panel_figure(out, dpi=SAVE_DPI)
                plt.close("all")
                images.append(("Stochastic phase portraits", out_png))
                p(f"  Saved phase portraits ({len(top)} genes)")
    except Exception as exc:
        p("  phase portraits failed:", exc)

    # 1d. Top velocity genes scatter per cluster ---------------------------
    p("\n--- 1d. Top velocity genes per cluster ---")
    try:
        rank_groupby = "clusters" if "clusters" in adata.obs.columns else groupby
        if rank_groupby and "rank_velocity_genes" in adata.uns:
            df_rv = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"])
            outline_labels = [
                x for x in ["Ngn3 high EP", "Pre-endocrine", "Beta"]
                if x in set(adata.obs[rank_groupby].astype(str))
            ]
            scatter_kwargs = dict(frameon=False, size=10, linewidth=1.5)
            if outline_labels:
                scatter_kwargs["add_outline"] = ", ".join(outline_labels)
            for cl in df_rv.columns[:10]:
                gl = [g for g in df_rv[cl][:5] if g in adata.var_names]
                if not gl:
                    continue
                cl_token = _sanitize_filename_token(cl)
                out = str(fig_stoch / f"stochastic_top_genes_{cl_token}.png")
                plt.close("all")
                scv.pl.scatter(
                    adata, gl, ylabel=str(cl), color=rank_groupby,
                    ncols=3, figsize=(20, 10), show=False, **scatter_kwargs
                )
                out_png = _save_gene_panel_figure(out, dpi=SAVE_DPI)
                plt.close("all")
                images.append((f"Top velocity genes -- {cl}", out_png))
            p("  Saved per-cluster gene scatters")
    except Exception as exc:
        p("  top gene scatter failed:", exc)

    # 1e. Cell cycle scoring -----------------------------------------------
    p("\n--- 1e. Cell cycle scoring ---")
    try:
        scv.tl.score_genes_cell_cycle(adata)
        for bk, bn in bases:
            for sk in ["S_score", "G2M_score"]:
                out = str(fig_stoch / f"cell_cycle_{sk}_{bk}.png")
                fig, ax = plt.subplots(figsize=(12, 9))
                scv.pl.scatter(adata, color=sk, color_map="coolwarm",
                               basis=bk, size=50, perc=[5, 95],
                               ax=ax, show=False)
                ax.set_title(f"{sk} ({bn})", fontsize=20, fontweight="bold")
                out_png = _save_figure(
                    fig, out, dpi=SAVE_DPI, bbox_inches="tight"
                )
                plt.close(fig)
                images.append((f"{sk} ({bn})", out_png))
        p("  Saved cell cycle plots")
    except Exception as exc:
        p("  cell cycle failed:", exc)

    # 1f. Speed & coherence ------------------------------------------------
    p("\n--- 1f. Speed & coherence ---")
    try:
        scv.tl.velocity_confidence(adata)
        for bk, bn in bases:
            for vk in ["velocity_length", "velocity_confidence"]:
                out = str(fig_stoch / f"{vk}_{bk}.png")
                fig, ax = plt.subplots(figsize=(12, 9))
                scv.pl.scatter(adata, color=vk, color_map="coolwarm",
                               basis=bk, size=50, perc=[5, 95],
                               ax=ax, show=False)
                ax.set_title(f"{vk} ({bn})", fontsize=20, fontweight="bold")
                out_png = _save_figure(
                    fig, out, dpi=SAVE_DPI, bbox_inches="tight"
                )
                plt.close(fig)
                images.append((f"{vk} ({bn})", out_png))
        if groupby:
            vkeys = [k for k in ["velocity_length", "velocity_confidence"]
                     if k in adata.obs]
            if vkeys:
                cdf = adata.obs.groupby(groupby)[vkeys].mean().T
                styled_html = cdf.style.background_gradient(cmap="coolwarm", axis=1).to_html()
                tp = str(tab_dir / "velocity_confidence_by_cluster.html")
                Path(tp).write_text(styled_html)
                tables.append(("Velocity confidence by cluster", tp))
        keys = [k for k in ["velocity_length", "velocity_confidence"]
                if k in adata.obs]
        if keys:
            out = str(fig_stoch / "velocity_confidence_combined.png")
            plt.close("all")
            scv.pl.scatter(adata, c=keys, cmap="coolwarm",
                           perc=[5, 95], ncols=1, show=False)
            out_png = _save_current_figure(
                out, dpi=SAVE_DPI, bbox_inches="tight"
            )
            plt.close("all")
            images.append(("Velocity length & confidence", out_png))
        p("  Saved speed & coherence")
    except Exception as exc:
        p("  speed/coherence failed:", exc)

    # 1g. Velocity graph visualisation -------------------------------------
    p("\n--- 1g. Velocity graph visualisation ---")
    try:
        for bk, bn in bases:
            out = str(fig_stoch / f"velocity_graph_{bk}.png")
            fig, ax = plt.subplots(figsize=(12, 9))
            scv.pl.velocity_graph(
                adata, basis=bk, color=groupby, threshold=0.1,
                ax=ax, show=False
            )
            out_png = _save_figure(
                fig, out, dpi=SAVE_DPI, bbox_inches="tight"
            )
            plt.close(fig)
            images.append((f"Velocity graph ({bn})", out_png))
        p("  Saved velocity graph")
    except Exception as exc:
        p("  velocity graph viz failed:", exc)
    # 1g2. Cell transitions ------------------------------------------------
    p("\n--- 1g2. Cell transitions ---")
    try:
        for bk, bn in bases:
            x, y = scv.utils.get_cell_transitions(adata, basis=bk, starting_cell=70)
            out = str(fig_stoch / f"cell_transitions_{bk}.png")
            fig, ax = plt.subplots(figsize=(12, 9))
            ax = scv.pl.velocity_graph(
                adata, basis=bk, c="lightgrey", edge_width=0.05,
                ax=ax, show=False
            )
            ax = scv.pl.scatter(adata, x=x, y=y, s=120, c="ascending", cmap="gnuplot", ax=ax, show=False)
            out_png = _save_figure(
                fig, out, dpi=SAVE_DPI, bbox_inches="tight"
            )
            plt.close(fig)
            images.append((f"Cell transitions ({bn})", out_png))
        p("  Saved cell transitions")
    except Exception as exc:
        p("  cell transitions failed:", exc)

    # 1h. Velocity pseudotime ----------------------------------------------
    p("\n--- 1h. Velocity pseudotime ---")
    try:
        scv.tl.velocity_pseudotime(adata)
        for bk, bn in bases:
            out = str(fig_stoch / f"velocity_pseudotime_{bk}.png")
            fig, ax = plt.subplots(figsize=(12, 9))
            scv.pl.scatter(adata, color="velocity_pseudotime",
                           color_map="gnuplot", basis=bk, size=50,
                           ax=ax, show=False)
            ax.set_title(f"Velocity pseudotime ({bn})", fontsize=20, fontweight="bold")
            out_png = _save_figure(
                fig, out, dpi=SAVE_DPI, bbox_inches="tight"
            )
            plt.close(fig)
            images.append((f"Velocity pseudotime ({bn})", out_png))
        p("  Saved velocity pseudotime")
    except Exception as exc:
        p("  velocity pseudotime failed:", exc)

    # 1i. PAGA velocity graph ----------------------------------------------
    p("\n--- 1i. PAGA velocity graph ---")
    try:
        if groupby:
            # workaround for older scvelo/scanpy compat
            if "neighbors" in adata.uns:
                adata.uns["neighbors"]["distances"] = \
                    adata.obsp["distances"]
                adata.uns["neighbors"]["connectivities"] = \
                    adata.obsp["connectivities"]
            scv.tl.paga(adata, groups=groupby)
            try:
                pdf = scv.get_df(
                    adata, "paga/transitions_confidence", precision=2).T
                styled_html = pdf.style.background_gradient(cmap="Blues").format("{:.2g}").to_html()
                tp = str(tab_dir / "paga_transitions.html")
                Path(tp).write_text(styled_html)
                tables.append(("PAGA transitions", tp))
            except Exception:
                pass
            for bk, bn in bases:
                out = str(fig_stoch / f"paga_{bk}.png")
                plt.close("all")
                scv.pl.paga(adata, basis=bk, size=50, alpha=0.1,
                            min_edge_width=2, node_size_scale=1.5,
                            legend_loc="right margin",
                            figsize=(12, 9), show=False)
                out_png = _save_current_figure(
                    out, dpi=SAVE_DPI, bbox_inches="tight"
                )
                plt.close("all")
                images.append((f"PAGA velocity ({bn})", out_png))
            p("  Saved PAGA")
    except Exception as exc:
        p("  PAGA failed:", exc)

    # Save stochastic checkpoint -------------------------------------------
    stoch_path = str(data_dir / "adata_scvelo_stochastic.h5ad")
    adata.write(stoch_path)
    p(f"\n  Saved stochastic h5ad: {stoch_path}")

    # ================================================================== #
    #  PART 2 -- Dynamical Modeling                                       #
    # ================================================================== #
    p("\n" + "=" * 60)
    p("  PART 2 -- Dynamical Modeling")
    p("=" * 60)

    # 2a. Recover dynamics -------------------------------------------------
    p("\n--- 2a. Recovering dynamics ---")
    scv.tl.recover_dynamics(adata, n_jobs=args.n_jobs)

    # 2b. Dynamical velocity -----------------------------------------------
    p("\n--- 2b. Dynamical velocity ---")
    scv.tl.velocity(adata, mode="dynamical")
    scv.tl.velocity_graph(adata)

    # 2c. Stream & arrow plots (dynamical) ---------------------------------
    p("\n--- 2c. Dynamical stream & arrow plots ---")
    _stream_arrows(adata, bases, color_vars, fig_dyn, "dynamical", images)

    # 2d. Kinetic rate parameter histograms --------------------------------
    p("\n--- 2d. Kinetic rate parameters ---")
    try:
        dv = adata.var
        mask = dv.get("fit_likelihood", pd.Series(dtype=float)) > 0.1
        if "velocity_genes" in dv.columns:
            mask = mask & dv["velocity_genes"]
        dvf = dv[mask]
        if len(dvf) > 0:
            fig, axes = plt.subplots(1, 3, figsize=(20, 6))
            for ax, col, lab in zip(
                axes,
                ["fit_alpha", "fit_beta", "fit_gamma"],
                ["Transcription rate", "Splicing rate",
                 "Degradation rate"],
            ):
                v = dvf[col].dropna()
                if col == "fit_beta" and "fit_scaling" in dvf.columns:
                    v = v * dvf.loc[v.index, "fit_scaling"]
                v = v[v > 0]
                if len(v):
                    ax.hist(v, bins=50, color="#474787", edgecolor="white")
                    ax.set_xscale("log")
                ax.set_xlabel(lab, fontsize=16)
                ax.set_ylabel("Count", fontsize=16)
                ax.tick_params(labelsize=13)
            fig.suptitle("Kinetic rate parameters", fontsize=22, fontweight="bold")
            fig.tight_layout(rect=[0, 0, 1, 0.95])
            out = str(fig_dyn / "kinetic_rates.png")
            out_png = _save_figure(
                fig, out, dpi=SAVE_DPI, bbox_inches="tight"
            )
            plt.close(fig)
            images.append(("Kinetic rate parameters", out_png))
            p("  Saved kinetic rate histograms")
    except Exception as exc:
        p("  kinetic rates failed:", exc)

    # 2e. Latent time ------------------------------------------------------
    p("\n--- 2e. Latent time ---")
    try:
        scv.tl.latent_time(adata)
        for bk, bn in bases:
            # continuous colour map
            out = str(fig_dyn / f"latent_time_{bk}.png")
            fig, ax = plt.subplots(figsize=(12, 9))
            scv.pl.scatter(adata, color="latent_time", color_map="gnuplot",
                           basis=bk, size=50, ax=ax, show=False)
            ax.set_title(f"Latent time ({bn})", fontsize=20, fontweight="bold")
            out_png = _save_figure(
                fig, out, dpi=SAVE_DPI, bbox_inches="tight"
            )
            plt.close(fig)
            images.append((f"Latent time ({bn})", out_png))
            # categorical overlays
            for cv in color_vars:
                out = str(fig_dyn / f"latent_time_{bk}_{cv}.png")
                fig, ax = plt.subplots(figsize=(12, 9))
                leg = "right margin" if cv == "TagIDs" else "on data"
                scv.pl.scatter(adata, color=cv, basis=bk, size=50,
                               legend_loc=leg, ax=ax, show=False)
                ax.set_title(f"Clusters ({bn}) -- {cv}", fontsize=20, fontweight="bold")
                out_png = _save_figure(
                    fig, out, dpi=SAVE_DPI, bbox_inches="tight"
                )
                plt.close(fig)
                images.append((f"Latent time view ({bn}) -- {cv}", out_png))
        p("  Saved latent time")
    except Exception as exc:
        p("  latent time failed:", exc)

    # 2f. Heatmap of top genes by latent time ------------------------------
    try:
        dyn_groupby = "clusters" if "clusters" in adata.obs.columns else groupby
        if "fit_likelihood" in adata.var.columns:
            # A smaller panel count is easier to read in the HTML report.
            t150 = (adata.var["fit_likelihood"]
                    .sort_values(ascending=False).index[:150])
            out = str(fig_dyn / "heatmap_top_genes.png")
            plt.close("all")
            scv.pl.heatmap(adata, var_names=t150, sortby="latent_time",
                           col_color=dyn_groupby, n_convolve=100,
                           figsize=(18, 12), show=False)
            _format_heatmap_figure()
            out_png = _save_current_figure(out, dpi=SAVE_DPI, bbox_inches=None)
            plt.close("all")
            images.append(("Top genes heatmap by latent time", out_png))
            p("  Saved heatmap")
    except Exception as exc:
        p("  heatmap failed:", exc)

    # 2g. Top-likelihood gene phase portraits ------------------------------
    p("\n--- 2g. Top-likelihood gene phase portraits ---")
    try:
        dyn_plot_group = "clusters" if "clusters" in adata.obs.columns else groupby
        if "fit_likelihood" in adata.var.columns:
            t15 = [g for g in adata.var["fit_likelihood"]
                   .sort_values(ascending=False).index[:15]
                   if g in adata.var_names]
            if t15 and dyn_plot_group:
                out = str(fig_dyn / "top_likelihood_genes.png")
                plt.close("all")
                scv.pl.scatter(adata, basis=t15,
                               color=dyn_plot_group, frameon=False,
                               ncols=3, figsize=(24, 18), show=False)
                out_png = _save_gene_panel_figure(out, dpi=SAVE_DPI)
                plt.close("all")
                images.append(
                    ("Top-likelihood gene phase portraits", out_png))
                p(f"  Saved {len(t15)} phase portraits")
            elif t15 and not dyn_plot_group:
                p("  Skipped top-likelihood portraits: no cluster/groupby key")
    except Exception as exc:
        p("  top-likelihood portraits failed:", exc)

    # 2h. Gene expression over latent time ---------------------------------
    p("\n--- 2h. Gene expression vs latent time ---")
    try:
        dyn_plot_group = "clusters" if "clusters" in adata.obs.columns else groupby
        if ("fit_likelihood" in adata.var.columns
                and "latent_time" in adata.obs.columns):
            t4 = [g for g in adata.var["fit_likelihood"]
                  .sort_values(ascending=False).index[:4]
                  if g in adata.var_names]
            if t4 and dyn_plot_group:
                out = str(fig_dyn / "gene_expression_phase.png")
                plt.close("all")
                scv.pl.scatter(adata, t4, color=dyn_plot_group,
                               frameon=False, ncols=2,
                               figsize=(20, 12), show=False)
                out_png = _save_gene_panel_figure(out, dpi=SAVE_DPI)
                plt.close("all")
                images.append(
                    ("Top gene phase portraits (dynamical)", out_png))
                out = str(fig_dyn / "gene_expression_latent_time.png")
                plt.close("all")
                scv.pl.scatter(adata, x="latent_time", y=t4,
                               color=dyn_plot_group, frameon=False,
                               ncols=2, figsize=(20, 12), show=False)
                out_png = _save_gene_panel_figure(out, dpi=SAVE_DPI)
                plt.close("all")
                images.append(
                    ("Gene expression over latent time", out_png))
                p("  Saved expression vs latent time")
            elif t4 and not dyn_plot_group:
                p("  Skipped expression-vs-time portraits: no cluster/groupby key")
    except Exception as exc:
        p("  expression vs latent time failed:", exc)

    # 2i. Rank dynamical genes + per-cluster scatter -----------------------
    p("\n--- 2i. Cluster-specific dynamical genes ---")
    try:
        dyn_groupby = "clusters" if "clusters" in adata.obs.columns else groupby
        if dyn_groupby:
            scv.tl.rank_dynamical_genes(
                adata, groupby=dyn_groupby, n_genes=200)
            ddf = scv.get_df(adata, "rank_dynamical_genes/names")
            tp = str(tab_dir / "dynamical_genes_by_cluster.tsv")
            ddf.to_csv(tp, sep="\t", index=False)
            tables.append(("Dynamical genes by cluster", tp))
            for cl in ddf.columns[:10]:
                gl = [g for g in ddf[cl][:5] if g in adata.var_names]
                if not gl:
                    continue
                cl_token = _sanitize_filename_token(cl)
                out = str(fig_dyn / f"dynamical_top_genes_{cl_token}.png")
                plt.close("all")
                scv.pl.scatter(adata, gl, ylabel=str(cl),
                               color=dyn_groupby, frameon=False,
                               ncols=3, figsize=(20, 10), show=False)
                out_png = _save_gene_panel_figure(out, dpi=SAVE_DPI)
                plt.close("all")
                images.append((f"Dynamical top genes -- {cl}", out_png))
            p("  Saved dynamical gene ranking + scatters")
    except Exception as exc:
        p("  dynamical gene ranking failed:", exc)

    # -- Save final h5ad + report ------------------------------------------
    cn = Path(args.case_dir).name if args.case_dir else "output"
    fp = str(data_dir / f"{cn}_dynamical.h5ad")
    adata.write(fp)
    p(f"\nSaved final h5ad: {fp}")
    generate_html(outdir, images, tables)

    # -- Rename outputs to report-matching names ---------------------------
    p("\n--- Renaming outputs to report-matching names ---")
    rename_outputs(outdir)

    # -- Render branded Quarto report --------------------------------------
    p("\n--- Rendering Quarto report ---")
    pipeline_dir = Path(__file__).resolve().parent.parent
    cn = Path(args.case_dir).name if args.case_dir else "output"
    render_quarto_report(pipeline_dir, results_dir=outdir, case_name=cn,
                          n_cells=adata.n_obs)

    p("\n=== All done. Outputs in", outdir, "===")


if __name__ == "__main__":
    main()
