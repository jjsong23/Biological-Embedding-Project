#!/usr/bin/env python3
"""Generate publication-quality figures for the manuscript results section."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib.patches import FancyArrowPatch

FONT_SM  = 8    # minor annotations, significance notes
FONT_MD  = 9    # axis labels, tick labels, legend
FONT_LG  = 10   # panel titles, suptitles

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.size":        FONT_MD,
    "axes.labelsize":   FONT_MD,
    "axes.titlesize":   FONT_LG,
    "xtick.labelsize":  FONT_MD,
    "ytick.labelsize":  FONT_MD,
    "legend.fontsize":  FONT_MD,
    "figure.dpi":       150,
    "axes.spines.top":  False,
    "axes.spines.right":False,
})

OUT    = Path("/home/jsong205/atac_out/figures_paper")
RES    = Path("/home/jsong205/atac_out/week3_results")
RES_RP = Path("/home/jsong205/atac_out/week3_results_rp")

CONDITIONS_ORDER = [
    "atacformer", "rp",
    "pubmedbert", "sapbert", "medcpt", "general",
    "cto_pubmedbert", "cto_sapbert", "cto_medcpt", "cto_general",
    "tfidf", "random",
]
DISPLAY = {
    "atacformer":    "Atacformer",
    "rp":            "RP",
    "pubmedbert":    "PubMedBERT\n(full)",
    "sapbert":       "SapBERT\n(full)",
    "medcpt":        "MedCPT\n(full)",
    "general":       "General\n(full)",
    "cto_pubmedbert":"PubMedBERT\n(name)",
    "cto_sapbert":   "SapBERT\n(name)",
    "cto_medcpt":    "MedCPT\n(name)",
    "cto_general":   "General\n(name)",
    "tfidf":         "TF-IDF",
    "random":        "Random",
}
COLORS = {
    "atacformer":    "#1D3557",
    "rp":            "#2A9D8F",
    "pubmedbert":    "#E63946",
    "sapbert":       "#F4A261",
    "medcpt":        "#457B9D",
    "general":       "#9B5DE5",
    "cto_pubmedbert":"#FF6B6B",
    "cto_sapbert":   "#FFB347",
    "cto_medcpt":    "#87CEEB",
    "cto_general":   "#C8A2C8",
    "tfidf":         "#6B8F71",
    "random":        "#AAAAAA",
}
K_VALUES = [5, 10, 20, 50]


# ── Figure 1: Study design schematic ──────────────────────────────────────────

def fig1_schematic():
    fig, ax = plt.subplots(figsize=(13, 2.8))
    ax.set_xlim(0, 13); ax.set_ylim(0, 2.8)
    ax.axis("off")

    stages = [
        ("Data",
         "32,245 ATAC-seq experiments\n1,108 cell types · 22 classes",
         "#DCEEFB", "#1A3A5C"),
        ("Embeddings",
         "Atacformer · RP\nPubMedBERT · SapBERT · MedCPT · General\n(full text & name-only) · TF-IDF · Random",
         "#C8E6C9", "#1A3A1A"),
        ("Similarity Matrices",
         "Corpus-centered cosine similarity\nCell-type agg. (1,108 × 1,108)\nFull dataset (32,245 × 32,245)",
         "#FFE8D6", "#5C2A00"),
        ("Comparison Metrics",
         "Spearman rank correlation\nMantel permutation test\nTop-k neighborhood overlap",
         "#EDE7F6", "#2A0060"),
    ]

    box_w = 2.8
    gap   = 0.4
    total = len(stages) * box_w + (len(stages) - 1) * gap
    x0    = (13 - total) / 2

    for i, (title, body, bg, fg) in enumerate(stages):
        x = x0 + i * (box_w + gap)
        rect = mpatches.FancyBboxPatch(
            (x, 0.25), box_w, 2.3,
            boxstyle="round,pad=0.1",
            facecolor=bg, edgecolor="#aaaaaa", linewidth=0.8)
        ax.add_patch(rect)
        ax.text(x + box_w/2, 2.25, title,
                ha="center", va="center", fontsize=FONT_MD,
                fontweight="bold", color=fg)
        ax.plot([x + 0.15, x + box_w - 0.15], [2.0, 2.0],
                color="#aaaaaa", lw=0.6)
        ax.text(x + box_w/2, 1.1, body,
                ha="center", va="center", fontsize=FONT_SM,
                color="#333333", linespacing=1.5)
        if i < len(stages) - 1:
            ax.annotate("",
                xy=(x + box_w + gap, 1.4),
                xytext=(x + box_w, 1.4),
                arrowprops=dict(arrowstyle="-|>", color="#555555", lw=1.4))

    plt.tight_layout(pad=0.2)
    plt.savefig(OUT / "fig1.pdf", bbox_inches="tight")
    plt.savefig(OUT / "fig1.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: fig1")


# ── Figure 2: Spearman ρ vs Atacformer ────────────────────────────────────────

def fig2_spearman():
    df_ct   = pd.read_csv(RES / "results_celltype.csv")
    df_full = pd.read_csv(RES / "results_full.csv")

    keys = [k for k in CONDITIONS_ORDER if k in df_ct["condition"].values]
    keys = [k for k in keys if k not in ("atacformer",)]

    x = np.arange(len(keys))
    w = 0.36

    fig, ax = plt.subplots(figsize=(13, 5))
    ct_idx   = df_ct.set_index("condition").reindex(keys)
    full_idx = df_full.set_index("condition").reindex(keys)
    bars_ct   = ct_idx["spearman_rho"].values
    bars_full = full_idx["spearman_rho"].values
    colors    = [COLORS[k] for k in keys]

    err_ct   = np.vstack([bars_ct   - ct_idx["rho_ci_low"].values,
                           ct_idx["rho_ci_high"].values   - bars_ct]).clip(min=0)
    err_full = np.vstack([bars_full - full_idx["rho_ci_low"].values,
                           full_idx["rho_ci_high"].values - bars_full]).clip(min=0)

    ax.bar(x - w/2, bars_ct,   w, color=colors, alpha=0.95, label="Cell-type aggregated",
           yerr=err_ct,   error_kw=dict(ecolor="#333333", elinewidth=1, capsize=3))
    ax.bar(x + w/2, bars_full, w, color=colors, alpha=0.5,  label="Full 32K experiments", hatch="//",
           yerr=err_full, error_kw=dict(ecolor="#333333", elinewidth=1, capsize=3))

    ax.set_xticks(x)
    ax.set_xticklabels([DISPLAY[k] for k in keys], fontsize=FONT_MD)
    ax.set_ylabel("Spearman ρ  (vs. Atacformer)")
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_ylim(-0.05, 0.55)
    ax.legend(frameon=False)

    sig_ct = ct_idx["mantel_pval"].values
    for i, (h_ct, p_ct, e_ct) in enumerate(zip(bars_ct, sig_ct, err_ct[1])):
        star = "***" if p_ct < 0.001 else "**" if p_ct < 0.01 else "*" if p_ct < 0.05 else ""
        if star:
            ax.text(x[i] - w/2, h_ct + e_ct + 0.01, star,
                    ha="center", va="bottom", fontsize=FONT_SM, color="#333333")

    ax.text(0.99, 0.98,
            "*** p < 0.001  ** p < 0.01  * p < 0.05\n"
            "Cell-type: Mantel permutation test\n"
            "Error bars: 95% bootstrap CI",
            transform=ax.transAxes, fontsize=FONT_SM, va="top", ha="right", color="#555555")

    plt.tight_layout()
    plt.savefig(OUT / "fig2.pdf", bbox_inches="tight")
    plt.savefig(OUT / "fig2.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: fig2")


# ── Figure 3: Top-k overlap — Atacformer vs RP side by side ───────────────────

def fig3_topk_comparison():
    df_ct_atac   = pd.read_csv(RES    / "results_celltype.csv")
    df_full_atac = pd.read_csv(RES    / "results_full.csv")
    df_ct_rp     = pd.read_csv(RES_RP / "results_celltype.csv")
    df_full_rp   = pd.read_csv(RES_RP / "results_full.csv")

    keys = [k for k in CONDITIONS_ORDER
            if k in df_full_atac["condition"].values and k != "atacformer"]

    all_vals = []
    for df in [df_full_atac, df_full_rp]:
        df_idx = df.set_index("condition")
        for key in keys:
            if key not in df_idx.index:
                continue
            for k in K_VALUES:
                col = f"topk_{k}"
                if col in df_idx.columns:
                    all_vals.append(df_idx.loc[key, col])
    ymax = max(all_vals) * 1.08

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)
    panels = [
        (axes[0], df_full_atac, "vs. Atacformer  (full 32K)"),
        (axes[1], df_full_rp,   "vs. Regulatory Potential  (full 32K)"),
    ]
    for ax, df, title in panels:
        df_idx = df.set_index("condition")
        for key in keys:
            if key not in df_idx.index:
                continue
            ys = [df_idx.loc[key, f"topk_{k}"] for k in K_VALUES
                  if f"topk_{k}" in df_idx.columns]
            ks = [k for k in K_VALUES if f"topk_{k}" in df_idx.columns]
            lw    = 2.2 if key in ("pubmedbert","sapbert","medcpt","general","rp") else 1.2
            alpha = 1.0 if key in ("pubmedbert","sapbert","medcpt","general","rp") else 0.5
            ax.plot(ks, ys, marker="o", label=DISPLAY[key],
                    color=COLORS[key], linewidth=lw, markersize=5, alpha=alpha)
        ax.set_xlabel("k")
        ax.set_ylabel("Top-k neighborhood overlap")
        ax.set_title(title, pad=8)
        ax.set_xticks(K_VALUES)
        ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
        ax.grid(axis="y", alpha=0.25, linestyle="--")
        ax.set_ylim(0, ymax)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=FONT_MD, loc="lower center",
               ncol=6, bbox_to_anchor=(0.5, -0.08), frameon=False)
    plt.tight_layout()
    plt.savefig(OUT / "fig3.pdf", bbox_inches="tight")
    plt.savefig(OUT / "fig3.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: fig3")

    # Individual panels for single-column placement
    for suffix, df, title in [
        ("a", df_full_atac, "vs. Atacformer  (full 32K)"),
        ("b", df_full_rp,   "vs. Regulatory Potential  (full 32K)"),
    ]:
        fig_p, ax_p = plt.subplots(figsize=(7, 5.5))
        df_idx = df.set_index("condition")
        for key in keys:
            if key not in df_idx.index:
                continue
            ys = [df_idx.loc[key, f"topk_{k}"] for k in K_VALUES
                  if f"topk_{k}" in df_idx.columns]
            ks = [k for k in K_VALUES if f"topk_{k}" in df_idx.columns]
            lw    = 2.2 if key in ("pubmedbert","sapbert","medcpt","general","rp") else 1.2
            alpha = 1.0 if key in ("pubmedbert","sapbert","medcpt","general","rp") else 0.5
            ax_p.plot(ks, ys, marker="o", label=DISPLAY[key],
                      color=COLORS[key], linewidth=lw, markersize=5, alpha=alpha)
        ax_p.set_xlabel("k")
        ax_p.set_ylabel("Top-k neighborhood overlap")
        ax_p.set_title(title, pad=8)
        ax_p.set_xticks(K_VALUES)
        ax_p.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
        ax_p.grid(axis="y", alpha=0.25, linestyle="--")
        ax_p.set_ylim(0, ymax)
        ax_p.legend(fontsize=FONT_MD, frameon=False, ncol=2)
        plt.tight_layout()
        plt.savefig(OUT / f"fig3{suffix}.pdf", bbox_inches="tight")
        plt.savefig(OUT / f"fig3{suffix}.png", dpi=300, bbox_inches="tight")
        plt.close()
        print(f"Saved: fig3{suffix}")


# ── Figure S1: Sanity check grid (UMAP + heatmap, all 12 conditions) ──────────

_SC_BASE = Path("/home/jsong205/atac_out")
_SC_CELL_TYPES = [
    "CD4+ T cells", "CD8+ T cells", "T cells", "CAR-T cells", "B cells",
    "Natural killer cells", "Monocytes", "Macrophages", "Neutrophils", "PBMC",
    "CD34+", "CD34 Hematopoietic stem cells", "Bone Marrow Cells",
    "ES cells", "iPS cells", "hESC H9",
    "Brain", "iPSC derived neural cells", "hESC derived neural cells",
    "Prefrontal Cortex", "Glioblastoma stem cells",
    "K-562", "Acute myeloid leukemia", "B-ALL", "Jurkat", "GM12878",
    "MCF-7", "A549", "HeLa",
    "HUVEC", "Fibroblasts", "Skin fibroblasts",
    "Pancreatic islets", "iPSC derived cardiac cells",
]
_SC_LINEAGE = {
    "CD4+ T cells": "immune adaptive", "CD8+ T cells": "immune adaptive",
    "T cells": "immune adaptive", "CAR-T cells": "immune adaptive",
    "B cells": "immune adaptive",
    "Natural killer cells": "immune innate", "Monocytes": "immune innate",
    "Macrophages": "immune innate", "Neutrophils": "immune innate",
    "PBMC": "immune innate",
    "CD34+": "hematopoietic", "CD34 Hematopoietic stem cells": "hematopoietic",
    "Bone Marrow Cells": "hematopoietic",
    "ES cells": "pluripotent", "iPS cells": "pluripotent", "hESC H9": "pluripotent",
    "Brain": "neural", "iPSC derived neural cells": "neural",
    "hESC derived neural cells": "neural", "Prefrontal Cortex": "neural",
    "Glioblastoma stem cells": "neural",
    "K-562": "blood cancer", "Acute myeloid leukemia": "blood cancer",
    "B-ALL": "blood cancer", "Jurkat": "blood cancer", "GM12878": "blood cancer",
    "MCF-7": "solid cancer", "A549": "solid cancer", "HeLa": "solid cancer",
    "HUVEC": "somatic", "Fibroblasts": "somatic", "Skin fibroblasts": "somatic",
    "Pancreatic islets": "somatic", "iPSC derived cardiac cells": "somatic",
}
_SC_PALETTE = {
    "immune adaptive": "#E63946",
    "immune innate":   "#F4A261",
    "hematopoietic":   "#2A9D8F",
    "pluripotent":     "#9B5DE5",
    "neural":          "#457B9D",
    "blood cancer":    "#C1121F",
    "solid cancer":    "#FF6B6B",
    "somatic":         "#6B8F71",
}
_SC_ABBREV = {
    "immune adaptive": "Im. Ad.",
    "immune innate":   "Im. In.",
    "hematopoietic":   "Hemato.",
    "pluripotent":     "Plurip.",
    "neural":          "Neural",
    "blood cancer":    "Bl. Ca.",
    "solid cancer":    "So. Ca.",
    "somatic":         "Somat.",
}
_SC_CONDITIONS = [
    ("Atacformer",         _SC_BASE / "atacformer_embeddings"),
    ("RP",                 _SC_BASE / "rp_embeddings"),
    ("PubMedBERT\n(full)", _SC_BASE / "text_embeddings/full/pubmedbert"),
    ("SapBERT\n(full)",    _SC_BASE / "text_embeddings/full/sapbert"),
    ("MedCPT\n(full)",     _SC_BASE / "text_embeddings/full/medcpt"),
    ("General\n(full)",    _SC_BASE / "text_embeddings/full/general"),
    ("PubMedBERT\n(name)", _SC_BASE / "text_embeddings/cell_type_only/pubmedbert"),
    ("SapBERT\n(name)",    _SC_BASE / "text_embeddings/cell_type_only/sapbert"),
    ("MedCPT\n(name)",     _SC_BASE / "text_embeddings/cell_type_only/medcpt"),
    ("General\n(name)",    _SC_BASE / "text_embeddings/cell_type_only/general"),
    ("TF-IDF",             _SC_BASE / "text_embeddings/tfidf"),
    ("Random",             _SC_BASE / "text_embeddings/random"),
]


def _sc_load(df, emb_dir):
    """Return (centered_embs, labels) for sentinel cell types, or None."""
    means, labels = [], []
    for ct in _SC_CELL_TYPES:
        ids = df.loc[df["cell_type"] == ct, "experiment_id"].tolist()
        vecs = [np.load(emb_dir / f"{eid}.npy")
                for eid in ids if (emb_dir / f"{eid}.npy").exists()]
        if vecs:
            means.append(np.mean(vecs, axis=0))
            labels.append(ct)
    if not means:
        return None
    embs = np.array(means, dtype=np.float32)
    paths = list(emb_dir.glob("*.npy"))[:5000]
    if not paths:
        return None
    cm = np.mean([np.load(p) for p in paths], axis=0).astype(np.float32)
    return embs - cm, labels


def figS1_sanity_grid():
    from sklearn.metrics.pairwise import cosine_similarity
    import umap as _umap

    df = pd.read_csv(_SC_BASE / "atac_metadata_final.csv")

    # Compute UMAP + similarity for all 12 conditions
    cond_data = {}
    for name, emb_dir in _SC_CONDITIONS:
        result = _sc_load(df, emb_dir)
        if result is None:
            continue
        centered, labels = result
        lineages = [_SC_LINEAGE[lb] for lb in labels]
        colors   = [_SC_PALETTE[ln] for ln in lineages]
        sim      = cosine_similarity(centered).astype(np.float32)
        reducer  = _umap.UMAP(n_neighbors=8, min_dist=0.3, metric="cosine",
                              random_state=42)
        coords   = reducer.fit_transform(centered)
        cond_data[name] = (coords, sim, labels, colors, lineages)

    # One shared legend used by both grids
    shared_patches = [mpatches.Patch(color=c, label=lb)
                      for lb, c in _SC_PALETTE.items()]

    ncols, nrows = 4, 3

    # ── UMAP grid: clean dots, no per-panel legend, no cell-type annotations ──
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 2.5))
    axes = axes.flatten()
    for idx, (name, (coords, sim, labels, colors, lineages)) in enumerate(cond_data.items()):
        ax = axes[idx]
        ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=20,
                   edgecolors="white", linewidths=0.3)
        ax.set_title(name, fontsize=FONT_MD, fontweight="bold", pad=3)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP 1", fontsize=FONT_SM)
        ax.set_ylabel("UMAP 2", fontsize=FONT_SM)
    for idx in range(len(cond_data), len(axes)):
        axes[idx].set_visible(False)
    fig.suptitle("Sanity check: UMAP (34 cell types, corpus-centered)",
                 fontsize=FONT_LG, fontweight="bold")
    fig.legend(handles=shared_patches, fontsize=FONT_SM, loc="lower center",
               ncol=4, bbox_to_anchor=(0.5, 0), frameon=False)
    plt.tight_layout(rect=[0, 0.10, 1, 0.97])
    plt.savefig(OUT / "figS1a.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved: figS1a_sanity_umaps.png")

    # ── Heatmap grid: tissue-color strip on y-axis, no cell-type tick labels ──
    strip_w = 3.0  # data-unit width of the tissue color strip

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3, nrows * 2.7))
    axes = axes.flatten()
    for idx, (name, (coords, sim, labels, colors, lineages)) in enumerate(cond_data.items()):
        ax = axes[idx]
        n = len(labels)
        ax.imshow(sim, cmap="RdYlBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_title(name, fontsize=FONT_MD, fontweight="bold", pad=3)
        ax.set_xticks([]); ax.set_yticks([])
        ax.tick_params(length=0)

        # Extend xlim left to fit the tissue-color strip
        ax.set_xlim(-strip_w - 0.5, n - 0.5)

        # Colored strip: one rectangle per row
        for i, color in enumerate(colors):
            ax.add_patch(mpatches.Rectangle(
                (-strip_w - 0.2, i - 0.5), strip_w, 1.0, color=color))

        # Tissue group label at midpoint of each group
        grp_idx: dict = {}
        for i, ln in enumerate(lineages):
            grp_idx.setdefault(ln, []).append(i)
        for ln, idxs in grp_idx.items():
            mid = np.mean(idxs)
            abbrev = _SC_ABBREV.get(ln, ln.split()[0].title())
            ax.text(-strip_w - 0.5, mid, abbrev,
                    ha="right", va="center",
                    fontsize=FONT_SM - 2, color=_SC_PALETTE[ln],
                    fontweight="bold", clip_on=False)

    for idx in range(len(cond_data), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(
        "Sanity check: cosine similarity heatmap (34 cell types, corpus-centered)",
        fontsize=FONT_LG, fontweight="bold")
    plt.tight_layout(rect=[0, 0.02, 1, 0.97])
    plt.savefig(OUT / "figS1b.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved: figS1b_sanity_heatmaps.png")


# ── Figure S2: Spearman ρ vs. RP ──────────────────────────────────────────────

def figS2_spearman_rp():
    df_ct   = pd.read_csv(RES_RP / "results_celltype.csv")
    df_full = pd.read_csv(RES_RP / "results_full.csv")

    keys = [k for k in CONDITIONS_ORDER if k in df_ct["condition"].values and k != "rp"]

    x = np.arange(len(keys))
    w = 0.36

    fig, ax = plt.subplots(figsize=(13, 5))
    ct_idx   = df_ct.set_index("condition").reindex(keys)
    full_idx = df_full.set_index("condition").reindex(keys)
    bars_ct   = ct_idx["spearman_rho"].values
    bars_full = full_idx["spearman_rho"].values
    colors    = [COLORS[k] for k in keys]

    err_ct   = np.vstack([bars_ct   - ct_idx["rho_ci_low"].values,
                           ct_idx["rho_ci_high"].values   - bars_ct]).clip(min=0)
    err_full = np.vstack([bars_full - full_idx["rho_ci_low"].values,
                           full_idx["rho_ci_high"].values - bars_full]).clip(min=0)

    ax.bar(x - w/2, bars_ct,   w, color=colors, alpha=0.95, label="Cell-type aggregated",
           yerr=err_ct,   error_kw=dict(ecolor="#333333", elinewidth=1, capsize=3))
    ax.bar(x + w/2, bars_full, w, color=colors, alpha=0.5,  label="Full 32K experiments", hatch="//",
           yerr=err_full, error_kw=dict(ecolor="#333333", elinewidth=1, capsize=3))

    ax.set_xticks(x)
    ax.set_xticklabels([DISPLAY[k] for k in keys], fontsize=FONT_MD)
    ax.set_ylabel("Spearman ρ  (vs. Regulatory Potential)")
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.set_ylim(-0.05, 0.55)
    ax.legend(frameon=False)

    sig_ct = ct_idx["mantel_pval"].values
    for i, (h_ct, p_ct, e_ct) in enumerate(zip(bars_ct, sig_ct, err_ct[1])):
        star = "***" if p_ct < 0.001 else "**" if p_ct < 0.01 else "*" if p_ct < 0.05 else ""
        if star:
            ax.text(x[i] - w/2, h_ct + e_ct + 0.01, star,
                    ha="center", va="bottom", fontsize=FONT_SM, color="#333333")

    ax.text(0.99, 0.98,
            "*** p < 0.001  ** p < 0.01  * p < 0.05\n"
            "Cell-type: Mantel permutation test\n"
            "Error bars: 95% bootstrap CI",
            transform=ax.transAxes, fontsize=FONT_SM, va="top", ha="right", color="#555555")

    plt.tight_layout()
    plt.savefig(OUT / "figS2.pdf", bbox_inches="tight")
    plt.savefig(OUT / "figS2.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: figS2")

    # Narrow single-column version: key conditions only (no name-only variants)
    keys_n = [k for k in keys if not k.startswith("cto_")]
    x_n = np.arange(len(keys_n))
    ct_idx_n   = df_ct.set_index("condition").reindex(keys_n)
    full_idx_n = df_full.set_index("condition").reindex(keys_n)
    bars_ct_n   = ct_idx_n["spearman_rho"].values
    bars_full_n = full_idx_n["spearman_rho"].values
    colors_n    = [COLORS[k] for k in keys_n]
    err_ct_n   = np.vstack([bars_ct_n   - ct_idx_n["rho_ci_low"].values,
                             ct_idx_n["rho_ci_high"].values   - bars_ct_n]).clip(min=0)
    err_full_n = np.vstack([bars_full_n - full_idx_n["rho_ci_low"].values,
                             full_idx_n["rho_ci_high"].values - bars_full_n]).clip(min=0)

    fig_n, ax_n = plt.subplots(figsize=(7, 5))
    ax_n.bar(x_n - w/2, bars_ct_n,   w, color=colors_n, alpha=0.95, label="Cell-type aggregated",
             yerr=err_ct_n,   error_kw=dict(ecolor="#333333", elinewidth=1, capsize=3))
    ax_n.bar(x_n + w/2, bars_full_n, w, color=colors_n, alpha=0.5,  label="Full 32K experiments", hatch="//",
             yerr=err_full_n, error_kw=dict(ecolor="#333333", elinewidth=1, capsize=3))
    ax_n.set_xticks(x_n)
    ax_n.set_xticklabels([DISPLAY[k] for k in keys_n], fontsize=FONT_MD)
    ax_n.set_ylabel("Spearman ρ  (vs. Regulatory Potential)")
    ax_n.axhline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax_n.set_ylim(-0.05, 0.55)
    ax_n.legend(frameon=False)
    sig_ct_n = ct_idx_n["mantel_pval"].values
    for i, (h_ct, p_ct, e_ct) in enumerate(zip(bars_ct_n, sig_ct_n, err_ct_n[1])):
        star = "***" if p_ct < 0.001 else "**" if p_ct < 0.01 else "*" if p_ct < 0.05 else ""
        if star:
            ax_n.text(x_n[i] - w/2, h_ct + e_ct + 0.01, star,
                      ha="center", va="bottom", fontsize=FONT_SM, color="#333333")
    ax_n.text(0.99, 0.98,
              "*** p < 0.001  ** p < 0.01\nMantel permutation test",
              transform=ax_n.transAxes, fontsize=FONT_SM, va="top", ha="right", color="#555555")
    plt.tight_layout()
    plt.savefig(OUT / "figS2n.pdf", bbox_inches="tight")
    plt.savefig(OUT / "figS2n.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: figS2n")


# ── Figure 4: LLM decomposition (partial correlation + cross top-k) ───────────

LLM_KEYS   = ["pubmedbert", "sapbert", "medcpt", "general"]
LLM_DISPLAY = {
    "pubmedbert": "PubMedBERT", "sapbert": "SapBERT",
    "medcpt": "MedCPT",        "general": "General",
}
LLM_COLORS = {
    "pubmedbert": "#E63946", "sapbert": "#F4A261",
    "medcpt": "#457B9D",     "general": "#9B5DE5",
}


def _pval_to_star(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


def _partial_pval(r, n):
    """Two-tailed p-value for partial correlation via t-distribution."""
    from scipy.stats import t as tdist
    df = n - 3
    t_stat = r * np.sqrt(df) / np.sqrt(max(1 - r**2, 1e-12))
    return float(2 * tdist.sf(abs(t_stat), df))


def fig4_llm_decomposition():
    df_part = pd.read_csv(OUT / "partial_corr_full.csv")
    df_topk = pd.read_csv(OUT / "cross_topk_full.csv")

    df_full = pd.read_csv(RES / "results_full.csv").set_index("condition")
    raw_pvals = {k: float(df_full.loc[k, "spearman_pval"]) for k in LLM_KEYS}

    N_SUBSAMPLE = 2_000_000
    partial_pvals = {
        row["model"]: _partial_pval(row["partial"], N_SUBSAMPLE)
        for _, row in df_part.iterrows()
    }

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8),
                             gridspec_kw={"width_ratios": [1.3, 1]})

    # ── Panel A: partial correlation ─────────────────────────────────────────
    ax = axes[0]
    x = np.arange(len(LLM_KEYS))
    w = 0.35
    colors = [LLM_COLORS[k] for k in LLM_KEYS]

    ax.bar(x - w/2, df_part["raw"],     w, color=colors, alpha=0.9,
           label="Raw ρ (vs. Atacformer)")
    ax.bar(x + w/2, df_part["partial"], w, color=colors, alpha=0.45,
           hatch="//", label="Partial ρ (controlling for TF-IDF)")
    ax.axhline(0, color="black", lw=0.9, ls="--", alpha=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([LLM_DISPLAY[k] for k in LLM_KEYS])
    ax.set_ylabel("Spearman ρ  (full 32K level)")
    ax.set_title("(A)  Does LLM–Atacformer correlation\nsurvive controlling for TF-IDF?")
    ax.legend(fontsize=FONT_MD, frameon=False, loc="upper right")
    ax.set_ylim(-0.08, 0.42)

    for i, (key, row) in enumerate(zip(LLM_KEYS, df_part.itertuples())):
        raw_star  = _pval_to_star(raw_pvals[key])
        part_star = _pval_to_star(partial_pvals[LLM_DISPLAY[key]])
        ax.text(i - w/2, row.raw + 0.005, raw_star,
                ha="center", va="bottom", fontsize=FONT_MD, color="#333333")
        h_part = row.partial
        pad = 0.005 if h_part >= 0 else -0.018
        va  = "bottom" if h_part >= 0 else "top"
        ax.text(i + w/2, h_part + pad, part_star,
                ha="center", va=va, fontsize=FONT_MD, color="#333333")

    ax.text(0.5, -0.18,
            "*** p < 0.001   ** p < 0.01   * p < 0.05\n"
            "Solid bars: Spearman p-value  |  Hatched bars: t-test (df = n−3, n = 2M pairs)",
            transform=ax.transAxes, fontsize=FONT_SM, va="top", ha="center", color="#666666")

    # ── Panel B: cross top-k ─────────────────────────────────────────────────
    ax = axes[1]
    for key in LLM_KEYS:
        sub = df_topk[df_topk["model"] == LLM_DISPLAY[key]]
        ax.plot(sub["k"], sub["overlap"], marker="o",
                label=LLM_DISPLAY[key], color=LLM_COLORS[key],
                linewidth=2, markersize=6)
    ax.set_xlabel("k")
    ax.set_ylabel("Fraction of top-k neighbors shared with TF-IDF")
    ax.set_title("(B)  How much do LLM and TF-IDF\nneighbor lists overlap?")
    ax.set_xticks(K_VALUES)
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
    ax.grid(axis="y", alpha=0.2, linestyle="--")
    ax.legend(frameon=False)
    ax.set_ylim(0, 0.6)

    plt.tight_layout(w_pad=3)
    plt.subplots_adjust(bottom=0.18)
    plt.savefig(OUT / "fig4.pdf", bbox_inches="tight")
    plt.savefig(OUT / "fig4.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: fig4")

    # Panel A only — 6.5" wide for single-column (scale ≈ 0.54)
    fig_a, ax_a = plt.subplots(figsize=(6.5, 5.8))
    ax_a.bar(x - w/2, df_part["raw"],     w, color=colors, alpha=0.9,
             label="Raw ρ (vs. Atacformer)")
    ax_a.bar(x + w/2, df_part["partial"], w, color=colors, alpha=0.45,
             hatch="//", label="Partial ρ (controlling for TF-IDF)")
    ax_a.axhline(0, color="black", lw=0.9, ls="--", alpha=0.6)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels([LLM_DISPLAY[k] for k in LLM_KEYS])
    ax_a.set_ylabel("Spearman ρ  (full 32K level)")
    ax_a.legend(fontsize=FONT_MD, frameon=False, loc="upper right")
    ax_a.set_ylim(-0.08, 0.42)
    for i, (key, row) in enumerate(zip(LLM_KEYS, df_part.itertuples())):
        raw_star  = _pval_to_star(raw_pvals[key])
        part_star = _pval_to_star(partial_pvals[LLM_DISPLAY[key]])
        ax_a.text(i - w/2, row.raw + 0.005, raw_star,
                  ha="center", va="bottom", fontsize=FONT_MD, color="#333333")
        h_part = row.partial
        pad_p = 0.005 if h_part >= 0 else -0.018
        va_p  = "bottom" if h_part >= 0 else "top"
        ax_a.text(i + w/2, h_part + pad_p, part_star,
                  ha="center", va=va_p, fontsize=FONT_MD, color="#333333")
    ax_a.text(0.5, -0.15,
              "*** p < 0.001   ** p < 0.01   * p < 0.05\n"
              "Solid bars: Spearman p-value  |  Hatched bars: t-test (df = n−3, n = 2M pairs)",
              transform=ax_a.transAxes, fontsize=FONT_SM, va="top", ha="center", color="#666666")
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)
    plt.savefig(OUT / "fig4a.pdf", bbox_inches="tight")
    plt.savefig(OUT / "fig4a.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: fig4a")

    # Panel B only
    fig_b, ax_b = plt.subplots(figsize=(6, 5.5))
    for key in LLM_KEYS:
        sub = df_topk[df_topk["model"] == LLM_DISPLAY[key]]
        ax_b.plot(sub["k"], sub["overlap"], marker="o",
                  label=LLM_DISPLAY[key], color=LLM_COLORS[key],
                  linewidth=2, markersize=6)
    ax_b.set_xlabel("k")
    ax_b.set_ylabel("Fraction of top-k neighbors shared with TF-IDF")
    ax_b.set_xticks(K_VALUES)
    ax_b.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1))
    ax_b.grid(axis="y", alpha=0.2, linestyle="--")
    ax_b.legend(frameon=False)
    ax_b.set_ylim(0.2, 0.5)
    plt.tight_layout()
    plt.savefig(OUT / "fig4b.pdf", bbox_inches="tight")
    plt.savefig(OUT / "fig4b.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: fig4b")


# ── Figure S3: Similarity scatter (TF-IDF vs LLM, colored by Atacformer) ─────

def figS3_similarity_scatter():
    from scipy.stats import spearmanr as _spearmanr
    from scipy.stats import ttest_ind as _ttest_ind

    SIM2 = Path("/home/jsong205/atac_out/similarity_matrices")
    rng2 = np.random.default_rng(42)

    print("  Loading full matrices for scatter...")
    atac  = np.load(SIM2 / "full_atacformer.npy")
    tfidf = np.load(SIM2 / "full_tfidf.npy")
    idx_all = np.triu_indices(atac.shape[0], k=1)
    sample  = rng2.choice(len(idx_all[0]), 25_000, replace=False)
    a_s = atac[idx_all][sample]
    t_s = tfidf[idx_all][sample]
    del atac, tfidf

    fig, axes = plt.subplots(2, 4, figsize=(16, 8.5), gridspec_kw={"height_ratios": [3, 2]})
    sc = None
    stats_rows = []
    stored_l_s = {}
    stored_rho = {}
    stored_above = {}
    stored_below = {}
    stored_means = {}

    for col, key in enumerate(LLM_KEYS):
        ax = axes[0, col]
        llm = np.load(SIM2 / f"full_{key}.npy")
        l_s = llm[idx_all][sample]
        del llm
        rho = float(_spearmanr(t_s, l_s)[0])
        stored_l_s[key] = l_s.copy()
        stored_rho[key] = rho
        sc = ax.scatter(t_s, l_s, c=a_s, cmap="RdYlBu_r", s=2, alpha=0.3,
                        vmin=np.percentile(a_s, 2), vmax=np.percentile(a_s, 98))
        lim = [min(t_s.min(), l_s.min()), max(t_s.max(), l_s.max())]
        ax.plot(lim, lim, color="black", lw=0.8, ls="--", alpha=0.5)
        ax.set_xlabel("TF-IDF similarity")
        ax.set_title(LLM_DISPLAY[key])
        ax.text(0.05, 0.95, f"ρ = {rho:.3f}",
                transform=ax.transAxes, fontsize=FONT_MD, va="top")
        if col == 0:
            ax.set_ylabel("LLM similarity")

        above = a_s[l_s > t_s]
        below = a_s[l_s <= t_s]
        bx = axes[1, col]
        parts = bx.violinplot([below, above], widths=0.7,
                               showmeans=True, showmedians=False, showextrema=True)
        for body, c in zip(parts["bodies"], ["#BBBBBB", LLM_COLORS[key]]):
            body.set_facecolor(c); body.set_alpha(0.7)
            body.set_edgecolor("black"); body.set_linewidth(0.8)
        for part_name in ("cbars", "cmins", "cmaxes", "cmeans"):
            parts[part_name].set_edgecolor("black"); parts[part_name].set_linewidth(0.8)
        bx.set_xticks([1, 2])
        bx.set_xticklabels(["Below", "Above"])

        t_stat, p_val = _ttest_ind(above, below, equal_var=False)
        pooled_sd = np.sqrt((np.var(above, ddof=1) + np.var(below, ddof=1)) / 2)
        cohens_d  = (np.mean(above) - np.mean(below)) / pooled_sd
        stats_rows.append({
            "model": LLM_DISPLAY[key],
            "n_below": len(below), "n_above": len(above),
            "mean_below": float(np.mean(below)), "mean_above": float(np.mean(above)),
            "t_stat": float(t_stat), "p_value": float(p_val), "cohens_d": float(cohens_d),
        })
        sig = "***" if p_val < 1e-3 else "**" if p_val < 1e-2 else "*" if p_val < 0.05 else "n.s."
        stored_above[key] = above.copy()
        stored_below[key] = below.copy()
        stored_means[key] = (float(np.mean(below)), float(np.mean(above)), sig)
        bx.set_title(f"mean: {np.mean(below):.3f} vs {np.mean(above):.3f}  ({sig})")
        if col == 0:
            bx.set_ylabel("Atacformer similarity")

    cbar = fig.colorbar(sc, ax=axes[0, -1], shrink=0.85)
    cbar.set_label("Atacformer similarity")
    fig.suptitle(
        "TF-IDF vs. LLM pairwise similarity (top, colored by Atacformer similarity) and\n"
        "Atacformer similarity for above- vs. below-diagonal pairs (bottom), full 32K level",
        fontsize=FONT_LG, y=1.02)
    plt.tight_layout()
    plt.savefig(OUT / "figS3.pdf", bbox_inches="tight")
    plt.savefig(OUT / "figS3.png", dpi=300, bbox_inches="tight")
    plt.close()

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(OUT / "figS3_above_below_ttest.csv", index=False)
    print(stats_df.to_string(index=False))
    print("Saved: figS3")

    # figS3a: scatter plots only — resized to 6.5" wide for single-column (scale ≈ 0.54)
    from matplotlib.gridspec import GridSpec as _GS
    fig_a = plt.figure(figsize=(6.5, 5.5))
    gs_a = _GS(2, 3, figure=fig_a, width_ratios=[1, 1, 0.08],
               left=0.10, right=0.88, top=0.90, bottom=0.10,
               wspace=0.45, hspace=0.50)
    sc_ref = None
    for i, key in enumerate(LLM_KEYS):
        ax = fig_a.add_subplot(gs_a[i // 2, i % 2])
        l_si  = stored_l_s[key]
        rho_i = stored_rho[key]
        sc_ref = ax.scatter(t_s, l_si, c=a_s, cmap="RdYlBu_r", s=1.5, alpha=0.25,
                            vmin=np.percentile(a_s, 2), vmax=np.percentile(a_s, 98))
        lim = [min(t_s.min(), l_si.min()), max(t_s.max(), l_si.max())]
        ax.plot(lim, lim, color="black", lw=0.8, ls="--", alpha=0.5)
        ax.set_xlabel("TF-IDF similarity")
        ax.set_ylabel("LLM similarity")
        ax.set_title(f"{LLM_DISPLAY[key]}  (r={rho_i:.3f})")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    cbar_ax  = fig_a.add_subplot(gs_a[:, 2])
    cbar_obj = fig_a.colorbar(sc_ref, cax=cbar_ax)
    cbar_obj.set_label("Atacformer similarity", fontsize=FONT_SM)
    cbar_ax.tick_params(labelsize=FONT_SM)
    plt.savefig(OUT / "figS3a.pdf", bbox_inches="tight")
    plt.savefig(OUT / "figS3a.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: figS3a")

    # figS3b: violin plots only, 2×2 single-column layout
    fig_b, axes_b = plt.subplots(2, 2, figsize=(7, 6.5))
    for i, key in enumerate(LLM_KEYS):
        ax = axes_b[i // 2, i % 2]
        below_i, above_i = stored_below[key], stored_above[key]
        mean_b, mean_a, sig_i = stored_means[key]
        parts = ax.violinplot([below_i, above_i], widths=0.7,
                               showmeans=True, showmedians=False, showextrema=True)
        for body, c in zip(parts["bodies"], ["#BBBBBB", LLM_COLORS[key]]):
            body.set_facecolor(c); body.set_alpha(0.7)
            body.set_edgecolor("black"); body.set_linewidth(0.8)
        for pname in ("cbars", "cmins", "cmaxes", "cmeans"):
            parts[pname].set_edgecolor("black"); parts[pname].set_linewidth(0.8)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Below", "Above"])
        ax.set_title(f"{LLM_DISPLAY[key]}: {mean_b:.3f} vs {mean_a:.3f} ({sig_i})")
        if i % 2 == 0:
            ax.set_ylabel("Atacformer similarity")
    plt.tight_layout()
    plt.savefig(OUT / "figS3b.pdf", bbox_inches="tight")
    plt.savefig(OUT / "figS3b.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved: figS3b")


if __name__ == "__main__":
    # fig1 is the BioRender PNG — not regenerated from code
    fig2_spearman()
    fig3_topk_comparison()
    fig4_llm_decomposition()
    figS1_sanity_grid()
    figS2_spearman_rp()
    figS3_similarity_scatter()
    print("All paper figures saved to", OUT)
