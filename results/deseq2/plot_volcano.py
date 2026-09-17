#!/usr/bin/env python3
"""Volcano plots (log2FC_shrunk vs padj) for DESeq2 result TSVs."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(DIR, "volcano_plots")
os.makedirs(OUT_DIR, exist_ok=True)

LFC_THRESH = 0.5
PADJ_THRESH = 0.05

COL_UP = "#e34948"      # categorical slot 8 (red)
COL_DOWN = "#2a78d6"    # categorical slot 1 (blue)
COL_NS = "#c3c2b7"      # muted ink (ns points)
COL_LINE = "#898781"    # muted ink for threshold lines
COL_TEXT = "#0b0b0b"
COL_SEC_TEXT = "#52514e"
COL_GRID = "#e1e0d9"
COL_SURFACE = "#fcfcfb"

FILES = [
    "GSE135251_main_F4_vs_F0.tsv",
    "GSE135251_sens1_F4_vs_F0F1.tsv",
    "GSE135251_sens2_F3F4_vs_F0F1.tsv",
    "GSE135251_sens3_F4_vs_F0ctrl.tsv",
    "GSE162694_main_F4_vs_F0.tsv",
    "GSE162694_main_F4_vs_ctrl.tsv",
    "GSE162694_sens1_F4_vs_F0F1.tsv",
    "GSE162694_sens2_F3F4_vs_F0F1.tsv",
    "GSE162694_sens3_F4_vs_F0ctrl.tsv",
]

for fname in FILES:
    path = os.path.join(DIR, fname)
    if not os.path.exists(path):
        print(f"SKIP {fname} (not found -- run the DESeq2 script first)")
        continue
    df = pd.read_csv(path, sep="\t")
    df = df.dropna(subset=["log2FC_shrunk", "padj"]).copy()

    neglog10_padj = -np.log10(df["padj"].clip(lower=1e-300))
    lfc = df["log2FC_shrunk"]

    is_up = (lfc > LFC_THRESH) & (df["padj"] < PADJ_THRESH)
    is_down = (lfc < -LFC_THRESH) & (df["padj"] < PADJ_THRESH)
    is_ns = ~(is_up | is_down)

    n_up, n_down, n_ns = is_up.sum(), is_down.sum(), is_ns.sum()

    fig, ax = plt.subplots(figsize=(7, 6), dpi=150)
    fig.patch.set_facecolor(COL_SURFACE)
    ax.set_facecolor(COL_SURFACE)

    ax.scatter(lfc[is_ns], neglog10_padj[is_ns], s=6, c=COL_NS,
               alpha=0.5, linewidths=0, label=f"NS (n={n_ns:,})")
    ax.scatter(lfc[is_down], neglog10_padj[is_down], s=6, c=COL_DOWN,
               alpha=0.7, linewidths=0, label=f"Down (n={n_down:,})")
    ax.scatter(lfc[is_up], neglog10_padj[is_up], s=6, c=COL_UP,
               alpha=0.7, linewidths=0, label=f"Up (n={n_up:,})")

    # threshold lines (dashed)
    ax.axvline(LFC_THRESH, color=COL_LINE, linestyle="--", linewidth=1)
    ax.axvline(-LFC_THRESH, color=COL_LINE, linestyle="--", linewidth=1)
    ax.axhline(-np.log10(PADJ_THRESH), color=COL_LINE, linestyle="--", linewidth=1)

    ax.set_xlabel("log2 Fold Change (shrunk)", color=COL_TEXT)
    ax.set_ylabel("-log10(padj)", color=COL_TEXT)
    title = fname.replace(".tsv", "")
    ax.set_title(f"{title}\n|log2FC_shrunk| > {LFC_THRESH} & padj < {PADJ_THRESH}",
                 color=COL_TEXT, fontsize=11)

    ax.grid(True, color=COL_GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color(COL_GRID)
    ax.tick_params(colors=COL_SEC_TEXT)

    legend = ax.legend(loc="upper left", frameon=False, fontsize=9)
    for text in legend.get_texts():
        text.set_color(COL_SEC_TEXT)

    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, f"{title}_volcano.png")
    fig.savefig(out_path, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"Saved {out_path}  (up={n_up}, down={n_down}, ns={n_ns})")
