"""
Differential Bins QC Plots
===========================
Input : df_dge_count.csv
Output: diff_bins_pval_per_tissue.svg/png
        diff_bins_fdr_per_tissue.svg/png
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ── Config ────────────────────────────────────────────────────────────────────
INPUT  = "/mnt/user-uploads/df_dge_count.csv"
OUT    = "/mnt/results/dge"
os.makedirs(OUT, exist_ok=True)

UP_COLOR   = "#E8735A"   # warm red  — up
DOWN_COLOR = "#5B8DB8"   # steel blue — down

sns.set_theme(style="ticks", font_scale=1.0)

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(INPUT, index_col=0)
df["tissue"]    = df["tissue_ct"].str.split("-").str[0]
df["cell_type"] = df["tissue_ct"].str.split("-", n=1).str[1]

tissue_order = sorted(df["tissue"].unique())

# ── Core plot function ────────────────────────────────────────────────────────
def plot_diff_bins(up_col, down_col, sort_col, xlabel, title, filename):
    """
    2×6 horizontal grouped bar plot.
    Each tissue = one panel; each cell type = two offset bars (up / down).

    Parameters
    ----------
    up_col   : column name for up counts   (e.g. 'up_pval' or 'up_fdr')
    down_col : column name for down counts (e.g. 'down_pval' or 'down_fdr')
    sort_col : column used to sort cell types within each panel
    xlabel   : x-axis label
    title    : figure suptitle
    filename : output filename (without extension)
    """
    bar_h, gap = 0.35, 0.05
    n_cols, n_rows = 6, 2

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.8 * n_cols, 6.5 * n_rows))

    for idx, tissue in enumerate(tissue_order):
        ax = axes.flat[idx]
        sub = df[df["tissue"] == tissue].copy().sort_values(sort_col, ascending=True)
        n_ct = len(sub)

        y    = np.arange(n_ct)
        y_up = y + gap / 2 + bar_h / 2   # up bar: above centre
        y_dn = y - gap / 2 - bar_h / 2   # down bar: below centre

        ax.barh(y_up, sub[up_col],   height=bar_h, color=UP_COLOR,   edgecolor="white", linewidth=0.3)
        ax.barh(y_dn, sub[down_col], height=bar_h, color=DOWN_COLOR, edgecolor="white", linewidth=0.3)

        # Value labels on every bar
        max_val = max(sub[[up_col, down_col]].max().max(), 1)
        for i, (u, d) in enumerate(zip(sub[up_col], sub[down_col])):
            ax.text(u + max_val * 0.015, y_up[i], str(int(u)), va="center", ha="left", fontsize=7)
            ax.text(d + max_val * 0.015, y_dn[i], str(int(d)), va="center", ha="left", fontsize=7)

        ax.set_yticks(y)
        ax.set_yticklabels(sub["cell_type"].tolist(), fontsize=9)
        ax.set_ylim(-0.7, n_ct - 0.3)
        ax.set_title(tissue, fontsize=13, fontweight="bold", pad=5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", labelsize=9)
        ax.set_xlim(0, max_val * 1.25)

    for col in range(n_cols):
        axes[1, col].set_xlabel(xlabel, fontsize=11)
    for idx in range(len(tissue_order), n_cols * n_rows):
        axes.flat[idx].set_visible(False)

    legend_handles = [
        mpatches.Patch(color=UP_COLOR,   label=up_col),
        mpatches.Patch(color=DOWN_COLOR, label=down_col),
    ]
    fig.legend(handles=legend_handles,
               loc="lower center", bbox_to_anchor=(0.5, -0.02),
               ncol=2, fontsize=12, frameon=False)

    fig.suptitle(title, fontsize=16, fontweight="bold", y=1.01)
    plt.tight_layout(h_pad=4, w_pad=3)

    for ext in ("svg", "png"):
        fig.savefig(f"{OUT}/{filename}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Saved: {filename}")


# ── Generate plots ────────────────────────────────────────────────────────────

# p-value threshold
plot_diff_bins(
    up_col   = "up_pval",
    down_col = "down_pval",
    sort_col = "total_pval",
    xlabel   = "Number of diff bins (p-val)",
    title    = "Differential bins per cell type (H3K9me3, p-val threshold)",
    filename = "diff_bins_pval_per_tissue",
)

# FDR threshold
df["total_fdr"] = df["up_fdr"] + df["down_fdr"]
plot_diff_bins(
    up_col   = "up_fdr",
    down_col = "down_fdr",
    sort_col = "total_fdr",
    xlabel   = "Number of diff bins (FDR)",
    title    = "Differential bins per cell type (H3K9me3, FDR threshold)",
    filename = "diff_bins_fdr_per_tissue",
)
