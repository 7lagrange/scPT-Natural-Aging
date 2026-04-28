"""
Cell Counts & Sequencing Depth QC
==================================
Input : cell-counts.csv
Output: /mnt/results/cell_type_qc/

Figures generated
-----------------
Tissue level (total across all cell types):
  total_reads_qc_per_tissue.svg          -- total reads per mouse, 2x6 bar
  cell_counts_qc_per_tissue.svg          -- total cell counts per mouse, 2x6 bar
  total_depth_per_tissue.svg             -- total reads  age3M+age27M, horizontal bar
  total_cellcount_per_tissue.svg         -- total counts age3M+age27M, horizontal bar
  rpc_tissue_level.svg                   -- RNA reads/cell per mouse, 2x6 bar

Cell type level (per tissue, 2x6 subplots):
  stacked_bar_reads_all_tissues.svg      -- total reads stacked by cell type (absolute)
  stacked_bar_reads_pct_all_tissues.svg  -- total reads stacked by cell type (%)
  stacked_bar_pct_all_tissues.svg        -- cell counts stacked by cell type (%)
  cellcount_by_celltype_per_tissue.svg   -- cell counts age3M+age27M, horizontal bar + 10k line
  totalreads_by_celltype_per_tissue.svg  -- total reads age3M+age27M, horizontal bar + 10M line
  rpc_celltype_level.svg                 -- RNA reads/cell age3M vs age27M, offset bar + 1500 line
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# ── Output directory ──────────────────────────────────────────────────────────
OUT = "/mnt/results/cell_type_qc"
os.makedirs(OUT, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv("/mnt/user-uploads/cell-counts.csv", index_col=0)
df["tissue"]    = df["tissue_ct"].str.split("-").str[0]
df["cell_type"] = df["tissue_ct"].str.split("-", n=1).str[1]

# ── Color maps ────────────────────────────────────────────────────────────────
df_meta = df[["mouse", "age"]].drop_duplicates()

n_3M  = df_meta[df_meta["age"] == "age3M"]["mouse"].nunique()
n_27M = df_meta[df_meta["age"] == "age27M"]["mouse"].nunique()

mouse_colors = {}
for age, palette, n in [("age3M", "Blues", n_3M), ("age27M", "Reds", n_27M)]:
    cmap = sns.color_palette(palette, n_colors=n + 2)[2:]
    for i, m in enumerate(sorted(df_meta[df_meta["age"] == age]["mouse"].unique())):
        mouse_colors[m] = cmap[i]

mice_3M      = sorted(df_meta[df_meta["age"] == "age3M"]["mouse"].unique())
mice_27M     = sorted(df_meta[df_meta["age"] == "age27M"]["mouse"].unique())
ordered_mice = mice_3M + mice_27M
tissue_order = sorted(df["tissue"].unique())

# Cell-type color palette per tissue (tab20 + tab20b)
tissue_ct_colors = {}
for tissue in tissue_order:
    cts = sorted(df[df["tissue"] == tissue]["cell_type"].unique())
    palette = (sns.color_palette("tab20", 20) + sns.color_palette("tab20b", 20))[: len(cts)]
    tissue_ct_colors[tissue] = dict(zip(cts, palette))

# ── Shared helpers ────────────────────────────────────────────────────────────
AGE3M_COLOR  = "#4393C3"
AGE27M_COLOR = "#D6604D"
REF_COLOR    = "#FF69B4"

sns.set_theme(style="ticks", font_scale=1.0)


def mouse_legend_handles():
    h = [mpatches.Patch(color="none", label="age3M")]
    h += [mpatches.Patch(color=mouse_colors[m], label=m) for m in mice_3M]
    h += [mpatches.Patch(color="none", label="")]
    h += [mpatches.Patch(color="none", label="age27M")]
    h += [mpatches.Patch(color=mouse_colors[m], label=m) for m in mice_27M]
    return h


def fmt_reads(x, _):
    return f"{x/1e6:.0f}M" if x >= 1e6 else (f"{x/1000:.0f}k" if x >= 1000 else str(int(x)))


def fmt_counts(x, _):
    return f"{x/1000:.0f}k" if x >= 1000 else str(int(x))


def save(fig, name):
    for ext in ("svg", "png"):
        fig.savefig(f"{OUT}/{name}.{ext}", bbox_inches="tight", dpi=150)
    plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Per-mouse, 2×6 bar plots
# ═════════════════════════════════════════════════════════════════════════════

def plot_per_mouse_2x6(metric, ylabel, title, filename, formatter):
    """Bar plot: one panel per tissue, X = mice (ordered), Y = metric."""
    agg = df.groupby(["tissue", "mouse"])[metric].sum().reset_index()

    fig, axes = plt.subplots(2, 6, figsize=(3.0 * 6, 6.0 * 2), sharey=False)
    for idx, tissue in enumerate(tissue_order):
        ax = axes.flat[idx]
        sub = agg[agg["tissue"] == tissue].set_index("mouse").reindex(
            [m for m in ordered_mice if m in agg[agg["tissue"] == tissue]["mouse"].values]
        ).reset_index()
        ax.bar(range(len(sub)), sub[metric],
               color=[mouse_colors[m] for m in sub["mouse"]],
               width=0.7, edgecolor="white", linewidth=0.5)
        ax.set_title(tissue, fontsize=13, fontweight="bold", pad=6)
        ax.set_xticks([])
        ax.yaxis.set_major_formatter(plt.FuncFormatter(formatter))
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", labelsize=11)
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4, integer=True))
    for row in range(2):
        axes[row, 0].set_ylabel(ylabel, fontsize=13)
    for idx in range(len(tissue_order), 12):
        axes.flat[idx].set_visible(False)
    fig.legend(handles=mouse_legend_handles(), loc="center right",
               bbox_to_anchor=(1.07, 0.5), fontsize=11, frameon=False,
               handlelength=1.4, handleheight=1.2)
    fig.suptitle(title, fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    save(fig, filename)


plot_per_mouse_2x6("total_reads",  "Total reads",      "Total reads per tissue",       "total_reads_qc_per_tissue",  fmt_reads)
plot_per_mouse_2x6("cell_count",   "Total cell count", "Total cell counts per tissue",  "cell_counts_qc_per_tissue",  fmt_counts)


# ── Tissue-level weighted RNA reads per cell ──────────────────────────────────
def plot_rpc_tissue():
    agg = (
        df.groupby(["tissue", "mouse"], group_keys=False)
        .apply(lambda x: x["total_reads"].sum() / x["cell_count"].sum())
        .reset_index(name="rpc")
    )
    fig, axes = plt.subplots(2, 6, figsize=(3.0 * 6, 6.0 * 2), sharey=False)
    for idx, tissue in enumerate(tissue_order):
        ax = axes.flat[idx]
        sub = agg[agg["tissue"] == tissue].set_index("mouse").reindex(
            [m for m in ordered_mice if m in agg[agg["tissue"] == tissue]["mouse"].values]
        ).reset_index()
        ax.bar(range(len(sub)), sub["rpc"],
               color=[mouse_colors[m] for m in sub["mouse"]],
               width=0.7, edgecolor="white", linewidth=0.5)
        ax.set_title(tissue, fontsize=13, fontweight="bold", pad=6)
        ax.set_xticks([])
        ax.yaxis.set_major_formatter(plt.FuncFormatter(fmt_counts))
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", labelsize=11)
        ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4))
    for row in range(2):
        axes[row, 0].set_ylabel("RNA reads per cell", fontsize=13)
    for idx in range(len(tissue_order), 12):
        axes.flat[idx].set_visible(False)
    fig.legend(handles=mouse_legend_handles(), loc="center right",
               bbox_to_anchor=(1.07, 0.5), fontsize=11, frameon=False,
               handlelength=1.4, handleheight=1.2)
    fig.suptitle("RNA reads per cell — tissue level", fontsize=15, fontweight="bold", y=1.02)
    plt.tight_layout()
    save(fig, "rpc_tissue_level")

plot_rpc_tissue()


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 2 — age3M + age27M combined, single horizontal bar per tissue
# ═════════════════════════════════════════════════════════════════════════════

def plot_combined_hbar(metric, xlabel, title, filename, formatter):
    """Horizontal stacked bar: one bar per tissue, blue=3M, red=27M."""
    agg = df.groupby(["tissue", "age"])[metric].sum().reset_index()
    pivot = agg.pivot(index="tissue", columns="age", values=metric).fillna(0)
    for col in ["age3M", "age27M"]:
        if col not in pivot.columns:
            pivot[col] = 0
    pivot["total"] = pivot["age3M"] + pivot["age27M"]
    pivot = pivot.sort_values("total", ascending=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    y = np.arange(len(pivot))
    ax.barh(y, pivot["age3M"],  height=0.55, color=AGE3M_COLOR,  label="age3M",  edgecolor="white", linewidth=0.4)
    ax.barh(y, pivot["age27M"], height=0.55, left=pivot["age3M"], color=AGE27M_COLOR, label="age27M", edgecolor="white", linewidth=0.4)
    for i, total in enumerate(pivot["total"]):
        ax.text(total + pivot["total"].max() * 0.01, i,
                formatter(total, None), va="center", ha="left", fontsize=10)
    ax.set_yticks(y)
    ax.set_yticklabels(pivot.index, fontsize=12)
    ax.set_xlabel(xlabel, fontsize=13)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(formatter))
    ax.tick_params(axis="x", labelsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=10)
    ax.legend(fontsize=11, frameon=False, loc="lower right")
    ax.set_xlim(0, pivot["total"].max() * 1.15)
    plt.tight_layout()
    save(fig, filename)


plot_combined_hbar("total_reads", "Total reads",      "Total sequencing depth per tissue\n(age3M + age27M)", "total_depth_per_tissue",     fmt_reads)
plot_combined_hbar("cell_count",  "Total cell count", "Total cell counts per tissue\n(age3M + age27M)",      "total_cellcount_per_tissue",  fmt_counts)


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Stacked bar by cell type, 2×6 subplots
# ═════════════════════════════════════════════════════════════════════════════

def plot_stacked_2x6(metric, ylabel, title, filename, formatter, normalize=False):
    """2×6 stacked bar: X = mice, stacked by cell type."""
    fig, axes = plt.subplots(2, 6, figsize=(3.8 * 6, 6.5 * 2))
    for idx, tissue in enumerate(tissue_order):
        ax = axes.flat[idx]
        sub = df[df["tissue"] == tissue]
        cell_types = sorted(sub["cell_type"].unique())
        ct_colors  = tissue_ct_colors[tissue]

        pivot = sub.pivot_table(index="mouse", columns="cell_type",
                                values=metric, aggfunc="sum").fillna(0)
        pivot = pivot.reindex([m for m in ordered_mice if m in pivot.index])
        if normalize:
            pivot = pivot.div(pivot.sum(axis=1), axis=0) * 100

        bottom = np.zeros(len(pivot))
        for ct in cell_types:
            vals = pivot[ct].values if ct in pivot.columns else np.zeros(len(pivot))
            ax.bar(range(len(pivot)), vals, bottom=bottom,
                   color=ct_colors[ct], label=ct, width=0.7,
                   edgecolor="white", linewidth=0.3)
            bottom += vals

        ax.set_xticks(range(len(pivot)))
        ax.set_xticklabels(pivot.index, rotation=45, ha="right",
                           fontsize=9, fontweight="bold")
        for tick, mouse in zip(ax.get_xticklabels(), pivot.index):
            tick.set_color(mouse_colors.get(mouse, "black"))

        if normalize:
            ax.set_ylim(0, 100)
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x)}%"))
            ax.yaxis.set_major_locator(plt.MultipleLocator(25))
        else:
            ax.yaxis.set_major_formatter(plt.FuncFormatter(formatter))
            ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=4))

        ax.set_title(tissue, fontsize=13, fontweight="bold", pad=5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", labelsize=10)

        ct_handles = [mpatches.Patch(color=ct_colors[ct], label=ct) for ct in cell_types]
        ax.legend(handles=ct_handles, fontsize=6.5, frameon=False,
                  loc="upper left", bbox_to_anchor=(1.01, 1),
                  title="Cell type", title_fontsize=7.5)

    for row in range(2):
        axes[row, 0].set_ylabel(ylabel, fontsize=12)
    for idx in range(len(tissue_order), 12):
        axes.flat[idx].set_visible(False)
    fig.suptitle(title, fontsize=16, fontweight="bold", y=1.01)
    plt.tight_layout(h_pad=4, w_pad=3)
    save(fig, filename)


plot_stacked_2x6("total_reads", "Total reads",               "Total reads by cell type per tissue",                    "stacked_bar_reads_all_tissues",     fmt_reads,  normalize=False)
plot_stacked_2x6("total_reads", "Total reads proportion (%)", "Total reads by cell type per tissue (% of total reads)", "stacked_bar_reads_pct_all_tissues",  fmt_reads,  normalize=True)
plot_stacked_2x6("cell_count",  "Cell type proportion (%)",   "Cell type composition per tissue (% of total cells)",    "stacked_bar_pct_all_tissues",        fmt_counts, normalize=True)


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 4 — age3M + age27M horizontal bar by cell type, 2×6 subplots
# ═════════════════════════════════════════════════════════════════════════════

def plot_celltype_hbar_2x6(metric, xlabel, title, filename, formatter, ref_line, ref_label):
    """2×6 horizontal stacked bar: Y = cell type, X = metric, blue=3M red=27M."""
    agg = df.groupby(["tissue", "cell_type", "age"])[metric].sum().reset_index()

    fig, axes = plt.subplots(2, 6, figsize=(4.5 * 6, 6.5 * 2))
    for idx, tissue in enumerate(tissue_order):
        ax = axes.flat[idx]
        sub = agg[agg["tissue"] == tissue]
        pivot = sub.pivot(index="cell_type", columns="age", values=metric).fillna(0)
        for col in ["age3M", "age27M"]:
            if col not in pivot.columns:
                pivot[col] = 0
        pivot["total"] = pivot["age3M"] + pivot["age27M"]
        pivot = pivot.sort_values("total", ascending=True)

        y = np.arange(len(pivot))
        ax.barh(y, pivot["age3M"],  height=0.6, color=AGE3M_COLOR,  edgecolor="white", linewidth=0.3)
        ax.barh(y, pivot["age27M"], height=0.6, left=pivot["age3M"], color=AGE27M_COLOR, edgecolor="white", linewidth=0.3)
        ax.axvline(x=ref_line, color=REF_COLOR, linewidth=1.5, linestyle="--", zorder=5)

        max_val = pivot["total"].max()
        for i, total in enumerate(pivot["total"]):
            ax.text(total + max_val * 0.02, i, formatter(total, None),
                    va="center", ha="left", fontsize=8)

        ax.set_yticks(y)
        ax.set_yticklabels(pivot.index, fontsize=9)
        ax.set_title(tissue, fontsize=13, fontweight="bold", pad=5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(formatter))
        ax.tick_params(axis="x", labelsize=9)
        ax.set_xlim(0, max_val * 1.2)

    for col in range(6):
        axes[1, col].set_xlabel(xlabel, fontsize=11)
    for idx in range(len(tissue_order), 12):
        axes.flat[idx].set_visible(False)

    legend_handles = [
        mpatches.Patch(color=AGE3M_COLOR,  label="age3M"),
        mpatches.Patch(color=AGE27M_COLOR, label="age27M"),
        plt.Line2D([0], [0], color=REF_COLOR, linewidth=1.5, linestyle="--", label=ref_label),
    ]
    fig.legend(handles=legend_handles, loc="lower right",
               bbox_to_anchor=(0.98, 0.02), fontsize=12, frameon=False)
    fig.suptitle(title, fontsize=16, fontweight="bold", y=1.01)
    plt.tight_layout(h_pad=4, w_pad=3)
    save(fig, filename)


plot_celltype_hbar_2x6("cell_count",  "Cell count",   "Cell counts by cell type per tissue (age3M + age27M)",   "cellcount_by_celltype_per_tissue",   fmt_counts, 10_000, "10k threshold")
plot_celltype_hbar_2x6("total_reads", "Total reads",  "Total reads by cell type per tissue (age3M + age27M)",   "totalreads_by_celltype_per_tissue",  fmt_reads,  10e6,   "10M threshold")


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 5 — RNA reads per cell, cell type level, offset grouped bar
# ═════════════════════════════════════════════════════════════════════════════

def plot_rpc_celltype():
    ct_agg = (
        df.groupby(["tissue", "cell_type", "age"], group_keys=False)
        .apply(lambda x: x["total_reads"].sum() / x["cell_count"].sum())
        .reset_index(name="rpc")
    )

    bar_h, gap = 0.35, 0.05
    fig, axes = plt.subplots(2, 6, figsize=(4.5 * 6, 6.5 * 2))

    for idx, tissue in enumerate(tissue_order):
        ax = axes.flat[idx]
        sub = ct_agg[ct_agg["tissue"] == tissue]
        pivot = sub.pivot(index="cell_type", columns="age", values="rpc").fillna(0)
        for col in ["age3M", "age27M"]:
            if col not in pivot.columns:
                pivot[col] = 0
        pivot["mean_rpc"] = (pivot["age3M"] + pivot["age27M"]) / (
            (pivot["age3M"] > 0).astype(int) + (pivot["age27M"] > 0).astype(int)
        )
        pivot = pivot.sort_values("mean_rpc", ascending=True)

        y     = np.arange(len(pivot))
        y_3M  = y + gap / 2 + bar_h / 2
        y_27M = y - gap / 2 - bar_h / 2

        ax.barh(y_3M,  pivot["age3M"],  height=bar_h, color=AGE3M_COLOR,  edgecolor="white", linewidth=0.3)
        ax.barh(y_27M, pivot["age27M"], height=bar_h, color=AGE27M_COLOR, edgecolor="white", linewidth=0.3)
        ax.axvline(x=1500, color=REF_COLOR, linewidth=1.5, linestyle="--", zorder=5)

        max_val = max(pivot[["age3M", "age27M"]].max().max(), 1)
        for i, (r3, r27) in enumerate(zip(pivot["age3M"], pivot["age27M"])):
            for val, ypos in [(r3, y_3M[i]), (r27, y_27M[i])]:
                if val > 0:
                    label = f"{val/1000:.1f}k" if val >= 1000 else f"{int(val)}"
                    ax.text(val + max_val * 0.015, ypos, label,
                            va="center", ha="left", fontsize=7)

        ax.set_yticks(y)
        ax.set_yticklabels(pivot.index, fontsize=9)
        ax.set_ylim(-0.7, len(pivot) - 0.3)
        ax.set_title(tissue, fontsize=13, fontweight="bold", pad=5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(fmt_counts))
        ax.tick_params(axis="x", labelsize=9)
        ax.set_xlim(0, max_val * 1.25)

    for col in range(6):
        axes[1, col].set_xlabel("RNA reads per cell", fontsize=11)
    for idx in range(len(tissue_order), 12):
        axes.flat[idx].set_visible(False)

    legend_handles = [
        mpatches.Patch(color=AGE3M_COLOR,  label="age3M"),
        mpatches.Patch(color=AGE27M_COLOR, label="age27M"),
        plt.Line2D([0], [0], color=REF_COLOR, linewidth=1.5, linestyle="--", label="1500 threshold"),
    ]
    fig.legend(handles=legend_handles, loc="lower right",
               bbox_to_anchor=(0.98, 0.02), fontsize=12, frameon=False)
    fig.suptitle("RNA reads per cell — cell type level (age3M vs age27M)",
                 fontsize=16, fontweight="bold", y=1.01)
    plt.tight_layout(h_pad=4, w_pad=3)
    save(fig, "rpc_celltype_level")

plot_rpc_celltype()

print("All figures saved to", OUT)
