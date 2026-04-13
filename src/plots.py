"""
plots.py
Generate Figures 4–7 (and a bonus figure) matching Liu & Zhu (2018).

Fig 4 : Sub-matrix size (number of nodes) vs traffic load — GORA vs AORA
Fig 5 : Computational complexity (time ms) vs traffic load — GORA vs AORA
Fig 6 : Routing capacity (% demands routed) vs traffic load
Fig 7 : Transmission delay (ms) vs traffic load
Bonus : Distribution of AORA sub-matrix sizes
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless / no display required
import matplotlib.pyplot as plt

# Shared style
_COLORS   = {"gora": "#1f77b4", "aora": "#d62728"}
_MARKERS  = {"gora": "o",        "aora": "s"}
_LABELS   = {"gora": "GORA",     "aora": "AORA"}
_FIG_SIZE = (6.5, 4.5)
_DPI      = 150


def _load_pct(load_levels):
    """Convert normalised load [0, 1.2] to percentage strings for x-axis."""
    return load_levels * 100


def plot_matrix_size_comparison(results: dict, output_dir: str = "figures"):
    """
    Fig 4 — Matrix node count vs traffic load.
    GORA is constant at M+N = 766; AORA varies by sub-graph size.
    """
    os.makedirs(output_dir, exist_ok=True)

    load_pct  = _load_pct(results["load_levels"])
    gora_size = results["gora_matrix_size"]          # constant 766
    aora_size = results["aora_avg_submatrix_size"]

    fig, ax = plt.subplots(figsize=_FIG_SIZE)
    ax.plot(load_pct, gora_size, color=_COLORS["gora"],
            marker=_MARKERS["gora"], label=_LABELS["gora"], linewidth=2)
    ax.plot(load_pct, aora_size, color=_COLORS["aora"],
            marker=_MARKERS["aora"], label=_LABELS["aora"], linewidth=2)

    ax.set_xlabel("Traffic Load (%)", fontsize=12)
    ax.set_ylabel("Matrix Size (number of nodes)", fontsize=12)
    ax.set_title("Fig. 4 — Correlation Matrix Size vs Traffic Load", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(left=0)

    path = os.path.join(output_dir, "fig4_matrix_size.png")
    fig.tight_layout()
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


def plot_computational_complexity(results: dict, output_dir: str = "figures"):
    """
    Fig 5 — Average Dijkstra computation time (ms) vs traffic load.
    """
    os.makedirs(output_dir, exist_ok=True)

    load_pct   = _load_pct(results["load_levels"])
    gora_time  = results["gora_avg_comp_time_ms"]
    aora_time  = results["aora_avg_comp_time_ms"]

    fig, ax = plt.subplots(figsize=_FIG_SIZE)
    ax.plot(load_pct, gora_time, color=_COLORS["gora"],
            marker=_MARKERS["gora"], label=_LABELS["gora"], linewidth=2)
    ax.plot(load_pct, aora_time, color=_COLORS["aora"],
            marker=_MARKERS["aora"], label=_LABELS["aora"], linewidth=2)

    ax.set_xlabel("Traffic Load (%)", fontsize=12)
    ax.set_ylabel("Average Computation Time (ms)", fontsize=12)
    ax.set_title("Fig. 5 — Computational Complexity vs Traffic Load", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(left=0)

    path = os.path.join(output_dir, "fig5_comp_time.png")
    fig.tight_layout()
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


def plot_routing_capacity(results: dict, output_dir: str = "figures"):
    """
    Fig 6 — Routing capacity (fraction of demands successfully routed) vs load.
    """
    os.makedirs(output_dir, exist_ok=True)

    load_pct   = _load_pct(results["load_levels"])
    gora_cap   = results["gora_routing_capacity"] * 100   # → %
    aora_cap   = results["aora_routing_capacity"] * 100

    fig, ax = plt.subplots(figsize=_FIG_SIZE)
    ax.plot(load_pct, gora_cap, color=_COLORS["gora"],
            marker=_MARKERS["gora"], label=_LABELS["gora"], linewidth=2)
    ax.plot(load_pct, aora_cap, color=_COLORS["aora"],
            marker=_MARKERS["aora"], label=_LABELS["aora"], linewidth=2)

    ax.set_xlabel("Traffic Load (%)", fontsize=12)
    ax.set_ylabel("Routing Capacity (%)", fontsize=12)
    ax.set_title("Fig. 6 — Routing Capacity vs Traffic Load", fontsize=13)
    ax.set_ylim([0, 105])
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(left=0)

    path = os.path.join(output_dir, "fig6_routing_capacity.png")
    fig.tight_layout()
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


def plot_transmission_delay(results: dict, output_dir: str = "figures"):
    """
    Fig 7 — Average end-to-end propagation delay (ms) vs traffic load.
    """
    os.makedirs(output_dir, exist_ok=True)

    load_pct    = _load_pct(results["load_levels"])
    gora_delay  = results["gora_avg_delay_ms"]
    aora_delay  = results["aora_avg_delay_ms"]

    # Replace inf with NaN so matplotlib skips missing points
    gora_delay = np.where(np.isinf(gora_delay), np.nan, gora_delay)
    aora_delay = np.where(np.isinf(aora_delay), np.nan, aora_delay)

    fig, ax = plt.subplots(figsize=_FIG_SIZE)
    ax.plot(load_pct, gora_delay, color=_COLORS["gora"],
            marker=_MARKERS["gora"], label=_LABELS["gora"], linewidth=2)
    ax.plot(load_pct, aora_delay, color=_COLORS["aora"],
            marker=_MARKERS["aora"], label=_LABELS["aora"], linewidth=2)

    ax.set_xlabel("Traffic Load (%)", fontsize=12)
    ax.set_ylabel("Average Transmission Delay (ms)", fontsize=12)
    ax.set_title("Fig. 7 — Transmission Delay vs Traffic Load", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(left=0)

    path = os.path.join(output_dir, "fig7_delay.png")
    fig.tight_layout()
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


def plot_submatrix_size_distribution(results: dict, output_dir: str = "figures"):
    """
    Bonus figure — Histogram of AORA sub-matrix sizes across all load levels.
    Illustrates the reduction achieved compared to the full 766-node matrix.
    """
    os.makedirs(output_dir, exist_ok=True)

    all_sizes = []
    for sizes in results["aora_submatrix_sizes"]:
        all_sizes.extend(sizes)

    full_size = results["full_matrix_size"]

    fig, ax = plt.subplots(figsize=_FIG_SIZE)
    ax.hist(all_sizes, bins=30, color=_COLORS["aora"], edgecolor="white",
            alpha=0.8, label="AORA sub-matrix")
    ax.axvline(full_size, color=_COLORS["gora"], linestyle="--",
               linewidth=2, label=f"GORA full ({full_size} nodes)")

    ax.set_xlabel("Sub-matrix Node Count", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title("AORA Sub-matrix Size Distribution", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)

    path = os.path.join(output_dir, "bonus_submatrix_distribution.png")
    fig.tight_layout()
    fig.savefig(path, dpi=_DPI)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


def plot_all(results: dict, output_dir: str = "figures"):
    """Generate and save all figures."""
    paths = [
        plot_matrix_size_comparison(results, output_dir),
        plot_computational_complexity(results, output_dir),
        plot_routing_capacity(results, output_dir),
        plot_transmission_delay(results, output_dir),
        plot_submatrix_size_distribution(results, output_dir),
    ]
    print(f"\nAll figures saved to '{output_dir}/'")
    return paths
