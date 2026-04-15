import os
import numpy as np
import matplotlib
matplotlib.use("Agg") # no display required to create the plots
import matplotlib.pyplot as plt

# shared style constants for all figures
colors = {"gora": "#1f77b4", "aora": "#d62728"}
markers = {"gora": "o", "aora": "s"}
labels = {"gora": "GORA", "aora": "AORA"}
figSize = (6.5, 4.5)
dpi = 150


# converts normalised load [0, 1.2] to percentage for x-axis
def loadToPct(loadLevels):
    return loadLevels * 100


# fig 4 - matrix node count vs traffic load
# GORA is constant at 766, AORA varies by sub-matrix size
def plot_matrix_size_comparison(results: dict, outputDir: str = "figures"):
    os.makedirs(outputDir, exist_ok=True)

    loadPct = loadToPct(results["load_levels"])
    goraSize = results["gora_matrix_size"]
    aoraSize = results["aora_avg_submatrix_size"]

    fig, ax = plt.subplots(figsize=figSize)
    ax.plot(loadPct, goraSize, color=colors["gora"], marker=markers["gora"], label=labels["gora"], linewidth=2)
    ax.plot(loadPct, aoraSize, color=colors["aora"], marker=markers["aora"], label=labels["aora"], linewidth=2)

    ax.set_xlabel("Traffic Load (%)", fontsize=12)
    ax.set_ylabel("Matrix Size (number of nodes)", fontsize=12)
    ax.set_title("Fig. 4 - Correlation Matrix Size vs Traffic Load", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(left=0)

    path = os.path.join(outputDir, "fig4_matrix_size.png")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


# fig 5 - average dijkstra computation time (ms) vs traffic load
def plot_computational_complexity(results: dict, outputDir: str = "figures"):
    os.makedirs(outputDir, exist_ok=True)

    loadPct = loadToPct(results["load_levels"])
    goraTime = results["gora_avg_comp_time_ms"]
    aoraTime = results["aora_avg_comp_time_ms"]

    fig, ax = plt.subplots(figsize=figSize)
    ax.plot(loadPct, goraTime, color=colors["gora"], marker=markers["gora"], label=labels["gora"], linewidth=2)
    ax.plot(loadPct, aoraTime, color=colors["aora"], marker=markers["aora"], label=labels["aora"], linewidth=2)

    ax.set_xlabel("Traffic Load (%)", fontsize=12)
    ax.set_ylabel("Average Computation Time (ms)", fontsize=12)
    ax.set_title("Fig. 5 - Computational Complexity vs Traffic Load", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(left=0)

    path = os.path.join(outputDir, "fig5_comp_time.png")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


# fig 6 - fraction of demands successfully routed vs traffic load
def plot_routing_capacity(results: dict, outputDir: str = "figures"):
    os.makedirs(outputDir, exist_ok=True)

    loadPct = loadToPct(results["load_levels"])
    goraCap = results["gora_routing_capacity"] * 100 # convert to %
    aoraCap = results["aora_routing_capacity"] * 100

    fig, ax = plt.subplots(figsize=figSize)
    ax.plot(loadPct, goraCap, color=colors["gora"], marker=markers["gora"], label=labels["gora"], linewidth=2)
    ax.plot(loadPct, aoraCap, color=colors["aora"], marker=markers["aora"], label=labels["aora"], linewidth=2)

    ax.set_xlabel("Traffic Load (%)", fontsize=12)
    ax.set_ylabel("Routing Capacity (%)", fontsize=12)
    ax.set_title("Fig. 6 - Routing Capacity vs Traffic Load", fontsize=13)
    ax.set_ylim([0, 105])
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(left=0)

    path = os.path.join(outputDir, "fig6_routing_capacity.png")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


# fig 7 - average end-to-end propagation delay (ms) vs traffic load
def plot_transmission_delay(results: dict, outputDir: str = "figures"):
    os.makedirs(outputDir, exist_ok=True)

    loadPct = loadToPct(results["load_levels"])
    goraDelay = results["gora_avg_delay_ms"]
    aoraDelay = results["aora_avg_delay_ms"]

    # replace inf with NaN so matplotlib skips missing points
    goraDelay = np.where(np.isinf(goraDelay), np.nan, goraDelay)
    aoraDelay = np.where(np.isinf(aoraDelay), np.nan, aoraDelay)

    fig, ax = plt.subplots(figsize=figSize)
    ax.plot(loadPct, goraDelay, color=colors["gora"], marker=markers["gora"], label=labels["gora"], linewidth=2)
    ax.plot(loadPct, aoraDelay, color=colors["aora"], marker=markers["aora"], label=labels["aora"], linewidth=2)

    ax.set_xlabel("Traffic Load (%)", fontsize=12)
    ax.set_ylabel("Average Transmission Delay (ms)", fontsize=12)
    ax.set_title("Fig. 7 - Transmission Delay vs Traffic Load", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.set_xlim(left=0)

    path = os.path.join(outputDir, "fig7_delay.png")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


# bonus - histogram of AORA sub-matrix sizes across all load levels
def plot_submatrix_size_distribution(results: dict, outputDir: str = "figures"):
    os.makedirs(outputDir, exist_ok=True)

    allSizes = []
    for sizes in results["aora_submatrix_sizes"]:
        allSizes.extend(sizes)

    fullSize = results["full_matrix_size"]

    fig, ax = plt.subplots(figsize=figSize)
    ax.hist(allSizes, bins=30, color=colors["aora"], edgecolor="white", alpha=0.8, label="AORA sub-matrix")
    ax.axvline(fullSize, color=colors["gora"], linestyle="--", linewidth=2, label=f"GORA full ({fullSize} nodes)")

    ax.set_xlabel("Sub-matrix Node Count", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title("AORA Sub-matrix Size Distribution", fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, linestyle="--", alpha=0.5)

    path = os.path.join(outputDir, "bonus_submatrix_distribution.png")
    fig.tight_layout()
    fig.savefig(path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {path}")
    return path


# generates and saves all figures
def plot_all(results: dict, outputDir: str = "figures"):
    paths = [
        plot_matrix_size_comparison(results, outputDir),
        plot_computational_complexity(results, outputDir),
        plot_routing_capacity(results, outputDir),
        plot_transmission_delay(results, outputDir),
        plot_submatrix_size_distribution(results, outputDir),
    ]
    print(f"\nAll figures saved to '{outputDir}/'")
    return paths
