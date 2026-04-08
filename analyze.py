#!/usr/bin/env python3
"""
analyze.py - Analyze sandpile simulation output.
Produces plots of avalanche size/duration distributions, power-law fits,
and connectivity analysis.

Usage: python3 analyze.py [data_directory]
"""

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter
from scipy import stats

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
OUT_DIR = DATA_DIR

def load_csv(fname, col):
    data = []
    path = os.path.join(DATA_DIR, fname)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split(",")
            if len(parts) > col:
                data.append(int(parts[col]))
    return np.array(data)

def power_law_fit(x, y):
    """Fit log-log linear regression, return slope and R^2."""
    mask = (x > 0) & (y > 0)
    lx, ly = np.log10(x[mask]), np.log10(y[mask])
    slope, intercept, r, p, se = stats.linregress(lx, ly)
    return slope, intercept, r**2

def plot_distribution(data, xlabel, title, filename, expected_exp=None):
    counts = Counter(data)
    sizes = np.array(sorted(counts.keys()))
    freqs = np.array([counts[s] for s in sizes])
    prob = freqs / freqs.sum()

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Log-log histogram
    ax = axes[0]
    ax.scatter(sizes, prob, s=8, alpha=0.6, color='steelblue')
    # Fit power law to tail (top 90% of range)
    mask = sizes >= np.percentile(sizes, 10)
    s_fit, p_fit = sizes[mask], prob[mask]
    if len(s_fit) > 5:
        slope, intercept, r2 = power_law_fit(s_fit, p_fit)
        x_line = np.logspace(np.log10(s_fit.min()), np.log10(s_fit.max()), 50)
        y_line = 10**intercept * x_line**slope
        ax.plot(x_line, y_line, 'r-', lw=2,
                label=f'$\\tau = {-slope:.2f}$ ($R^2={r2:.3f}$)')
        if expected_exp:
            ax.set_title(f'{title}\nFitted exponent: {-slope:.2f} (expected ~{expected_exp})')
        else:
            ax.set_title(f'{title}\nFitted exponent: {-slope:.2f}')
        ax.legend()
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(xlabel); ax.set_ylabel('P(' + xlabel + ')')
    ax.grid(True, alpha=0.3)

    # CCDF
    ax2 = axes[1]
    sorted_data = np.sort(data)
    ccdf = 1 - np.arange(1, len(sorted_data)+1) / len(sorted_data)
    ax2.plot(sorted_data, ccdf, color='darkgreen', lw=0.5)
    ax2.set_xscale('log'); ax2.set_yscale('log')
    ax2.set_xlabel(xlabel); ax2.set_ylabel('CCDF')
    ax2.set_title(f'Complementary CDF of {xlabel}')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, filename), dpi=150)
    plt.close()
    print(f"  Saved {filename}")
    return slope if len(s_fit) > 5 else None

def plot_connectivity(fname):
    sizes, durations, sites, radii = [], [], [], []
    path = os.path.join(DATA_DIR, fname)
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.strip().split(",")
            if len(parts) >= 5:
                sizes.append(int(parts[1]))
                durations.append(int(parts[2]))
                sites.append(int(parts[3]))
                radii.append(int(parts[4]))
    sizes = np.array(sizes)
    durations = np.array(durations)
    sites = np.array(sites)
    radii = np.array(radii)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Size vs unique sites
    ax = axes[0]
    ax.scatter(sizes, sites, s=2, alpha=0.1, color='steelblue')
    ax.set_xlabel('Avalanche Size (topplings)')
    ax.set_ylabel('Unique Sites Toppled')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title('Avalanche Size vs Spatial Extent')
    ax.grid(True, alpha=0.3)

    # Size vs duration
    ax = axes[1]
    ax.scatter(sizes, durations, s=2, alpha=0.1, color='darkorange')
    ax.set_xlabel('Avalanche Size'); ax.set_ylabel('Duration (rounds)')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title('Avalanche Size vs Duration')
    ax.grid(True, alpha=0.3)

    # Size vs max radius
    ax = axes[2]
    ax.scatter(sizes, radii, s=2, alpha=0.1, color='forestgreen')
    ax.set_xlabel('Avalanche Size'); ax.set_ylabel('Max Radius')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_title('Avalanche Size vs Max Radius (Connectivity)')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "connectivity_analysis.png"), dpi=150)
    plt.close()
    print("  Saved connectivity_analysis.png")

def plot_grid_snapshot(fname):
    """Plot the last grid snapshot as a heatmap."""
    path = os.path.join(DATA_DIR, fname)
    grid = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                grid = []  # reset, keep last snapshot
                continue
            parts = line.strip().split(",")
            if parts and parts[0]:
                grid.append([int(x) for x in parts])
    if not grid:
        return
    grid = np.array(grid)
    fig, ax = plt.subplots(figsize=(6, 6))
    im = ax.imshow(grid, cmap='YlOrRd', interpolation='nearest')
    plt.colorbar(im, ax=ax, label='Stress level')
    ax.set_title('Final Grid State (Stress Distribution)')
    ax.set_xlabel('x'); ax.set_ylabel('y')
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "grid_snapshot.png"), dpi=150)
    plt.close()
    print("  Saved grid_snapshot.png")

if __name__ == "__main__":
    print("Analyzing simulation data...")

    # Avalanche size distribution (expected tau ~ 1.1 for 2D BTW)
    sizes = load_csv("avalanche_sizes.csv", 1)
    print(f"  Total avalanches: {len(sizes)}, max size: {sizes.max()}")
    tau = plot_distribution(sizes, "Avalanche Size", "Avalanche Size Distribution",
                           "size_distribution.png", expected_exp="1.1")

    # Duration distribution
    durations = load_csv("avalanche_durations.csv", 1)
    plot_distribution(durations, "Duration", "Avalanche Duration Distribution",
                      "duration_distribution.png", expected_exp="1.25")

    # Connectivity analysis
    plot_connectivity("connectivity_stats.csv")

    # Grid heatmap
    plot_grid_snapshot("grid_snapshots.csv")

    print("Analysis complete.")
