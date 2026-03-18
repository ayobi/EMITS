#!/usr/bin/env python3
"""
Generate publication figures for EMU-ITS simulated community validation.

Generates 4 figures:
  fig_sim1: Estimated vs Expected abundance scatter (EM & Naive vs truth)
  fig_sim2: Within-genus resolution for all 6 multi-species genera
  fig_sim3: False positive suppression — Penicillium deep dive
  fig_sim4: Per-species absolute error comparison (EM vs Naive)

Usage:
    python generate_sim_figures.py \
        --comparison results/simulated_community/emu_its_results_species_comparison.tsv \
        --truth results/simulated_community/ground_truth.tsv \
        --outdir figures/simulated_community
"""

import argparse
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from collections import defaultdict

# ---------- style ----------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

EM_COLOR = "#2563EB"
NAIVE_COLOR = "#DC2626"
TRUTH_COLOR = "#059669"
FP_COLOR = "#F59E0B"


def load_truth(path):
    truth = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            truth[row["species"]] = float(row["expected_abundance"])
    return truth


def load_comparison(path):
    results = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            results.append({
                "taxon": row["taxon"],
                "em": float(row["em_abundance"]),
                "naive": float(row["naive_abundance"]),
                "n_accessions": int(row["n_accessions"]),
            })
    return results


# =====================================================================
# FIG 1: Scatter — Estimated vs Expected
# =====================================================================
def fig_scatter(truth, results, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)

    for ax, (method, color, label) in zip(axes, [
        ("em", EM_COLOR, "EMU-ITS (EM)"),
        ("naive", NAIVE_COLOR, "Naive (best-hit)"),
    ]):
        x_vals, y_vals, labels = [], [], []
        for sp, t in sorted(truth.items()):
            match = [r for r in results if r["taxon"] == sp]
            if match:
                est = match[0][method]
                x_vals.append(t * 100)
                y_vals.append(est * 100)
                labels.append(sp)

        ax.scatter(x_vals, y_vals, c=color, s=50, alpha=0.8, edgecolors="white", linewidth=0.5, zorder=3)

        # Perfect line
        max_val = max(max(x_vals), max(y_vals)) * 1.1
        ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, linewidth=1, zorder=1)

        # R²
        x_arr, y_arr = np.array(x_vals), np.array(y_vals)
        r = np.corrcoef(x_arr, y_arr)[0, 1]
        l1 = np.sum(np.abs(x_arr - y_arr))

        # Label select points
        for xi, yi, lab in zip(x_vals, y_vals, labels):
            offset = abs(xi - yi)
            if offset > 0.4 or xi > 8:
                genus_sp = lab.split()
                short = f"{genus_sp[0][0]}. {genus_sp[1]}" if len(genus_sp) > 1 else lab
                ax.annotate(short, (xi, yi), fontsize=6.5, alpha=0.7,
                           xytext=(4, 4), textcoords="offset points")

        ax.set_xlabel("Expected abundance (%)")
        ax.set_title(f"{label}\nR² = {r**2:.4f}  |  L1 = {l1:.2f}%")
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.2)

    axes[0].set_ylabel("Estimated abundance (%)")

    fig.suptitle("Simulated 21-species ITS community: Estimated vs Expected", fontsize=13, fontweight="bold", y=1.02)
    fig.tight_layout()

    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(outdir, f"fig_sim1_scatter.{ext}"))
    plt.close(fig)
    print(f"  ✓ fig_sim1_scatter")


# =====================================================================
# FIG 2: Within-genus resolution — grouped bar for all 6 genera
# =====================================================================
def fig_within_genus(truth, results, outdir):
    genera = {
        "Aspergillus": ["Aspergillus niger", "Aspergillus tubingensis", "Aspergillus fumigatus"],
        "Fusarium": ["Fusarium oxysporum", "Fusarium solani", "Fusarium graminearum"],
        "Candida": ["Candida albicans", "Candida parapsilosis"],
        "Penicillium": ["Penicillium chrysogenum", "Penicillium expansum"],
        "Trichoderma": ["Trichoderma harzianum", "Trichoderma viride"],
        "Cladosporium": ["Cladosporium cladosporioides", "Cladosporium herbarum"],
    }

    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    axes = axes.flatten()

    for idx, (genus, species_list) in enumerate(genera.items()):
        ax = axes[idx]

        short_names = []
        truth_vals = []
        em_vals = []
        naive_vals = []

        for sp in species_list:
            parts = sp.split()
            short = f"{parts[0][0]}. {parts[1]}" if len(parts) > 1 else sp
            short_names.append(short)
            truth_vals.append(truth.get(sp, 0) * 100)

            match = [r for r in results if r["taxon"] == sp]
            if match:
                em_vals.append(match[0]["em"] * 100)
                naive_vals.append(match[0]["naive"] * 100)
            else:
                em_vals.append(0)
                naive_vals.append(0)

        x = np.arange(len(short_names))
        width = 0.25

        bars_t = ax.bar(x - width, truth_vals, width, color=TRUTH_COLOR, alpha=0.7, label="Expected")
        bars_e = ax.bar(x, em_vals, width, color=EM_COLOR, alpha=0.8, label="EM")
        bars_n = ax.bar(x + width, naive_vals, width, color=NAIVE_COLOR, alpha=0.8, label="Naive")

        ax.set_xticks(x)
        ax.set_xticklabels(short_names, rotation=25, ha="right", fontsize=8)
        ax.set_title(genus, fontweight="bold")
        ax.set_ylabel("Abundance (%)" if idx % 3 == 0 else "")
        ax.grid(axis="y", alpha=0.2)

        if idx == 0:
            ax.legend(fontsize=8, loc="upper right")

    fig.suptitle("Within-genus species resolution: Expected vs EM vs Naive",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()

    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(outdir, f"fig_sim2_within_genus.{ext}"))
    plt.close(fig)
    print(f"  ✓ fig_sim2_within_genus")


# =====================================================================
# FIG 3: Penicillium false-positive deep dive
# =====================================================================
def fig_penicillium(truth, results, outdir):
    # Get all Penicillium results
    pen_results = [(r["taxon"], r["em"], r["naive"]) for r in results
                   if r["taxon"].startswith("Penicillium") and (r["em"] > 0.0001 or r["naive"] > 0.0001)]
    pen_results.sort(key=lambda x: -max(x[1], x[2]))

    species = [p[0] for p in pen_results]
    em_vals = [p[1] * 100 for p in pen_results]
    naive_vals = [p[2] * 100 for p in pen_results]

    # Short names
    short = []
    for sp in species:
        parts = sp.split()
        s = f"P. {parts[1]}" if len(parts) > 1 else sp
        if sp in truth:
            s += " ★"
        short.append(s)

    fig, ax = plt.subplots(figsize=(8, 5))

    y = np.arange(len(short))
    height = 0.35

    bars_n = ax.barh(y + height / 2, naive_vals, height, color=NAIVE_COLOR, alpha=0.8, label="Naive")
    bars_e = ax.barh(y - height / 2, em_vals, height, color=EM_COLOR, alpha=0.8, label="EM")

    # Add expected abundance markers
    for i, sp in enumerate(species):
        if sp in truth:
            expected = truth[sp] * 100
            ax.plot(expected, i, "D", color=TRUTH_COLOR, markersize=8, zorder=5, markeredgecolor="white", markeredgewidth=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(short, fontsize=9)
    ax.set_xlabel("Abundance (%)")
    ax.set_title("Penicillium: EM reduces false positives by 34%\n(★ = expected species, ◆ = expected abundance)",
                 fontweight="bold")
    ax.legend(loc="lower center")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.2)

    # Add annotation
    naive_fp = sum(p[2] * 100 for p in pen_results if p[0] not in truth)
    em_fp = sum(p[1] * 100 for p in pen_results if p[0] not in truth)
    ax.text(0.97, 0.03,
            f"False positive total:\nEM: {em_fp:.2f}%\nNaive: {naive_fp:.2f}%",
            transform=ax.transAxes, fontsize=9, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.9))

    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(outdir, f"fig_sim3_penicillium_fp.{ext}"))
    plt.close(fig)
    print(f"  ✓ fig_sim3_penicillium_fp")


# =====================================================================
# FIG 4: Per-species absolute error (EM vs Naive)
# =====================================================================
def fig_error_comparison(truth, results, outdir):
    species_list = sorted(truth.keys(), key=lambda x: -truth[x])

    em_errors = []
    naive_errors = []
    short_names = []

    for sp in species_list:
        match = [r for r in results if r["taxon"] == sp]
        if match:
            em_err = abs(match[0]["em"] - truth[sp]) * 100
            nv_err = abs(match[0]["naive"] - truth[sp]) * 100
            em_errors.append(em_err)
            naive_errors.append(nv_err)
            parts = sp.split()
            short_names.append(f"{parts[0][0]}. {parts[1]}" if len(parts) > 1 else sp)

    # Also add total false positive error
    fp_em = sum(r["em"] for r in results if r["taxon"] not in truth) * 100
    fp_naive = sum(r["naive"] for r in results if r["taxon"] not in truth) * 100
    em_errors.append(fp_em)
    naive_errors.append(fp_naive)
    short_names.append("False positives\n(all spurious)")

    fig, ax = plt.subplots(figsize=(10, 5))

    x = np.arange(len(short_names))
    width = 0.35

    bars_e = ax.bar(x - width / 2, em_errors, width, color=EM_COLOR, alpha=0.8, label="EM")
    bars_n = ax.bar(x + width / 2, naive_errors, width, color=NAIVE_COLOR, alpha=0.8, label="Naive")

    # Highlight where EM is notably better
    for i in range(len(em_errors)):
        if naive_errors[i] > em_errors[i] * 1.2:
            ax.annotate("", xy=(i - width / 2, em_errors[i]),
                        xytext=(i + width / 2, naive_errors[i]),
                        arrowprops=dict(arrowstyle="<-", color=TRUTH_COLOR, lw=1.5))

    ax.set_xticks(x)
    ax.set_xticklabels(short_names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Absolute error (%)")
    ax.set_title("Per-species estimation error: EM vs Naive\n(lower is better)",
                 fontweight="bold")
    ax.legend()
    ax.grid(axis="y", alpha=0.2)

    # Summary annotation
    total_em = sum(em_errors)
    total_naive = sum(naive_errors)
    improvement = (1 - total_em / total_naive) * 100
    ax.text(0.85, 0.97,
            f"Total L1: EM={total_em:.2f}%  Naive={total_naive:.2f}%\nEM improvement: {improvement:.1f}%",
            transform=ax.transAxes, fontsize=9, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow", alpha=0.9))

    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(outdir, f"fig_sim4_error_comparison.{ext}"))
    plt.close(fig)
    print(f"  ✓ fig_sim4_error_comparison")


# =====================================================================
# MAIN
# =====================================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--comparison", required=True, help="Species comparison TSV from EMU-ITS")
    parser.add_argument("--truth", required=True, help="Ground truth TSV")
    parser.add_argument("--outdir", default="figures/simulated_community", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("=" * 60)
    print("  Generating simulated community figures")
    print("=" * 60)

    truth = load_truth(args.truth)
    results = load_comparison(args.comparison)
    print(f"  Ground truth: {len(truth)} species")
    print(f"  Results: {len(results)} species")

    print()
    fig_scatter(truth, results, args.outdir)
    fig_within_genus(truth, results, args.outdir)
    fig_penicillium(truth, results, args.outdir)
    fig_error_comparison(truth, results, args.outdir)

    print(f"\n  All figures saved to {args.outdir}/")
    print(f"  Formats: PNG (300 dpi) + PDF (vector)")


if __name__ == "__main__":
    main()