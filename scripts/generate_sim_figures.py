#!/usr/bin/env python3
"""
Generate publication figures for EMITS simulated community validation.

Generates 4 figures:
  fig_sim1: Estimated vs Expected abundance scatter (EMITS, Naive, optionally EMU)
  fig_sim2: Within-genus resolution for all 6 multi-species genera
  fig_sim3: False positive suppression — Penicillium deep dive
  fig_sim4: Per-species absolute error comparison

Usage:
    # Two-way (EMITS vs Naive) — original behaviour
    python generate_sim_figures.py \\
        --comparison results/simulated_community/emu_its_results_species_comparison.tsv \\
        --truth results/simulated_community/ground_truth.tsv \\
        --outdir figures/simulated_community

    # Three-way (EMITS vs Naive vs EMU) — new
    python generate_sim_figures.py \\
        --comparison results/simulated_community/emu_its_results_species_comparison.tsv \\
        --truth results/simulated_community/ground_truth.tsv \\
        --emu results/emu_synthetic/synthetic_emu_rel-abundance.tsv \\
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

EM_COLOR = "#2563EB"        # EMITS-EM, blue
NAIVE_COLOR = "#DC2626"     # Naive, red
EMU_COLOR = "#F59E0B"       # EMU, orange
TRUTH_COLOR = "#059669"     # Truth, green


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


def load_emu(path):
    """Load EMU rel-abundance.tsv into a dict keyed on species name (with spaces)."""
    emu = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tax_id = row["tax_id"]
            if tax_id in ("unmapped", "mapped_filtered", "mapped_unclassified"):
                continue
            sp_field = row.get("species", "").strip()
            if not sp_field:
                continue
            sp = sp_field.replace("_", " ").strip()
            try:
                ab = float(row["abundance"])
            except (ValueError, KeyError):
                continue
            emu[sp] = emu.get(sp, 0.0) + ab
    return emu


# =====================================================================
# FIG 1: Scatter — Estimated vs Expected
# =====================================================================
def fig_scatter(truth, results, outdir, emu=None):
    n_panels = 3 if emu is not None else 2
    fig_width = 14 if emu is not None else 10
    fig, axes = plt.subplots(1, n_panels, figsize=(fig_width, 4.5), sharey=True)
    if n_panels == 1:
        axes = [axes]

    panel_specs = [
        ("em", EM_COLOR, "EMITS (EM)"),
        ("naive", NAIVE_COLOR, "Naive (best-hit)"),
    ]
    if emu is not None:
        panel_specs.append(("emu", EMU_COLOR, "EMU"))

    for ax, (method, color, label) in zip(axes, panel_specs):
        x_vals, y_vals, labels = [], [], []
        for sp, t in sorted(truth.items()):
            if method == "emu":
                est = emu.get(sp, 0.0)
                x_vals.append(t * 100)
                y_vals.append(est * 100)
                labels.append(sp)
            else:
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

    # R1 (revision, applied for consistency): suptitle removed.
    fig.tight_layout()

    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(outdir, f"fig_sim1_scatter.{ext}"))
    plt.close(fig)
    print(f"  ✓ fig_sim1_scatter")


# =====================================================================
# FIG 2: Within-genus resolution — grouped bar for all 6 genera
# =====================================================================
def fig_within_genus(truth, results, outdir, emu=None):
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
        emu_vals = []

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

            if emu is not None:
                emu_vals.append(emu.get(sp, 0.0) * 100)

        x = np.arange(len(short_names))

        if emu is not None:
            # 4 bars per group: Truth, EMITS, Naive, EMU
            width = 0.20
            ax.bar(x - 1.5 * width, truth_vals, width, color=TRUTH_COLOR, alpha=0.7, label="Expected")
            ax.bar(x - 0.5 * width, em_vals, width, color=EM_COLOR, alpha=0.85, label="EMITS")
            ax.bar(x + 0.5 * width, naive_vals, width, color=NAIVE_COLOR, alpha=0.85, label="Naive")
            ax.bar(x + 1.5 * width, emu_vals, width, color=EMU_COLOR, alpha=0.85, label="EMU")
        else:
            # Original 3-bar grouping
            width = 0.25
            ax.bar(x - width, truth_vals, width, color=TRUTH_COLOR, alpha=0.7, label="Expected")
            ax.bar(x, em_vals, width, color=EM_COLOR, alpha=0.85, label="EMITS")
            ax.bar(x + width, naive_vals, width, color=NAIVE_COLOR, alpha=0.85, label="Naive")

        ax.set_xticks(x)
        ax.set_xticklabels(short_names, rotation=25, ha="right", fontsize=8)
        ax.set_title(genus, fontweight="bold")
        ax.set_ylabel("Abundance (%)" if idx % 3 == 0 else "")
        ax.grid(axis="y", alpha=0.2)

        if idx == 0:
            ax.legend(fontsize=8, loc="upper right")

    # R1 (revision, applied for consistency): suptitle removed.
    fig.tight_layout()

    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(outdir, f"fig_sim2_within_genus.{ext}"))
    plt.close(fig)
    print(f"  ✓ fig_sim2_within_genus")


# =====================================================================
# FIG 3: Penicillium false-positive deep dive
# =====================================================================
def fig_penicillium(truth, results, outdir, emu=None):
    # Build a unified species list with EMITS and (optionally) EMU
    pen_data = {}
    for r in results:
        if r["taxon"].startswith("Penicillium ") and (r["em"] > 0.0001 or r["naive"] > 0.0001):
            pen_data[r["taxon"]] = {"em": r["em"], "naive": r["naive"], "emu": 0.0}
    if emu is not None:
        for sp, ab in emu.items():
            if sp.startswith("Penicillium ") and ab > 0.0001:
                if sp not in pen_data:
                    pen_data[sp] = {"em": 0.0, "naive": 0.0, "emu": ab}
                else:
                    pen_data[sp]["emu"] = ab

    # Sort by max abundance across all available methods
    pen_sorted = sorted(
        pen_data.keys(),
        key=lambda t: -max(pen_data[t]["em"], pen_data[t]["naive"], pen_data[t]["emu"])
    )

    species = pen_sorted
    em_vals = [pen_data[t]["em"] * 100 for t in species]
    naive_vals = [pen_data[t]["naive"] * 100 for t in species]
    emu_vals = [pen_data[t]["emu"] * 100 for t in species] if emu is not None else None

    # Short names with star for expected species
    short = []
    for sp in species:
        parts = sp.split()
        s = f"P. {parts[1]}" if len(parts) > 1 else sp
        if sp in truth:
            s += " ★"
        short.append(s)

    fig_height = max(5, len(short) * 0.45 + 1.5)
    fig, ax = plt.subplots(figsize=(9, fig_height))

    y = np.arange(len(short))

    if emu is not None:
        height = 0.27
        ax.barh(y - height, naive_vals, height, color=NAIVE_COLOR, alpha=0.85, label="Naive")
        ax.barh(y, em_vals, height, color=EM_COLOR, alpha=0.85, label="EMITS")
        ax.barh(y + height, emu_vals, height, color=EMU_COLOR, alpha=0.85, label="EMU")
    else:
        height = 0.35
        ax.barh(y + height / 2, naive_vals, height, color=NAIVE_COLOR, alpha=0.85, label="Naive")
        ax.barh(y - height / 2, em_vals, height, color=EM_COLOR, alpha=0.85, label="EMITS")

    # Add expected abundance markers (diamond)
    for i, sp in enumerate(species):
        if sp in truth:
            expected = truth[sp] * 100
            ax.plot(expected, i, "D", color=TRUTH_COLOR, markersize=8, zorder=5,
                    markeredgecolor="white", markeredgewidth=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(short, fontsize=9)
    ax.set_xlabel("Abundance (%)")

    em_fp = sum(pen_data[t]["em"] * 100 for t in species if t not in truth)
    naive_fp = sum(pen_data[t]["naive"] * 100 for t in species if t not in truth)
    fp_red = (1 - em_fp / naive_fp) * 100 if naive_fp > 0 else 0

    # R1 (revision): title removed; the marker key and the FP reduction now
    # live in the LaTeX caption. Legend moved off "upper right", where it
    # obscured the expected-abundance diamond, into the lower-left slot vacated
    # by the removed annotation box.
    ax.legend(loc="lower left")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.2)

    # R1 (revision): in-panel annotation box removed. Totals are printed so they
    # can be quoted in the caption. Both the plotted subset and the whole-genus
    # totals are reported, because the caption must state which it is using.
    emu_fp = (sum(pen_data[t]["emu"] * 100 for t in species if t not in truth)
              if emu is not None else None)
    print(f"[Fig4] PLOTTED FP species (n={sum(1 for t in species if t not in truth)}): "
          f"EMITS={em_fp:.2f}%  Naive={naive_fp:.2f}%"
          + (f"  EMU={emu_fp:.2f}%" if emu_fp is not None else "")
          + f"   -> EMITS reduction vs naive {fp_red:.1f}%")
    all_em = sum(v["em"] * 100 for t, v in pen_data.items() if t not in truth)
    all_nv = sum(v["naive"] * 100 for t, v in pen_data.items() if t not in truth)
    all_emu = (sum(v["emu"] * 100 for t, v in pen_data.items() if t not in truth)
               if emu is not None else None)
    all_red = (1 - all_em / all_nv) * 100 if all_nv > 0 else 0
    print(f"[Fig4] ALL unexpected Penicillium (n={sum(1 for t in pen_data if t not in truth)}): "
          f"EMITS={all_em:.2f}%  Naive={all_nv:.2f}%"
          + (f"  EMU={all_emu:.2f}%" if all_emu is not None else "")
          + f"   -> EMITS reduction vs naive {all_red:.1f}%")

    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(outdir, f"fig_sim3_penicillium_fp.{ext}"))
    plt.close(fig)
    print(f"  ✓ fig_sim3_penicillium_fp")


# =====================================================================
# FIG 4: Per-species absolute error
# =====================================================================
def fig_error_comparison(truth, results, outdir, emu=None):
    species_list = sorted(truth.keys(), key=lambda x: -truth[x])

    em_errors = []
    naive_errors = []
    emu_errors = []
    short_names = []

    for sp in species_list:
        match = [r for r in results if r["taxon"] == sp]
        if match:
            em_err = abs(match[0]["em"] - truth[sp]) * 100
            nv_err = abs(match[0]["naive"] - truth[sp]) * 100
            em_errors.append(em_err)
            naive_errors.append(nv_err)
            if emu is not None:
                emu_err = abs(emu.get(sp, 0.0) - truth[sp]) * 100
                emu_errors.append(emu_err)
            parts = sp.split()
            short_names.append(f"{parts[0][0]}. {parts[1]}" if len(parts) > 1 else sp)

    # Total false positive error
    fp_em = sum(r["em"] for r in results if r["taxon"] not in truth) * 100
    fp_naive = sum(r["naive"] for r in results if r["taxon"] not in truth) * 100
    em_errors.append(fp_em)
    naive_errors.append(fp_naive)
    if emu is not None:
        fp_emu = sum(ab for sp, ab in emu.items() if sp not in truth) * 100
        emu_errors.append(fp_emu)
    short_names.append("False positives\n(all spurious)")

    fig, ax = plt.subplots(figsize=(11, 5))

    x = np.arange(len(short_names))

    if emu is not None:
        width = 0.27
        ax.bar(x - width, em_errors, width, color=EM_COLOR, alpha=0.85, label="EMITS")
        ax.bar(x, naive_errors, width, color=NAIVE_COLOR, alpha=0.85, label="Naive")
        ax.bar(x + width, emu_errors, width, color=EMU_COLOR, alpha=0.85, label="EMU")
    else:
        width = 0.35
        ax.bar(x - width / 2, em_errors, width, color=EM_COLOR, alpha=0.85, label="EMITS")
        ax.bar(x + width / 2, naive_errors, width, color=NAIVE_COLOR, alpha=0.85, label="Naive")

    ax.set_xticks(x)
    ax.set_xticklabels(short_names, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Absolute error (%)")
    # R1 (revision): title removed; content moved to the LaTeX caption.
    ax.legend()
    ax.grid(axis="y", alpha=0.2)

    total_em = sum(em_errors)
    total_naive = sum(naive_errors)
    # R1 (revision): in-panel annotation box removed; values printed for the caption.
    summary = f"[Fig5] Total L1 -- EMITS={total_em:.2f}%  Naive={total_naive:.2f}%"
    if emu is not None:
        total_emu = sum(emu_errors)
        summary += f"  EMU={total_emu:.2f}%"
    summary += f"   (EMITS improvement vs naive: {(1 - total_em / total_naive) * 100:.1f}%)"
    print(summary)

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
    parser.add_argument("--comparison", required=True, help="Species comparison TSV from EMITS")
    parser.add_argument("--truth", required=True, help="Ground truth TSV")
    parser.add_argument("--emu", default=None, help="(Optional) EMU rel-abundance.tsv")
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

    emu = None
    if args.emu:
        emu = load_emu(args.emu)
        print(f"  EMU: {len(emu)} species")

    print()
    fig_scatter(truth, results, args.outdir, emu=emu)
    fig_within_genus(truth, results, args.outdir, emu=emu)
    fig_penicillium(truth, results, args.outdir, emu=emu)
    fig_error_comparison(truth, results, args.outdir, emu=emu)

    print(f"\n  All figures saved to {args.outdir}/")


if __name__ == "__main__":
    main()
