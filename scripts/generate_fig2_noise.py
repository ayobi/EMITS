#!/usr/bin/env python3
"""
Generate Figure 2: L1 error vs alignment score noise level.

This script can work in two modes:
  1. Parse the output of `emu-its simulate --all` (pipe or file)
  2. Use hardcoded values from the Rust simulation results

Usage:
  # Mode 1: Parse simulation output
  cargo run --release -- simulate --all 2>&1 | python3 scripts/generate_fig2_noise.py --stdin

  # Mode 2: Use saved values (if you already ran the sim and recorded them)
  python3 scripts/generate_fig2_noise.py

  # Custom output path
  python3 scripts/generate_fig2_noise.py --outdir figures/

Output: fig2_noise_curve.pdf
"""

import argparse
import re
import sys
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np


def parse_simulation_output(text):
    """Parse Experiment 6 noise sensitivity table from emu-its simulate --all output."""
    noise_levels = []
    naive_l1 = []
    em_l1 = []

    # Look for lines like: +/-30           0.069200   0.028710   ...
    # The format from Rust is:
    #   +/-0            0.023200   0.021138   ...   no
    #   +/-10           0.034800   0.012661   ...   YES
    pattern = re.compile(
        r'\+/-(\d+)\s+'           # noise level
        r'([\d.]+)\s+'            # naive L1
        r'([\d.]+)\s+'            # EM L1
        r'([\d.]+)\s+'            # naive BC
        r'([\d.]+)\s+'            # EM BC
        r'(\w+)'                  # YES/no marker
    )

    for line in text.splitlines():
        m = pattern.search(line)
        if m:
            noise_levels.append(int(m.group(1)))
            naive_l1.append(float(m.group(2)))
            em_l1.append(float(m.group(3)))

    return noise_levels, naive_l1, em_l1


# Hardcoded values from the validated Rust simulation (Experiment 6)
# These are the results from `cargo run --release -- simulate --all`
# Community: 3 species (F. oxysporum 60%, F. solani 20%, T. harzianum 20%)
# Cross-mapping: 50%, true_score=1900, cross_score=1880
DEFAULT_NOISE =    [0,      10,     20,     30,     40,     50,     60,     80,     100]
DEFAULT_NAIVE_L1 = [0.0232, 0.0348, 0.0560, 0.0692, 0.1040, 0.1420, 0.1720, 0.2200, 0.2753]
DEFAULT_EM_L1 =    [0.0211, 0.0127, 0.0140, 0.0140, 0.0142, 0.0138, 0.0141, 0.0145, 0.0140]

# Approximate % of best-hits that are wrong at each noise level
# (from the Python v2 validation results)
DEFAULT_WRONG_PCT = [0, 1.5, 3.2, 5.7, 12.0, 22.0, 28.4, 38.0, 45.0]


def generate_figure(noise, naive_l1, em_l1, wrong_pct=None, outdir='.'):
    """Generate the publication-quality noise curve figure."""

    fig, ax1 = plt.subplots(1, 1, figsize=(7, 4.5))

    # Main plot: L1 error vs noise
    ax1.plot(noise, naive_l1, 'o-', color='#E24B4A', linewidth=2, markersize=6,
             label='Naive (best-hit)', zorder=3)
    ax1.plot(noise, em_l1, 's-', color='#378ADD', linewidth=2, markersize=6,
             label='EMU-ITS (EM)', zorder=3)

    # Fill between to emphasize the gap
    ax1.fill_between(noise, em_l1, naive_l1, alpha=0.08, color='#E24B4A', zorder=1)

    ax1.set_xlabel('Alignment score noise (±)', fontsize=12)
    ax1.set_ylabel('L1 error', fontsize=12)
    ax1.set_xlim(-2, max(noise) + 2)
    ax1.set_ylim(0, max(naive_l1) * 1.12)

    # Add annotations for key points — positioned above the shaded region
    # ±30 point
    idx_30 = noise.index(30) if 30 in noise else None
    if idx_30 is not None:
        improvement = (1 - em_l1[idx_30] / naive_l1[idx_30]) * 100
        ax1.annotate(f'EM {improvement:.0f}% better',
                     xy=(noise[idx_30], naive_l1[idx_30]),
                     xytext=(noise[idx_30] - 18, naive_l1[idx_30] + 0.06),
                     fontsize=9, color='#5F5E5A',
                     arrowprops=dict(arrowstyle='->', color='#5F5E5A', lw=0.8))

    # ±60 point
    idx_60 = noise.index(60) if 60 in noise else None
    if idx_60 is not None:
        improvement = (1 - em_l1[idx_60] / naive_l1[idx_60]) * 100
        ax1.annotate(f'EM {improvement:.0f}% better',
                     xy=(noise[idx_60], naive_l1[idx_60]),
                     xytext=(noise[idx_60] - 18, naive_l1[idx_60] + 0.06),
                     fontsize=9, color='#5F5E5A',
                     arrowprops=dict(arrowstyle='->', color='#5F5E5A', lw=0.8))

    # Annotate the EM flatness
    em_mean = np.mean(em_l1[1:])  # exclude noise=0
    ax1.axhline(y=em_mean, color='#378ADD', linestyle=':', alpha=0.4, linewidth=1)
    ax1.text(max(noise) * 0.95, em_mean + 0.008,
             f'EM ≈ {em_mean:.3f}', fontsize=9, color='#378ADD',
             ha='right', style='italic')

    # Add secondary x-axis on top showing % wrong best-hits
    if wrong_pct is not None and len(wrong_pct) == len(noise):
        ax2 = ax1.twiny()
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(noise)
        ax2.set_xticklabels([f'{w:.0f}%' if w > 0 else '0%' for w in wrong_pct],
                            fontsize=8, color='#888780')
        ax2.set_xlabel('Best-hits assigned to wrong species', fontsize=10, color='#888780')
        ax2.tick_params(axis='x', colors='#888780')

    ax1.legend(fontsize=11, loc='upper left', framealpha=0.9)
    ax1.grid(True, alpha=0.15)

    # Title
    ax1.set_title('EM robustness to alignment score noise\n'
                   '(simulated 3-species community, 50% cross-mapping)',
                   fontsize=12, fontweight='bold', pad=30 if wrong_pct else 10)

    plt.tight_layout()

    # Save
    os.makedirs(outdir, exist_ok=True)
    pdf_path = os.path.join(outdir, 'fig2_noise_curve.pdf')
    png_path = os.path.join(outdir, 'fig2_noise_curve.png')
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")
    return pdf_path


def main():
    parser = argparse.ArgumentParser(description='Generate Fig 2: L1 error vs noise level')
    parser.add_argument('--stdin', action='store_true',
                        help='Parse simulation output from stdin')
    parser.add_argument('--input', type=str, default=None,
                        help='Parse simulation output from file')
    parser.add_argument('--outdir', type=str, default='.',
                        help='Output directory for figures')
    parser.add_argument('--no-wrong-pct', action='store_true',
                        help='Skip the secondary axis showing %% wrong best-hits')
    args = parser.parse_args()

    if args.stdin:
        text = sys.stdin.read()
        noise, naive_l1, em_l1 = parse_simulation_output(text)
        if not noise:
            print("WARNING: Could not parse noise data from stdin. Using defaults.", file=sys.stderr)
            noise, naive_l1, em_l1 = DEFAULT_NOISE, DEFAULT_NAIVE_L1, DEFAULT_EM_L1
    elif args.input:
        with open(args.input) as f:
            text = f.read()
        noise, naive_l1, em_l1 = parse_simulation_output(text)
        if not noise:
            print("WARNING: Could not parse noise data from file. Using defaults.", file=sys.stderr)
            noise, naive_l1, em_l1 = DEFAULT_NOISE, DEFAULT_NAIVE_L1, DEFAULT_EM_L1
    else:
        print("Using hardcoded simulation values (from validated Rust output)")
        noise, naive_l1, em_l1 = DEFAULT_NOISE, DEFAULT_NAIVE_L1, DEFAULT_EM_L1

    wrong_pct = None if args.no_wrong_pct else DEFAULT_WRONG_PCT

    print(f"\nNoise levels: {noise}")
    print(f"Naive L1:     {naive_l1}")
    print(f"EM L1:        {em_l1}")
    print()

    for n, nl, el in zip(noise, naive_l1, em_l1):
        imp = (1 - el / nl) * 100 if nl > 0 else 0
        marker = "✓ EM wins" if el < nl else "✗ naive wins"
        print(f"  ±{n:>3d}  naive={nl:.4f}  EM={el:.4f}  improvement={imp:+.1f}%  {marker}")

    generate_figure(noise, naive_l1, em_l1, wrong_pct, args.outdir)


if __name__ == '__main__':
    main()