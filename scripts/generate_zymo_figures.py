#!/usr/bin/env python3
"""
Generate Zymo mock community figures from EMU-ITS comparison TSV.

Produces:
  - fig3_mock_within_genus.pdf  (Main Fig 3 — Trichophyton, Penicillium, Aspergillus)
  - figS2_mock_within_genus_all.pdf  (Supp — if needed, broader view)
  - figS3_nakaseomyces_consolidation.pdf  (Supp — accession consolidation)

Usage:
  python3 generate_zymo_figures.py --comparison emu_its_abundances_comparison.tsv --outdir figures/
"""

import argparse
import csv
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def parse_comparison_tsv(filepath):
    """Parse the EMU-ITS comparison TSV and extract species info."""
    rows = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Parse UNITE header: Species_name|accession|SH...|type|taxonomy
            taxon = row['taxon']
            parts = taxon.split('|')
            species_name = parts[0].replace('_', ' ') if parts else taxon

            # Extract genus
            genus = species_name.split()[0] if ' ' in species_name else species_name

            # Extract taxonomy string if present
            taxonomy = parts[4] if len(parts) > 4 else ''

            # Parse species from taxonomy (s__ field)
            tax_species = ''
            if 's__' in taxonomy:
                s_match = re.search(r's__([^;]+)', taxonomy)
                if s_match:
                    tax_species = s_match.group(1).replace('_', ' ')

            rows.append({
                'taxon_full': taxon,
                'species_name': species_name,
                'genus': genus,
                'tax_species': tax_species,
                'accession': parts[1] if len(parts) > 1 else '',
                'em': float(row['em_abundance']),
                'naive': float(row['naive_abundance']),
                'em_error': float(row['em_error']),
                'naive_error': float(row['naive_error']),
            })
    return rows


def get_genus_data(rows, genus_name, top_n=5):
    """Get top N species for a genus, sorted by max(EM, naive)."""
    genus_rows = [r for r in rows if r['genus'] == genus_name]
    genus_rows.sort(key=lambda r: max(r['em'], r['naive']), reverse=True)
    return genus_rows[:top_n]


def shorten_species(name, genus):
    """Shorten 'Genus species' to 'G. species'."""
    if name.startswith(genus):
        return genus[0] + '.' + name[len(genus):]
    return name


# ═══════════════════════════════════════════════════════════════
# FIGURE 3: Within-genus species resolution (3 panels)
# ═══════════════════════════════════════════════════════════════
def generate_fig3(rows, outdir):
    """Main Figure 3: Trichophyton, Penicillium, Aspergillus panels."""

    genera = [
        ('Trichophyton', 4),
        ('Penicillium', 5),
        ('Aspergillus', 4),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    for idx, (genus, top_n) in enumerate(genera):
        ax = axes[idx]
        data = get_genus_data(rows, genus, top_n)

        if not data:
            ax.text(0.5, 0.5, f'No {genus} data', transform=ax.transAxes, ha='center')
            continue

        species = [shorten_species(r['species_name'], genus) for r in data]
        em_vals = [r['em'] * 100 for r in data]
        naive_vals = [r['naive'] * 100 for r in data]

        x = np.arange(len(species))
        width = 0.35

        bars_em = ax.bar(x - width/2, em_vals, width, label='EM',
                         color='#378ADD', alpha=0.85, edgecolor='white', linewidth=0.5)
        bars_naive = ax.bar(x + width/2, naive_vals, width, label='Naive',
                            color='#E24B4A', alpha=0.85, edgecolor='white', linewidth=0.5)

        ax.set_xlabel('')
        ax.set_ylabel('Abundance (%)' if idx == 0 else '')
        ax.set_title(f'$\\it{{{genus}}}$', fontsize=13, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f'$\\it{{{s}}}$' for s in species],
            rotation=35, ha='right', fontsize=9
        )
        ax.tick_params(axis='y', labelsize=9)

        if idx == 0:
            ax.legend(fontsize=10, loc='upper right')

        ax.grid(axis='y', alpha=0.15)
        ax.set_axisbelow(True)

        # Set y-axis limit with headroom for labels
        max_val = max(max(em_vals), max(naive_vals))
        ax.set_ylim(0, max_val * 1.25)

        # Add value labels on top of bars (only if they fit within the axes)
        for bar_group in [bars_em, bars_naive]:
            for bar in bar_group:
                height = bar.get_height()
                if height > 0.05:
                    label_y = height + max_val * 0.02
                    if label_y < ax.get_ylim()[1] * 0.97:
                        ax.text(bar.get_x() + bar.get_width()/2., label_y,
                                f'{height:.1f}', ha='center', va='bottom', fontsize=7,
                                clip_on=True)

    fig.suptitle('Within-genus species resolution: ONT mock community (EM vs Naive)',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()

    pdf_path = os.path.join(outdir, 'fig3_mock_within_genus.pdf')
    png_path = os.path.join(outdir, 'fig3_mock_within_genus.png')
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Fig 3 saved: {pdf_path}")


# ═══════════════════════════════════════════════════════════════
# FIGURE S3: Nakaseomyces accession consolidation
# ═══════════════════════════════════════════════════════════════
def generate_figS3(rows, outdir):
    """Supplementary: Nakaseomyces accession consolidation."""

    # Get all Nakaseomyces entries with any abundance
    nak_rows = [r for r in rows if r['genus'] == 'Nakaseomyces'
                and (r['em'] > 1e-5 or r['naive'] > 1e-5)]
    nak_rows.sort(key=lambda r: r['em'], reverse=True)

    if not nak_rows:
        print("WARNING: No Nakaseomyces data found")
        return

    # Label accessions
    labels = []
    for i, r in enumerate(nak_rows):
        acc = r['accession'][:10] if r['accession'] else f'acc_{i}'
        sp = shorten_species(r['species_name'], 'Nakaseomyces')
        labels.append(f'{sp}\n({acc})')

    em_vals = [r['em'] * 100 for r in nak_rows]
    naive_vals = [r['naive'] * 100 for r in nak_rows]

    fig, ax = plt.subplots(figsize=(9, max(4, len(labels) * 0.5 + 1.5)))

    y = np.arange(len(labels))
    height = 0.35

    ax.barh(y - height/2, em_vals, height, label='EM',
            color='#378ADD', alpha=0.85, edgecolor='white', linewidth=0.5)
    ax.barh(y + height/2, naive_vals, height, label='Naive',
            color='#E24B4A', alpha=0.85, edgecolor='white', linewidth=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Abundance (%)', fontsize=11)
    ax.set_title('$\\it{Nakaseomyces}$ $\\it{glabratus}$: EM consolidates across UNITE accessions',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(axis='x', alpha=0.15)
    ax.set_axisbelow(True)
    ax.invert_yaxis()

    # Add value labels
    for i, (em, naive) in enumerate(zip(em_vals, naive_vals)):
        if em > 0.1:
            ax.text(em + 0.15, i - height/2, f'{em:.1f}%', va='center', fontsize=8, color='#185FA5')
        if naive > 0.1:
            ax.text(naive + 0.15, i + height/2, f'{naive:.1f}%', va='center', fontsize=8, color='#A32D2D')

    plt.tight_layout()

    pdf_path = os.path.join(outdir, 'figS3_nakaseomyces_consolidation.pdf')
    png_path = os.path.join(outdir, 'figS3_nakaseomyces_consolidation.png')
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Fig S3 saved: {pdf_path}")


# ═══════════════════════════════════════════════════════════════
# BONUS: Print genus-level summary for the paper
# ═══════════════════════════════════════════════════════════════
def print_genus_summary(rows):
    """Print genus-level summary matching Table 2 in the paper."""

    target_genera = [
        'Aspergillus', 'Cryptococcus', 'Trichophyton', 'Penicillium',
        'Fusarium', 'Candida', 'Nakaseomyces', 'Malassezia',
        'Saccharomyces', 'Cutaneotrichosporon'
    ]

    print("\n" + "=" * 80)
    print("  Genus-level summary (ONT mock community)")
    print("=" * 80)

    for genus in target_genera:
        genus_rows = [r for r in rows if r['genus'] == genus]
        em_total = sum(r['em'] for r in genus_rows) * 100
        naive_total = sum(r['naive'] for r in genus_rows) * 100
        n_refs = len(genus_rows)

        # Top species
        top = sorted(genus_rows, key=lambda r: r['em'], reverse=True)

        print(f"\n  {genus} (expected ~10%, {n_refs} refs in DB)")
        print(f"    Genus total: EM={em_total:.2f}%  Naive={naive_total:.2f}%")
        for r in top[:3]:
            if r['em'] > 0.0001 or r['naive'] > 0.0001:
                sp = shorten_species(r['species_name'], genus)
                print(f"    {sp:<30s}  EM={r['em']*100:.2f}%  Naive={r['naive']*100:.2f}%")


def main():
    parser = argparse.ArgumentParser(description='Generate Zymo mock community figures')
    parser.add_argument('--comparison', type=str,
                        default='emu_its_abundances_comparison.tsv',
                        help='Path to EMU-ITS comparison TSV')
    parser.add_argument('--outdir', type=str, default='.',
                        help='Output directory for figures')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Reading: {args.comparison}")
    rows = parse_comparison_tsv(args.comparison)
    print(f"Loaded {len(rows)} taxa")

    # Generate figures
    generate_fig3(rows, args.outdir)
    generate_figS3(rows, args.outdir)

    # Print summary
    print_genus_summary(rows)

    print(f"\nAll figures saved to {args.outdir}/")


if __name__ == '__main__':
    main()