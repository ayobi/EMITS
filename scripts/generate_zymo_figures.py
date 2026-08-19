#!/usr/bin/env python3
"""
Generate ATCC mock community figures from EMITS comparison TSV (and optionally EMU output).

Produces:
  - fig3_mock_within_genus.pdf  (Main Fig 3 — Trichophyton, Penicillium, Aspergillus)
  - figS3_nakaseomyces_consolidation.pdf  (Supp — accession consolidation)

Usage:
  # Two-way (EMITS vs Naive) — original behaviour
  python3 generate_zymo_figures.py \\
      --comparison results/emu_its_v2_species_comparison.tsv \\
      --outdir figures/

  # Three-way (EMITS vs Naive vs EMU) — new
  python3 generate_zymo_figures.py \\
      --comparison results/emu_its_v2_species_comparison.tsv \\
      --emu results/emu_atcc/atcc_emu_rel-abundance.tsv \\
      --outdir figures/

NOTE: This script expects a SPECIES-LEVEL comparison TSV with columns
  taxon (species name, e.g. 'Trichophyton mentagrophytes'), em_abundance,
  naive_abundance, n_accessions
  i.e. the same format as `emu_its_v2_species_comparison.tsv`.
"""

import argparse
import csv
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


EM_COLOR = "#378ADD"      # EMITS, blue (matches existing fig style)
NAIVE_COLOR = "#E24B4A"   # Naive, red (matches existing fig style)
EMU_COLOR = "#F59E0B"     # EMU, orange


def parse_comparison_tsv(filepath):
    """Parse the SPECIES-LEVEL EMITS comparison TSV.

    Expects columns: taxon (e.g. 'Trichophyton mentagrophytes'),
    em_abundance, naive_abundance, n_accessions
    """
    rows = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            taxon = row['taxon'].strip()
            species_name = taxon.replace('_', ' ').strip()
            genus = species_name.split()[0] if ' ' in species_name else species_name

            rows.append({
                'taxon_full': taxon,
                'species_name': species_name,
                'genus': genus,
                'em': float(row['em_abundance']),
                'naive': float(row['naive_abundance']),
                'n_accessions': int(row.get('n_accessions', 0) or 0),
            })
    return rows


def parse_emu_tsv(filepath):
    """Parse EMU rel-abundance.tsv into species-keyed dict."""
    emu = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            tax_id = row['tax_id']
            if tax_id in ('unmapped', 'mapped_filtered', 'mapped_unclassified'):
                continue
            sp_field = row.get('species', '').strip()
            if not sp_field:
                continue
            sp = sp_field.replace('_', ' ').strip()
            try:
                ab = float(row['abundance'])
            except (ValueError, KeyError):
                continue
            emu[sp] = emu.get(sp, 0.0) + ab
    return emu


# R1 (revision): species known to be present in the ATCC mock community.
# Pinned to the left-most position of each panel so the correct species reads
# first; congeners follow in descending abundance.
EXPECTED_SPECIES = {
    'Trichophyton': 'Trichophyton mentagrophytes',
    'Penicillium': 'Penicillium flavigenum',   # UNITE synonym of P. chrysogenum
    'Aspergillus': 'Aspergillus fumigatus',
    'Nakaseomyces': 'Nakaseomyces glabratus',
}


def get_genus_data(rows, genus_name, top_n=5, emu=None):
    """Get top N species for a genus, sorted by max abundance across methods.

    R1 (revision): the expected species is moved to the front. Only its single
    highest-ranked entry is moved -- sorting the whole list on an "is expected"
    key would pull every entry of that species forward and push the congeners
    out of the panel.
    """
    genus_rows = {r['species_name']: r for r in rows if r['genus'] == genus_name}

    if emu is not None:
        for sp, ab in emu.items():
            if sp.startswith(genus_name + ' ') and ab > 0.0001 and sp not in genus_rows:
                genus_rows[sp] = {
                    'species_name': sp,
                    'genus': genus_name,
                    'em': 0.0,
                    'naive': 0.0,
                    'n_accessions': 0,
                }

    out = []
    for sp, r in genus_rows.items():
        rcopy = dict(r)
        rcopy['emu'] = emu.get(sp, 0.0) if emu is not None else 0.0
        out.append(rcopy)

    out.sort(key=lambda r: -max(r['em'], r['naive'], r.get('emu', 0)))

    # R1 (revision): move the expected species to the left-most position.
    expected = EXPECTED_SPECIES.get(genus_name)
    if expected:
        idx = next((i for i, r in enumerate(out) if r['species_name'] == expected), None)
        if idx:  # found, and not already first
            out.insert(0, out.pop(idx))
    return out[:top_n]


def shorten_species(name, genus):
    if name.startswith(genus + ' '):
        return genus[0] + '.' + name[len(genus):]
    return name


# ═══════════════════════════════════════════════════════════════
# FIGURE 3: Within-genus species resolution (3 panels)
# ═══════════════════════════════════════════════════════════════
def generate_fig3(rows, outdir, emu=None):
    genera = [
        ('Trichophyton', 4),
        ('Penicillium', 5),
        ('Aspergillus', 4),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))

    for idx, (genus, top_n) in enumerate(genera):
        ax = axes[idx]
        data = get_genus_data(rows, genus, top_n, emu=emu)

        if not data:
            ax.text(0.5, 0.5, f'No {genus} data', transform=ax.transAxes, ha='center')
            continue

        species = [shorten_species(r['species_name'], genus) for r in data]
        em_vals = [r['em'] * 100 for r in data]
        naive_vals = [r['naive'] * 100 for r in data]
        emu_vals = [r.get('emu', 0) * 100 for r in data] if emu is not None else None

        x = np.arange(len(species))

        if emu is not None:
            width = 0.27
            bars_em = ax.bar(x - width, em_vals, width, label='EMITS',
                             color=EM_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
            bars_naive = ax.bar(x, naive_vals, width, label='Naive',
                                color=NAIVE_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
            bars_emu = ax.bar(x + width, emu_vals, width, label='EMU',
                              color=EMU_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
            bar_groups = [bars_em, bars_naive, bars_emu]
        else:
            width = 0.35
            bars_em = ax.bar(x - width / 2, em_vals, width, label='EMITS',
                             color=EM_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
            bars_naive = ax.bar(x + width / 2, naive_vals, width, label='Naive',
                                color=NAIVE_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
            bar_groups = [bars_em, bars_naive]

        ax.set_xlabel('')
        ax.set_ylabel('Abundance (%)' if idx == 0 else '')
        ax.set_title(f'$\\it{{{genus}}}$', fontsize=13, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(
            [f'$\\it{{{s}}}$' for s in species],
            rotation=35, ha='right', fontsize=9
        )
        ax.tick_params(axis='y', labelsize=9)

        # R1 (revision): legend moved to the top-left of the left-most panel so
        # it no longer sits over the bars.
        if idx == 0:
            ax.legend(fontsize=10, loc='upper left')

        ax.grid(axis='y', alpha=0.15)
        ax.set_axisbelow(True)

        max_val = max(max(em_vals), max(naive_vals),
                      max(emu_vals) if emu_vals else 0)
        ax.set_ylim(0, max_val * 1.30)

        for bar_group in bar_groups:
            for bar in bar_group:
                height = bar.get_height()
                if height > 0.05:
                    label_y = height + max_val * 0.02
                    if label_y < ax.get_ylim()[1] * 0.97:
                        ax.text(bar.get_x() + bar.get_width() / 2., label_y,
                                f'{height:.1f}', ha='center', va='bottom', fontsize=7,
                                clip_on=True)

    # R1 (revision): figure suptitle removed; that text now opens the LaTeX
    # caption. Panel titles (genus names) are retained, as R1 requested.
    plt.tight_layout()

    pdf_path = os.path.join(outdir, 'fig3_mock_within_genus.pdf')
    png_path = os.path.join(outdir, 'fig3_mock_within_genus.png')
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Fig 3 saved: {pdf_path}")


# ═══════════════════════════════════════════════════════════════
# FIGURE S3: Nakaseomyces species-level consolidation
# ═══════════════════════════════════════════════════════════════
def generate_figS3(rows, outdir, emu=None):
    nak_rows = [r for r in rows if r['genus'] == 'Nakaseomyces'
                and (r['em'] > 1e-5 or r['naive'] > 1e-5)]
    nak_rows.sort(key=lambda r: r['em'], reverse=True)

    if not nak_rows and (emu is None or not any(sp.startswith('Nakaseomyces ') for sp in emu)):
        print("WARNING: No Nakaseomyces data found")
        return

    labels = [shorten_species(r['species_name'], 'Nakaseomyces') for r in nak_rows]
    em_vals = [r['em'] * 100 for r in nak_rows]
    naive_vals = [r['naive'] * 100 for r in nak_rows]

    emu_vals = None
    if emu is not None:
        emu_vals = [emu.get(r['species_name'], 0.0) * 100 for r in nak_rows]
        seen = {r['species_name'] for r in nak_rows}
        for sp, ab in emu.items():
            if sp.startswith('Nakaseomyces ') and sp not in seen and ab > 1e-5:
                labels.append(shorten_species(sp, 'Nakaseomyces'))
                em_vals.append(0)
                naive_vals.append(0)
                emu_vals.append(ab * 100)

    fig, ax = plt.subplots(figsize=(9, max(4, len(labels) * 0.6 + 1.5)))

    y = np.arange(len(labels))

    if emu_vals is not None:
        height = 0.27
        ax.barh(y - height, em_vals, height, label='EMITS',
                color=EM_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
        ax.barh(y, naive_vals, height, label='Naive',
                color=NAIVE_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
        ax.barh(y + height, emu_vals, height, label='EMU',
                color=EMU_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
    else:
        height = 0.35
        ax.barh(y - height / 2, em_vals, height, label='EMITS',
                color=EM_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)
        ax.barh(y + height / 2, naive_vals, height, label='Naive',
                color=NAIVE_COLOR, alpha=0.9, edgecolor='white', linewidth=0.5)

    ax.set_yticks(y)
    ax.set_yticklabels([f'$\\it{{{s}}}$' for s in labels], fontsize=10)
    ax.set_xlabel('Abundance (%)', fontsize=11)
    title = ('$\\it{Nakaseomyces}$: EMITS vs Naive vs EMU consolidation'
             if emu_vals is not None else
             '$\\it{Nakaseomyces}$: EMITS consolidates across UNITE accessions')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(axis='x', alpha=0.15)
    ax.set_axisbelow(True)
    ax.invert_yaxis()

    for i in range(len(labels)):
        if em_vals[i] > 0.1:
            y_offset = -height if emu_vals is not None else -height / 2
            ax.text(em_vals[i] + 0.15, i + y_offset, f'{em_vals[i]:.2f}%',
                    va='center', fontsize=8, color='#185FA5')
        if naive_vals[i] > 0.1:
            y_offset = 0 if emu_vals is not None else height / 2
            ax.text(naive_vals[i] + 0.15, i + y_offset, f'{naive_vals[i]:.2f}%',
                    va='center', fontsize=8, color='#A32D2D')
        if emu_vals is not None and emu_vals[i] > 0.1:
            ax.text(emu_vals[i] + 0.15, i + height, f'{emu_vals[i]:.2f}%',
                    va='center', fontsize=8, color='#A66708')

    plt.tight_layout()

    pdf_path = os.path.join(outdir, 'figS3_nakaseomyces_consolidation.pdf')
    png_path = os.path.join(outdir, 'figS3_nakaseomyces_consolidation.png')
    fig.savefig(pdf_path, dpi=300, bbox_inches='tight')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Fig S3 saved: {pdf_path}")


# ═══════════════════════════════════════════════════════════════
# Genus-level summary
# ═══════════════════════════════════════════════════════════════
def print_genus_summary(rows, emu=None):
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
        emu_total = 0
        if emu is not None:
            emu_total = sum(ab for sp, ab in emu.items() if sp.startswith(genus + ' ')) * 100

        n_refs = len(genus_rows)
        top = sorted(genus_rows, key=lambda r: r['em'], reverse=True)

        print(f"\n  {genus} (expected ~10%, {n_refs} EMITS-detected species)")
        if emu is not None:
            print(f"    Genus total: EMITS={em_total:.2f}%  Naive={naive_total:.2f}%  EMU={emu_total:.2f}%")
        else:
            print(f"    Genus total: EMITS={em_total:.2f}%  Naive={naive_total:.2f}%")
        for r in top[:3]:
            if r['em'] > 0.0001 or r['naive'] > 0.0001:
                sp = shorten_species(r['species_name'], genus)
                line = f"    {sp:<30s}  EMITS={r['em']*100:.2f}%  Naive={r['naive']*100:.2f}%"
                if emu is not None:
                    line += f"  EMU={emu.get(r['species_name'], 0)*100:.2f}%"
                print(line)


def main():
    parser = argparse.ArgumentParser(description='Generate ATCC mock community figures')
    parser.add_argument('--comparison', type=str, required=True,
                        help='Path to species-level EMITS comparison TSV')
    parser.add_argument('--emu', type=str, default=None,
                        help='(Optional) Path to EMU rel-abundance.tsv for three-way comparison')
    parser.add_argument('--outdir', type=str, default='.',
                        help='Output directory for figures')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Reading EMITS comparison: {args.comparison}")
    rows = parse_comparison_tsv(args.comparison)
    print(f"  Loaded {len(rows)} taxa")

    emu = None
    if args.emu:
        print(f"Reading EMU output: {args.emu}")
        emu = parse_emu_tsv(args.emu)
        print(f"  Loaded {len(emu)} species")

    generate_fig3(rows, args.outdir, emu=emu)
    generate_figS3(rows, args.outdir, emu=emu)
    print_genus_summary(rows, emu=emu)

    print(f"\nAll figures saved to {args.outdir}/")


if __name__ == '__main__':
    main()
