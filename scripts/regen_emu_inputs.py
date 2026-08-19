#!/usr/bin/env python3
"""
Generate a fresh seq2tax.map and taxonomy.tsv directly from a UNITE FASTA,
using Biopython to parse headers exactly as EMU will.

This bypasses the original build_emu_database.py output entirely, fixing
any parsing inconsistencies between the original raw-text parser and EMU's
Biopython-based reader.

Usage:
    python regen_emu_inputs.py \\
        --fasta data/unite/sh_general_release_dynamic_19.02.2025.fasta \\
        --outdir data/unite_emu_inputs
"""

import argparse
import os
import re
import sys
from Bio import SeqIO


RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def parse_unite_taxonomy(taxonomy_str):
    """Parse 'k__Fungi;p__Ascomycota;...;s__Species_name' into a dict."""
    ranks = {r: "unidentified" for r in RANK_ORDER}
    rank_letters = {"k": "kingdom", "p": "phylum", "c": "class",
                    "o": "order", "f": "family", "g": "genus", "s": "species"}
    for field in taxonomy_str.split(";"):
        field = field.strip()
        m = re.match(r"^([kpcofgs])__(.*)$", field)
        if m:
            letter, name = m.group(1), m.group(2).strip() or "unidentified"
            ranks[rank_letters[letter]] = name
    return ranks


def parse_record_header(record):
    """
    Parse the UNITE header from a Biopython SeqRecord.
    Returns (record_id, taxonomy_dict) or (None, None) if malformed.

    UNITE headers don't have spaces, so record.id == record.description.
    Format: Species|accession|SH_id|reps_or_refs|k__;p__;...;s__
    """
    header = record.description  # full header
    parts = header.split("|")
    if len(parts) < 5:
        return None, None
    taxonomy = parse_unite_taxonomy(parts[4])
    return header, taxonomy


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Pass 1: parse all records, assign hierarchical taxids
    print(f"[1] Parsing FASTA: {args.fasta}")
    record_taxonomies = []  # list of (header, taxonomy_dict)
    skipped = 0

    for record in SeqIO.parse(args.fasta, "fasta"):
        header, taxonomy = parse_record_header(record)
        if header is None:
            skipped += 1
            continue
        record_taxonomies.append((header, taxonomy))

    print(f"  Parsed {len(record_taxonomies):,} records")
    if skipped:
        print(f"  Skipped {skipped:,} malformed headers")

    # Assign hierarchical taxids — same lineage = same taxid at every rank
    print(f"\n[2] Assigning taxids")
    lineage_to_taxid = {}
    next_id = 2  # 1 reserved for root
    record_species_taxids = []

    for header, tax in record_taxonomies:
        species_taxid = None
        for depth in range(1, len(RANK_ORDER) + 1):
            lineage = tuple(tax[r] for r in RANK_ORDER[:depth])
            if lineage not in lineage_to_taxid:
                lineage_to_taxid[lineage] = next_id
                next_id += 1
            if depth == len(RANK_ORDER):
                species_taxid = lineage_to_taxid[lineage]
        record_species_taxids.append(species_taxid)

    n_species = len({t for t in record_species_taxids})
    print(f"  Assigned {n_species:,} unique species-level taxids")

    # Write seq2tax.map (headerless, two columns)
    seq2tax_path = os.path.join(args.outdir, "seq2tax.map")
    print(f"\n[3] Writing {seq2tax_path}")
    with open(seq2tax_path, "w") as f:
        for (header, _), taxid in zip(record_taxonomies, record_species_taxids):
            f.write(f"{header}\t{taxid}\n")

    # Write taxonomy.tsv (tax_id + ranks; one row per species-level taxid)
    tax_path = os.path.join(args.outdir, "taxonomy.tsv")
    print(f"\n[4] Writing {tax_path}")

    # Build species_taxid -> lineage map (deduplicated)
    species_lineages = {}
    for (header, tax), tid in zip(record_taxonomies, record_species_taxids):
        if tid not in species_lineages:
            species_lineages[tid] = tax

    with open(tax_path, "w") as f:
        # Header row — first col MUST be tax_id
        f.write("tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\n")
        for tid in sorted(species_lineages.keys()):
            tax = species_lineages[tid]
            f.write(f"{tid}\t{tax['species']}\t{tax['genus']}\t{tax['family']}\t"
                    f"{tax['order']}\t{tax['class']}\t{tax['phylum']}\t{tax['kingdom']}\n")

    print(f"\n[5] Done.")
    print(f"  Records written:        {len(record_taxonomies):,}")
    print(f"  Unique species taxids:  {n_species:,}")
    print(f"\n  Now run EMU's database builder:")
    print(f"\n  emu build-database data/unite_emu_db \\")
    print(f"      --sequences {args.fasta} \\")
    print(f"      --seq2tax {seq2tax_path} \\")
    print(f"      --taxonomy-list {tax_path}")


if __name__ == "__main__":
    main()
