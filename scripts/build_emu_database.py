#!/usr/bin/env python3
"""
Build an EMU-compatible reference database from a UNITE general FASTA release.

EMU (Curry et al. 2022, Nature Methods) requires a database directory containing:
  - species_taxid.fasta  : sanitized FASTA (no spaces in headers, taxid as ID)
  - names_df.tsv         : taxonomic names for each taxid
  - nodes_df.tsv         : parent-child taxonomic relationships
  - unique_taxids.tsv    : list of unique species-level taxids
  - seq2tax.map          : maps each sequence header → species taxid

This script parses UNITE FASTA headers of the form:
    Species_name|accession|SH_id|reps_or_refs|k__Fungi;p__...;s__Species_name

and produces all five files in a single output directory.

Usage:
    python build_emu_database.py \\
        --unite data/unite/sh_general_release_dynamic_fungi.fasta \\
        --outdir data/unite_emu_db

Requirements:
    Python 3.7+, pandas

Notes on database compatibility:
    - We use sequential integer taxids since UNITE does not provide NCBI taxids
      consistently. Multiple UNITE accessions for the same species share a taxid.
    - Species-level taxids are at rank 7 (kingdom=1, phylum=2, ..., species=7).
    - This follows the Nagpal et al. 2024 / Erlandson et al. 2024 approach for
      adapting EMU to ITS data.

Citation:
    If you use this script, please cite EMU and the UNITE prior-art adaptations:
      Curry et al. 2022, Nature Methods 19:845-853
      Nagpal et al. 2024 / Erlandson et al. 2024 (ITS adaptations)
"""

import argparse
import os
import re
import sys
from collections import defaultdict


def parse_unite_header(header):
    """
    Parse a UNITE FASTA header into structured fields.

    Expected format:
        Species_name|accession|SH_id|reps_or_refs|k__;p__;c__;o__;f__;g__;s__

    Returns a dict with keys: species, accession, sh_id, type, taxonomy
    or None if the header is malformed.
    """
    parts = header.split("|")
    if len(parts) < 5:
        return None

    species_field = parts[0].strip()
    accession = parts[1].strip()
    sh_id = parts[2].strip()
    seq_type = parts[3].strip()  # 'reps' or 'refs'
    taxonomy_str = parts[4].strip()

    # Parse the taxonomy string: k__Fungi;p__Ascomycota;...;s__Species_name
    ranks = {}
    for rank_field in taxonomy_str.split(";"):
        rank_field = rank_field.strip()
        m = re.match(r"^([kpcofgs])__(.*)$", rank_field)
        if m:
            rank_letter, name = m.group(1), m.group(2).strip()
            ranks[rank_letter] = name if name else "unidentified"

    return {
        "species_field": species_field,
        "accession": accession,
        "sh_id": sh_id,
        "type": seq_type,
        "kingdom": ranks.get("k", "unidentified"),
        "phylum": ranks.get("p", "unidentified"),
        "class": ranks.get("c", "unidentified"),
        "order": ranks.get("o", "unidentified"),
        "family": ranks.get("f", "unidentified"),
        "genus": ranks.get("g", "unidentified"),
        "species": ranks.get("s", "unidentified"),
    }


def parse_unite_fasta(fasta_path):
    """Yield (header, sequence) tuples from a UNITE FASTA file."""
    current_header = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header and current_seq:
                    yield current_header, "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header and current_seq:
            yield current_header, "".join(current_seq)


def assign_taxids(parsed_records):
    """
    Walk the parsed records and assign integer taxids per unique taxonomic
    lineage at each rank. Uses a hierarchical scheme so children point at parents.

    Returns:
        taxid_map  : dict[(rank_index, lineage_tuple)] -> taxid
        record_taxids : list aligned with parsed_records of species-level taxids
    """
    rank_order = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxid_map = {}
    next_taxid = [2]  # 1 is reserved for the root

    # Root node
    taxid_map[(0, ())] = 1

    record_taxids = []

    for rec in parsed_records:
        # Build the lineage tuple at each rank, register each unique one
        species_taxid = None
        for depth, rank in enumerate(rank_order, start=1):
            lineage = tuple(rec[r] for r in rank_order[:depth])
            key = (depth, lineage)
            if key not in taxid_map:
                taxid_map[key] = next_taxid[0]
                next_taxid[0] += 1
            if depth == len(rank_order):
                species_taxid = taxid_map[key]
        record_taxids.append(species_taxid)

    return taxid_map, record_taxids


def write_emu_database(parsed_records, sequences, taxid_map, record_taxids,
                       original_headers, outdir):
    """
    Write the five files EMU needs:
      species_taxid.fasta, names_df.tsv, nodes_df.tsv, unique_taxids.tsv, seq2tax.map
    """
    os.makedirs(outdir, exist_ok=True)
    rank_order = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    rank_names = {
        1: "superkingdom",
        2: "phylum",
        3: "class",
        4: "order",
        5: "family",
        6: "genus",
        7: "species",
    }

    # ---- species_taxid.fasta ----
    # EMU expects sequence IDs to be unique. We use the original UNITE accession
    # as the visible header (sanitized) and rely on seq2tax.map for the taxid link.
    fasta_out = os.path.join(outdir, "species_taxid.fasta")
    seq2tax_out = os.path.join(outdir, "seq2tax.map")

    print(f"  Writing {fasta_out}")
    print(f"  Writing {seq2tax_out}")

    with open(fasta_out, "w") as fasta_f, open(seq2tax_out, "w") as map_f:
        # seq2tax.map header (tab-separated: seq_id, taxid)
        for orig_header, seq, taxid in zip(original_headers, sequences, record_taxids):
            # Sanitize: replace whitespace with underscores; EMU forbids spaces
            sanitized = orig_header.replace(" ", "_").replace("\t", "_")
            # Use accession as the unique ID — this is what EMU's minimap2 will see
            # We keep the full sanitized header to preserve traceability
            seq_id = sanitized
            fasta_f.write(f">{seq_id}\n{seq}\n")
            map_f.write(f"{seq_id}\t{taxid}\n")

    # ---- names_df.tsv ----
    # EMU's names_df.tsv: tab-separated, columns include tax_id and name fields
    # at every rank. Format follows EMU's expected schema:
    #   tax_id, species, genus, family, order, class, phylum, superkingdom
    names_out = os.path.join(outdir, "names_df.tsv")
    print(f"  Writing {names_out}")

    # Build a reverse lookup: taxid -> lineage at each rank
    # We only need entries for species-level taxids since that's what seq2tax.map points at
    taxid_to_lineage = {}
    for (depth, lineage), tid in taxid_map.items():
        if depth == 7:  # species level
            taxid_to_lineage[tid] = lineage

    with open(names_out, "w") as f:
        # Header line — EMU uses these column names
        f.write("tax_id\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom\n")
        for tid in sorted(taxid_to_lineage.keys()):
            lineage = taxid_to_lineage[tid]
            # lineage is (kingdom, phylum, class, order, family, genus, species)
            kingdom, phylum, cls, order, family, genus, species = lineage
            # Order in EMU's TSV: species, genus, family, order, class, phylum, superkingdom
            f.write(f"{tid}\t{species}\t{genus}\t{family}\t{order}\t{cls}\t{phylum}\t{kingdom}\n")

    # ---- nodes_df.tsv ----
    # Parent-child relationships and rank labels
    nodes_out = os.path.join(outdir, "nodes_df.tsv")
    print(f"  Writing {nodes_out}")

    with open(nodes_out, "w") as f:
        f.write("tax_id\tparent_tax_id\trank\n")
        # Root
        f.write("1\t1\troot\n")
        # All other taxids
        for (depth, lineage), tid in sorted(taxid_map.items(), key=lambda x: x[1]):
            if tid == 1:
                continue
            # Parent is the same lineage truncated by one rank
            if depth == 1:
                parent_tid = 1  # root
            else:
                parent_lineage = lineage[:-1]
                parent_key = (depth - 1, parent_lineage)
                parent_tid = taxid_map.get(parent_key, 1)
            rank_label = rank_names.get(depth, "no_rank")
            f.write(f"{tid}\t{parent_tid}\t{rank_label}\n")

    # ---- unique_taxids.tsv ----
    # List of all species-level taxids (one per line)
    unique_out = os.path.join(outdir, "unique_taxids.tsv")
    print(f"  Writing {unique_out}")

    species_taxids = sorted(set(record_taxids))
    with open(unique_out, "w") as f:
        f.write("tax_id\n")
        for tid in species_taxids:
            f.write(f"{tid}\n")

    print(f"\n  Total sequences:        {len(sequences):,}")
    print(f"  Total taxa (all ranks): {len(taxid_map):,}")
    print(f"  Unique species taxids:  {len(species_taxids):,}")


def main():
    parser = argparse.ArgumentParser(
        description="Build an EMU-compatible database from a UNITE FASTA release",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--unite", required=True, help="Path to UNITE FASTA file")
    parser.add_argument("--outdir", required=True, help="Output directory for EMU database files")
    parser.add_argument("--skip-malformed", action="store_true",
                        help="Skip sequences with unparseable headers (default: error)")
    args = parser.parse_args()

    if not os.path.isfile(args.unite):
        print(f"ERROR: UNITE FASTA not found: {args.unite}", file=sys.stderr)
        sys.exit(1)

    print(f"=" * 65)
    print(f"  UNITE → EMU database builder")
    print(f"=" * 65)
    print(f"\n[1] Parsing UNITE FASTA: {args.unite}")

    sequences = []
    parsed_records = []
    original_headers = []
    skipped = 0

    for header, seq in parse_unite_fasta(args.unite):
        rec = parse_unite_header(header)
        if rec is None:
            if args.skip_malformed:
                skipped += 1
                continue
            print(f"ERROR: Malformed UNITE header: {header[:80]}", file=sys.stderr)
            print(f"       Use --skip-malformed to ignore these.", file=sys.stderr)
            sys.exit(1)
        sequences.append(seq)
        parsed_records.append(rec)
        original_headers.append(header)

    print(f"  Parsed {len(parsed_records):,} sequences")
    if skipped:
        print(f"  Skipped {skipped:,} malformed entries")

    print(f"\n[2] Assigning taxids (hierarchical, sequential integers)")
    taxid_map, record_taxids = assign_taxids(parsed_records)

    print(f"\n[3] Writing EMU database files to {args.outdir}")
    write_emu_database(parsed_records, sequences, taxid_map, record_taxids,
                       original_headers, args.outdir)

    print(f"\n{'=' * 65}")
    print(f"  Database ready! Use with EMU as:")
    print(f"{'=' * 65}")
    print(f"""
  export EMU_DATABASE_DIR={os.path.abspath(args.outdir)}

  emu abundance \\
      --type map-ont \\
      --N 10 \\
      --keep-counts \\
      --threads 8 \\
      --output-dir results/emu_synthetic \\
      --output-basename synthetic_emu \\
      data/simulated_community/simulated_community.fastq

  # For the ATCC mock community:
  emu abundance \\
      --type map-ont \\
      --N 10 \\
      --keep-counts \\
      --threads 8 \\
      --output-dir results/emu_atcc \\
      --output-basename atcc_emu \\
      data/combined_its_reads.fastq
    """)


if __name__ == "__main__":
    main()
