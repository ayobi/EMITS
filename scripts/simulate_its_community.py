#!/usr/bin/env python3
"""
Simulate a synthetic ITS mock community from UNITE reference sequences.

Workflow:
1. Search UNITE FASTA for target species, pick one representative per species
2. Use badread to simulate ONT reads with realistic error profiles
3. Combine into a single FASTQ with known ground truth
4. Ready to align with minimap2 → run EMU-ITS

Usage:
    python simulate_its_community.py \
        --unite /path/to/unite_refs.fasta \
        --outdir /path/to/output \
        --total-reads 50000 \
        --seed 42
"""

import argparse
import os
import sys
import random
import subprocess
import re
from collections import defaultdict

# === COMMUNITY DESIGN ===
COMMUNITY = {
    # Multi-species genera (6 genera, 15 species) — where EM shines
    "Aspergillus niger":            0.10,
    "Aspergillus tubingensis":      0.04,
    "Aspergillus fumigatus":        0.06,
    "Fusarium oxysporum":           0.08,
    "Fusarium solani":              0.05,
    "Fusarium graminearum":         0.03,
    "Candida albicans":             0.07,
    "Candida tropicalis":           0.04,
    "Candida parapsilosis":         0.02,
    "Penicillium chrysogenum":      0.06,
    "Penicillium expansum":         0.03,
    "Trichoderma harzianum":        0.05,
    "Trichoderma viride":           0.03,
    "Cladosporium cladosporioides": 0.04,
    "Cladosporium herbarum":        0.02,
    # Singleton genera — baseline
    "Saccharomyces cerevisiae":     0.08,
    "Botrytis cinerea":             0.06,
    "Alternaria alternata":         0.05,
    "Rhizopus oryzae":              0.03,
    "Mucor racemosus":              0.02,
    # Rare species — sensitivity test
    "Cryptococcus neoformans":      0.02,
    "Malassezia globosa":           0.01,
    "Trichophyton rubrum":          0.01,
}


def parse_unite_fasta(fasta_path):
    """Parse UNITE FASTA, return dict of header -> sequence."""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_seq:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    if current_header and current_seq:
        sequences[current_header] = "".join(current_seq)
    
    return sequences


def find_species_in_unite(sequences, target_species):
    """
    For each target species, find matching UNITE entries.
    Returns dict: species_name -> [(header, sequence), ...]
    """
    matches = defaultdict(list)
    
    for header, seq in sequences.items():
        # UNITE header format: Species|accession|SH|type|k__;p__;...;s__Species_name
        # Extract species from the taxonomy string
        parts = header.split("|")
        if len(parts) >= 5:
            taxonomy = parts[4]
            # Extract s__Species_name
            s_match = re.search(r's__(\S+)', taxonomy)
            if s_match:
                db_species = s_match.group(1).replace("_", " ")
                for target in target_species:
                    if target.lower() == db_species.lower():
                        matches[target].append((header, seq))
    
    return matches


def select_representatives(matches, seed=42):
    """
    For each species, select one representative sequence.
    Prefer: longer sequences, 'reps' type, no ambiguous bases.
    """
    random.seed(seed)
    representatives = {}
    
    for species, entries in matches.items():
        # Score each entry
        scored = []
        for header, seq in entries:
            score = 0
            # Prefer longer sequences
            score += len(seq)
            # Prefer 'reps' (representative) over 'refs'
            if "|reps|" in header:
                score += 10000
            # Penalize ambiguous bases
            ambig = sum(1 for c in seq.upper() if c not in "ACGT")
            score -= ambig * 100
            # Penalize very short sequences
            if len(seq) < 400:
                score -= 50000
            scored.append((score, header, seq))
        
        scored.sort(key=lambda x: -x[0])
        if scored:
            representatives[species] = (scored[0][1], scored[0][2])
    
    return representatives


def write_reference_fasta(representatives, outpath):
    """Write selected references to a FASTA file (one per species)."""
    with open(outpath, "w") as f:
        for species, (header, seq) in sorted(representatives.items()):
            f.write(f">{header}\n{seq}\n")
    print(f"  Wrote {len(representatives)} reference sequences to {outpath}")


def write_ground_truth(community, total_reads, outpath):
    """Write ground truth TSV."""
    with open(outpath, "w") as f:
        f.write("species\texpected_abundance\texpected_reads\n")
        for sp, ab in sorted(community.items(), key=lambda x: -x[1]):
            f.write(f"{sp}\t{ab:.4f}\t{int(ab * total_reads)}\n")
    print(f"  Wrote ground truth to {outpath}")


def simulate_reads_badread(representatives, community, total_reads, outdir, seed=42):
    """
    Simulate ONT reads using badread for each species proportionally.
    """
    combined_fastq = os.path.join(outdir, "simulated_community.fastq")
    read_origin = os.path.join(outdir, "read_origins.tsv")
    
    all_reads = []
    origins = []
    
    for species, (header, seq) in representatives.items():
        if species not in community:
            continue
        
        n_reads = int(community[species] * total_reads)
        if n_reads == 0:
            continue
        
        # Write temporary single-sequence FASTA for badread
        tmp_fasta = os.path.join(outdir, "tmp_ref.fasta")
        with open(tmp_fasta, "w") as f:
            f.write(f">{species.replace(' ', '_')}\n{seq}\n")
        
        # Calculate quantity: we want n_reads reads of roughly len(seq) each
        # badread quantity = total bases = n_reads * len(seq)
        quantity = n_reads * len(seq)
        
        # Run badread simulate
        # Using ONT R10 error model approximation
        cmd = [
            "badread", "simulate",
            "--reference", tmp_fasta,
            "--quantity", str(quantity),
            "--identity", "92,98,4",  # ONT R10.4.1 typical: mean 95%, sd ~3%
            "--length", f"{len(seq)},{int(len(seq)*0.2)}",  # mean=reflen, sd=20%
            "--error_model", "random",  # random errors (no specific model needed)
            "--qscore_model", "nanopore2020",
            "--seed", str(seed),
            "--start_adapter_seq", "",
            "--end_adapter_seq", "",
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            # Parse output FASTQ from stdout
            lines = result.stdout.strip().split("\n")
            read_count = 0
            i = 0
            while i < len(lines) - 3:
                if lines[i].startswith("@"):
                    # Rename read to include species origin
                    read_id = f"@sim_{species.replace(' ', '_')}_{read_count}"
                    all_reads.append(f"{read_id}\n{lines[i+1]}\n+\n{lines[i+3]}\n")
                    origins.append(f"{read_id[1:]}\t{species}\n")
                    read_count += 1
                    i += 4
                else:
                    i += 1
            
            print(f"  {species}: {read_count} reads simulated (target: {n_reads})")
        
        except FileNotFoundError:
            print(f"\n  ERROR: badread not found! Install with: pip install badread")
            print(f"  Alternatively, use a simple simulation (see below)")
            return None
        except subprocess.CalledProcessError as e:
            print(f"  WARNING: badread failed for {species}: {e.stderr[:200]}")
            # Fall back to simple simulation
            print(f"  Using simple simulation fallback for {species}")
            for j in range(n_reads):
                mutated = simple_simulate_read(seq, identity=0.95, seed=seed+j)
                read_id = f"@sim_{species.replace(' ', '_')}_{j}"
                qual = "I" * len(mutated)  # placeholder quality
                all_reads.append(f"{read_id}\n{mutated}\n+\n{qual}\n")
                origins.append(f"{read_id[1:]}\t{species}\n")
            print(f"  {species}: {n_reads} reads (simple simulation)")
    
    # Clean up
    tmp = os.path.join(outdir, "tmp_ref.fasta")
    if os.path.exists(tmp):
        os.remove(tmp)
    
    # Shuffle reads
    random.seed(seed)
    combined = list(zip(all_reads, origins))
    random.shuffle(combined)
    all_reads, origins = zip(*combined) if combined else ([], [])
    
    # Write combined FASTQ
    with open(combined_fastq, "w") as f:
        f.writelines(all_reads)
    
    # Write origins
    with open(read_origin, "w") as f:
        f.write("read_id\tspecies\n")
        f.writelines(origins)
    
    print(f"\n  Total reads: {len(all_reads)}")
    print(f"  Output: {combined_fastq}")
    print(f"  Origins: {read_origin}")
    
    return combined_fastq


def simple_simulate_read(seq, identity=0.95, seed=0):
    """Simple read simulation: introduce random substitutions/indels."""
    rng = random.Random(seed)
    result = []
    error_rate = 1.0 - identity
    
    for base in seq:
        r = rng.random()
        if r < error_rate * 0.6:  # substitution
            result.append(rng.choice([b for b in "ACGT" if b != base.upper()]))
        elif r < error_rate * 0.8:  # deletion
            pass  # skip this base
        elif r < error_rate:  # insertion
            result.append(rng.choice("ACGT"))
            result.append(base)
        else:
            result.append(base)
    
    return "".join(result)


def simulate_simple(representatives, community, total_reads, outdir, seed=42):
    """
    Simple simulation without badread — introduce random errors directly.
    Use this if badread is not installed.
    """
    combined_fastq = os.path.join(outdir, "simulated_community.fastq")
    read_origin = os.path.join(outdir, "read_origins.tsv")
    
    all_reads = []
    origins = []
    rng = random.Random(seed)
    
    for species, (header, seq) in sorted(representatives.items()):
        if species not in community:
            continue
        
        n_reads = int(community[species] * total_reads)
        
        for j in range(n_reads):
            # Vary read length slightly (±15%)
            trim_start = rng.randint(0, max(0, int(len(seq) * 0.05)))
            trim_end = rng.randint(0, max(0, int(len(seq) * 0.05)))
            subseq = seq[trim_start:len(seq)-trim_end] if trim_end > 0 else seq[trim_start:]
            
            # Simulate with ~92-98% identity (ONT R10 range)
            identity = rng.gauss(0.95, 0.02)
            identity = max(0.85, min(0.99, identity))
            
            mutated = simple_simulate_read(subseq, identity=identity, seed=seed+j*1000+hash(species)%10000)
            
            read_id = f"@sim_{species.replace(' ', '_')}_{j}"
            qual = "I" * len(mutated)
            all_reads.append(f"{read_id}\n{mutated}\n+\n{qual}\n")
            origins.append(f"{read_id[1:]}\t{species}\n")
        
        print(f"  {species}: {n_reads} reads (simple simulation, ~95% identity)")
    
    # Shuffle
    combined = list(zip(all_reads, origins))
    rng.shuffle(combined)
    all_reads, origins = zip(*combined) if combined else ([], [])
    
    with open(combined_fastq, "w") as f:
        f.writelines(all_reads)
    
    with open(read_origin, "w") as f:
        f.write("read_id\tspecies\n")
        f.writelines(origins)
    
    print(f"\n  Total reads: {len(all_reads)}")
    print(f"  Output: {combined_fastq}")
    
    return combined_fastq


def main():
    parser = argparse.ArgumentParser(description="Simulate ITS mock community from UNITE")
    parser.add_argument("--unite", required=True, help="Path to UNITE FASTA database")
    parser.add_argument("--outdir", default="simulated_community", help="Output directory")
    parser.add_argument("--total-reads", type=int, default=50000, help="Total reads to simulate")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--simple", action="store_true", help="Use simple simulation (no badread needed)")
    args = parser.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print("=" * 65)
    print("  EMU-ITS Synthetic Community Simulation")
    print("=" * 65)
    
    # 1. Parse UNITE
    print(f"\n[1] Parsing UNITE database: {args.unite}")
    sequences = parse_unite_fasta(args.unite)
    print(f"  Loaded {len(sequences)} sequences")
    
    # 2. Find target species
    print(f"\n[2] Searching for {len(COMMUNITY)} target species...")
    matches = find_species_in_unite(sequences, list(COMMUNITY.keys()))
    
    found = []
    missing = []
    for sp in sorted(COMMUNITY.keys()):
        if sp in matches and len(matches[sp]) > 0:
            found.append(sp)
            print(f"  ✓ {sp}: {len(matches[sp])} entries")
        else:
            missing.append(sp)
            print(f"  ✗ {sp}: NOT FOUND")
    
    if missing:
        print(f"\n  WARNING: {len(missing)} species not found in UNITE!")
        print(f"  Consider substituting or checking species names.")
        
        # Renormalize abundances for found species only
        found_total = sum(COMMUNITY[sp] for sp in found)
        adjusted = {sp: COMMUNITY[sp] / found_total for sp in found}
        print(f"  Renormalized {len(found)} species (was {found_total:.2f}, now 1.00)")
    else:
        adjusted = COMMUNITY.copy()
    
    # 3. Select representatives
    print(f"\n[3] Selecting representative sequences...")
    representatives = select_representatives(matches, seed=args.seed)
    print(f"  Selected {len(representatives)} representatives")
    
    # Write reference FASTA
    ref_path = os.path.join(args.outdir, "reference_sequences.fasta")
    write_reference_fasta(representatives, ref_path)
    
    # 4. Write ground truth
    print(f"\n[4] Writing ground truth...")
    truth_path = os.path.join(args.outdir, "ground_truth.tsv")
    write_ground_truth(adjusted, args.total_reads, truth_path)
    
    # 5. Simulate reads
    print(f"\n[5] Simulating {args.total_reads} reads...")
    if args.simple:
        fastq = simulate_simple(representatives, adjusted, args.total_reads, args.outdir, args.seed)
    else:
        fastq = simulate_reads_badread(representatives, adjusted, args.total_reads, args.outdir, args.seed)
        if fastq is None:
            print("\n  Falling back to simple simulation...")
            fastq = simulate_simple(representatives, adjusted, args.total_reads, args.outdir, args.seed)
    
    # 6. Print next steps
    print(f"\n{'=' * 65}")
    print(f"  DONE! Next steps:")
    print(f"{'=' * 65}")
    print(f"""
  # Align against UNITE
  minimap2 -c --secondary=yes -N 10 -p 0.9 -x map-ont \\
      {args.unite} \\
      {fastq} \\
      > {args.outdir}/sim_vs_unite.paf

  # Run EMU-ITS
  cargo run --release -- run \\
      --preset ont-r10 \\
      --input {args.outdir}/sim_vs_unite.paf \\
      --output {args.outdir}/emu_its_results.tsv \\
      --compare --rank species
    """)


if __name__ == "__main__":
    main()