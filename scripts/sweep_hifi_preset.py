#!/usr/bin/env python3
"""
Retune the pacbio-hifi preset filters against the simulated HiFi community.

R1 (revision) asked for PacBio benchmarking. The first benchmark showed that the
shipped `pacbio-hifi` preset is OUTPERFORMED by `ont-r10` on the same HiFi reads
(L1 6.43% vs 6.04%), because min_identity=0.95 and min_mapq=5 discard 25.9% of
alignments and collapse the candidate set from 187 taxa to 31 -- leaving EM with
nothing to weigh. This sweeps those two filters to find a defensible setting.

Usage:
    python3 scripts/sweep_hifi_preset.py \
        --paf results/sim_hifi/sim_vs_unite_hifi.paf \
        --truth results/sim_hifi/ground_truth.tsv \
        --outdir results/sim_hifi/sweep
"""

import argparse
import csv
import os
import subprocess
import sys

IDENTITIES = [0.80, 0.85, 0.90, 0.95]
MAPQS = [0, 5]


def read_truth(path):
    with open(path) as f:
        return {r["species"]: float(r["expected_abundance"])
                for r in csv.DictReader(f, delimiter="\t")}


def score(comparison_tsv, truth):
    """Total L1 error for EM and naive, plus false-positive count."""
    with open(comparison_tsv) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    em = {r["taxon"]: float(r["em_abundance"]) for r in rows}
    nv = {r["taxon"]: float(r["naive_abundance"]) for r in rows}
    allsp = set(truth) | set(em)
    l1_em = sum(abs(em.get(s, 0) - truth.get(s, 0)) for s in allsp) * 100
    l1_nv = sum(abs(nv.get(s, 0) - truth.get(s, 0)) for s in allsp) * 100
    fp = sum(v for s, v in em.items() if s not in truth) * 100
    fn = sum(1 for s in truth if em.get(s, 0) == 0)
    return l1_em, l1_nv, fp, fn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--outdir", default="results/sim_hifi/sweep")
    ap.add_argument("--binary", default="./target/release/emits")
    ap.add_argument("--preset", default="pacbio-hifi")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    truth = read_truth(args.truth)

    print(f"\nSweeping {args.preset}: min_identity x min_mapq")
    print(f"Ground truth: {len(truth)} species\n")
    header = (f"{'min_id':>7}{'mapq':>6}{'reads':>9}{'taxa':>6}{'iters':>7}"
              f"{'EM L1':>9}{'naive L1':>10}{'EM gain':>9}{'FP%':>8}{'FN':>4}")
    print(header)
    print("-" * len(header))

    results = []
    for ident in IDENTITIES:
        for mapq in MAPQS:
            tag = f"id{ident}_mq{mapq}"
            out = os.path.join(args.outdir, f"{tag}.tsv")
            cmd = [args.binary, "run",
                   "--preset", args.preset,
                   "--min-identity", str(ident),
                   "--min-mapq", str(mapq),
                   "-i", args.paf, "-o", out,
                   "--compare", "--rank", "species"]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                print(f"  FAILED {tag}: {proc.stderr.strip()[:160]}")
                continue

            blob = proc.stdout + proc.stderr
            reads = taxa = iters = "?"
            for line in blob.splitlines():
                if "unique reads" in line:
                    reads = line.split()[-3]
                if "candidate taxa" in line:
                    taxa = line.split()[-3]
                if "converged after" in line:
                    iters = line.split("converged after")[1].split()[0]
                elif "did not converge after" in line:
                    iters = line.split("after")[1].split()[0] + "*"

            comp = out.replace(".tsv", "_species_comparison.tsv")
            if not os.path.exists(comp):
                print(f"  MISSING {comp}")
                continue
            l1e, l1n, fp, fn = score(comp, truth)
            gain = (1 - l1e / l1n) * 100 if l1n else 0.0
            results.append((l1e, ident, mapq, gain))
            print(f"{ident:>7.2f}{mapq:>6}{reads:>9}{taxa:>6}{iters:>7}"
                  f"{l1e:>8.2f}%{l1n:>9.2f}%{gain:>8.1f}%{fp:>7.3f}%{fn:>4}")

    if results:
        results.sort()
        best_l1, best_id, best_mq, best_gain = results[0]
        print("\n" + "=" * 70)
        print(f"Lowest EM L1: min_identity={best_id}, min_mapq={best_mq} "
              f"-> {best_l1:.2f}% (EM gain {best_gain:.1f}%)")
        print("Shipped preset is min_identity=0.95, min_mapq=5.")
        print("ont-r10 on the same reads gave EM L1 6.04% (EM gain 6.5%).")
        print("=" * 70)
        print("\n* = EM did not converge within max_iterations.")
        print("Pick on EM L1, but check FN and convergence before changing the preset:")
        print("a setting that wins on L1 while dropping true species is not an improvement.")


if __name__ == "__main__":
    sys.exit(main())
