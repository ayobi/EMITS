# EMITS
[![bioconda downloads](https://img.shields.io/conda/dn/bioconda/emits.svg?label=bioconda%20downloads)](https://anaconda.org/bioconda/emits)
[![CI](https://github.com/ayobi/emits/actions/workflows/ci.yml/badge.svg)](https://github.com/ayobi/emits/actions/workflows/ci.yml)
[![Release](https://github.com/ayobi/emits/actions/workflows/release.yml/badge.svg)](https://github.com/ayobi/emits/actions/workflows/release.yml)
[![Docker](https://github.com/ayobi/emits/actions/workflows/docker.yml/badge.svg)](https://github.com/ayobi/emits/actions/workflows/docker.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![GitHub release](https://img.shields.io/github/v/release/ayobi/emits)](https://github.com/ayobi/emits/releases)
[![Docker Image](https://img.shields.io/badge/Docker-ghcr.io%2Fayobi%2Femits-blue?logo=docker)](https://ghcr.io/ayobi/emits)

Expectation-Maximization abundance estimation for fungal ITS communities from long-read sequencing

EMITS applies expectation-maximization (EM) to resolve ambiguous read-to-reference mappings in fungal ITS amplicon sequencing, producing probabilistic species-level abundance estimates from minimap2 alignments against the [UNITE](https://unite.ut.ee) database.

## Why EMITS?

Naive best-hit classification assigns each read entirely to its top-scoring reference. This fails when closely related species share similar ITS sequences, the "best hit" isn't always the correct one, especially with ONT sequencing noise. EMITS iteratively refines abundance estimates by considering *all* candidate alignments for each read, using community-level priors to resolve ambiguity.

**Key results:**
- **77–91% reduction** in L1 error under realistic alignment noise (simulations)
- **Correct within-genus species resolution** on the ONT ATCC fungal mock community (*Trichophyton*, *Penicillium*, *Aspergillus*)
- **54% reduction** in false positive species calls on a 21-species synthetic community
- **Accession consolidation**, resolves UNITE database redundancy (e.g., 13 entries for *N. glabratus*)

## Installation

### Bioconda (recommended)

```bash
conda install -c bioconda emits
```

### From source

```bash
git clone https://github.com/ayobi/emits.git
cd emits
cargo build --release
# Binary at target/release/emits
```

## Quick start

```bash
# 1. Align reads against UNITE with secondary alignments enabled
minimap2 --secondary=yes -N 10 -p 0.9 -c -t 8 \
    unite_database.fasta reads.fastq > alignments.paf

# 2. Run EMITS
emits run --input alignments.paf --output abundances.tsv --preset ont-r10

# 3. Compare EM vs naive (optional)
emits run --input alignments.paf --output abundances.tsv --preset ont-r10 --compare
```

## Usage

### `emits run` — Abundance estimation

```
emits run [OPTIONS] --input <PAF_FILE>

Options:
  -i, --input <PAF_FILE>       Input PAF file from minimap2 (use - for stdin)
  -o, --output <FILE>          Output TSV file [default: abundances.tsv]
      --preset <PRESET>        Platform preset: ont-r10, ont-r9, pacbio-hifi, ont-duplex
      --min-identity <FLOAT>   Minimum alignment identity (0.0-1.0) [default: 0.8]
      --min-mapq <INT>         Minimum mapping quality [default: 0]
      --max-iter <INT>         Maximum EM iterations [default: 100]
      --threshold <FLOAT>      EM convergence threshold [default: 1e-6]
      --rank <RANK>            Taxonomic rank for aggregation: species, genus [default: species]
      --compare                Also output naive best-hit estimates for comparison
```

### `emits preset` — Show platform parameters

```bash
# Display preset parameters and suggested minimap2 command
emits preset ont-r10
emits preset pacbio-hifi
```

### `emits simulate` — Built-in validation

```bash
# Run all simulation experiments
emits simulate --all

# Custom simulation
emits simulate --n-reads 10000 --cross-rate 0.5
```

## Platform presets

| Preset | Temperature (τ) | Min. identity | Use case |
|--------|:-:|:-:|----------|
| `ont-r10` | 0.50 | 0.80 | ONT R10.4.1 (default) |
| `ont-r9` | 0.80 | 0.70 | ONT R9.4.1 |
| `pacbio-hifi` | 0.15 | 0.95 | PacBio HiFi / CCS |
| `ont-duplex` | 0.20 | 0.90 | ONT Duplex |

The temperature parameter (τ) controls sensitivity to alignment score differences. Lower values make the model more decisive; higher values allow more ambiguity. Presets are empirically tuned for each platform's error profile.

## Output format

### Abundance TSV

| Column | Description |
|--------|-------------|
| `taxon` | UNITE reference identifier |
| `abundance` | EM-estimated relative abundance |
| `reads` | Estimated read count |

### Species-level aggregated TSV (with `--rank species`)

| Column | Description |
|--------|-------------|
| `species` | Species name (parsed from UNITE header) |
| `abundance` | Summed abundance across all accessions |
| `n_accessions` | Number of UNITE accessions contributing |

### Comparison TSV (with `--compare`)

| Column | Description |
|--------|-------------|
| `taxon` | UNITE reference identifier |
| `em_abundance` | EM-estimated abundance |
| `naive_abundance` | Naive best-hit abundance |
| `em_error` | Absolute error (EM) |
| `naive_error` | Absolute error (naive) |

## How it works

1. **Parse PAF** — Read minimap2 alignments, retaining all secondary alignments above identity/quality thresholds.
2. **Score-to-likelihood** — Convert normalized alignment scores to likelihoods via temperature-scaled exponential: `L(read, taxon) = exp(score / (τ × query_length))`.
3. **EM iteration** — E-step computes posterior assignment probabilities for each read; M-step updates taxon abundances by summing fractional assignments. Repeat until convergence.
4. **Taxonomic aggregation** — Parse UNITE headers and collapse abundances across accessions belonging to the same species.

## Complete pipeline with ITSxRust

EMITS pairs with [ITSxRust](https://github.com/ayobi/itsxrust) for a complete ONT fungal amplicon pipeline:

```bash
# 1. Extract ITS region
itsxrust --input reads.fastq --region full --preset ont --output its_reads.fastq

# 2. Align to UNITE
minimap2 --secondary=yes -N 10 -p 0.9 -c -t 8 \
    unite_database.fasta its_reads.fastq > alignments.paf

# 3. Estimate abundances
emits run --input alignments.paf --output abundances.tsv --preset ont-r10 --compare
```

## Requirements

- **Rust** ≥ 1.70 (for building from source)
- **minimap2** (for alignment; not bundled)
- **UNITE database** — download from [unite.ut.ee](https://unite.ut.ee/repository.php)

## Benchmarks

### Simulated communities

Under controlled simulations with tunable alignment noise, EM error remains flat (~0.014 L1) regardless of noise level, while naive best-hit error climbs to 0.275:

| Noise (±) | % wrong best-hits | Naive L1 | EM L1 | EM improvement |
|:-:|:-:|:-:|:-:|:-:|
| 0 | 0% | 0.023 | 0.021 | +9% |
| 30 | 6% | 0.069 | 0.014 | +80% |
| 60 | 28% | 0.172 | 0.014 | +92% |
| 100 | 45% | 0.275 | 0.014 | +95% |

### ONT mock community (ATCC Mycobiome, 10 species)

All 10 expected genera detected (99.95% abundance in expected taxa). EM resolves within-genus species where naive fails:

| Genus | Correct species | EM (%) | Naive (%) |
|-------|----------------|:-:|:-:|
| *Trichophyton* | *T. mentagrophytes* | **2.21** | 0.38 |
| | *T. simii* (wrong) | 1.24 | **3.12** |
| *Penicillium* | *P. flavigenum* | **2.85** | 0.78 |
| | *P. paneum* (wrong) | 0.32 | **2.06** |
| *Nakaseomyces* | *N. glabratus* (primary acc.) | **11.91** | 1.27 |
| | *N. glabratus* (2nd acc.) | 0.29 | **11.01** |

### Synthetic community (21 species)

- Overall L1: EM 7.48% vs naive 8.64% (**13.4% improvement**)
- False positives: EM 0.46% vs naive 1.01% (**54% reduction**)
- R²: EM 0.977 vs naive 0.970

## License

MIT License. See [LICENSE](LICENSE) for details.
