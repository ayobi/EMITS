#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# EMU-ITS Real Data Validation Pipeline
# ============================================================================
# Downloads ONT fungal ITS mock community data and UNITE database,
# runs minimap2 alignment, then tests EMU-ITS against known composition.
#
# Requirements:
#   - minimap2 (conda install -c bioconda minimap2)
#   - aws CLI (conda install -c conda-forge awscli) OR just use curl/wget
#   - emu-its (cargo build --release, binary at target/release/emu-its)
#   - ~10 GB free disk space
#
# Usage:
#   chmod +x scripts/run_real_data_validation.sh
#   ./scripts/run_real_data_validation.sh
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="${PROJECT_DIR}/data"
RESULTS_DIR="${PROJECT_DIR}/results"
EMU_ITS="${PROJECT_DIR}/target/release/emu-its"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log() { echo -e "${GREEN}[$(date '+%H:%M:%S')]${NC} $1"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] WARNING:${NC} $1"; }
info() { echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} $1"; }

# ── Check dependencies ──────────────────────────────────────────────
log "Checking dependencies..."

if ! command -v minimap2 &> /dev/null; then
    echo "minimap2 not found. Install with:"
    echo "  conda install -c bioconda minimap2"
    exit 1
fi

if [ ! -f "$EMU_ITS" ]; then
    log "Building emu-its..."
    cd "$PROJECT_DIR" && cargo build --release
fi

mkdir -p "$DATA_DIR" "$RESULTS_DIR"

# ── Step 1: Download ONT Fungal ITS Mock Community ──────────────────
# Source: ONT EPI2ME Open Data - ATCC Mycobiome Genomic DNA Mix
# 10 fungal species, even mix, R10.4.1 chemistry, GridION
# https://epi2me.nanoporetech.com/fungal_its_2025.09/
#
# We only need the basecalls (FASTQ), not the raw POD5 or analysis.
# The basecalls/ folder contains demultiplexed FASTQ files per barcode.
# Each barcode is a biological replicate.

FASTQ_DIR="${DATA_DIR}/ont_fungal_its"

if [ ! -d "$FASTQ_DIR" ]; then
    log "Downloading ONT fungal ITS mock community FASTQ files..."
    log "This is ~7.2 GB — downloading basecalls only."

    if command -v aws &> /dev/null; then
        aws s3 sync --no-sign-request \
            s3://ont-open-data/fungal_ITS_2025.09/basecalls \
            "$FASTQ_DIR"
    else
        warn "AWS CLI not found. You can install it with:"
        warn "  conda install -c conda-forge awscli"
        warn ""
        warn "Alternatively, browse and download manually from:"
        warn "  https://42basepairs.com/browse/s3/ont-open-data/fungal_ITS_2025.09/basecalls"
        warn ""
        warn "Place the FASTQ files in: ${FASTQ_DIR}/"
        exit 1
    fi
else
    log "ONT fungal ITS data already present at ${FASTQ_DIR}"
fi

# ── Step 2: Download UNITE Database ─────────────────────────────────
# General FASTA release - all fungi, dynamic species hypotheses
# https://unite.ut.ee/repository.php

UNITE_DIR="${DATA_DIR}/unite"
UNITE_FASTA="${UNITE_DIR}/sh_general_release_dynamic_fungi.fasta"

if [ ! -f "$UNITE_FASTA" ]; then
    log "Downloading UNITE general FASTA release..."
    mkdir -p "$UNITE_DIR"

    info "Please download the UNITE general FASTA release manually:"
    info "  1. Go to https://unite.ut.ee/repository.php"
    info "  2. Under 'General FASTA release', download the latest version"
    info "     (look for 'Fungi' — the dynamic version is recommended)"
    info "  3. Extract and place the .fasta file at:"
    info "     ${UNITE_FASTA}"
    info ""
    info "  Alternatively, if you have a direct URL from the UNITE DOI page:"
    info "    curl -L <URL> -o ${UNITE_DIR}/unite.tar.gz"
    info "    tar -xzf ${UNITE_DIR}/unite.tar.gz -C ${UNITE_DIR}"
    info "    # Then find and rename/symlink the .fasta file"
    info ""
    warn "UNITE requires you to agree to their license (CC BY-SA 4.0)."
    warn "Automated download is not provided to respect their terms."

    # Check if user has already downloaded it under a different name
    EXISTING=$(find "$UNITE_DIR" -name "*.fasta" -o -name "*.fa" 2>/dev/null | head -1)
    if [ -n "$EXISTING" ]; then
        log "Found existing FASTA: $EXISTING"
        UNITE_FASTA="$EXISTING"
    else
        exit 1
    fi
else
    log "UNITE database already present at ${UNITE_FASTA}"
fi

# ── Step 3: Prepare input FASTQ ─────────────────────────────────────
# Concatenate all FASTQ files from one barcode (one biological replicate)
# or all barcodes for a combined analysis

COMBINED_FASTQ="${DATA_DIR}/combined_its_reads.fastq"

if [ ! -f "$COMBINED_FASTQ" ]; then
    log "Combining FASTQ files..."

    # Find all .fastq.gz files in the basecalls directory
    FASTQ_FILES=$(find "$FASTQ_DIR" -name "*.fastq" -type f | sort)
    FASTQ_COUNT=$(echo "$FASTQ_FILES" | wc -l | tr -d ' ')

    if [ "$FASTQ_COUNT" -eq 0 ]; then
        warn "No FASTQ files found in ${FASTQ_DIR}"
        warn "Check the download and directory structure."
        exit 1
    fi

    log "Found ${FASTQ_COUNT} FASTQ files"
    cat $FASTQ_FILES > "$COMBINED_FASTQ"
    log "Combined into ${COMBINED_FASTQ}"
else
    log "Combined FASTQ already present"
fi

# ── Step 4: Build minimap2 index for UNITE ──────────────────────────
UNITE_INDEX="${UNITE_FASTA}.mmi"

if [ ! -f "$UNITE_INDEX" ]; then
    log "Building minimap2 index for UNITE..."
    minimap2 -d "$UNITE_INDEX" "$UNITE_FASTA"
else
    log "UNITE minimap2 index already present"
fi

# ── Step 5: Align reads to UNITE ────────────────────────────────────
# Key flags:
#   -c          : output CIGAR and alignment details
#   --secondary=yes : keep secondary alignments (needed for multi-mapping!)
#   -N 10       : output up to 10 secondary alignments per read
#   -p 0.9      : keep secondary if score >= 90% of primary
#   -K 500M     : batch size for memory efficiency

PAF_FILE="${RESULTS_DIR}/its_vs_unite.paf"

if [ ! -f "$PAF_FILE" ]; then
    log "Aligning ITS reads to UNITE with minimap2..."
    log "  Keeping secondary alignments for EM multi-mapping resolution"

    minimap2 \
        -c \
        --secondary=yes \
        -N 10 \
        -p 0.9 \
        -t $(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4) \
        "$UNITE_INDEX" \
        "$COMBINED_FASTQ" \
        > "$PAF_FILE"

    TOTAL_ALIGNMENTS=$(wc -l < "$PAF_FILE" | tr -d ' ')
    log "Generated ${TOTAL_ALIGNMENTS} alignments"
else
    log "PAF file already present at ${PAF_FILE}"
fi

# ── Step 6: Run EMU-ITS ─────────────────────────────────────────────
EM_OUTPUT="${RESULTS_DIR}/emu_its_abundances.tsv"
NAIVE_OUTPUT="${RESULTS_DIR}/naive_abundances.tsv"

log "Running EMU-ITS..."

# Run with comparison mode
"$EMU_ITS" run \
    --input "$PAF_FILE" \
    --output "$EM_OUTPUT" \
    --min-identity 0.8 \
    --compare

log "EMU-ITS complete!"

# ── Step 7: Also run naive counting for comparison ──────────────────
# (The --compare flag already generates this, but let's be explicit)

# ── Step 8: Print summary ───────────────────────────────────────────
echo ""
echo "============================================================================"
echo "  EMU-ITS Real Data Validation — Results Summary"
echo "============================================================================"
echo ""
echo "  Dataset:  ONT Fungal ITS Mock Community (ATCC Mycobiome gDNA Mix)"
echo "  Species:  10 fungal species, even mix (~10% each)"
echo "  Database: UNITE general FASTA release"
echo ""
echo "  Expected species:"
echo "    - Aspergillus fumigatus"
echo "    - Cryptococcus neoformans"
echo "    - Trichophyton interdigitale"
echo "    - Penicillium chrysogenum"
echo "    - Fusarium keratoplasticum"
echo "    - Candida albicans"
echo "    - Nakaseomyces glabratus"
echo "    - Malassezia globosa"
echo "    - Saccharomyces cerevisiae"
echo "    - Cutaneotrichosporon dermatis"
echo ""
echo "  Results:"
echo "    EM abundances:    ${EM_OUTPUT}"
echo "    Comparison:       ${RESULTS_DIR}/emu_its_abundances_comparison.tsv"
echo ""
echo "  Next: compare EM vs naive estimates against the known 10% composition."
echo "============================================================================"
