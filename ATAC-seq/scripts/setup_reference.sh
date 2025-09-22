#!/bin/bash
# setup_reference.sh - Download and index GRCh38 reference genome

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

REF_DIR="reference_genome"
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
REF_FILE="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
INDEX_NAME="GRCh38_index"

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

echo "Setting up reference genome..."

# Create reference directory
mkdir -p "$REF_DIR"
cd "$REF_DIR"

# Download reference genome if not already present
if [[ ! -f "$REF_FILE" ]]; then
    echo "Downloading GRCh38 reference genome..."
    wget "$REF_URL"
else
    echo "Reference genome already downloaded"
fi

# Extract if needed
if [[ ! -f "${REF_FILE%.gz}" ]]; then
    echo "Extracting reference genome..."
    gunzip -k "$REF_FILE"
fi

# Build Bowtie2 index if not already present
if [[ ! -f "${INDEX_NAME}.1.bt2" ]]; then
    echo "Building Bowtie2 index..."
    bowtie2-build "${REF_FILE%.gz}" "$INDEX_NAME"
else
    echo "Bowtie2 index already exists"
fi

echo "Reference setup complete!"