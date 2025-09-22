#!/bin/bash
# setup_blacklist.sh - Download ENCODE blacklist for hg38

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

BLACKLIST_DIR="reference_genome"  # Directory to store blacklist files
BLACKLIST_VERSION="v2"  # Version of the blacklist
BLACKLIST_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
BLACKLIST_FILE="hg38-blacklist.v2.bed.gz"  # Downloaded file name
EXTRACTED_FILE="hg38-blacklist.v2.bed"  # Extracted file name

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

echo "Setting up hg38 blacklist..."

# Create blacklist directory
mkdir -p "$BLACKLIST_DIR"
cd "$BLACKLIST_DIR"

# Download blacklist if not already present
if [[ ! -f "$BLACKLIST_FILE" ]]; then
    echo "Downloading hg38 blacklist $BLACKLIST_VERSION..."
    wget "$BLACKLIST_URL"
else
    echo "Blacklist file already exists"
fi

# Extract if needed
if [[ ! -f "$EXTRACTED_FILE" ]]; then
    echo "Extracting blacklist file..."
    gunzip -k "$BLACKLIST_FILE"
else
    echo "Blacklist already extracted"
fi

echo "Blacklist setup complete!"