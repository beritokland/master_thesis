#!/bin/bash

# Generalized Trimmomatic script for ATAC-seq data processing
# This script processes paired-end sequencing data in a directory structure

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Set base directory containing sample folders
BASE_DIR="01.RawData"  # Change this to your raw data directory

# Set adapter file path
ADAPTERS="NexteraPE-PE.fa"  # Local adapter file or full path

# Set number of threads (adjust based on your system)
THREADS=4  # Adjust based on your CPU cores

# Output directory for trimmed reads
OUTPUT_DIR="trimmed_reads"

# =============================================================================
# SCRIPT EXECUTION - Generally no need to modify below this line
# =============================================================================

echo "Starting Trimmomatic processing..."
echo "Base directory: $BASE_DIR"
echo "Using $THREADS threads"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if adapter file exists
if [[ ! -f "$ADAPTERS" ]]; then
    echo "Error: Adapter file '$ADAPTERS' not found!"
    echo "Please ensure the adapter file exists or update the ADAPTERS variable."
    exit 1
fi

# Loop through each sample directory
for SAMPLE_DIR in "$BASE_DIR"/*; do
    if [ -d "$SAMPLE_DIR" ]; then  # Ensure it's a directory
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        echo "Processing sample: $SAMPLE_NAME"
        
        # Create sample-specific output directory
        SAMPLE_OUTPUT="$OUTPUT_DIR/$SAMPLE_NAME"
        mkdir -p "$SAMPLE_OUTPUT"
        
        # Find paired-end files
        R1_FILE=$(ls "$SAMPLE_DIR"/*_1.fq.gz 2>/dev/null | head -n 1)
        R2_FILE=$(ls "$SAMPLE_DIR"/*_2.fq.gz 2>/dev/null | head -n 1)
        
        if [[ -z "$R1_FILE" || -z "$R2_FILE" ]]; then
            echo "Warning: Missing R1 or R2 file in $SAMPLE_DIR. Skipping..."
            continue
        fi
        
        # Extract the common prefix
        PREFIX=$(basename "$R1_FILE" | sed 's/_1\.fq\.gz//')
        echo "  Processing files: $PREFIX"
        
        # Define output files
        R1_PAIRED="$SAMPLE_OUTPUT/${PREFIX}_R1_paired.fastq.gz"
        R1_UNPAIRED="$SAMPLE_OUTPUT/${PREFIX}_R1_unpaired.fastq.gz"
        R2_PAIRED="$SAMPLE_OUTPUT/${PREFIX}_R2_paired.fastq.gz"
        R2_UNPAIRED="$SAMPLE_OUTPUT/${PREFIX}_R2_unpaired.fastq.gz"
        
        # Run Trimmomatic
        echo "  Running Trimmomatic..."
        trimmomatic PE -threads $THREADS -phred33 \
            "$R1_FILE" "$R2_FILE" \
            "$R1_PAIRED" "$R1_UNPAIRED" "$R2_PAIRED" "$R2_UNPAIRED" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True \
            LEADING:3 TRAILING:3 MINLEN:36
        
        if [ $? -eq 0 ]; then
            echo "  ✓ Successfully processed $SAMPLE_NAME"
        else
            echo "  ✗ Error processing $SAMPLE_NAME"
        fi
    fi
done

echo "All samples processed!"
echo "Trimmed reads saved in: $OUTPUT_DIR/"
