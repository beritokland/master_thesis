#!/bin/bash
# trimmomatic.sh - Trimmomatic adapter trimming for ATAC-seq data

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
BASE_DIR="01.RawData"  # Raw data directory
OUTPUT_DIR="trimmed_reads"  # Output directory for trimmed reads

# Processing parameters
ADAPTERS="NexteraPE-PE.fa"  # Local adapter file or full path
THREADS=24  # Number of threads (adjust based on your CPU)

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

echo "Starting Trimmomatic processing..."
mkdir -p "$OUTPUT_DIR"

for SAMPLE_DIR in "$BASE_DIR"/*; do
    if [ -d "$SAMPLE_DIR" ]; then
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        echo "Processing: $SAMPLE_NAME"
        
        # Create sample output directory
        SAMPLE_OUTPUT="$OUTPUT_DIR/$SAMPLE_NAME"
        mkdir -p "$SAMPLE_OUTPUT"
        
        # Find paired-end files
        R1_FILE=$(ls "$SAMPLE_DIR"/*_1.fq.gz 2>/dev/null | head -n 1)
        R2_FILE=$(ls "$SAMPLE_DIR"/*_2.fq.gz 2>/dev/null | head -n 1)
        
        # Extract sample prefix
        PREFIX=$(basename "$R1_FILE" | sed 's/_1\.fq\.gz//')
        
        # Define output files
        R1_PAIRED="$SAMPLE_OUTPUT/${PREFIX}_R1_paired.fastq.gz"
        R1_UNPAIRED="$SAMPLE_OUTPUT/${PREFIX}_R1_unpaired.fastq.gz"
        R2_PAIRED="$SAMPLE_OUTPUT/${PREFIX}_R2_paired.fastq.gz"
        R2_UNPAIRED="$SAMPLE_OUTPUT/${PREFIX}_R2_unpaired.fastq.gz"
        
        # Run Trimmomatic
        trimmomatic PE -threads $THREADS -phred33 \
            "$R1_FILE" "$R2_FILE" \
            "$R1_PAIRED" "$R1_UNPAIRED" "$R2_PAIRED" "$R2_UNPAIRED" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True \
            LEADING:3 TRAILING:3 MINLEN:36
    fi
done

echo "Trimmomatic processing complete!"