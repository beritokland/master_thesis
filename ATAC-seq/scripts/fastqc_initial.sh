#!/bin/bash
# fastqc.sh - Quality control analysis of raw sequencing reads

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
RAW_DATA_DIR="01.RawData"  # Directory containing raw FASTQ files
OUTPUT_DIR="fastqc_initial"  # Output directory for FastQC reports

# Processing parameters
THREADS=4  # Number of threads for FastQC (adjust based on your CPU)

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run FastQC on all samples
echo "Running FastQC on all samples..."
for sample_dir in "$RAW_DATA_DIR"/*/; do
    echo "Processing $(basename "$sample_dir")..."
    fastqc -t $THREADS -o "$OUTPUT_DIR" "$sample_dir"*.fq.gz
done

echo "FastQC analysis complete. Results saved in $OUTPUT_DIR/"