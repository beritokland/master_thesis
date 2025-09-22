#!/bin/bash
# fastqc_deduplication.sh - Run FastQC on deduplicated BAM files

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
INPUT_DIR="deduplicated_reads"  # Directory containing deduplicated BAM files
OUTPUT_DIR="fastqc_deduplicated"  # Output directory for FastQC results

# Processing parameters
THREADS=4  # Number of threads for FastQC
BAM_PATTERN="*_deduplicated.bam"  # Pattern to match deduplicated BAM files

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

echo "Running FastQC on deduplicated BAM files..."

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Find deduplicated BAM files
bam_files=($(find "$INPUT_DIR" -name "$BAM_PATTERN" -type f))

# Run FastQC analysis on all BAM files
fastqc -t $THREADS "${bam_files[@]}" -o "$OUTPUT_DIR"

# Run MultiQC to aggregate results if available
if command -v multiqc &> /dev/null; then
    cd "$OUTPUT_DIR"
    multiqc .
    cd - > /dev/null
fi

echo "Results saved in: $OUTPUT_DIR/"