#!/bin/bash
# Run FastQC on deduplicated BAM files

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
INPUT_DIR="deduplicated_reads"
OUTPUT_DIR="fastqc_deduplicated"

# Processing parameters
THREADS=4
BAM_PATTERN="*_deduplicated.bam"

# =============================================================================
# SCRIPT EXECUTION - Generally no need to modify below this line
# =============================================================================

echo "Running FastQC on deduplicated BAM files..."

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Find BAM files
bam_files=($(find "$INPUT_DIR" -name "$BAM_PATTERN" -type f))

if [ ${#bam_files[@]} -eq 0 ]; then
    echo "Error: No BAM files found"
    exit 1
fi

# Run FastQC
fastqc -t $THREADS "${bam_files[@]}" -o "$OUTPUT_DIR"

# Run MultiQC if available
if command -v multiqc &> /dev/null; then
    cd "$OUTPUT_DIR"
    multiqc .
    cd - > /dev/null
fi

echo "Results saved in: $OUTPUT_DIR/"