#!/bin/bash
# remove_blacklist.sh - Remove ENCODE blacklisted regions from BAM files

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
INPUT_DIR="filtered_chromosomes"  # Directory containing chromosome-filtered BAM files
OUTPUT_DIR="blacklist_removed"  # Output directory for blacklist-removed BAM files

# Processing parameters
THREADS=8  # Number of threads for samtools operations
BLACKLIST_FILE="reference_genome/hg38-blacklist.v2.bed"  # Path to blacklist BED file
INPUT_PATTERN="*_filtered.bam"  # Pattern to match input BAM files

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

echo "Removing blacklisted regions from BAM files..."
mkdir -p "$OUTPUT_DIR"

for bam_file in "$INPUT_DIR"/$INPUT_PATTERN; do
    sample_name=$(basename "$bam_file" _filtered.bam)
    echo "Processing: $sample_name"
    
    output_bam="$OUTPUT_DIR/${sample_name}_clean.bam"
    
    # Remove blacklisted regions using bedtools intersect
    bedtools intersect -abam "$bam_file" -b "$BLACKLIST_FILE" -v | \
    samtools sort -@ $THREADS -o "$output_bam" -
    
    if [ $? -eq 0 ]; then
        # Index the cleaned BAM file
        samtools index -@ $THREADS "$output_bam"
        echo "Successfully processed $sample_name"
    fi
done

echo "Blacklist removal complete!"