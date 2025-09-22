#!/bin/bash
# filter_chromosomes.sh - Filter BAM files for standard chromosomes

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
INPUT_DIR="deduplicated_reads"  # Directory containing deduplicated BAM files
OUTPUT_DIR="filtered_chromosomes"  # Output directory for chromosome-filtered BAM files

# Processing parameters
THREADS=4  # Number of threads for samtools operations
INPUT_PATTERN="*_deduplicated.bam"  # Pattern to match input BAM files

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

echo "Filtering BAM files for standard chromosomes..."
mkdir -p "$OUTPUT_DIR"

for bam_file in "$INPUT_DIR"/$INPUT_PATTERN; do
    sample_name=$(basename "$bam_file" .bam)
    output_bam="$OUTPUT_DIR/${sample_name}_filtered.bam"
    
    echo "Processing: $sample_name"
    
    # Filter out non-standard chromosomes and sort
    samtools view -@ $THREADS -h "$bam_file" | \
        grep -v chrM | \
        grep -v chrUn | \
        grep -v random | \
        grep -v chrEBV | \
        samtools sort -@ $THREADS -o "$output_bam" -
    
    if [ $? -eq 0 ]; then
        # Index the filtered BAM file
        samtools index -@ $THREADS "$output_bam"
        echo "Successfully processed $sample_name"
    fi
done

echo "Chromosome filtering complete!"