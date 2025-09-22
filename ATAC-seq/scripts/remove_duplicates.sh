#!/bin/bash
# remove_duplicates.sh - Remove PCR duplicates using Picard

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
INPUT_DIR="aligned_reads"  # Directory containing filtered BAM files
OUTPUT_DIR="deduplicated_reads"  # Output directory for deduplicated BAM files

# Software paths
PICARD_JAR="tools/picard.jar"  # Path to Picard JAR file

# Processing parameters
THREADS=4  # Number of threads for samtools operations
JAVA_MEMORY="4g"  # Memory allocation for Java/Picard (adjust based on your system)
INPUT_PATTERN="*_filtered.bam"  # Pattern to match input BAM files

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

echo "Starting duplicate removal with Picard..."
mkdir -p "$OUTPUT_DIR"

for bam_file in "$INPUT_DIR"/$INPUT_PATTERN; do
    sample_name=$(basename "$bam_file" _filtered.bam)
    echo "Processing: $sample_name"
    
    # Define output file paths
    temp_bam="$OUTPUT_DIR/${sample_name}_temp.bam"  # Temporary unsorted BAM
    final_bam="$OUTPUT_DIR/${sample_name}_deduplicated.bam"  # Final sorted BAM
    metrics_file="$OUTPUT_DIR/${sample_name}_dup_metrics.txt"  # Duplication metrics
    
    # Remove duplicates with Picard MarkDuplicates
    java -Xmx${JAVA_MEMORY} -jar "$PICARD_JAR" MarkDuplicates \
        INPUT="$bam_file" \
        OUTPUT="$temp_bam" \
        METRICS_FILE="$metrics_file" \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY=LENIENT \
        CREATE_INDEX=false
    
    if [ $? -eq 0 ]; then
        # Sort and index the deduplicated BAM file
        samtools sort --threads $THREADS -o "$final_bam" "$temp_bam"
        samtools index "$final_bam"
        rm "$temp_bam"  # Remove temporary file to save space
        echo "Successfully processed $sample_name"
    fi
done

echo "Duplicate removal complete!"