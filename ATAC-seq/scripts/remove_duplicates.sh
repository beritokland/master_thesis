#!/bin/bash
# Generalized script for removing PCR duplicates using Picard

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

# File patterns
INPUT_PATTERN="*_filtered.bam"  # Pattern to match input BAM files
VALIDATION_STRINGENCY="LENIENT"  # Picard validation stringency (STRICT, LENIENT, SILENT)

# =============================================================================
# SCRIPT EXECUTION - Generally no need to modify below this line
# =============================================================================

echo "Starting duplicate removal with Picard..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Using $THREADS threads"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Validate Picard installation
if [[ ! -f "$PICARD_JAR" ]]; then
    echo "Error: Picard JAR not found at $PICARD_JAR"
    echo "Please run setup_picard.sh first or update the PICARD_JAR variable."
    exit 1
fi

# Check Java installation
if ! command -v java &> /dev/null; then
    echo "Error: Java not found. Please install Java to run Picard."
    exit 1
fi

# Find input BAM files
input_files=($(find "$INPUT_DIR" -name "$INPUT_PATTERN" -type f))

if [ ${#input_files[@]} -eq 0 ]; then
    echo "Error: No BAM files found matching pattern '$INPUT_PATTERN' in $INPUT_DIR"
    exit 1
fi

echo "Found ${#input_files[@]} BAM files to process"

# Process each BAM file
successful=0
failed=0

for bam_file in "${input_files[@]}"; do
    # Extract sample information
    relative_path=$(realpath --relative-to="$INPUT_DIR" "$bam_file")
    sample_dir=$(dirname "$relative_path")
    filename=$(basename "$bam_file")
    
    # Create output directory structure
    output_sample_dir="$OUTPUT_DIR/$sample_dir"
    mkdir -p "$output_sample_dir"
    
    # Extract sample name (remove file extensions and suffixes)
    sample_name=$(echo "$filename" | sed 's/_aligned.*\.bam$//' | sed 's/_filtered\.bam$//' | sed 's/\.bam$//')
    
    echo "Processing sample: $sample_name"
    echo "  Input: $bam_file"
    
    # Define output files
    temp_bam="$output_sample_dir/${sample_name}_temp.bam"
    final_bam="$output_sample_dir/${sample_name}_deduplicated.bam"
    metrics_file="$output_sample_dir/${sample_name}_dup_metrics.txt"
    
    echo "  Output: $final_bam"
    
    # Step 1: Remove duplicates with Picard
    echo "    Removing duplicates with Picard..."
    java -Xmx${JAVA_MEMORY} -jar "$PICARD_JAR" MarkDuplicates \
        INPUT="$bam_file" \
        OUTPUT="$temp_bam" \
        METRICS_FILE="$metrics_file" \
        REMOVE_DUPLICATES=true \
        VALIDATION_STRINGENCY="$VALIDATION_STRINGENCY" \
        CREATE_INDEX=false
    
    if [ $? -eq 0 ]; then
        # Step 2: Sort the deduplicated BAM file
        echo "    Sorting BAM file..."
        samtools sort --threads $THREADS -o "$final_bam" "$temp_bam"
        
        if [ $? -eq 0 ]; then
            # Step 3: Index the final BAM file
            echo "    Indexing BAM file..."
            samtools index "$final_bam"
            
            if [ $? -eq 0 ]; then
                # Clean up temporary file
                rm "$temp_bam"
                echo "  ✓ Successfully processed $sample_name"
                ((successful++))
            else
                echo "  ✗ Error indexing BAM file for $sample_name"
                ((failed++))
            fi
        else
            echo "  ✗ Error sorting BAM file for $sample_name"
            rm -f "$temp_bam"  # Clean up on failure
            ((failed++))
        fi
    else
        echo "  ✗ Error removing duplicates for $sample_name"
        rm -f "$temp_bam"  # Clean up on failure
        ((failed++))
    fi
    
    echo "  ----------------------------------------"
done

# Summary
echo ""
echo "==============================================="
echo "DUPLICATE REMOVAL SUMMARY"
echo "==============================================="
echo "Successfully processed: $successful files"
echo "Failed: $failed files"
echo "Results saved in: $OUTPUT_DIR/"
echo ""

if [ $successful -gt 0 ]; then
    echo "Deduplication metrics saved as *_dup_metrics.txt files"
    echo "You can review these files to assess duplication rates"
fi

echo "Duplicate removal complete!"