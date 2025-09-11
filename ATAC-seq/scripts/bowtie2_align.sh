#!/bin/bash
# Generalized Bowtie2 alignment for ATAC-seq data

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
INPUT_DIR="trimmed_reads"  # Directory containing trimmed reads
OUTPUT_DIR="aligned_reads"  # Output directory for alignment results

# Reference genome index
BOWTIE2_INDEX="reference_genome/GRCh38_index"  # Path to Bowtie2 index

# Processing parameters
THREADS=4  # Number of threads (adjust based on your CPU)
MAX_FRAGMENT_SIZE=2000  # Maximum fragment length for ATAC-seq
MAPQ_THRESHOLD=30  # Minimum mapping quality threshold

# =============================================================================
# SCRIPT EXECUTION - Generally no need to modify below this line
# =============================================================================

echo "Starting Bowtie2 alignment..."
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Using $THREADS threads"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Check if Bowtie2 index exists
if [[ ! -f "${BOWTIE2_INDEX}.1.bt2" ]]; then
    echo "Error: Bowtie2 index not found at $BOWTIE2_INDEX"
    echo "Please run setup_reference.sh first or update the BOWTIE2_INDEX variable."
    exit 1
fi

# Process each sample directory
for SAMPLE_DIR in "$INPUT_DIR"/*; do
    if [ -d "$SAMPLE_DIR" ]; then
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        echo "Processing sample: $SAMPLE_NAME"
        
        # Create sample output directory
        SAMPLE_OUTPUT="$OUTPUT_DIR/$SAMPLE_NAME"
        mkdir -p "$SAMPLE_OUTPUT"
        
        # Find paired-end files (supports multiple naming conventions)
        R1_FILE=$(ls "$SAMPLE_DIR"/*_R1_paired.fastq.gz "$SAMPLE_DIR"/*_1.fq.gz "$SAMPLE_DIR"/*_R1_001.fastq.gz 2>/dev/null | head -n 1)
        R2_FILE=$(ls "$SAMPLE_DIR"/*_R2_paired.fastq.gz "$SAMPLE_DIR"/*_2.fq.gz "$SAMPLE_DIR"/*_R2_001.fastq.gz 2>/dev/null | head -n 1)
        
        if [[ -z "$R1_FILE" || -z "$R2_FILE" ]]; then
            echo "  Warning: Missing R1 or R2 file in $SAMPLE_DIR. Skipping..."
            continue
        fi
        
        # Extract sample prefix
        if [[ "$R1_FILE" == *_R1_paired.fastq.gz ]]; then
            PREFIX=$(basename "$R1_FILE" | sed 's/_R1_paired\.fastq\.gz//')
        elif [[ "$R1_FILE" == *_R1_001.fastq.gz ]]; then
            PREFIX=$(basename "$R1_FILE" | sed 's/_R1_001\.fastq\.gz//')
        else
            PREFIX=$(basename "$R1_FILE" | sed 's/_1\.fq\.gz//')
        fi
        
        echo "  Aligning files: $PREFIX"
        
        # Define output files
        SAM_FILE="$SAMPLE_OUTPUT/${PREFIX}_aligned.sam"
        BAM_SORTED="$SAMPLE_OUTPUT/${PREFIX}_sorted.bam"
        BAM_FILTERED="$SAMPLE_OUTPUT/${PREFIX}_filtered.bam"
        
        # Step 1: Bowtie2 alignment
        echo "    Running Bowtie2 alignment..."
        bowtie2 --threads $THREADS --very-sensitive --no-discordant -X $MAX_FRAGMENT_SIZE \
            -x "$BOWTIE2_INDEX" \
            -1 "$R1_FILE" -2 "$R2_FILE" \
            -S "$SAM_FILE"
        
        if [ $? -ne 0 ]; then
            echo "    ✗ Error during alignment for $SAMPLE_NAME"
            continue
        fi
        
        # Step 2: Convert SAM to BAM and sort
        echo "    Converting SAM to BAM and sorting..."
        samtools view -bS -@ $THREADS "$SAM_FILE" | \
        samtools sort -@ $THREADS -o "$BAM_SORTED"
        
        # Step 3: Filter BAM (remove low-quality reads, keep only properly paired)
        echo "    Filtering BAM file..."
        samtools view -bh -q $MAPQ_THRESHOLD -f 3 -@ $THREADS "$BAM_SORTED" > "$BAM_FILTERED"
        
        # Step 4: Index the filtered BAM file
        echo "    Indexing BAM file..."
        samtools index -@ $THREADS "$BAM_FILTERED"
        
        # Clean up intermediate files
        rm "$SAM_FILE" "$BAM_SORTED"
        
        echo "  ✓ Successfully processed $SAMPLE_NAME"
    fi
done

echo "All samples aligned, filtered, and indexed!"
echo "Results saved in: $OUTPUT_DIR/"