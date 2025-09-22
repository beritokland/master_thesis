#!/bin/bash
# bowtie2_align.sh - Bowtie2 alignment for ATAC-seq data

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
# SCRIPT EXECUTION
# =============================================================================

echo "Starting Bowtie2 alignment..."
mkdir -p "$OUTPUT_DIR"

for SAMPLE_DIR in "$INPUT_DIR"/*; do
    if [ -d "$SAMPLE_DIR" ]; then
        SAMPLE_NAME=$(basename "$SAMPLE_DIR")
        echo "Processing: $SAMPLE_NAME"
        
        SAMPLE_OUTPUT="$OUTPUT_DIR/$SAMPLE_NAME"
        mkdir -p "$SAMPLE_OUTPUT"
        
        # Find paired-end files (supports multiple naming conventions)
        R1_FILE=$(ls "$SAMPLE_DIR"/*_R1_paired.fastq.gz "$SAMPLE_DIR"/*_1.fq.gz "$SAMPLE_DIR"/*_R1_001.fastq.gz 2>/dev/null | head -n 1)
        R2_FILE=$(ls "$SAMPLE_DIR"/*_R2_paired.fastq.gz "$SAMPLE_DIR"/*_2.fq.gz "$SAMPLE_DIR"/*_R2_001.fastq.gz 2>/dev/null | head -n 1)
        
        # Extract sample prefix based on file naming convention
        if [[ "$R1_FILE" == *_R1_paired.fastq.gz ]]; then
            PREFIX=$(basename "$R1_FILE" | sed 's/_R1_paired\.fastq\.gz//')
        elif [[ "$R1_FILE" == *_R1_001.fastq.gz ]]; then
            PREFIX=$(basename "$R1_FILE" | sed 's/_R1_001\.fastq\.gz//')
        else
            PREFIX=$(basename "$R1_FILE" | sed 's/_1\.fq\.gz//')
        fi
        
        # Define output file paths
        SAM_FILE="$SAMPLE_OUTPUT/${PREFIX}_aligned.sam"  # Raw alignment output
        BAM_SORTED="$SAMPLE_OUTPUT/${PREFIX}_sorted.bam"  # Sorted BAM file
        BAM_FILTERED="$SAMPLE_OUTPUT/${PREFIX}_filtered.bam"  # Quality filtered BAM file
        
        # Bowtie2 alignment with ATAC-seq optimized parameters
        bowtie2 --threads $THREADS --very-sensitive --no-discordant -X $MAX_FRAGMENT_SIZE \
            -x "$BOWTIE2_INDEX" \
            -1 "$R1_FILE" -2 "$R2_FILE" \
            -S "$SAM_FILE"
        
        # Convert SAM to BAM and sort by genomic coordinates
        samtools view -bS -@ $THREADS "$SAM_FILE" | \
        samtools sort -@ $THREADS -o "$BAM_SORTED"
        
        # Filter for high-quality, properly paired reads
        samtools view -bh -q $MAPQ_THRESHOLD -f 3 -@ $THREADS "$BAM_SORTED" > "$BAM_FILTERED"
        
        # Create index for fast random access
        samtools index -@ $THREADS "$BAM_FILTERED"
        
        # Remove intermediate files to save space
        rm "$SAM_FILE" "$BAM_SORTED"
    fi
done

echo "Alignment complete!"