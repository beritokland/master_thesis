# Install trimmomatic
sudo apt install trimmomatic

# Manually download the adapter file 
wget https://raw.githubusercontent.com/usadellab/Trimmomatic/master/adapters/NexteraPE-PE.fa -O NexteraPE-PE.fa

# Find path to the adapter fasta file provided with the Trimmomatic program
locate NexteraPE-PE.fa
# /usr/local/share/trimmomatic/adapters/NexteraPE-PE.fa # Path to fasta file

# Thomas performed trimmomatic for me on genlab7, due to storage issues

# Define base directory containing sample folders
BASE_DIR="/open/work/Thomas/EMT/ATAC-seq/X204SC24100951-Z01-F001/01.RawData"
ADAPTERS="/usr/local/share/trimmomatic/adapters/NexteraPE-PE.fa"
THREADS=24

# Loop through each sample directory
for SAMPLE_DIR in "$BASE_DIR"/*; do
    if [ -d "$SAMPLE_DIR" ]; then  # Ensure it's a directory
        echo "Processing directory: $SAMPLE_DIR"

        # Find the correct file prefix (removing _1.fq.gz or _2.fq.gz)
        R1_FILE=$(ls "$SAMPLE_DIR"/*_1.fq.gz 2>/dev/null | head -n 1)
        R2_FILE=$(ls "$SAMPLE_DIR"/*_2.fq.gz 2>/dev/null | head -n 1)

        if [[ -z "$R1_FILE" || -z "$R2_FILE" ]]; then
            echo "Warning: Missing R1 or R2 file in $SAMPLE_DIR. Skipping..."
            continue
        fi

        # Extract the common prefix (remove _1.fq.gz or _2.fq.gz)
        PREFIX=$(basename "$R1_FILE" | sed 's/_1\.fq\.gz//')

        echo "Detected sample: $PREFIX"

        # Define output files
        R1_PAIRED="$SAMPLE_DIR/${PREFIX}_R1_paired.fastq.gz"
        R1_UNPAIRED="$SAMPLE_DIR/${PREFIX}_R1_unpaired.fastq.gz"
        R2_PAIRED="$SAMPLE_DIR/${PREFIX}_R2_paired.fastq.gz"
        R2_UNPAIRED="$SAMPLE_DIR/${PREFIX}_R2_unpaired.fastq.gz"

        # Run Trimmomatic
        trimmomatic PE -threads $THREADS -phred33 \
            "$R1_FILE" "$R2_FILE" \
            "$R1_PAIRED" "$R1_UNPAIRED" "$R2_PAIRED" "$R2_UNPAIRED" \
            ILLUMINACLIP:"$ADAPTERS":2:30:10:2:True \
            LEADING:3 TRAILING:3 MINLEN:36

        echo "Finished processing $PREFIX"
    fi
done

echo "All samples processed!"

# Input: fastq read files (forward and reverse) and fasta file with adapter and primer sequences (overreprested sequences in fastqc file?)
# Output: four fastq files with trimmed reads
# PE: paired-end mode
# threads: number used cores
# phred 33 quality scores 
# trimlog: creates a log of all read trimmings
# ILLUMINACLIP: PAHT/TO/FASTAFILE:
