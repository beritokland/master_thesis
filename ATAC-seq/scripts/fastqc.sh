#!/bin/bash

# Set variables for easy customization
RAW_DATA_DIR="01.RawData"
OUTPUT_DIR="fastqc_initial"
THREADS=4

# Create output directory
mkdir -p ${OUTPUT_DIR}

# Run FastQC on all samples
echo "Running FastQC on all samples..."
for sample_dir in ${RAW_DATA_DIR}/*/; do
    echo "Processing $(basename ${sample_dir})..."
    fastqc -t ${THREADS} -o ${OUTPUT_DIR} ${sample_dir}*.fq.gz
done

echo "FastQC analysis complete. Results saved in ${OUTPUT_DIR}/"