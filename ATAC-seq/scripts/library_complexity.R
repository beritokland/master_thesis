#!/usr/bin/env Rscript
# library_complexity.R - Library complexity analysis for ATAC-seq data

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
BAM_FOLDER <- "aligned_reads"  # Directory containing BAM files
OUTPUT_DIR <- "library_complexity_results"  # Output directory for PDF reports

# Processing parameters
NUM_CORES <- 4  # Number of cores (adjust based on your system)
BAM_PATTERN <- ".*\\.filtered\\.bam$"  # Pattern to match BAM files

# PDF output settings
PDF_WIDTH <- 11
PDF_HEIGHT <- 8.5

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

# Install and load required packages
if (!require("ATACseqQC", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("ATACseqQC")
}

library(ATACseqQC)
library(parallel)

cat("Starting library complexity analysis...\n")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE)

# Find BAM files
bam_files <- list.files(BAM_FOLDER, pattern = BAM_PATTERN, full.names = TRUE, recursive = TRUE)

# Process each BAM file
for (bam_file in bam_files) {
  sample_name <- gsub("\\.filtered\\.bam$|\\.bam$", "", basename(bam_file))
  cat("Processing:", sample_name, "\n")
  
  # Create PDF output
  pdf_file <- file.path(OUTPUT_DIR, paste0(sample_name, "_library_complexity.pdf"))
  
  # Read duplicate frequencies and generate plot
  dup_freq <- readsDupFreq(bam_file)
  pdf(pdf_file, width = PDF_WIDTH, height = PDF_HEIGHT)
  estimateLibComplexity(dup_freq)
  title(main = paste("Library Complexity Analysis:", sample_name))
  dev.off()
}

cat("Library complexity analysis complete!\n")