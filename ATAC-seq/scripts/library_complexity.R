#!/usr/bin/env Rscript
# Generalized script for ATAC-seq library complexity analysis

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
BAM_FOLDER <- "aligned_reads"  # Directory containing BAM files
OUTPUT_DIR <- "library_complexity_results"  # Output directory for PDF reports

# Processing parameters
NUM_CORES <- 4  # Number of cores to use (adjust based on your system)
BAM_PATTERN <- ".*\\.filtered\\.bam$"  # Pattern to match BAM files

# PDF output settings
PDF_WIDTH <- 11
PDF_HEIGHT <- 8.5

# =============================================================================
# PACKAGE INSTALLATION AND LOADING
# =============================================================================

cat("Setting up R environment...\n")

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      if (pkg %in% c("ATACseqQC", "GenomicAlignments")) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org/")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org/")
      }
    }
  }
}

# Install and load required packages
required_packages <- c("ATACseqQC", "parallel")
install_if_missing(required_packages)

# Load libraries
library(ATACseqQC)
library(parallel)

# =============================================================================
# MAIN PROCESSING FUNCTIONS
# =============================================================================

# Function to process a single BAM file
process_bam_file <- function(bam_file, output_dir) {
  # Extract sample name from file path
  sample_name <- gsub("\\.filtered\\.bam$|\\.bam$", "", basename(bam_file))
  
  # Create PDF output path
  pdf_file <- file.path(output_dir, paste0(sample_name, "_library_complexity.pdf"))
  
  cat("Processing:", sample_name, "\n")
  
  # Error handling wrapper
  tryCatch({
    # Validate input file
    if (!file.exists(bam_file)) {
      stop("File does not exist: ", bam_file)
    }
    
    # Check if BAM index exists
    bai_file <- paste0(bam_file, ".bai")
    if (!file.exists(bai_file)) {
      cat("  Warning: BAM index not found, creating...\n")
      # Note: This requires samtools to be installed and in PATH
      system(paste("samtools index", bam_file))
    }
    
    # Read duplicate frequencies
    cat("  Reading duplicate frequencies...\n")
    dup_freq <- readsDupFreq(bam_file)
    
    # Generate library complexity plot
    cat("  Estimating library complexity...\n")
    pdf(pdf_file, width = PDF_WIDTH, height = PDF_HEIGHT)
    estimateLibComplexity(dup_freq)
    title(main = paste("Library Complexity Analysis:", sample_name))
    dev.off()
    
    cat("  ✓ Successfully processed", sample_name, "\n")
    return(list(success = TRUE, sample = sample_name))
    
  }, error = function(e) {
    cat("  ✗ Error processing", sample_name, ":", e$message, "\n")
    if (dev.cur() > 1) dev.off()  # Close any open PDF device
    return(list(success = FALSE, sample = sample_name, error = e$message))
  })
}

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

cat("Starting library complexity analysis...\n")
cat("BAM folder:", BAM_FOLDER, "\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("Using", NUM_CORES, "cores for processing\n")

# Validate input directory
if (!dir.exists(BAM_FOLDER)) {
  stop("Input directory does not exist: ", BAM_FOLDER)
}

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Find BAM files
bam_files <- list.files(
  path = BAM_FOLDER, 
  pattern = BAM_PATTERN, 
  full.names = TRUE, 
  recursive = TRUE
)

if (length(bam_files) == 0) {
  stop("No BAM files found matching pattern: ", BAM_PATTERN)
}

cat("Found", length(bam_files), "BAM files\n")

# Process files
start_time <- Sys.time()

if (NUM_CORES > 1 && length(bam_files) > 1) {
  cat("Processing files in parallel...\n")
  # Parallel processing
  results <- mclapply(
    bam_files, 
    function(x) process_bam_file(x, OUTPUT_DIR), 
    mc.cores = min(NUM_CORES, length(bam_files))
  )
} else {
  cat("Processing files sequentially...\n")
  # Sequential processing
  results <- lapply(bam_files, function(x) process_bam_file(x, OUTPUT_DIR))
}

# Summarize results
end_time <- Sys.time()
successful <- sum(sapply(results, function(x) x$success))
failed <- length(results) - successful

cat("\n" + "="*50 + "\n")
cat("SUMMARY\n")
cat("="*50 + "\n")
cat("Processing completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes\n")
cat("Successfully processed:", successful, "files\n")
cat("Failed:", failed, "files\n")
cat("Results saved in:", OUTPUT_DIR, "\n")

# Report failed files if any
if (failed > 0) {
  cat("\nFailed files:\n")
  failed_files <- results[!sapply(results, function(x) x$success)]
  for (fail in failed_files) {
    cat("  -", fail$sample, "(", fail$error, ")\n")
  }
}

cat("\nLibrary complexity analysis complete!\n")