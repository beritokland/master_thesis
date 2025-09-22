#!/usr/bin/env Rscript
# atacseq_qc.R - Quality control analysis for ATAC-seq data

# =============================================================================
# CONFIGURATION - Modify these variables for your setup
# =============================================================================

# Input/Output directories
BAM_DIR <- "blacklist_removed"  # Directory containing cleaned BAM files
OUTPUT_DIR <- "atacseq_qc_plots"  # Output directory for QC plots and scores

# Processing parameters
BAM_PATTERN <- "*_clean.bam$"  # Pattern to match input BAM files

# PDF output settings
PDF_WIDTH <- 11
PDF_HEIGHT <- 8.5

# =============================================================================
# SCRIPT EXECUTION
# =============================================================================

# Install and load required packages
required_packages <- c("ATACseqQC", "TxDb.Hsapiens.UCSC.hg38.knownGene", "GenomicAlignments")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install(pkg)
  }
}

library(ATACseqQC)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicAlignments)

cat("Starting ATAC-seq quality control analysis...\n")

# Create output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Find BAM files
bam_files <- list.files(path = BAM_DIR, pattern = BAM_PATTERN, full.names = TRUE)

# Get transcripts for TSS enrichment analysis
txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Initialize dataframe to store TSSE scores
tsse_scores <- data.frame(sample_name = character(), tsse_score = numeric(), stringsAsFactors = FALSE)

# Process each BAM file
for (bam_file in bam_files) {
  sample_name <- gsub("_clean.bam$", "", basename(bam_file))
  cat("Processing:", sample_name, "\n")
  
  # Fragment size distribution
  pdf_path <- file.path(OUTPUT_DIR, paste0(sample_name, ".fragment_sizes.pdf"))
  pdf(pdf_path, width = PDF_WIDTH, height = PDF_HEIGHT)
  fragSize <- fragSizeDist(bam_file, sample_name)
  dev.off()
  
  # TSS enrichment analysis
  tss_pdf_path <- file.path(OUTPUT_DIR, paste0(sample_name, ".TSS_enrichment.pdf"))
  pdf(tss_pdf_path, width = PDF_WIDTH, height = PDF_HEIGHT)
  
  # Read alignments and subset to standard chromosomes
  alignments <- readGAlignments(bam_file)
  standard_chroms <- c(paste0("chr", 1:22), "chrX", "chrY")
  alignments_subset <- alignments[seqnames(alignments) %in% standard_chroms]
  
  # Calculate and plot TSS enrichment
  tsse <- TSSEscore(alignments_subset, txs)
  tsse_scores <- rbind(tsse_scores, data.frame(sample_name = sample_name, tsse_score = tsse$TSSEscore))
  
  plot(100*(-9:10-.5), tsse$values, type = "b", 
       xlab = "Distance to TSS (bp)", ylab = "Aggregate TSS score",
       main = paste0(sample_name, " - TSS Enrichment (Score: ", round(tsse$TSSEscore, 2), ")"))
  dev.off()
}

# Write TSSE scores summary
write.table(tsse_scores, file = file.path(OUTPUT_DIR, "all_tsse_scores.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("ATAC-seq quality control analysis complete!\n")