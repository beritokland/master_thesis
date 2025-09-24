# Differential Accessibility Analysis 

#### Requirements

- R (version 4.0 or higher recommended)
- Required R packages (install via `install.packages()` or `BiocManager::install()`)

#### Required Libraries

```r
library(GenomicRanges)
library(csaw)
library(limma)
library(edgeR)
library(ggplot2)
library(stringr)
```

#### Customization Instructions

**To customize for your experiment:**
1. Modify `cell_lines` to include your cell lines
2. Change `treatment_labels` and `treatment_codes` for your conditions  
3. Adjust `samples_per_treatment` if you have different sample sizes
4. Update file paths (`peak_dir`, `bam_dir`) as needed

#### User Parameters

```r
# Experimental design parameters
cell_lines <- c("d492", "hmle")                    # Cell line names (lowercase for file paths)
treatment_labels <- c("epithelial", "mesenchymal") # Full treatment names for analysis
treatment_codes <- c("E", "M")                     # Short codes used in BAM file names
samples_per_treatment <- 3                         # Number of samples per treatment group

# File paths
peak_dir <- "merged_peaks"
bam_dir <- "/persistent01/shifted_bams"

# File naming patterns
peak_file <- "All_Samples.fwp.filter.non_overlapping.bed"
```

#### Read merged peak files and add column names

```r
# Read peak files for each cell line
peak_data <- list()
for(cell_line in cell_lines) {
  # sep="\t" means the file is tab separated; select the first 3 columns
  peak_data[[cell_line]] <- read.table(file.path(peak_dir, cell_line, peak_file), sep="\t")[,1:3]
  colnames(peak_data[[cell_line]]) <- c("chrom", "start", "end")
}
```

#### Convert to GRanges objects

```r
# Peak files are text files --> need to be converted to be used with genomic packages
granges_list <- list()
for(cell_line in cell_lines) {
  granges_list[[cell_line]] <- GRanges(peak_data[[cell_line]])
}
```

#### Specify paired-end BAMs

```r
# Get BAM files sorted with treatment groups in specified order
# BAM files need to be coordinate sorted with index files (.bai) in same directory
bam_files <- list()
for(cell_line in cell_lines) {
  cell_line_upper <- toupper(cell_line)
  pattern <- paste0("^", cell_line_upper, "_.*_shifted_sorted\\.bam$")
  bams <- list.files(bam_dir, pattern = pattern, full.names = TRUE)
  # Sort by treatment groups (second treatment code comes last)
  bam_files[[cell_line]] <- bams[order(grepl(paste0("_", treatment_codes[2], "_"), bams), bams)]
}
```

#### Define read parameters for read counting in the consensus peaks

```r
standard.chr <- paste0("chr", c(1:22, "X", "Y")) # only use standard chromosomes 
param <- readParam(max.frag=1000, pe="both", restrict=standard.chr) 
# max.frag: maximum fragment length
# pe="both": include paired-end reads
```

#### Count the reads in windows specified by our consensus peak set

```r
# Using the regionCounts function from the csaw package
count_objects <- list()
for(cell_line in cell_lines) {
  count_objects[[cell_line]] <- regionCounts(bam_files[[cell_line]], 
                                            granges_list[[cell_line]], 
                                            param=param)
}
```

#### Filter out low abundance peaks

```r
# Remove peaks with low counts (logCPM < -3) to improve statistical power
for(cell_line in cell_lines) {
  abundances <- aveLogCPM(asDGEList(count_objects[[cell_line]]))
  count_objects[[cell_line]] <- count_objects[[cell_line]][abundances > -3, ]
}
```

#### Normalization method: csaw loess normalization

```r
# Apply loess normalization to account for library size and bias differences
for(cell_line in cell_lines) {
  count_objects[[cell_line]] <- normOffsets(count_objects[[cell_line]], se.out=TRUE)
}
```

#### Setup design matrices for differential analysis

```r
# Create DGEList objects with sample metadata for each cell line
dgelist_objects <- list()
design_matrices <- list()

for(cell_line in cell_lines) {
  # Convert to DGEList object for edgeR analysis
  dgelist_objects[[cell_line]] <- asDGEList(count_objects[[cell_line]])
  
  # Set sample names — must match BAM column order
  sample_names <- paste0(rep(treatment_labels, each=samples_per_treatment), "_", 1:samples_per_treatment)
  colnames(dgelist_objects[[cell_line]]$counts) <- sample_names
  rownames(dgelist_objects[[cell_line]]$samples) <- sample_names
  
  # Assign group labels in the same order
  dgelist_objects[[cell_line]]$samples$group <- rep(treatment_labels, each=samples_per_treatment)
  
  # Build design matrix with no intercept
  design_matrices[[cell_line]] <- model.matrix(~0 + group, data = dgelist_objects[[cell_line]]$samples)
  colnames(design_matrices[[cell_line]]) <- treatment_labels
}
```

#### PCA Analysis

```r
# Create log2CPM matrices for PCA (applies normalization offsets)
pca_data <- list()
for(cell_line in cell_lines) {
  pca_data[[cell_line]] <- cpm(dgelist_objects[[cell_line]], log=TRUE, prior.count=2)
}

# Perform PCA for each cell line
pca_results <- list()
for(cell_line in cell_lines) {
  pca_results[[cell_line]] <- prcomp(t(pca_data[[cell_line]]))  # transpose so samples are rows
}
```

```r
# Create PCA plots for each cell line
for(cell_line in cell_lines) {
  # Calculate variance explained
  variance_explained <- summary(pca_results[[cell_line]])$importance[2,] * 100
  
  # Prepare data frame for plotting
  pca_df <- data.frame(
    PC1 = pca_results[[cell_line]]$x[,1],
    PC2 = pca_results[[cell_line]]$x[,2],
    sample = colnames(pca_data[[cell_line]]),
    group = dgelist_objects[[cell_line]]$samples$group
  )
  
  # Create plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size=4) +
    labs(title = paste("PCA of", toupper(cell_line), "samples"), 
         x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
         color = "Group")  
  
  # Display plot
  print(p)
  
  # Save plot
  ggsave(paste0("DAA_plots/", cell_line, "_PCA.png"), plot = p, width = 7, height = 6, dpi = 300)
}
```

