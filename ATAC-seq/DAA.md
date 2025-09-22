# Differential Accessibility Analysis

## Required libraries 

```r
library(GenomicRanges)
library(csaw)
library(limma)
library(edgeR)
library(ggplot2)
library(stringr)
```

## User Parameters

```r
# File paths
peak_dir <- "merged_peaks"
bam_dir <- "/persistent01/shifted_bams"

# File naming patterns
peak_file <- "All_Samples.fwp.filter.non_overlapping.bed"

# Treatment groups (modify for your experimental design)
treatment_groups <- c("E", "M")  # E = Epithelial, M = Mesenchymal (E samples will come first)
```

## Read merged peak files and add column names

```r
# sep="\t" means the file is tab separated; select the first 3 columns
d492 <- read.table(file.path(peak_dir, "d492", peak_file), sep="\t")[,1:3]
hmle <- read.table(file.path(peak_dir, "hmle", peak_file), sep="\t")[,1:3]

colnames(d492) <- c("chrom", "start", "end")
colnames(hmle) <- c("chrom", "start", "end")
```

## Convert to GRanges object

```r
d492.GR <- GRanges(d492)
hmle.GR <- GRanges(hmle)
```

## Specify paired-end BAMs

```r
# Get D492 BAM files sorted with treatment groups in specified order
d492.bams <- list.files(bam_dir, pattern = "^D492_.*_shifted_sorted\\.bam$", full.names = TRUE)
d492.bams <- d492.bams[order(grepl(paste0("_", treatment_groups[2], "_"), d492.bams), d492.bams)]

# Get HMLE BAM files sorted with treatment groups in specified order  
hmle.bams <- list.files(bam_dir, pattern = "^HMLE_.*_shifted_sorted\\.bam$", full.names = TRUE)
hmle.bams <- hmle.bams[order(grepl(paste0("_", treatment_groups[2], "_"), hmle.bams), hmle.bams)]
```

## Define read parameters for read counting in the consensus peaks

```r
standard.chr <- paste0("chr", c(1:22, "X", "Y")) # only use standard chromosomes 
param <- readParam(max.frag=1000, pe="both", restrict=standard.chr) 
```

## Count the reads in windows specified by our consensus peak set

```r
d492.counts <- regionCounts(d492.bams, d492.GR, param=param)
hmle.counts <- regionCounts(hmle.bams, hmle.GR, param=param)
```

## Filter out low abundance peaks

```r
d492.abundances <- aveLogCPM(asDGEList(d492.counts))
d492.counts <- d492.counts[d492.abundances > -3, ] # only use peaks logCPM > -3

hmle.abundances <- aveLogCPM(asDGEList(hmle.counts))
hmle.counts <- hmle.counts[hmle.abundances > -3, ] # only use peaks logCPM > -3
```

## Normalization method: csaw loess normalization

```r
d492.counts <- normOffsets(d492.counts, se.out=TRUE) # type="loess" is now default
hmle.counts <- normOffsets(hmle.counts, se.out=TRUE) # type="loess" is now default
```

## Summary

Objects created:
- `d492.counts`, `hmle.counts` - normalized count objects ready for differential analysis
- `d492.GR`, `hmle.GR` - GRanges objects with peak coordinates  
- `d492.bams`, `hmle.bams` - BAM file paths for each cell line

Your data is now ready for differential accessibility analysis!