## Preprocessing of ATAC-seq data

#### EVALUATION OF READ QUALITY USING FastQC

```{bash}
# Create and enter new directory for inital quality control
mkdir fastqc_initial
cd fastqc_initial

# Run FastQC on all fastq files containing the sequencing reads at the same time
find 01.RawData -name "*.fq.gz" | xargs fastqc -t 4 -o fastqc_initial
# -t 4 specifies the number of cores to be used 

# Use MultiQC to summarize the results for all files in one html report
cd ./fastqc_initial/
multiqc .
```

#### REMOVAL OF ADAPTER SEQUENCES USING Trimmomatic

Run [trimmomatic.sh](trimmomatic.sh)
