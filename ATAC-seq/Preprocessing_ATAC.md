## Preprocessing of ATAC-seq Data

#### Setup 

1. Download the `scripts/` directory to your working directory

2. Ensure your raw data is in `01.RawData/` with the expected structure

3. Make all scripts executable: 

```bash
chmod +x scripts/*.sh scripts/*.R
```

#### Quality Control

To evaluate the quality of the raw sequencing reads, run [fastqc_initial.sh](scripts/fastqc_initial.sh):

```bash
./scripts/fastqc_initial.sh
```

Use `MultiQC` to summarize the results for all files in one html report:

```bash
cd ./fastqc_initial/
multiqc .
```

#### Removal of Adapter Sequences

ATAC-seq fragments can be shorter than the read length, which can cause sequencing read-through into adapters. To remove this, run [trimmomatic.sh](scripts/trimmomatic.sh):

```bash
./scripts/trimmomatic.sh
```

#### Alignment of Sequencing Reads to Reference Genome (GRCh38)

Run [setup_reference.sh](scripts/setup_reference.sh) to download and index the reference genome (run once):

```bash
./scripts/setup_reference.sh
```

To align the trimmed sequencing reads, run [bowtie2_align.sh](scripts/bowtie2_align.sh):

```bash
./scripts/bowtie2_align.sh
```

#### Evalutation of the Library Complexity 

To assess library complexity to determine if sequencing depth was sufficient, run [library_complexity.R](scripts/library_complexity.R):

```bash
Rscript scripts/library_complexity.R
```

#### Removal of Duplicate Reads

Run [setup_picard.sh](scripts/setup_picard.sh) to download `Picard` (run once):

```bash
./scripts/setup_picard.sh
```

To remove technical duplicates from PCR, run [remove_duplicates.sh](scripts/remove_duplicates.sh):

```bash
./scripts/remove_duplicates.sh
```

Validate the removal of duplicate reads using running [fastqc_deduplication.sh](scripts/fastqc_deduplication.sh):

```bash
./scripts/fastqc_deduplication.sh
```

#### Filtering for Standard Chromosomes

To remove mitochondrial reads (chrM), unplaced contigs (chrUn), viral contamination (chrEBV) and random chromosomes, keeping only standard autosomes and sex chromosomes, run [filter_chromosomes.sh](scripts/filter_chromosomes.sh):

```bash
./scripts/filter_chromosomes.sh
```

#### Removal of Blacklisted Regions 

Run [setup_blacklist.sh](scripts/setup_blacklist.sh) to download ENCODE blacklist for hg38 (run once):

```bash
./scripts/setup_blacklist.sh
```

To remove ENCODE blacklist regions, run [remove_blacklist.sh](scripts/remove_blacklist.sh):

```bash
./scripts/remove_blacklist.sh
```

### Fragment Size Distribution (FSD) and Transcriptional Start Site Enrichment (TSSE) 

To generate plots of FSD and TSSE, run [FSD_TSSE.R](scripts/FSD_TSSE.R): 

```bash
./scripts/FSD_TSSE.R
```








