## Preprocessing of ATAC-seq data

#### Setup 

1. Download the `scripts/` directory to your working directory

2. Ensure your raw data is in `01.RawData/` with the expected structure

3. Make all scripts executable: 

```bash
chmod +x scripts/*.sh scripts/*.R
```

#### Quality Control

Run the [fastqc_initial.sh](scripts/fastqc_initial.sh) analysis script:

```bash
./scripts/fastqc_initial.sh
```

Use `MultiQC` to summarize the results for all files in one html report:

```bash
cd ./fastqc_initial/
multiqc .
```

#### Removal of adapter sequences using Trimmomatic

Run the [trimmomatic.sh](scripts/trimmomatic.sh) script:

```bash
./scripts/trimmomatic.sh
```

#### Alignment of sequencing reads to GRCh38 

Run the [setup_reference.sh](scripts/setup_reference.sh) script for downloading and indexing the reference genome (run once):

```bash
./scripts/setup_reference.sh
```

Run the [bowtie2_align.sh](scripts/bowtie2_align.sh) script:

```bash
./scripts/bowtie2_align.sh
```

#### Evalutation of the library complexity 

Run the [library_complexity.R](scripts/library_complexity.R) script:

```bash
Rscript scripts/library_complexity.R
```

#### Removal of duplicate reads

Run the [setup_picard.sh](scripts/setup_picard.sh) script  for downloading `Picard` (run once):

```bash
./scripts/setup_picard.sh
```

Run the [remove_duplicates.sh](scripts/remove_duplicates.sh) script:

```bash
./scripts/remove_duplicates.sh
```

Validate the removal of duplicate reads using running the [fastqc_deduplication.sh](scripts/fastqc_deduplication.sh) script:

```bash
./scripts/fastqc_deduplication.sh
```





