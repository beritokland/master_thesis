## Preprocessing of ATAC-seq data

#### Setup 

1. Download the `scripts/` directory to your working directory

2. Ensure your raw data is in `01.RawData/` with the expected structure

3. Make all scripts executable: 

```bash
chmod +x scripts/*.sh
```

#### Quality Control

Run the FastQC analysis script [fastqc.sh](scripts/fastqc.sh):

```bash
./scripts/fastqc.sh
```

Use MultiQC to summarize the results for all files in one html report:

```bash
cd ./fastqc_initial/
multiqc .
```

#### Removal of adapter sequences using Trimmomatic

Run the Trimmomatic script [trimmomatic.sh](scripts/trimmomatic.sh):

```bash
./scripts/trimmomatic.sh
```

#### Alignment of sequencing reads to GRCh38 

Setup script [setup_reference.sh](scripts/setup_reference.sh) for downloading and indexing the reference genome (run once):

```bash
./scripts/setup_reference.sh
```

Run the Bowtie2 alignment script [bowtie2_align.sh](scripts/bowtie2_align.sh):

```bash
./scripts/bowtie2_align.sh
```








