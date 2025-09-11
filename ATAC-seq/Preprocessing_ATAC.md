## Preprocessing of ATAC-seq data

#### Setup 

1. Download the `scripts/` directory to your working directory

2. Ensure your raw data is in `01.RawData/` with the expected structure

3. Make all scripts executable: 

```{bash}
chmod +x scripts/*.sh
```

#### Quality Control

Run the FastQC analysis script [fastqc.sh](scripts/fastqc.sh):

```bash
./scripts/fastqc.sh
```

Use MultiQC to summarize the results for all files in one html report:

```{bash}
cd ./fastqc_initial/
multiqc .
```

#### REMOVAL OF ADAPTER SEQUENCES USING Trimmomatic

Run [trimmomatic.sh](trimmomatic.sh)

