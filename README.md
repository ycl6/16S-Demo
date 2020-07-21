# 16S rDNA amplicon sequencing analysis using R

### Demonstration of 16S rRNA-based microbiome analysis using dada2, phyloseq, LEfSe, picrust2 and other tools

-----

**16S-rDNA-V3-V4 repository:** https://github.com/ycl6/16S-rDNA-V3-V4

**Demo Dataset:** [PRJEB27564](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB27564) from *[Gut Microbiota in Parkinson's Disease: Temporal Stability and Relations to Disease Progression.](https://pubmed.ncbi.nlm.
nih.gov/31221587/) EBioMedicine. 2019;44:691-707*

**License:** GPL-3.0

## Introduction

In this tutorial, I will use the sequencing data from [PRJEB27564](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB27564) to demonstrate how to use `dada2`, `phyloseq`, `LEfSe`, `picrust2` and other tools to       process and analyse 16S rDNA amplicon sequencing data. We will perform analysis on fecal microbiome data obtained from 32 Parkinson's patients and 32 control subjects.

The universal bacterial primers used in the study are:

* 341F (5'-CCTACGGGNGGCWGCAG-3')
* 785R (5'-GACTACHVGGGTATCTAATCC-3')

## Prerequisites

Please refer to my GitHub repository [here](https://github.com/ycl6/16S-rDNA-V3-V4) to install the required software, R packages and databases.

## Example folder structure

```
/ngs/16S-rDNA-V3-V4-examples/run_trimming.pl                     # Perl script to run cutadapt
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564                          # project folder
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/raw/*_1.fastq.gz         # raw fastq files, read 1
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/raw/*_2.fastq.gz         # raw fastq files, read 2
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/trimmed/*_1.fastq.gz     # cutadapt trimmed files, read 1
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/trimmed/*_2.fastq.gz     # cutadapt trimmed files, read 2
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/filt/*_1.fastq.gz        # dada2 trimmed files, read 1
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/filt/*_1.fastq.gz        # dada2 trimmed files, read 2
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/images                   # output images/PDFs
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/outfiles                 # output files
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/lefse                    # LEfSe & GraPhlAn output files
/ngs/16S-rDNA-V3-V4-examples/PRJEB27564/picrust2_out_stratified  # picrust2 output files
```

## Download demo data

Download the demo data (252 compressed fastq files, total 3.5GB) from the European Nucleotide Archive (ENA) FTP server using the included Shell script `download.sh` in the `raw` folder

```
cd /ngs/16S-Demo/PRJEB27564/raw

./download.sh
```

# Workflow

### 1. dada2 tutorial

### 2. phyloseq tutorial

### 3. lefse tutorial

### 4. picrust2 ALDEx2 tutorial

