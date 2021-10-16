---
title: Playing with bioinformatics I
author: Lucas M. Taniguti
comments: true
date: 2021-10-13 11:33:00 +0800
categories: [Bioinformatics, File formats]
tags: [BLAST]
---

## Introduction

Bioinformatics is known to have plenty of file formats. Every time a discussion about file formats starts I remember [this comic](http://www.niso.org/sites/default/files/inline-images/How%20Standards%20Proliferate%20%281%29.png).

In this post I'll guide you throught some of these formats, always paying atention to our precious illumina reads contents. Here is our plan:

1. Download public data to play with. We'll use sars-cov-2 genome and Illumina reads from an whole genome amplicon experiment.
2. Convert FASTQ files to unaligned BAM (uBAM), just to show how it works.
3. Transform uBAM into a single interleave FASTQ file.
4. Align FASTQs files (interleaved and paired) on sars-cov-2 genome, creating a BAM file.
5. Convert BAM to CRAM and also to a lossy CRAM
6. Convert back from CRAM to pairs of FASTQs

## Tools

I'm assuming you are in a unix like system (ubuntu, fedora, WSL2, macOS, etc)

- [Picard](https://broadinstitute.github.io/picard/) - To convert FASTQ <-> uBAM and SAM <-> BAM <-> CRAM
- [Samtools](http://www.htslib.org/) - To create .fai index
- [BWA](https://github.com/lh3/bwa) - To align reads to genome
- [Crumble](https://www.sanger.ac.uk/tool/crumble/) - To try lossy compression on alignment file

> See [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Docker](https://www.docker.com/) to easily get your environment ready. I plan to make a post about it soon.

## Download public dataset

Firstly, go to a directory you like. For example `/tmp/bioinfo` if you don't want to persist files from this post after you reboot your system.

We'll use two repositories to retrieve some data to our experiment:
- [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra): stores raw data from NGS sequencing experiments.
- FTP from [NCBI/genomes](https://www.ncbi.nlm.nih.gov/genome): organizes information on genomes.

### Raw sequencing data

SRA organize it's data using unique identifiers, like this one under our `GENOME_SRA_ID` variable. If you want more information visit [ERR6335730](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR6335730) page to inspect.

The tool `fastq-dump` is used to download reads and store it in two FASTQ files. One for R1 and other for R2 reads.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L5-L7&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

### Reference sequence

Usually reference sequences (or genomes) are stored in the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).

In the code bellow, we download it with `wget` program, writing it as genome.fa.gz inside a directory. Then we decompress it with gunzip and finally we create a required index file (.fai) using samtools.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L9-L17&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

## Transforming reads

Here we'll transform our R1 and R2 files into two different formats, just to illustrate our diversity of file formats.

### The unaligned BAM format

In this step we'll convert R1 and R2 reads into a single [SAM-like](https://samtools.github.io/hts-specs/SAMv1.pdf) file named unaligned BAM (uBAM). The major advantage I see in this format, when comparing to the FASTQs, is that one can store metadata in its header, so you can have informations like what was the sequencing run name, the library, the date of the experiment, the facility it was sequenced and so on.

This is the format used by [Broad Institute in their workflows](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README/#input-requirements-and-expectations).

Here is how to have it:

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L20-L38&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

### The interleaved FASTQ format

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L40-L46&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>