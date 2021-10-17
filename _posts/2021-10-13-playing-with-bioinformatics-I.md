---
title: Playing with bioinformatics I
author: Lucas M. Taniguti
comments: true
date: 2021-10-15 11:33:00 +0800
categories: [Bioinformatics, File formats]
tags: [BLAST]
---

Bioinformatics is known to have plenty of file formats. Every time a discussion about file formats starts, I remember [this comic](http://www.niso.org/sites/default/files/inline-images/How%20Standards%20Proliferate%20%281%29.png).

I'll guide you through some of these formats in this post, always paying attention to our precious Illumina reads contents. Here is our plan:

1. Download public data to play with. We'll use the sars-cov-2 genome, and Illumina reads from a whole-genome amplicon experiment.
2. Convert FASTQ files to unaligned BAM (uBAM), just to show how it works.
3. Transform uBAM into a single interleave FASTQ file.
4. Align FASTQs files (interleaved and paired) on the sars-cov-2 genome, creating a BAM file.
5. Convert BAM to CRAM and also to a lossy CRAM
6. Convert back from CRAM to pairs of FASTQs

You can find the shell script we use here in this [GitHub repository](https://github.com/lmtani/lmtani.github.io/blob/wip-illumina-format-patterns/_code/playing-with-bioinformatics-files/playing-with-bioinformatics-files.sh).

## Tools

I'm assuming you are in a Unix-like system (ubuntu, fedora, WSL2, macOS, etc.)

- [Picard](https://broadinstitute.github.io/picard/) - To convert FASTQ ⇄ uBAM and SAM ⇄ BAM ⇄ CRAM
- [Samtools](http://www.htslib.org/) - To create .fai index
- [BWA](https://github.com/lh3/bwa) - To align reads to genome
- [Crumble](https://www.sanger.ac.uk/tool/crumble/) - To try lossy compression on alignment file

> See [miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Docker](https://www.docker.com/) to easily get your environment ready. I plan to make a post about it soon.

## Download public dataset

Firstly, go to a directory you like. For example, `/tmp/bioinfo` if you don't want to persist files from this post after reboot your system.

We'll use two repositories to retrieve some data for our experiment:
- [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra): stores raw data from NGS sequencing experiments.
- FTP from [NCBI/genomes](https://www.ncbi.nlm.nih.gov/genome): organizes information on genomes.

### Raw sequencing data

SRA organizes its data using unique identifiers, like this one under our `GENOME_SRA_ID` variable. If you want more information, visit [ERR6335730](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR6335730) page to inspect.

The tool `fastq-dump` is used to download reads and store them in two FASTQ files. One for R1 and the other for R2 reads.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L5-L7&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

### Reference sequence

Usually, reference sequences (or genomes) are stored in the [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).

In the code below, we download it with `wget` program, writing it as genome.fa.gz inside a directory. Then we decompress it with gunzip, and finally, we create a required index file (.fai) using samtools.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L9-L17&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

## Transforming reads

Here we'll transform our R1 and R2 files into two different formats, just to illustrate our diversity of file formats.

### The unaligned BAM format

In this step, we'll convert R1 and R2 reads into a single [SAM-like](https://samtools.github.io/hts-specs/SAMv1.pdf) file named unaligned BAM (uBAM). The significant advantage I see in this format is that one can store metadata in its header, so you can have pieces of information like the sequencing run name, the library, the date of the experiment, the facility it was sequenced so on.

This is the format used by [Broad Institute in their workflows](https://broadinstitute.github.io/warp/docs/Pipelines/Whole_Genome_Germline_Single_Sample_Pipeline/README/#input-requirements-and-expectations).

Here is how to have it:

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L20-L38&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

### The interleaved FASTQ format

This one is not as popular as the others, but I'll mention it here just in case you see it out there. The two reads (R1 and R2) are in the same file and marked as /1 and /2 to inform their origin.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L40-L46&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

## Align to reference sequence

This step usually is the most intense in genotyping workflows. This step is responsible for identifying where is the best region to align each of your NGS reads. There is a lot of complexity in this step, but it's out of the scope of this post. For now, you can consider the best region as the interval where your sequencing read and the reference sequence are more similar.

For comparison, we will align the pair of FASTQs (R1 and R2) and also the interleaved one.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L52-L61&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

It creates two alignment files named BAM. Remember that whenever you see SAM/BAM/CRAM files, they essentially store the same information: how your reads align to the reference sequence. SAM is not compressed, and CRAM is smaller in file size. Note that you'll need to use your reference sequence to restore information from CRAM files.

Now that we have our BAM files, we can calculate some statistics using Picard. The goal in using this is to see that both files contain essentially the same data. Minor differences are expected.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L63-L68&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

See that both stats have very similar metrics. From now on, we will use only the BAM obtained from the alignment of the interleaved FASTQ.

```sh
# zgrep to select only some lines and columns
zgrep -E "^PAIR|^CATEGORY" algn-from-pairs.bam.stats | cut -f 1,2,3,4,5,6,7,8 | column -t
CATEGORY  TOTAL_READS  PF_READS  PCT_PF_READS  PF_NOISE_READS  PF_READS_ALIGNED  PCT_PF_READS_ALIGNED  PF_ALIGNED_BASES
PAIR      200000       200000    1             0               199302            0.99651               21604991

zgrep -E "^PAIR|^CATEGORY" algn-from-interleaved.bam.stats | cut -f 1,2,3,4,5,6,7,8 | column -t
CATEGORY  TOTAL_READS  PF_READS  PCT_PF_READS  PF_NOISE_READS  PF_READS_ALIGNED  PCT_PF_READS_ALIGNED  PF_ALIGNED_BASES
PAIR      200000       200000    1             0               199302            0.99651               21604976
```

## Storing as CRAM file

CRAM is the most compressed alignment format, and because of this, it's usually the adopted form for long-term storage. From what I see in human genome/exome data, it's about 1/3 of the equivalent BAM file.

To create it:

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L70-L73&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

> Remember: you will always need the reference sequence (FASTA) to decompress this format.

### Example of lossy compression

Because of the amount of data that today's sequencers can produce, storage is a challenge. You can use local or cloud-based systems, but the problem remains the same.

Crumble is a tool to perform controlled loss of information to improve the compression of alignments. Here we'll execute it with the most "aggressive" parameters to show how much smaller an alignment could be if we throw away some of its information.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L75-L80&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

### Alignment file size comparison

- BAM: 16M
- CRAM: 9.1M
- Lossy CRAM: 3.4M

## Recover raw reads from alignment

Remember that you aligned your reads into a reference sequence, and this alignment is stored as SAM/BAM/CRAM. So, if you don't have used any destructive parameter, you could retrieve the same read information from the alignment. This is important because you could delete the original sequences or put them in a coldline storage system if you prefer.

> There are a lot of cloud storage categories where it's cheap to store, but you'll pay more if you need to retrieve the data.

Here we compare the content of one read. Headers are different, and because of this, we need to compare the sequence and quality only.

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Flmtani%2Flmtani.github.io%2Fblob%2Fwip-illumina-format-patterns%2F_code%2Fplaying-with-bioinformatics-files%2Fplaying-with-bioinformatics-files.sh%23L82-L90&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>

If you calculate the md5sum hash for these two files, you'll see that they have the same content. You can also inspect it manually by opening in any text editor or using the `cat` command.

```sh
# Original
zgrep -A 3 -w "@ERR6335730.9" ERR6335730_1.fastq.gz | grep -E -v "^\@|^\+"
GAACTCACTTTCCATCCAACTTTTGTTGTTTTTGTGGTTAGAAGTAACACCCAAAAATGGATCATTACAAAATTGAAATTCACAGACTTTAATAACAACATTAGTAGCTGTCTCTTATACACATCTCCGAGCCCACGAGACGAGAATGGT
>3AAAFFFFFFFFFGD4FFGGGHHGHHGHHHHGGHH2GFGFCHDGGFGHGGGGHEEGFFHHHHEGHBHHHHGHHHFGBGHHEFGHFHHGHHGHFHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGEGGGGGGEGGFGCDDHHFH

# Recovered from CRAM
zgrep -A 3 -w "@ERR6335730.9" STEP_5/restored_R1.fq.gz | grep -E -v "^\@|^\+"
GAACTCACTTTCCATCCAACTTTTGTTGTTTTTGTGGTTAGAAGTAACACCCAAAAATGGATCATTACAAAATTGAAATTCACAGACTTTAATAACAACATTAGTAGCTGTCTCTTATACACATCTCCGAGCCCACGAGACGAGAATGGT
>3AAAFFFFFFFFFGD4FFGGGHHGHHGHHHHGGHH2GFGFCHDGGFGHGGGGHEEGFFHHHHEGHBHHHHGHHHFGBGHHEFGHFHHGHHGHFHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGGGEGGGGGGEGGFGCDDHHFH

# If you want to know how it stands after crumble
# you could revert its CRAM to FASTQ and inspect.
zgrep -A 3 -w "@ERR6335730.9" from-crumble_1.fq.gz | grep -E -v "^\@|^\+"
GAACTCACTTTCCATCCAACTTTTGTTGTTTTTGTGGTTAGAAGTAACACCCAAAAATGGATCATTACAAAATTGAAATTCACAGACTTTAATAACAACATTAGTAGCTGTCTCTTATACACATCTCCGAGCCCACGAGACGAGAATGGT
88DDDDDDDDDDDDD<<GGGGGGGGGGGGGGGGGGG2EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
```

That's it. Hope you find it helpful to get familiar with some common bioinformatics steps and file formats. Let me know if you have any questions or suggestions in the comments below.
