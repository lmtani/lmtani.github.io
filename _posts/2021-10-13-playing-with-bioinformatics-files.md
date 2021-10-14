---
title: Playing with bioinformatics files
author: Lucas M. Taniguti
comments: true
date: 2021-10-13 11:33:00 +0800
categories: [Bioinformatics, File formats]
tags: [BLAST]
---

## Introduction

Bioinformatics is known to have plenty of file formats. Every time a discussion about file formats starts I remember [this comic](http://www.niso.org/sites/default/files/inline-images/How%20Standards%20Proliferate%20%281%29.png).

In this post I'll guide you throught three different file format, always paying atention to our precious illumina reads contents. Here is our plan:

1. Start with the standard paired-end FASTQ format.
2. Make an uBAM file from it. This format stands for unaligned BAM.
3. Transform it into a single interleaved FASTQ file.
4. Align it to a reference genome to create a BAM file.
5. Extract reads back to paired-end FASTQ format. It shoud have the same content from our starting files.

## Requirements

1. picard - convert FASTQ <-> uBAM
1. bwa - align reads to genome
1. samtools - convert SAM <-> BAM <-> CRAM
1. crumble

<script src="https://emgithub.com/embed.js?target=https%3A%2F%2Fgithub.com%2Fepi2me-labs%2Fwf-artic%2Fblob%2Ff7c1309f1802d18ccf57ad4b861b7e94bb461bcd%2Fenvironment.yaml%23L8-L9&style=zenburn&showLineNumbers=on&showFileMeta=on&showCopy=on"></script>