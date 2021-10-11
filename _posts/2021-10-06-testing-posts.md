---
title: BLAST using Sequence Server
author: Lucas M. Taniguti
date: 2021-10-06 11:33:00 +0800
categories: [Bioinformatics, Sequence Analysis]
tags: [BLAST]
---

## Introduction

BLAST is one of the most used bioinformatics tools. It allows researchers to quickly find regions of similarity between a sequence of interest and a database with thousands of other sequences.


<h2 data-toc-skip>Common scenarios</h2>

- You've sequenced an amplicon - e.g., using sanger - and want to know if the desired region was amplified. So you can align your sequence against the whole NCBI nucleotide database to infer it.

- You just received your FASTQ files with NGS data and want to know if it contains the desired specie. It's not the best strategy, but one could copy the content of one read and align it using BLAST against the NCBI database.

- You've assembled a new genome for an organism that still does not exists in public databases. You can create a database with that species so you and your collaborators can align known sequences against it. For example, you have a gene of interest from a similar species and want to know where it is in your whole-genome assembly.


Here we will configure a local instance of [Sequence Server](https://sequenceserver.com/) and configure it to store some genome and protein sequences.

At the end of this tutorial, you'll have a BLAST server [like this one](https://antgenomes.sequenceserver.com/)


## Requirements

- Docker or Podman

In the [oficial documentation](https://sequenceserver.com/) they show how to install it without Docker, but here we'll follow with the containerized solution.


## Running Sequence Server

Sequence Server developers distribute their container images through DockerHub. So in the following code snippet, we're going to:

1. Create a directory to store your FASTA sequences
2. Download a sample genome sequence just to get it running
3. Decompress it and move it to the above-created directory
4. Run `makeblastdb` to create the database for this sequence
5. Run SequenceServer.


```bash
# 1. Create directory
mkdir my-database

# 2. Download sample genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/010/845/GCA_001010845.1_ASM101084v1/GCA_001010845.1_ASM101084v1_genomic.fna.gz

# 3. Decompress
gunzip GCA_001010845.1_ASM101084v1_genomic.fna.gz

# Move to database directory
mv GCA_001010845.1_ASM101084v1_genomic.fna my-database/

# 4. Create BLAST dabase
# Assuming you are at $HOME (ex: /home/your-user/)
docker run -it \
    -v /tmp/blast/my-database:/db wurmlab/sequenceserver:2.0.0.rc8 \
    bash -c "
        makeblastdb -dbtype nucl \
                    -title 'Sporisorium scitamineum genome' \
                    -in /db/GCA_001010845.1_ASM101084v1_genomic.fna \
                    -parse_seqids"

# 5. Run the Sequence Server
docker run --rm -it -p 4567:4567 -v /tmp/blast/my-database:/db \
  wurmlab/sequenceserver:2.0.0.rc8
```

And that is it. Now you have a BLAST Sequence Server running under you [http://127.0.0.1:4567](http://127.0.0.1:4567) address.

Every time you want to add a new sequence to your database, you can execute step 4 for it.
