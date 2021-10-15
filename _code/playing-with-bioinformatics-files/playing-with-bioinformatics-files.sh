#!/bin/env bash

set -eou pipefail

# Data source
GENOME_SRA_ID=ERR6335730
GENOME_FASTA=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz
SAMPLE_NAME=ERS6670931
LIBRARY_NAME=MISEQ_WGA_SARSCOV2

# Obtain some illumina reads
fastq-dump --gzip -X 100000 --split-files "$GENOME_SRA_ID"

#1. Convert FASTQ to uBAM
mkdir STEP_1
R1="${GENOME_SRA_ID}_1.fastq.gz"
R2="${GENOME_SRA_ID}_2.fastq.gz"
OUTPUT_BASENAME=example

picard FastqToSam \
    FASTQ=$R1 \
    FASTQ2=$R2 \
    OUTPUT="STEP_1/$OUTPUT_BASENAME.unaligned.bam" \
    READ_GROUP_NAME=LIB_A \
    SAMPLE_NAME=$SAMPLE_NAME \
    LIBRARY_NAME=$LIBRARY_NAME \
    PLATFORM_UNIT=SERIAL_NUMBER_OF_SEQUENCER \
    PLATFORM=illumina \
    SEQUENCING_CENTER=UNKNOWN \
    RUN_DATE=2019-01-20T00:00:00-0400

#2. Convert uBAM to interleaved FASTQ
mkdir STEP_2
picard SamToFastq \
        INPUT="STEP_1/$OUTPUT_BASENAME.unaligned.bam" \
        FASTQ="STEP_2/$OUTPUT_BASENAME.interleave.fq.gz" \
        INTERLEAVE=true \
        NON_PF=true

#3. Align interleaved FASTQ into genome
#3.1 Download genome
mkdir reference_genome
wget "$GENOME_FASTA" \
    -O reference_genome/genome.fa.gz
gunzip reference_genome/genome.fa.gz

REFERENCE_GENOME=reference_genome/genome.fa

#3.2 Create reference index
bwa index "$REFERENCE_GENOME"

#3.3 Align illumina reads
mkdir STEP_3
bwa mem "$REFERENCE_GENOME" $R1 $R2 2>STEP_3/pair.stdout | samtools sort -o STEP_3/algn-from-pairs.sam
bwa mem -p "$REFERENCE_GENOME" "STEP_2/$OUTPUT_BASENAME.interleave.fq.gz" 2>STEP_3/interleaved.stdout | samtools sort -o STEP_3/algn-from-interleaved.sam

#3.4 Check you have exactly the same alignment statistics
samtools flagstats STEP_3/algn-from-pairs.sam > STEP_3/algn-from-pairs.stats
samtools flagstats STEP_3/algn-from-interleaved.sam > STEP_3/algn-from-interleaved.stats

md5sum STEP_3/algn-from-pairs.stats STEP_3/algn-from-interleaved.stats

#4. Convert SAM to BAM and BAM to CRAM
mkdir STEP_4
samtools view -O BAM -o STEP_4/algn-from-pairs.bam STEP_3/algn-from-pairs.sam
samtools index STEP_4/algn-from-pairs.bam

samtools view -O CRAM -T "$REFERENCE_GENOME" -o STEP_4/algn-from-pairs.cram STEP_4/algn-from-pairs.bam
samtools index STEP_4/algn-from-pairs.cram

#4.1 Do lossy compression on CRAM
crumble -z -v -9 \
    -O cram,reference="$REFERENCE_GENOME",lossy_names,seqs_per_slice=10000,nthreads=4 \
    STEP_4/algn-from-pairs.bam \
    STEP_4/alignment_nv9_lossnames.cram
