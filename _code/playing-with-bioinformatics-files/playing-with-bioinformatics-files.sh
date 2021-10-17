#!/bin/env bash

set -eoux pipefail

# Download Illumina paired-end reads
GENOME_SRA_ID=ERR6335730  # Sars-cov-2
fastq-dump --gzip -X 100000 --split-files "$GENOME_SRA_ID"

# Download reference sequence
mkdir reference_genome
GENOME_FASTA=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/858/895/GCA_009858895.3_ASM985889v3/GCA_009858895.3_ASM985889v3_genomic.fna.gz
REFERENCE_GENOME=reference_genome/genome.fa
wget "$GENOME_FASTA" \
    -O $REFERENCE_GENOME.gz

gunzip $REFERENCE_GENOME.gz
samtools faidx $REFERENCE_GENOME  # Create required index file


# 1. Convert FASTQ to uBAM
mkdir STEP_1
R1="${GENOME_SRA_ID}_1.fastq.gz"
R2="${GENOME_SRA_ID}_2.fastq.gz"
SAMPLE_NAME=ERS6670931
LIBRARY_NAME=MISEQ_WGA_SARSCOV2
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

# 2. Convert uBAM to interleaved FASTQ
mkdir STEP_2
picard SamToFastq \
        INPUT="STEP_1/$OUTPUT_BASENAME.unaligned.bam" \
        FASTQ="STEP_2/$OUTPUT_BASENAME.interleave.fq.gz" \
        INTERLEAVE=true \
        NON_PF=true

# 3. Align interleaved FASTQ into genome
# 3.1 Create BWA sequence index
bwa index "$REFERENCE_GENOME"

# 3.2 Align illumina reads. First from paired-fastqs and secondly from interleaved fastq
mkdir STEP_3
ALGN_BAM_PAIRS=STEP_3/algn-from-pairs.bam
ALGN_BAM_INTERLEAVE=STEP_3/algn-from-interleaved.bam

bwa mem "$REFERENCE_GENOME" $R1 $R2 \
    | picard SortSam O="$ALGN_BAM_PAIRS" I=/dev/stdin SORT_ORDER=coordinate

bwa mem -p "$REFERENCE_GENOME" "STEP_2/$OUTPUT_BASENAME.interleave.fq.gz" \
    | picard SortSam O="$ALGN_BAM_INTERLEAVE" I=/dev/stdin SORT_ORDER=coordinate

# 3.3 Check you have exactly the same alignment statistics
picard CollectAlignmentSummaryMetrics --INPUT "$ALGN_BAM_PAIRS" \
                                      --OUTPUT "$ALGN_BAM_PAIRS.stats"

picard CollectAlignmentSummaryMetrics --INPUT "$ALGN_BAM_INTERLEAVE" \
                                      --OUTPUT "$ALGN_BAM_INTERLEAVE.stats"

# 4. Convert BAM to CRAM
mkdir STEP_4
ALGN_CRAM_PAIRS=STEP_4/algn-from-pairs.cram
samtools view -O CRAM -T "$REFERENCE_GENOME" -o "$ALGN_CRAM_PAIRS" "$ALGN_BAM_PAIRS"

# 4.1 Do lossy compression on CRAM
ALGN_CRAM_LOSSY=STEP_4/alignment_nv9_lossnames.cram
crumble -z -v -9 \
    -O cram,reference="$REFERENCE_GENOME",lossy_names,seqs_per_slice=10000,nthreads=4 \
    "$ALGN_CRAM_PAIRS" \
    "$ALGN_CRAM_LOSSY"

# 5 Convert CRAM back to the original FASTQ
mkdir STEP_5
picard RevertSam I="$ALGN_CRAM_PAIRS" \
                 R="$REFERENCE_GENOME" \
                 O=/dev/stdout \
                 | picard SamToFastq I=/dev/stdin \
                                     F=STEP_5/restored_R1.fq.gz \
                                     F2=STEP_5/restored_R2.fq.gz \
                                     FU=unpaired.fq.gz

# 6 Check if recovered reads matchs the original ones
READ_NAME=@ERR6335730.9
zgrep -A 3 -w "$READ_NAME" $R1| grep -E -v "^\@|^\+" > original.txt
zgrep -A 3 -w "$READ_NAME" STEP_5/restored_R1.fq.gz| grep -E -v "^\@|^\+" > recovered.txt
